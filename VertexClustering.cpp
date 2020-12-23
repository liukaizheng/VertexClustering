#include "VertexClustering.h"
#include "Octree.h"
#include <queue>
#include <functional>
#include <unordered_set>
#include <iostream>

#define printMsg(x) \
    std::cout << #x << ": " << x << std::endl;

namespace sigma {

int numNodes(const Octant* octant) {
    if (!octant) {
        return 0;
    }
    int res = 1;
    for (unsigned i = 0; i < 8u; i++) {
        res += numNodes(octant->children_[i]);
    }
    return res;
}
VertexClustering::VertexClustering(const double &rate, const unsigned& min_n_vertices):
    rate_(rate), min_n_vertices_(min_n_vertices), mesh_(nullptr), tree_(nullptr) 
{
}

VertexClustering::~VertexClustering()
{
    if (tree_) {
        delete tree_;
    }
}

void VertexClustering::initialize(Mesh& mesh)
{
    mesh_ = &mesh;
    meshConnections(mesh, v2c_, next_corner_around_v_);
    std::vector<index_t> removed;
    for (index_t f = 0; f < mesh.facets_.size(); f++) {
        if (!connect(f)) {
            if (removed.empty()) {
                removed.assign(mesh.facets_.size(), 0);
            }
            removed[f] = 1;
        }
    }

    if (!removed.empty()) {
        mesh.facets_.deleteElements(removed, false);
    }

    meshConnections(mesh, v2c_, next_corner_around_v_);


    normals_.resize(mesh.facets_.size() * 3);
    for (index_t f = 0; f < mesh.facets_.size(); f++) {
        MeshFacetNormal(mesh, f, &normals_[f * 3]);
    }

    vertex_weights_.assign(mesh.vertices_.size(), 0.0);
    for (index_t v = 0; v < mesh.vertices_.size(); v++) {

        index_t c = v2c_[v];
        if (c == NO_CORNER) {
            continue;
        }

        if (vertexIsOnBorder(v)) {
            vertex_weights_[v] = 1.0;
            continue;
        }

        index_t cc = c;
        do {
            const index_t f = cc / 3;
            const auto n = Eigen::Map<const Eigen::Vector3d>(&normals_[f *3]);

            std::array<index_t, 2> adj_facets;
            adj_facets[0] = mesh.facet_corners_.adjacentFacet(cc);
            adj_facets[1] = mesh.facet_corners_.adjacentFacet(mesh.facets_.prevCornerAroundFacet(f, cc));

            for (index_t i = 0; i < adj_facets.size(); i++) {
                const auto cf = adj_facets[i];
                if (cf == NO_FACET) {
                    continue;
                }
                vertex_weights_[v] = std::max(vertex_weights_[v],
                    (1.0 - Eigen::Map<const Eigen::Vector3d>(&normals_[cf * 3]).dot(n)) * 0.5);
            }
        } while (cc != c);
    }

    /*for (index_t v = 0; v < mesh.vertices_.size(); v++) {
        assert(vertex_weights_[v] >= 0 && vertex_weights_[v] <= 1.0);
        printMsg(vertex_weights_[v]);
    }*/

    Eigen::Map<Eigen::Matrix<double, -1, 3, Eigen::RowMajor>> V(mesh.vertices_.getVertex(0), mesh.vertices_.size(), 3);
    const double extent = (V.colwise().maxCoeff() - V.colwise().minCoeff()).maxCoeff() * 0.5 / std::pow(double(mesh.vertices_.size()) * rate_, 0.33);


    if (tree_) {
        delete tree_;
    }
    tree_ = new Octree(1, extent);
    tree_->initialize(mesh.vertices_.getVertex(0), mesh.vertices_.size());
    Octant* root = tree_->root();

    std::queue<Octant*> Q;
    Q.emplace(root);

    std::priority_queue<Octant*, std::vector<Octant*>, std::function<bool (const Octant*, const Octant*)>> heap(
        [&](const Octant* a, const Octant* b) -> bool {
            return a->geometric_error_ < b->geometric_error_;
        });

    while (!Q.empty()) {
        Octant* cur = Q.front();
        Q.pop();
        if (cur->is_leaf_) {
            computeCellGeometricError(cur);
            heap.emplace(cur);
        } else {
            for (unsigned i = 0; i < cur->children_.size(); i++) {
                if (cur->children_[i]) {
                    Q.emplace(cur->children_[i]);
                }
            }
        }
    }

    const std::size_t n_max_vertices = std::max(min_n_vertices_, static_cast<unsigned>(mesh_->vertices_.size() * rate_));
    while (heap.size() < n_max_vertices) {
        Octant* cur = heap.top();
        if (cur->size_ == 1) {
            break;
        }

        heap.pop();


        tree_->split(cur);
        for (unsigned i = 0; i < cur->children_.size(); i++) {
            if (cur->children_[i]) {
                computeCellGeometricError(cur->children_[i]);
                heap.emplace(cur->children_[i]);
            }
        }
    }

    typedef Eigen::Matrix<index_t, -1, 1> VectorXu;
    auto vmap = VectorXu::LinSpaced(mesh.vertices_.size(), 0, mesh.vertices_.size() - 1).eval();
    while (!heap.empty()) {
        auto cur = heap.top();
        heap.pop();
        unsigned v_idx = cur->start_;
        for (unsigned _ = 0; _ < cur->size_; _++) {
            vmap[v_idx] = cur->represent_;
            v_idx = tree_->next(v_idx);
        }
    }

    for (index_t f = 0; f < mesh.facets_.size(); f++) {
        for (index_t c = mesh.facets_.cornerBegin(f); c < mesh.facets_.cornerEnd(f); c++) {
            mesh.facet_corners_.setVertex(c, vmap[mesh.facet_corners_.vertex(c)]);
        }
    }

    std::unordered_set<index_t> set;
    for (index_t f = 0; f < mesh.facets_.size(); f++) {
        for (index_t c = mesh.facets_.cornerBegin(f); c < mesh.facets_.cornerEnd(f); c++) {
            set.emplace(mesh.facet_corners_.vertex(c));
        }
    }
    printMsg(set.size());

    removed.clear();
    for (index_t f = 0; f < mesh.facets_.size(); f++) {
        const index_t* fv = mesh.facet_corners_.vertexPtr(f * 3);
        if (fv[0] == fv[1] || fv[1] == fv[2] || fv[2] == fv[0]) {
            if (removed.empty()) {
                removed.assign(mesh.facets_.size(), 0);
            }
            removed[f] = 1;
        }
    }
    mesh.facets_.deleteElements(removed);
}

bool VertexClustering::vertexIsOnBorder(const index_t v)
{
    const index_t c = v2c_[v];
    if (c == NO_CORNER) {
        return true;
    }

    index_t cc = c;

    do {
        if (mesh_->facet_corners_.adjacentFacet(cc) == NO_FACET) {
            return true;
        } else if (mesh_->facet_corners_
            .adjacentFacet(mesh_->facets_.prevCornerAroundFacet(cc / 3, cc)) == NO_FACET)
        {
            return true;
        }
        cc = next_corner_around_v_[cc];
    } while (c != cc);
    return false;
}

bool VertexClustering::getAdjacentCorner(index_t t1, index_t* adj_corner)
{
    for (index_t c1 = mesh_->facets_.cornerBegin(t1); c1 < mesh_->facets_.cornerEnd(t1); c1++) {
        index_t v2 = mesh_->facet_corners_.vertex(
            mesh_->facets_.nextCornerAroundFacet(t1, c1));
        *adj_corner = NO_CORNER;
        index_t c2 = next_corner_around_v_[c1];
        while (c2 != c1) {
            const index_t t2 = c2 / 3;
            index_t c3 = mesh_->facets_.prevCornerAroundFacet(t2, c2);
            index_t v3 = mesh_->facet_corners_.vertex(c3);
            if (v3 == v2) {
                if (*adj_corner == NO_CORNER) {
                    *adj_corner = c3;
                } else {
                    return false;
                }
            }
            c3 = mesh_->facets_.nextCornerAroundFacet(t2, c2);
            v3 = mesh_->facet_corners_.vertex(c3);
            if (v3 == v2) {
                if (*adj_corner == NO_CORNER) {
                *adj_corner = c2;
                } else {
                    return false;
                }
            }
            c2 = next_corner_around_v_[c2];
        }
        ++adj_corner;
    }
    return true;
}

void VertexClustering::connectAdjacentCorners(index_t t, index_t* adj_coners)
{
    for (index_t i = 0; i < 3; i++) {
        if (adj_coners[i] != NO_CORNER) {
            const index_t c =  mesh_->facets_.cornerBegin(t) + i;
            mesh_->facet_corners_.setAdjacentFacet(c, adj_coners[i] / 3);
            mesh_->facet_corners_.setAdjacentFacet(adj_coners[i], t);
        }
    }
}

bool VertexClustering::connect(const index_t& t) 
{
    index_t adj_corner[3];
    if (!getAdjacentCorner(t, adj_corner)) {
        return false;
    }
    connectAdjacentCorners(t, adj_corner);
    return true;
}

void VertexClustering::computeCellGeometricError(Octant* octant)
{
    octant->geometric_error_ = 0.0;
    octant->represent_ = octant->start_;
    if (octant->size_ == 1) {
        return;
    }

    unsigned v_idx = octant->start_;
    for (unsigned _ = 1; _ < octant->size_; _++) {
        v_idx = tree_->next(v_idx);
        if (vertex_weights_[v_idx] > vertex_weights_[octant->represent_]) {
            octant->represent_ = v_idx;
        }
    }

    auto getError = [&](const index_t* fv, index_t f) -> double {
        Eigen::Map<const Eigen::Vector3d> n(&normals_[f * 3]);
        Eigen::Map<const Eigen::Vector3d> v0(mesh_->vertices_.getVertex(fv[0]));
        Eigen::Map<const Eigen::Vector3d> v1(mesh_->vertices_.getVertex(fv[1]));
        Eigen::Map<const Eigen::Vector3d> v2(mesh_->vertices_.getVertex(fv[2]));
        return 1.0 - (v1 - v0).cross(v2 - v0).normalized().dot(n);
    };

    v_idx = octant->start_;
    for (unsigned _ = 0; _ < octant->size_; _++, v_idx = tree_->next(v_idx)) {
        if (v_idx == octant->represent_) {
            continue;
        }
        const index_t c = v2c_[v_idx];
        if (c == NO_CORNER) {
            continue;
        }


        index_t cc = c;
        int n_neighbor = 0;
        double sum_w = 0.0;
        do {
            const index_t f = cc / 3;
            std::array<index_t, 3> fv;
            unsigned ii = 0;
            for(index_t fc = mesh_->facets_.cornerBegin(f); fc < mesh_->facets_.cornerEnd(f); fc++) {
                if (fc != cc) {
                    fv[ii++] = mesh_->facet_corners_.vertex(fc);
                } else {
                    fv[ii++] = octant->represent_;
                }
            }
            n_neighbor += 1;
            sum_w += getError(fv.data(), f);

            cc = next_corner_around_v_[cc];
        } while (cc != c);

        if (n_neighbor != 0) {
            octant->geometric_error_ = std::max(octant->geometric_error_, sum_w / static_cast<double>(n_neighbor) * 0.5);
        }
    }

}

}
