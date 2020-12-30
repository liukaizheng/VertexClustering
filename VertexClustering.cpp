#include "VertexClustering.h"
#include "Octree.h"
#include <queue>
#include <functional>
#include <unordered_set>
#include <iostream>

#define printMsg(x) \
    std::cout << #x << ": " << x << std::endl;

namespace sigma {

bool facetIsOnBorder(const Mesh& mesh, index_t f)
{
    for (index_t c = mesh.facets_.cornerBegin(f); c < mesh.facets_.cornerEnd(f); c++) {
        if (mesh.facet_corners_.adjacentFacet(c) == NO_FACET) {
            return true;
        }
    }
    return false;
}

void computeBorderDistance(
    Mesh& mesh, std::vector<unsigned>& dist, unsigned max_iter)
{
    /*
     * todo: optimization(BFS)
     * */
    dist.assign(mesh.facets_.size(), max_iter);
    for (index_t i = 0; i < mesh.facets_.size(); i++) {
        if (facetIsOnBorder(mesh, i)) {
            dist[i] = 0;
        }
    }
    for (unsigned i = 1; i < max_iter; i++) {
        for (index_t f = 0; f < mesh.facets_.size(); f++) {
            if (dist[f] == max_iter) {
                for (index_t c = mesh.facets_.cornerBegin(f); c < mesh.facets_.cornerEnd(f); c++) {
                    const index_t f2 = mesh.facet_corners_.adjacentFacet(c);
                    if (f2 != NO_FACET && dist[f2] == i - 1) {
                        dist[f2] = i;
                    }
                }
            }
        }
    }
}

int getRelativeOrientation(const Mesh& mesh, index_t f1, index_t c11, index_t f2)
{
    const index_t c12 = mesh.facets_.nextCornerAroundFacet(f1, c11);
    const index_t v11 = mesh.facet_corners_.vertex(c11);
    const index_t v12 = mesh.facet_corners_.vertex(c12);
    for (index_t c21 = mesh.facets_.cornerBegin(f2); c21 < mesh.facets_.cornerEnd(f2); c21++) {
        const index_t c22 = mesh.facets_.nextCornerAroundFacet(f2, c21);
        const index_t v21 = mesh.facet_corners_.vertex(c21);
        const index_t v22 = mesh.facet_corners_.vertex(c22);
        if (v11 == v21 && v12 == v22) {
            return -1;
        } else if (v11 == v22 && v12 == v21) {
            return 1;
        }
    }
    return 0;
}

void repairDissociate(Mesh& mesh, index_t f1, index_t f2)
{
    for (index_t c = mesh.facets_.cornerBegin(f1); c < mesh.facets_.cornerEnd(f1); c++) {
        if (mesh.facet_corners_.adjacentFacet(c) == f2) {
            mesh.facet_corners_.setAdjacentFacet(c, NO_FACET);
        }
    }
    for (index_t c = mesh.facets_.cornerBegin(f2); c < mesh.facets_.cornerEnd(f2); c++) {
        if (mesh.facet_corners_.adjacentFacet(c) == f1) {
            mesh.facet_corners_.setAdjacentFacet(c, NO_FACET);
        }
    }
}

void repairPropagateOrientation(
    Mesh& mesh, index_t f, const std::vector<bool>& visited,
    std::size_t& n_moebius, std::vector<index_t>* moebius_facets = nullptr)
{
    std::size_t n_plus = 0;
    std::size_t n_minus = 0;
    for (index_t c = mesh.facets_.cornerBegin(f); c < mesh.facets_.cornerEnd(f); c++) {
        const index_t f2 = mesh.facet_corners_.adjacentFacet(c);
        if (f2 != NO_FACET && visited[f2]) {
            int ori = getRelativeOrientation(mesh, f, c, f2);
            switch(ori) {
                case 1:
                    n_plus += 1;
                    break;
                case -1:
                    n_minus += 1;
                    break;
                case 0:
                    break;
            }
        }
    }

    if (n_plus != 0 && n_minus != 0) {
        n_moebius += 1;
        if (moebius_facets) {
            moebius_facets->resize(mesh.facets_.size(), 0);
            (*moebius_facets)[f] = 1;
            for (index_t c = mesh.facets_.cornerBegin(f); c < mesh.facets_.cornerEnd(f); c++) {
                const index_t f2 = mesh.facet_corners_.adjacentFacet(c);
                if (f2 != NO_FACET) {
                    (*moebius_facets)[f2] = 1;
                }
            }
        }

        if (n_plus > n_minus) {
            n_minus = 0;
            for (index_t c = mesh.facets_.cornerBegin(f); c < mesh.facets_.cornerEnd(f); c++) {
                const index_t f2 = mesh.facet_corners_.adjacentFacet(c);
                if (f2 != NO_FACET && visited[f2] && getRelativeOrientation(mesh, f, c, f2) < 0) {
                    repairDissociate(mesh, f, f2);
                }
            }
        } else {
            n_plus = 0;
            for (index_t c = mesh.facets_.cornerBegin(f); c < mesh.facets_.cornerEnd(f); c++) {
                const index_t f2 = mesh.facet_corners_.adjacentFacet(c);
                if (f2 != NO_FACET && visited[f2] && getRelativeOrientation(mesh, f, c, f2) > 0) {
                    repairDissociate(mesh, f, f2);
                }
            }
        }
    }

    if (n_minus != 0) {
        mesh.facets_.flip(f);
        printMsg(f);
    }
}

void meshReorient(Mesh& mesh, std::vector<index_t>* moebius_facets = nullptr)
{
    const unsigned max_iter = 5;
    std::vector<unsigned> depth;
    computeBorderDistance(mesh, depth, max_iter);
    std::vector<bool> visited(mesh.facets_.size(), false);

    std::size_t n_moebius = 0;
    index_t n_visited = 0;
    std::queue<index_t> queue;
    for (unsigned i = 0; i <= max_iter; i++) {
        for (index_t f = 0; f < mesh.facets_.size(); f++) {
            if (!visited[f] && depth[f] == max_iter - i) {
                queue.emplace(f);
                visited[f] = true;
                n_visited += 1;
                while(!queue.empty()) {
                    const index_t cf = queue.front();
                    queue.pop();
                    for (index_t c = mesh.facets_.cornerBegin(cf); c < mesh.facets_.cornerEnd(cf); c++) {
                        const index_t f2 = mesh.facet_corners_.adjacentFacet(c);
                        if (f2 != NO_FACET && !visited[f2]) {
                            visited[f2] = true;
                            n_visited += 1;
                            repairPropagateOrientation(mesh, f2, visited, n_moebius, moebius_facets);
                            queue.emplace(f2);
                        }
                    }
                }
            }
            if (n_visited == visited.size()) {
                break;
            }
        }
    }
}

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
		meshConnections(mesh, v2c_, next_corner_around_v_);
    }

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
            cc = next_corner_around_v_[cc];
        } while (cc != c);
    }

    Eigen::Map<Eigen::Matrix<double, -1, 3, Eigen::RowMajor>> V(mesh.vertices_.getVertex(0), mesh.vertices_.size(), 3);
    const double extent = (V.colwise().maxCoeff() - V.colwise().minCoeff()).maxCoeff() * 0.5 / std::pow(double(mesh.vertices_.size()) * rate_, 0.33);


    if (tree_) {
        delete tree_;
    }
    tree_ = new Octree(1, extent);
    tree_->initialize(mesh.vertices_.getVertex(0), mesh.vertices_.size());
    Octant* root = tree_->root();
    tree_->setMinExtent(1e-3);

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
        if (!tree_->split(cur)) {
            break;
        }
        heap.pop();
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

    removed.clear();

    meshReorient(mesh, &removed);
    /*for (index_t f = 0; f < mesh.facets_.size(); f++) {
        if (!connect(f)) {
            if (removed.empty()) {
                removed.assign(mesh.facets_.size(), 0);
            }
            removed[f] = 1;
        }
    }*/
    if (!removed.empty()) {
        mesh.facets_.deleteElements(removed);
    }

    repairNonmanifoldByCopy();
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

index_t findCorner(Mesh& mesh, index_t f, index_t v) {
    for (index_t c = mesh.facets_.cornerBegin(f); c < mesh.facets_.cornerEnd(f); c++) {
        if (mesh.facet_corners_.vertex(c) == v) {
            return c;
        }
    }
    return NO_CORNER;
} 

void VertexClustering::repairNonmanifoldByCopy()
{
    std::vector<bool> visited(mesh_->facet_corners_.size(), false);
    std::vector<bool> is_used(mesh_->vertices_.size(), false);

    std::vector<double> new_vertices;
    index_t n_vertices = mesh_->vertices_.size();
    for (index_t f = 0; f < mesh_->facets_.size(); f++) {
        for (index_t c = mesh_->facets_.cornerBegin(f); c < mesh_->facets_.cornerEnd(f); c++) {
            if (visited[c]) {
                continue;
            }
            const index_t v = mesh_->facet_corners_.vertex(c);
            index_t nv = v;
            if (is_used[v]) {
                nv = n_vertices;
                n_vertices += 1;
                const double* v_ptr = mesh_->vertices_.getVertex(v);
                new_vertices.emplace_back(v_ptr[0]);
                new_vertices.emplace_back(v_ptr[1]);
                new_vertices.emplace_back(v_ptr[2]);
            } else {
                is_used[v] = true;
            }

            index_t cf = f;
            index_t cc = c;
            do {
                visited[cc] = true;
                mesh_->facet_corners_.setVertex(cc, nv);
                cf = mesh_->facet_corners_.adjacentFacet(cc);
                if (cf == NO_FACET || cf == f) {
                    break;
                }
                cc = findCorner(*mesh_, cf, v);
            } while (true);

            if (cf == NO_FACET) {
                cf = f;
                cc = c;
                do {
                    cc = mesh_->facets_.prevCornerAroundFacet(cf, cc);
                    cf = mesh_->facet_corners_.adjacentFacet(cc);
                    if (cf == NO_FACET) {
                        break;
                    }
                    cc = findCorner(*mesh_, cf, v);
                    visited[cc] = true;
                    mesh_->facet_corners_.setVertex(cc, nv);
                } while (true);
            }
        }
    }
    if (new_vertices.empty()) {
        return;
    }
    index_t first_v = mesh_->vertices_.createVertices(index_t(new_vertices.size() / 3));
    double* tar = mesh_->vertices_.getVertex(first_v);
    for (index_t i = 0; i < new_vertices.size(); i++) {
        tar[i] = new_vertices[i];
    }
}

}
