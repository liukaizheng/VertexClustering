#include "Mesh.h"
#include <stack>
#include <algorithm>

namespace sigma {

MeshSubElementsStore::MeshSubElementsStore(Mesh& mesh):
    mesh_(mesh), nb_(0)
{
}

index_t MeshSubElementsStore::createSubElement()
{
    nb_ += 1;
    return nb_ - 1;
}

index_t MeshSubElementsStore::createSubElements(const index_t& n)
{
    nb_ += n;
    return nb_ - n;
}

void MeshSubElementsStore::resizeStore(const index_t& new_size)
{
    nb_ = new_size;
}

MeshVertices::MeshVertices(Mesh& mesh):
    MeshSubElementsStore(mesh)
{
}

index_t MeshVertices::createVertex()
{
    index_t res = MeshSubElementsStore::createSubElement();
    points_.resize(nb_ * 3);
    return res;
}

index_t MeshVertices::createVertices(const index_t& n)
{
    index_t res = MeshSubElementsStore::createSubElements(n);
    points_.resize(nb_ * 3);
    return res;
}

index_t MeshVertices::createVertex(const double* coords)
{
    index_t res = this->createVertex();
    double* p = &points_[res * 3];
    for (unsigned i = 0; i < 3; i++) {
        p[i] = coords[i];
    }
    return res;
}

void MeshVertices::deleteElements(std::vector<index_t>& to_delete, bool remove_isolated_vertices)
{
    std::vector<index_t>& old2new = to_delete;
    if (hasNoZero(to_delete)) {
        index_t idx = 0;
        for (std::size_t i = 0; i < old2new.size(); i++) {
            if (old2new[i] == 0) {
                old2new[i] = idx;
                idx += 1;
            } else {
                old2new[i] = NO_VERTEX;
            }
        }
        for (index_t i = 0; i < old2new.size(); i++) {
            const index_t j = old2new[i];
            if (j != NO_FACET && i != j) {
                double* src = &points_[i * 3];
                double* tar = &points_[j * 3];
                tar[0] = src[0];
                tar[1] = src[1];
                tar[2] = src[2];
            }
        }

        resizeStore(idx);

        MeshEdges& edges = mesh_.edges_;
        for (index_t e = 0; e < edges.size(); e++) {
            for (index_t i = 0; i < 2; i++) {
                index_t v = edges.vertex(e, i);
                v = old2new[v];
                edges.setVertex(e, i, v);
            }
        }

        MeshFacetCornersStore& facet_coners = mesh_.facet_corners_;
        for (index_t c = 0; c < mesh_.facet_corners_.size(); c++) {
            index_t v = facet_coners.vertex(c);
            v = old2new[v];
            facet_coners.setVertex(c, v);
        }
    }
}

void MeshVertices::removeIsolated()
{
    std::vector<index_t> to_delete(size(), 1);
    MeshEdges& edges = mesh_.edges_;
    for (index_t e = 0; e < edges.size(); e++) {
        for (index_t i = 0; i < 2; i++) {
            index_t v = edges.vertex(e, i);
            to_delete[v] = 0;
        }
    }

    for (index_t f = 0; f < mesh_.facets_.size(); f++) {
        for (index_t c = mesh_.facets_.cornerBegin(f); c < mesh_.facets_.cornerEnd(f); c++) {
            index_t v = mesh_.facet_corners_.vertex(c);
            to_delete[v] = 0;
        }
    }
    deleteElements(to_delete);
}

void MeshVertices::pop()
{
    nb_ -= 1;
    points_.resize(nb_ * 3);
}

MeshEdges::MeshEdges(Mesh& mesh):
    MeshSubElementsStore(mesh)
{
}

index_t MeshEdges::createEdge()
{
    index_t res = MeshSubElementsStore::createSubElement();
    edge_vertex_.resize(nb_ * 2);
    return res;
}

index_t MeshEdges::createEdges(const index_t& n)
{
    index_t res = MeshSubElementsStore::createSubElements(n);
    edge_vertex_.resize(nb_ * 2);
    return res;
}

index_t MeshEdges::createEdge(index_t v1, index_t v2)
{
    const index_t res = createEdge();
    index_t* pe = &edge_vertex_[nb_ * 2];
    pe[0] = v1;
    pe[1] = v2;
    return res;
}

void MeshEdges::resizeStore(const index_t& new_size) 
{
    edge_vertex_.resize(new_size * 2, NO_VERTEX);
    MeshSubElementsStore::resizeStore(new_size);
}

void MeshEdges::deleteElements(std::vector<index_t>& to_delete, bool remove_isolated_vertices)
{
    if (!hasNoZero(to_delete)) {
        if (remove_isolated_vertices) {
            mesh_.vertices_.removeIsolated();
        }
        return;
    }

    std::vector<index_t>& old2new = to_delete;
    index_t n_edges = 0;
    for (index_t e = 0;  e < size(); e++) {
        if (old2new[e] != 0) {
            old2new[e] = NO_EDGE;
        } else {
            old2new[e] = n_edges;
            if (n_edges != e) {
                edge_vertex_[n_edges * 2] = edge_vertex_[e * 2];
                edge_vertex_[n_edges * 2 + 1] = edge_vertex_[e * 2 + 1];
            }
            n_edges += 1;
        }
    }

    resizeStore(n_edges);

    if (remove_isolated_vertices) {
        mesh_.vertices_.removeIsolated();
    }
}

void MeshEdges::pop()
{
    resizeStore(size() - 1);
}

MeshFacetsStore::MeshFacetsStore(Mesh& mesh):
    MeshSubElementsStore(mesh)
{
}

index_t MeshFacetsStore::createSubElement()
{
    facet_ptr_.push_back(NO_CORNER);
    return MeshSubElementsStore::createSubElement();
}

index_t MeshFacetsStore::createSubElements(const index_t& n)
{
    const index_t res = MeshSubElementsStore::createSubElements(n);
    facet_ptr_.resize(nb_, NO_CORNER);
    return res;
}

void MeshFacetsStore::resizeStore(const index_t& new_size)
{
    facet_ptr_.resize(new_size + 1);
    MeshSubElementsStore::resizeStore(new_size);
}

MeshFacetCornersStore::MeshFacetCornersStore(Mesh& mesh) :
    MeshSubElementsStore(mesh), vertices_(mesh.vertices_), facets_(mesh.facets_)
{
}

void MeshFacetCornersStore::resizeStore(const index_t& new_size)
{
    corner_vertex_.resize(new_size);
    corner_adjacent_facet_.resize(new_size);
    MeshSubElementsStore::resizeStore(new_size);
}

MeshFacets::MeshFacets(Mesh& mesh) :
    MeshFacetsStore(mesh), vertices_(mesh.vertices_), facet_corners_(mesh.facet_corners_)
{
}

void MeshFacets::deleteElements(std::vector<index_t>& to_delete, bool remove_isolated_vertices)
{
    if (!hasNoZero(to_delete)) {
        if (remove_isolated_vertices) {
            mesh_.vertices_.removeIsolated();
        }
        return;
    }

    std::vector<index_t>& old2new = to_delete;
    std::vector<index_t>& corner_vertex = facet_corners_.corner_vertex_;
    std::vector<index_t>& corner_adjacent_facet = facet_corners_.corner_adjacent_facet_;

    index_t n_facets = 0;
    index_t n_corners = 0;

    for (index_t f = 0; f < size(); f++) {
        if (old2new[f] != 0) {
            old2new[f] = NO_FACET;
        } else {
            old2new[f] = n_facets;
            facet_ptr_[n_facets] = n_corners;
            for (index_t c = cornerBegin(f); c < cornerEnd(f); c++) {
                if (c != n_corners) {
                    corner_vertex[n_corners] = corner_vertex[c];
                    corner_adjacent_facet[n_corners] = corner_adjacent_facet[c];
                }
                n_corners += 1;
            }
            n_facets += 1;
        }
    }
    facet_ptr_[n_facets] = n_corners;
    for (index_t c = 0; c < facet_corners_.size(); c++) {
        index_t f = corner_adjacent_facet[c];
        if (f != NO_FACET) {
            corner_adjacent_facet[c] = old2new[f];
        }
    }

    resizeStore(n_facets);
    facet_corners_.resizeStore(n_corners);

    if (remove_isolated_vertices) {
        mesh_.vertices_.removeIsolated();
    }
}

void MeshFacets::pop()
{
    index_t new_n_corner = facet_ptr_[size() - 1];
    resizeStore(size() - 1);
    facet_corners_.resizeStore(new_n_corner);
}

void MeshFacets::connect()
{
    std::vector<index_t> next_corner_around_vertex(facet_corners_.size(), NO_CORNER);
    std::vector<index_t> v2c(vertices_.size(), NO_CORNER);
    std::vector<index_t> c2f(facet_corners_.size(), NO_FACET);

    for (index_t f = 0; f < size(); f++) {
        for (index_t c = cornerBegin(f); c < cornerEnd(f); c++) {
            index_t v = facet_corners_.vertex(c);
            next_corner_around_vertex[c] = v2c[v];
            v2c[v] = c;
            c2f[c] = f;
        }
    }

    for (index_t f1 = 0; f1 < size(); f1++) {
        for (index_t c1 = cornerBegin(f1); c1 < cornerEnd(f1); c1++) {
            if (facet_corners_.adjacentFacet(c1) == NO_FACET) {
                index_t v2 = facet_corners_.vertex(nextCornerAroundFacet(f1, c1));
                for (index_t c2 = next_corner_around_vertex[c1]; c2 != NO_CORNER; c2 = next_corner_around_vertex[c2]) {
                    if (c2 != c1) {
                        index_t f2 = c2f[c2];
                        index_t c3 = prevCornerAroundFacet(f2, c2);
                        index_t v3 = facet_corners_.vertex(c3);
                        if (v3 == v2) {
                            facet_corners_.setAdjacentFacet(c1, f2);
                            facet_corners_.setAdjacentFacet(c3, f1);
                            break;
                        }
                    }
                }
            }
        }
    }
}

void MeshFacets::flip(const index_t& f) 
{
    index_t d = numVertices(f);
    std::vector<index_t> corner_vertex_index(d);
    std::vector<index_t> corner_adjacent_facet(d);

    index_t c0 = cornerBegin(f);
    for (index_t i = 0; i < d; i++) {
        corner_vertex_index[i] = facet_corners_.vertex(c0 + i);
        corner_adjacent_facet[i] = facet_corners_.adjacentFacet(c0 + i);
    }

    for (index_t i = 0; i < d; i++) {
        index_t i_v = d - 1 - i;
        index_t i_f = (i_v == 0) ? d - 1 : i_v - 1;
        facet_corners_.setVertex(c0 + i, corner_vertex_index[i_v]);
        facet_corners_.setAdjacentFacet(c0 + i, corner_adjacent_facet[i_f]);
    }
}

void MeshFacets::assignTriangleMesh(const std::vector<index_t>& triangles) 
{
    const index_t n_triangles = static_cast<index_t>(triangles.size() / 3);
    facet_ptr_.clear();
    resizeStore(n_triangles);
    for (index_t i = 0; i <= n_triangles; i++) {
        facet_ptr_[i] = i * 3;
    }
    facet_corners_.corner_vertex_ = triangles;
    facet_corners_.resizeStore(n_triangles * 3);
    facet_corners_.corner_adjacent_facet_.assign(n_triangles * 3, NO_FACET);
}

Mesh::Mesh():
    vertices_(*this),
    edges_(*this),
    facets_(*this),
    facet_corners_(*this)
{
}

Mesh::Mesh(const double* points, index_t n_points):
    vertices_(*this),
    edges_(*this),
    facets_(*this),
    facet_corners_(*this)
{
    vertices_.createVertices(n_points);
    std::copy(points, points + vertices_.points_.size(), vertices_.points_.begin());
}

index_t getConnectedComponents(const Mesh& mesh, std::vector<index_t>& comp)
{
    const static index_t NO_COMPONENT = index_t(-1);
    index_t n_comp = 0;
    comp.assign(mesh.facets_.size(), NO_COMPONENT);
    for (index_t f = 0; f < mesh.facets_.size(); f++) {
        if (comp[f] != NO_COMPONENT) {
            continue;
        }
        std::stack<index_t> S;
        S.emplace(f);
        comp[f] = n_comp;
        while (!S.empty()) {
            index_t cf = S.top();
            S.pop();
            for (index_t c = mesh.facets_.cornerBegin(cf); c < mesh.facets_.cornerEnd(cf); c++) {
                index_t adj_f = mesh.facet_corners_.adjacentFacet(c);
                if (adj_f != NO_FACET && comp[adj_f] == NO_COMPONENT) {
                    S.emplace(adj_f);
                    comp[adj_f] = n_comp;
                }
            }
        }
        n_comp += 1;
    }
    return n_comp;
}
double meshFacetArea(const Mesh& mesh, index_t f)
{
    const index_t* fi = mesh.facet_corners_.vertexPtr(f * 3);
    return (Eigen::Map<const Eigen::Vector3d>(mesh.vertices_.getVertex(fi[1])) -
            Eigen::Map<const Eigen::Vector3d>(mesh.vertices_.getVertex(fi[0])))
               .cross(Eigen::Map<const Eigen::Vector3d>(
                          mesh.vertices_.getVertex(fi[2])) -
                      Eigen::Map<const Eigen::Vector3d>(
                          mesh.vertices_.getVertex(fi[0])))
               .norm() *
           0.5;
}

void meshConnections(const Mesh& mesh, std::vector<index_t>& v2c, std::vector<index_t>& next_corner_around_vertex)
{
    v2c.assign(mesh.vertices_.size(), NO_CORNER);
    next_corner_around_vertex.resize(mesh.facet_corners_.size());
    for (index_t f = 0; f < mesh.facets_.size(); f++) {
        for (index_t c = mesh.facets_.cornerBegin(f); c < mesh.facets_.cornerEnd(f); c++) {
            const index_t v = mesh.facet_corners_.vertex(c);
            if (v2c[v] == NO_CORNER) {
                v2c[v] = c;
                next_corner_around_vertex[c] = c;
            } else {
                next_corner_around_vertex[c] = next_corner_around_vertex[v2c[v]];
                next_corner_around_vertex[v2c[v]] = c;
            }
        }
    }
}

void MeshFacetNormal(const Mesh& mesh, index_t f, double* normal)
{
    Eigen::Map<Eigen::Vector3d> n(normal);
    const index_t* fi = mesh.facet_corners_.vertexPtr(f * 3);
    n = (Eigen::Map<const Eigen::Vector3d>(mesh.vertices_.getVertex(fi[1])) -
            Eigen::Map<const Eigen::Vector3d>(mesh.vertices_.getVertex(fi[0])))
               .cross(Eigen::Map<const Eigen::Vector3d>(
                          mesh.vertices_.getVertex(fi[2])) -
                      Eigen::Map<const Eigen::Vector3d>(
                          mesh.vertices_.getVertex(fi[0]))).normalized();
    if (std::abs(n.squaredNorm() - 1.0) > 0.001) {
        n[0] = 0.0;
        n[1] = 0.0;
        n[2] = 1.0;
    }
}

bool MeshHalfedges::moveToNextAroundVertex(Halfedge& H) const 
{
    const index_t v = mesh_.facet_corners_.vertex(H.corner);
    const index_t f = mesh_.facet_corners_.adjacentFacet(H.corner);
    if (f == NO_FACET) {
        return false;
    }

    for (index_t c = mesh_.facets_.cornerBegin(f); c != mesh_.facets_.cornerEnd(f); c++) {
        index_t pc = mesh_.facets_.prevCornerAroundFacet(f, c);
        if (mesh_.facet_corners_.vertex(c) == v &&
            mesh_.facet_corners_.adjacentFacet(pc) == H.facet) 
        {
            H.corner = c;
            H.facet = f;
            return true;
        }
    }
    //not reached
    assert(false);
    return false;
}

bool MeshHalfedges::moveToPrevAroundVertex(Halfedge& H) const
{
    const index_t v = mesh_.facet_corners_.vertex(H.corner);
    const index_t pc = mesh_.facets_.prevCornerAroundFacet(H.facet, H.corner);
    const index_t f = mesh_.facet_corners_.adjacentFacet(pc);
    if (f == NO_FACET) {
        return false;
    }
    for (index_t c = mesh_.facets_.cornerBegin(f); c < mesh_.facets_.cornerEnd(f); c++) {
        if (mesh_.facet_corners_.vertex(c) == v &&
            mesh_.facet_corners_.adjacentFacet(c) == H.facet) 
        {
            H.corner = c;
            H.facet = f;
            return true;
        }
    }
    //not reached
    assert(false);
    return false;
}

void MeshHalfedges::moveToNextAroundBorder(Halfedge& H) const
{
    moveToNextAroundFacet(H);
    while (moveToNextAroundVertex(H)) {
    }
}

void MeshHalfedges::moveToPrevAroundBorder(Halfedge& H) const
{
    while (moveToPrevAroundVertex(H)) {
    }
    moveToPrevAroundFacet(H);
}

void MeshHalfedges::moveToOpposite(Halfedge& H) const
{
    const index_t v = mesh_.facet_corners_.vertex(
        mesh_.facets_.nextCornerAroundFacet(H.facet, H.corner));
    index_t f = mesh_.facet_corners_.adjacentFacet(H.corner);
    for (index_t c = mesh_.facets_.cornerBegin(f); c < mesh_.facets_.cornerEnd(f); c++) {
        if (mesh_.facet_corners_.vertex(c) == v) {
            H.facet = f;
            H.corner = c;
            return;
        }
    }
}
}
