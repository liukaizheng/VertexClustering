#ifndef MESH_H
#define MESH_H

#include <Eigen/Dense>
#include <algorithm>
#include <vector>

namespace sigma {

class Mesh;
class MeshFacets;

typedef unsigned index_t;

const index_t NO_VERTEX = index_t(-1);
const index_t NO_EDGE = index_t(-1);
const index_t NO_FACET = index_t(-1);
const index_t NO_CORNER = index_t(-1);

class MeshSubElementsStore {
public: 
    MeshSubElementsStore(Mesh& mesh);
    virtual ~MeshSubElementsStore() {};

    index_t size() const { return nb_; }
protected:
    virtual index_t createSubElement(); 
    virtual index_t createSubElements(const index_t& n);

    virtual void resizeStore(const index_t& new_size);

    Mesh& mesh_;
    index_t nb_;      //number
};

class MeshElements {
public:
    MeshElements() {};
    virtual ~MeshElements() {};

    virtual void deleteElements(std::vector<index_t>& to_delete, bool remove_isolated_vertices = true) = 0;

    virtual void pop() = 0;

protected:
    static bool hasNoZero(const std::vector<index_t>& indices) {
        for (std::size_t i = 0; i < indices.size(); i++) {
            if (indices[i] != 0) {
                return true;
            }
        }
        return false;
    }
};

class MeshVertices:
    public MeshSubElementsStore,
    public MeshElements {

public:
    MeshVertices(Mesh& mesh);
    virtual ~MeshVertices() {};

    index_t createVertex();
    index_t createVertex(const double* coords);
    index_t createVertices(const index_t& n);
    const double* getVertex(const index_t& i) const { return &points_[i * 3]; }
    double* getVertex(const index_t& i) { return &points_[i * 3]; }

    virtual void deleteElements(std::vector<index_t>& to_delete, bool remove_isolated_vertices = true);
    void removeIsolated();
    virtual void pop();

protected:
    virtual void resizeStore(const index_t& new_size) {
        points_.resize(new_size * 3);
        MeshSubElementsStore::resizeStore(new_size); 
    }

    std::vector<double> points_;

    friend class Mesh;
};

class MeshEdges:
    public MeshSubElementsStore,
    public MeshElements {

public:
    MeshEdges(Mesh& mesh);
    virtual ~MeshEdges() {};

    index_t vertex(index_t e, index_t i) const { return edge_vertex_[e * 2 + i]; }
    void setVertex(index_t e, index_t i, index_t v) { edge_vertex_[e * 2 + i] = v; }

    index_t createEdge();
    index_t createEdges(const index_t& n);
    index_t createEdge(index_t v1, index_t v2);

    virtual void deleteElements(std::vector<index_t>& to_delete, bool remove_isolated_vertices = true);
    virtual void pop();

protected:
    void resizeStore(const index_t& new_size);

    std::vector<index_t> edge_vertex_;

    friend class Mesh;
};

class MeshFacetsStore: public MeshSubElementsStore {
public:
    MeshFacetsStore(Mesh& mesh);
    virtual ~MeshFacetsStore() {};

    index_t cornerBegin(const index_t& f) const { return facet_ptr_[f]; }
    index_t cornerEnd(const index_t& f) const { return facet_ptr_[f + 1]; }
    index_t numCorners(const index_t& f) const { return facet_ptr_[f + 1] - facet_ptr_[f]; }
    index_t corner(index_t f, index_t i) const { return facet_ptr_[f] + i; }


protected:
    index_t createSubElement();
    index_t createSubElements(const index_t& n);

    virtual void resizeStore(const index_t &new_size);

    std::vector<index_t> facet_ptr_;
};

class MeshFacetCornersStore: public MeshSubElementsStore {
public:
    MeshFacetCornersStore(Mesh& mesh);
    virtual ~MeshFacetCornersStore() {};

    index_t vertex(const index_t& i) const { return corner_vertex_[i]; }

    index_t* vertexPtr(const index_t& i)  { return &corner_vertex_[i]; }
    const index_t* vertexPtr(const index_t& i) const  { return &corner_vertex_[i]; }

    index_t adjacentFacet(const index_t i) const { return corner_adjacent_facet_[i]; }

    void setVertex(index_t c, index_t v) { corner_vertex_[c] = v; }
    void setAdjacentFacet(index_t c, index_t f) { corner_adjacent_facet_[c] = f; }

protected:
    index_t createSubElement(index_t v, index_t f = NO_FACET) {
        corner_vertex_.push_back(v);
        corner_adjacent_facet_.push_back(f);
        return MeshSubElementsStore::createSubElement();
    }

    index_t createSubElements(const index_t& n) {
        const index_t res = MeshSubElementsStore::createSubElements(n);
        corner_vertex_.resize(res, NO_VERTEX);
        corner_adjacent_facet_.resize(res, NO_FACET);
        return res;
    }

    virtual void resizeStore(const index_t& new_size);

    MeshVertices& vertices_;
    MeshFacetsStore& facets_;
    std::vector<index_t> corner_vertex_;
    std::vector<index_t> corner_adjacent_facet_;

    friend class Mesh;
    friend class MeshFacets;
};

class MeshFacets:
    public MeshFacetsStore,
    public MeshElements {

public:
    MeshFacets(Mesh& mesh);
    virtual ~MeshFacets() {}

    index_t numVertices(const index_t& f) const { return numCorners(f); }
    index_t vertex(index_t f, index_t i) const { return facet_corners_.vertex(corner(f, i)); }
    void setVertex(index_t f, index_t i, index_t v) { facet_corners_.setVertex(corner(f, i), v); }
    index_t findVertex(index_t f, index_t v) const {
        for (index_t i = 0; i < numVertices(f); i++) {
            if (vertex(f, i) == v) {
                return i;
            }
        }
        return NO_VERTEX;
    }

    index_t adjacent(index_t f, index_t i) const { return facet_corners_.adjacentFacet(corner(f, i)); }
    index_t findAdjacent(index_t f, index_t f2) const {
        for (index_t i = 0; i < numVertices(f); i++) {
            if (adjacent(f, i) == f2) {
                return i;
            }
        }
        return NO_FACET;
    }

    void setAdjacent(index_t f, index_t i, index_t f2) { facet_corners_.setAdjacentFacet(corner(f, i), f2); }
    index_t nextCornerAroundFacet(index_t f, index_t i) const { return i + 1 == cornerEnd(f) ? cornerBegin(f) : i + 1; }
    index_t prevCornerAroundFacet(index_t f, index_t i) const { return i == cornerBegin(f) ? cornerEnd(f) - 1 : i - 1; }
    index_t createFacets(index_t n_facets, index_t n_vertices_per_polygon) {
        index_t first_facet = size();
        index_t corner = facet_corners_.size();
        facet_corners_.createSubElements(n_facets * n_vertices_per_polygon);
        index_t result = createSubElements(n_facets);

        for (index_t i = first_facet; i <= first_facet + n_facets; i++) {
            facet_ptr_[i] = corner;
            corner += n_vertices_per_polygon;
        }
        return result;
    }

    index_t createTriangle(index_t v1, index_t v2, index_t v3)
    {
        facet_corners_.createSubElement(v1);
        facet_corners_.createSubElement(v2);
        facet_corners_.createSubElement(v3);
        index_t result = createSubElement();
        facet_ptr_[result + 1] = facet_corners_.size();
        return result;
    }

    index_t createTriangles(const index_t& n_triangles) { return createFacets(n_triangles, 3); }


    virtual void deleteElements(std::vector<index_t>& to_delete, bool remove_isolated_vertices = true);
    virtual void pop();

    void connect();
    void flip(const index_t& f);
    void assignTriangleMesh(const std::vector<index_t>& triangles);

protected:
    MeshVertices& vertices_;
    MeshFacetCornersStore& facet_corners_;

    friend class Mesh;
};

class Mesh {
public:
    MeshVertices vertices_;
    MeshEdges edges_;
    MeshFacets facets_;
    MeshFacetCornersStore facet_corners_;
    Mesh();
    Mesh(const double* points, index_t n_points);
    virtual ~Mesh() {};
};

index_t getConnectedComponents(const Mesh& mesh, std::vector<index_t>& components);

double meshFacetArea(const Mesh& mesh, index_t f);
void MeshFacetNormal(const Mesh& mesh, index_t f, double* normal);
void meshConnections(const Mesh& mesh, std::vector<index_t>& v2c, std::vector<index_t>& next_corner_around_vertex);

class MeshHalfedges {
public:
    struct Halfedge {
        Halfedge():
            facet(NO_FACET), corner(NO_CORNER) {}

        Halfedge(index_t f, index_t c):
            facet(f), corner(c) {}

        bool operator== (const Halfedge& rhs) const {
            return facet == rhs.facet && corner == rhs.corner;
        }

        bool operator!= (const Halfedge& rhs) const {
            return !(rhs == *this);
        }

        index_t facet;
        index_t corner;
    };

    MeshHalfedges(Mesh& mesh): mesh_(mesh) {}

    Mesh& mesh() { return mesh_; }
    const Mesh& mesh() const { return mesh_; }

    bool halfedgeIsBorder(const Halfedge& H) const {
        return mesh_.facet_corners_.adjacentFacet(H.corner) == NO_FACET;
    }

    void moveToNextAroundFacet(Halfedge& H) const {
        H.corner = mesh_.facets_.nextCornerAroundFacet(H.facet, H.corner);
    }

    void moveToPrevAroundFacet(Halfedge& H) const {
        H.corner = mesh_.facets_.prevCornerAroundFacet(H.facet, H.corner);
    }

    bool moveToNextAroundVertex(Halfedge& H) const;

    bool moveToPrevAroundVertex(Halfedge& H) const;

    void moveToNextAroundBorder(Halfedge& H) const;

    void moveToPrevAroundBorder(Halfedge& H) const;

    void moveToOpposite(Halfedge& H) const;

private:

    Mesh& mesh_;
};

}

#endif
