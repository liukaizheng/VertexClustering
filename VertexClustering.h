#ifndef VERTEXCLUSTERING_H
#define VERTEXCLUSTERING_H

#include "Mesh.h"

namespace sigma
{
struct Octant;
class Octree;
class VertexClustering {
public:
    VertexClustering(const double& rate = 0.5, const unsigned& min_n_vertices = 20);
    ~VertexClustering();
    void initialize(Mesh& mesh);

private: 
    bool vertexIsOnBorder(const index_t v);
    bool getAdjacentCorner(index_t,index_t* adj_corner);
    void connectAdjacentCorners(index_t t, index_t* adj_corner);
    bool connect(const index_t& t);
    void computeCellGeometricError(Octant* octant);

    double rate_;
    unsigned min_n_vertices_;
    Mesh* mesh_;
    std::vector<index_t> v2c_;
    std::vector<index_t> next_corner_around_v_;
    std::vector<double> normals_;
    std::vector<double> vertex_weights_;
    Octree* tree_;
};
}

#endif
