#ifndef OCTREE_H
#define OCTREE_H

#include "MemoryPool.h"
#include <vector>
#include <array>

namespace sigma
{
struct Octant
{
    bool is_leaf_;
    std::array<double, 3> center_;
    double extent_;
    unsigned start_;
    unsigned end_;
    unsigned size_;
    std::array<Octant*, 8> children_;
    double geometric_error_;
    unsigned represent_;
};
class Octree
{
public:
    Octree(const unsigned& bucket_size_ = 1, const double& min_extent_ = 1e-6);
    ~Octree();

    void initialize(const double* pts, const unsigned& n_points);
    void radiusNeibors(const double* query, double sq_dist, std::vector<unsigned>& indices);

    void split(Octant* octant);

    Octant* root() { return root_; }
    const Octant* root() const { return root_; }

    unsigned next(const unsigned idx) {
        return successors_[idx];
    }

protected:
    Octree(Octree&);
    Octree& operator=(const Octree&);

    Octant* createOctant(const double* center, double extent, unsigned start_idx, unsigned end_idx, unsigned size);
    void radiusNeibors(const Octant* octant, const double* query, double radius, double sq_dist, std::vector<unsigned>& indices);

    unsigned bucket_size_;
    double min_extent_;
    const double* points_;
    std::vector<unsigned> successors_;
    Octant* root_;
    PooledAllocator pool_;

};


}

#endif
