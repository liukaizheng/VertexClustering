#include "Octree.h"
#include <Eigen/Dense>
#include <algorithm>

namespace sigma {
Octree::Octree(const unsigned& bucket_size, const double& min_extent):
    bucket_size_(bucket_size), min_extent_(min_extent), points_(nullptr), root_(nullptr)
{
}

Octree::~Octree() {
    pool_.FreeAll();
}

void Octree::initialize(const double* points, const unsigned& n_points) {
    points_ = points;
    successors_.resize(n_points);
    for (unsigned i = 0; i < n_points; i++) {
        successors_[i] = i + 1;
    }

    Eigen::Map<const Eigen::Matrix<double, -1, 3, Eigen::RowMajor>> V(points, n_points, 3);
    auto min_corner = V.colwise().minCoeff();
    auto max_corner = V.colwise().maxCoeff();

    auto center = ((min_corner + max_corner)* 0.5).eval();
    auto max_extent = (center - min_corner).maxCoeff();
    root_ = createOctant(center.data(), max_extent, 0, n_points - 1, n_points);
}

Octant* Octree::createOctant(const double* center, double extent, unsigned start_idx, unsigned end_idx, unsigned size) {
    Octant* octant = pool_.template Allocate<Octant>();
    octant->children_.fill(nullptr);
    octant->is_leaf_ = true;
    std::copy(center, center + 3, octant->center_.begin());
    octant->extent_ = extent;
    octant->start_ = start_idx;
    octant->end_ = end_idx;
    octant->size_ = size;
    static const double factor[] = {-0.5, 0.5};
    if (size > this->bucket_size_ && extent > 2 * min_extent_) {
        octant->is_leaf_ = false;
        std::array<unsigned, 8> child_starts;
        std::array<unsigned, 8> child_ends;
        std::array<unsigned, 8> child_size;
        child_size.fill(0);
        unsigned idx = start_idx;
        for (unsigned i = 0; i < size; i++) {
            const double* p = &this->points_[idx * 3];
            unsigned morton_code = 0;
            for (int j = 0; j < 3; j++) {
                if (p[j] > center[j]) {
                    morton_code |= 1 << j;
                }
            }

            if (child_size[morton_code] == 0) {
                child_starts[morton_code] = idx;
            } else {
                this->successors_[child_ends[morton_code]] = idx;
            }
            child_size[morton_code] += 1;
            child_ends[morton_code] = idx;
            idx = this->successors_[idx];
        }

        const double child_extent = 0.5 * extent;
        bool first_time = true;
        unsigned last_child_idx = 0;
        for (unsigned i = 0; i < 8u; i++) {
            if (child_size[i] == 0) {
                continue;
            }
            auto child_center = (Eigen::Map<const Eigen::Array3d>(center) +
                Eigen::Array3d(
                    factor[(i & 1)> 0] * extent,
                    factor[(i & 2)> 0] * extent,
                    factor[(i & 4)> 0] * extent)).eval();
            octant->children_[i] = createOctant(child_center.data(), child_extent, child_starts[i], child_ends[i], child_size[i]);
            if (first_time) {
                octant->start_ = octant->children_[i]->start_;
            } else {
                successors_[octant->children_[last_child_idx]->end_] = octant->children_[i]->start_;
            }
            last_child_idx = i;
            octant->end_ = octant->children_[i]->end_;
            first_time = false;
        }
    }
    return octant;
}

bool Octree::split(Octant* octant) {
    if (octant->size_ <= 1 || octant->extent_ < min_extent_) {
        return false;
    }
    octant->is_leaf_ = false;
    std::array<unsigned, 8> child_starts;
    std::array<unsigned, 8> child_ends;
    std::array<unsigned, 8> child_size;
    child_size.fill(0);
    unsigned idx = octant->start_;
    for (unsigned i = 0; i < octant->size_; i++) {
        const double* p = &this->points_[idx * 3];
        unsigned morton_code = 0;
        for (int j = 0; j < 3; j++) {
            if (p[j] > octant->center_[j]) {
                morton_code |= 1 << j;
            }
        }

        if (child_size[morton_code] == 0) {
            child_starts[morton_code] = idx;
        } else {
            this->successors_[child_ends[morton_code]] = idx;
        }
        child_size[morton_code] += 1;
        child_ends[morton_code] = idx;
        idx = this->successors_[idx];
    }

    const double child_extent = 0.5 * octant->extent_;
    bool first_time = true;
    unsigned last_child_idx = 0;
    static const double factor[] = {-0.5, 0.5};
    for (unsigned i = 0; i < 8u; i++) {
        if (child_size[i] == 0) {
            continue;
        }
        auto child_center = (Eigen::Map<const Eigen::Array3d>(octant->center_.data()) +
            Eigen::Array3d(
                factor[(i & 1)> 0] * octant->extent_,
                factor[(i & 2)> 0] * octant->extent_,
                factor[(i & 4)> 0] * octant->extent_)).eval();

        octant->children_[i] = pool_.template Allocate<Octant>();
        octant->children_[i]->is_leaf_ = true;
        octant->children_[i]->children_.fill(nullptr);
        std::copy(child_center.data(), child_center.data() + 3, octant->children_[i]->center_.begin());
        octant->children_[i]->extent_ = child_extent;
        octant->children_[i]->start_ = child_starts[i];
        octant->children_[i]->end_ = child_ends[i];
        octant->children_[i]->size_ = child_size[i];

        if (first_time) {
            octant->start_ = octant->children_[i]->start_;
        } else {
            successors_[octant->children_[last_child_idx]->end_] = octant->children_[i]->start_;
        }
        last_child_idx = i;
        octant->end_ = octant->children_[i]->end_;
        first_time = false;
    }
    return true;
}

void Octree::radiusNeibors(const double* query, double sq_dist, std::vector<unsigned>& indices) {
    if (!root_) {
        return;
    }

    indices.clear();
    this->radiusNeibors(root_, query, std::sqrt(sq_dist), sq_dist, indices);
}

void Octree::radiusNeibors(const Octant* octant, const double* query, double radius, double sq_dist, std::vector<unsigned>& indices) {
    auto q = Eigen::Map<const Eigen::Vector3d>(query);

    auto contains = [&]() -> bool {
        auto c = Eigen::Map<const Eigen::Vector3d>(octant->center_.data());
        return ((q - c).cwiseAbs().array() + octant->extent_).matrix().squaredNorm() < sq_dist;
    };

    auto overlaps = [&](const Octant* octant) -> bool {
        auto p = (q - Eigen::Map<const Eigen::Vector3d>(octant->center_.data())).cwiseAbs().eval();
        const double max_dist = radius + octant->extent_;
        if (p[0] > max_dist || p[1] > max_dist || p[2] > max_dist) {
            return false;
        }

        int n_less_extent = 
            (p[0] < octant->extent_) + 
            (p[1] < octant->extent_) + 
            (p[2] < octant->extent_);

        if (n_less_extent > 1) {
            return true;
        }

        return (p.array() - octant->extent_).max(Eigen::Array3d(0.0)).matrix().squaredNorm() < sq_dist;
    };

    if (contains()) {
        unsigned idx = octant->start_;
        for (unsigned _ = 0; _ < octant->size_; _++) {
            indices.emplace_back(idx);
            idx = this->successors_[idx];
        }
        return;
    }

    if (octant->is_leaf_) {
        unsigned idx = octant->start_;
        for (unsigned _ = 0; _ < octant->size_; _++) {
            auto p = Eigen::Map<const Eigen::Vector3d>(&points_[idx * 3]);
            if ((p - q).squaredNorm() < sq_dist) {
                indices.emplace_back(idx);
            }
            idx = this->successors_[idx];
        }
        return;
    }

    for (unsigned i = 0; i < 8; i++) {
        if (!octant->children_[i] || !overlaps(octant->children_[i])) {
            continue;
        }
        radiusNeibors(octant->children_[i], query, radius, sq_dist, indices);
    }
}

}
