#include <iostream>
#include "Octree.h"
#include "ReadXYZ.h"
#include <Eigen/Dense>
#include <fstream>

using namespace sigma;

void writeXYZ(const std::string& name, const double* points, const std::size_t& n)
{
    std::ofstream out(name);
    for (std::size_t i = 0; i < n; i++) {
        const double* src = &points[i * 3];
        out << src[0] << " " << src[1] << " " << src[2] << " " << "\n";
    }
    out.close();
}
/*void writeMesh(const std::string& name, const Mesh& mesh)
{
    Eigen::Map<const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>> V(mesh.vertices_.getVertex(0), mesh.vertices_.size(), 3);
    Eigen::Matrix<index_t, -1, -1> F(mesh.facets_.size() ,3);
    for (Eigen::Index i = 0; i < F.rows(); i++) {
        F(i, 0) = mesh.facet_corners_.vertex(i * 3 + 0);
        F(i, 1) = mesh.facet_corners_.vertex(i * 3 + 1);
        F(i, 2) = mesh.facet_corners_.vertex(i * 3 + 2);
    }
    writeOBJ(name, V, F);
}*/

int main()
{
    std::vector<Eigen::Vector3d> in_points;
    std::vector<Eigen::Vector3d> in_normal;
    ReadXYZ("test.xyz", in_points, in_normal);
    std::vector<double> raw_points(in_points.size() * 3);
    std::vector<double> raw_normals(in_normal.size() * 3);
    for (std::size_t i = 0; i < in_points.size(); i++) {
        const auto& p = in_points[i];
        double* target = &raw_points[i * 3];
        target[0] = p[0];
        target[1] = p[1];
        target[2] = p[2];

        auto& n = in_normal[i];
        n.normalize();
        target = &raw_normals[i * 3];
        target[0] = n[0];
        target[1] = n[1];
        target[2] = n[2];

    }

    sigma::Octree tree;
    tree.initialize(raw_points.data(), in_points.size());
    return 0;
}
