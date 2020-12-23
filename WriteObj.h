#ifndef WRITEOBJ_H
#define WRITEOBJ_H
#include<vector>
#include<string>
#include<Eigen/Dense>
#include<fstream>
#include<iostream>

template <typename DerivedV, typename DerivedF>
inline bool writeOBJ(
        const std::string &file_name,
        const Eigen::MatrixBase<DerivedV> &V,
        const Eigen::MatrixBase<DerivedF> &F
        )
{
    std::ofstream stream(file_name);
    if(!stream.is_open()) {
        std::cout << "Maybe there is an error on file name\n";
        stream.close();
        return false;
    }

    for(size_t i = 0; i < V.rows(); i++) {
        stream << "v ";
        for(int j = 0; j < V.cols()-1; j++)
            stream << V(i,j) << " ";
        stream << V(i, V.cols()-1) << "\n";
    }
    for(size_t i = 0; i < F.rows(); i++) {
        stream << "f ";
        for(int j = 0; j < F.cols()-1; j++)
            stream << F(i,j)+1 << " ";
        stream << F(i, F.cols()-1)+1 << "\n";
    }
    stream.close();
    return true;
}

#endif // WRITEOBJ_H
