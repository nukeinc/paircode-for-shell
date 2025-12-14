//
// Created by wang- on 24-12-15.
//

#ifndef GETBASIS_H
#define GETBASIS_H

#include <iostream>
#include <vector>
#include "WignerSymbol.hpp"
#include <cmath>
#include <algorithm>
#include "global.h"
#include <unordered_set>
#include <functional>
#include <set>
#include <tuple> // 用于 std::tie 比较
#include <map>

struct CoupledBasis {
    int ni; // neutronBasis 在 basisvector 中的索引
    int pi;  // protonBasis 在 basisvector 中的索引
    int J;            // 耦合后的总角动量

    // 按 J -> ni -> pi 的顺序排序
    bool operator<(const CoupledBasis& other) const {
        if (J != other.J) {
            return J < other.J;
        }
        if (ni != other.ni) {
            return ni < other.ni;
        }
        return pi < other.pi;
    }
};
// 重载 << 操作符以输出 CoupledBasis
std::ostream& operator<<(std::ostream& os, const CoupledBasis& basis) {
    os << "{ni=" << basis.ni << ", pi=" << basis.pi << ", J=" << basis.J << "}";
    return os;
}

// 重载 << 操作符以输出 vector<CoupledBasis>
std::ostream& operator<<(std::ostream& os, const std::vector<CoupledBasis>& vec) {
    os << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        if (i != 0) os << ", ";
        os << vec[i]<< ",\n";
    }
    os << "]";
    return os;
}

// 重载 << 操作符以输出 map<int, vector<CoupledBasis>>
std::ostream& operator<<(std::ostream& os, const std::map<int, std::vector<CoupledBasis>>& map) {
    os << "{\n";
    for (const auto& pair : map) {
        os << "  " << pair.first << ": " << pair.second << ",\n";
    }
    os << "}";
    return os;
}
#include"basismaker.cpp"
std::vector<std::pair<std::vector<int>, std::vector<int>>> generateValidPairs(int m, const std::vector<std::vector<int>>& r, int maxSum);
std::vector<basis> calculateBasisWithRange(int m,
                                        const std::vector<std::pair<std::vector<int>, std::vector<int>>>& validPairs) ;
void getystrbasis(
    std::vector<std::vector<std::vector<double>>>& ystr1,
    const std::vector<std::vector<double>>& ystrget,
    const std::vector<int>& ro1);
std::vector<std::vector<double>> gramschmidt(const std::vector<basis>& inputBases,
    std::vector<std::vector<std::vector<std::vector<double>>>>ystrall);
void calculateystr(
    std::vector<std::vector<std::vector<std::vector<double>>>>&ystrall,
    const std::vector<std::vector<double>>& ystrget,
    std::vector<basis> allbasis);
double coupleoverlap(const std::vector<std::vector<double>> overlapmatn,
    const std::vector<std::vector<double>> overlapmatp,
    CoupledBasis cpball1,CoupledBasis cpball2);
std::vector<std::vector<double>> gramSchmidtInOverlap(
std::vector<std::vector<double>>& cpomat,
std::vector<basis>& allbasisp,
std::vector<std::vector<std::vector<std::vector<double>>>>& ystrallp);
std::vector<std::vector<double>> gramSchmidtInOverlap(
    std::vector<std::vector<double>>& cpomat,
    std::vector<basis>& allbasisp,
    std::vector<std::vector<std::vector<std::vector<double>>>>& ystrallp,
    const std::vector<std::pair<int, int>>& blocks);
std::map<int, std::vector<CoupledBasis>> coupleBases(
    const std::vector<basis>& neutronBases,
    const std::vector<basis>& protonBases,std::vector<std::vector<double> > overlapp,std::vector<std::vector<double> > overlapn
);
std::vector<std::vector<double>> cpomatcal(const std::vector<std::vector<double>> overlapmatn,
    const std::vector<std::vector<double>> overlapmatp,std::vector<CoupledBasis>& inputBases);
std::vector<std::vector<double>> gramschmidt(const std::vector<CoupledBasis> inputBases,
    const std::vector<std::vector<double>> cpomat);
std::vector<std::vector<double>> invertMatrix(const std::vector<std::vector<double>>& matrix);
std::vector<std::vector<double>> multiplyMatrices(
    const std::vector<std::vector<double>>& A,
    const std::vector<std::vector<double>>& B);
std::vector<std::vector<double>> transposeMatrix(const std::vector<std::vector<double>>& matrix);
std::vector<std::vector<double>> transposeMatrix(const std::vector<std::vector<double>>& matrix);
std::vector<basism> calculateBasisform(int nucnum,std::vector<std::vector<int>> rorderrp);
Eigen::MatrixXd ComputeTransformationMatrix(
    const std::vector<basism>& allbasisc ,const std::vector<basis>& allbasis  // 核子对类型
);
std::vector<std::vector<Eigen::MatrixXd>> computeystrm(const std::vector<basism>& allbasisc,const std::vector<std::vector<double> >& ystrgetp);
SparseMatrix4D calq(const basism& bas1,const basism& bas2,const std::vector<Eigen::MatrixXd>&ystrm1,const std::vector<Eigen::MatrixXd>&ystrm2);
// Matrix4D calg(const basism& bas1,const basism& bas2,const std::vector<Eigen::MatrixXd>&ystrm1,const std::vector<Eigen::MatrixXd>&ystrm2);
std::vector<std::vector<Matrix4D>> calgall(const std::vector<basism>& allbasisc,const std::vector<std::vector<Eigen::MatrixXd>>& ystrm);

Eigen::Tensor<double, 4> matrix4d_to_tensor(const Matrix4D& mat4d);
void printNonZeroElements(const Matrix4D& tensor);
void printAsymmetricElements(const Matrix4D& tensor);
void checkInequality(const std::vector<std::vector<Matrix4D>>& galltest);
void checkInequalityt(const Matrix4D& fall) ;
Eigen::MatrixXd ham1(const std::vector<basism>& bas,const std::vector<std::vector<Eigen::MatrixXd>>&ystrm,
    std::vector<std::map<int, Matrix4D>> buildVValue1,
    std::vector<double>strength);
std::vector<std::vector<double>> readAndConvertFile(const std::string& filename);
std::vector<basis> removeNegativeParity(const std::vector<basis>& bases);
#endif //GETBASIS_H