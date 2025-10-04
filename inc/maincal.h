//
// Created by wang- on 25-7-29.
//

#ifndef MAINCAL_H
#define MAINCAL_H

#include <iostream>
#include <fstream>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include "mtcal.h"
#include "qcal.h"
#include <cstdlib>
#include "lambda_lanczos/lambda_lanczos.hpp"
#include "Eigen/Dense"
#include <vector>
#include"writeout.h"
#include"bcscal.h"
#include <sstream>
#include <string>
#include <map>
#include <utility>
#include "getbasis.h"

using MatrixIndex = std::pair<int, int>;
using std::setprecision;
using lambda_lanczos::LambdaLanczos;
struct DataRow {
    int t;       // 第1列
    int a, b, c, d; // 第2-5列
    int J;       // 第6列
    double value; // 第7列
};



#include"maincal.cpp"
std::vector<std::map<int, Matrix4D>> buildVValue(const std::vector<DataRow>& data,int num);

std::vector<std::map<int, Matrix4D>> processFile(const std::string& filename);

std::vector<std::map<int, Matrix4D>> buildVValuepn(const std::vector<DataRow>& data);

std::vector<std::map<int, Matrix4D>> processFilepn(const std::string& filename);

std::map<int,Matrix4D>getvpnval(const std::vector<std::map<int, Matrix4D>>& V_value,std::vector<double> strengthvec);

void decompose_right_tensor(const Eigen::Tensor<double, 4>& right_tensor,
                          int rank,
                          std::vector<Eigen::MatrixXd>& q_pi,
                          Eigen::MatrixXd& V_it,
                          std::vector<Eigen::MatrixXd>& q_nu);

void getsvdresult(const std::map<int,Matrix4D>& vpnval,
                          std::vector<std::vector<Eigen::MatrixXd>>& q_pi,  // 输出: q_pi[i](j_alpha, j_beta)
                          std::vector<Eigen::MatrixXd>& V_it,               // 输出: V_it(i, t)
                          std::vector<std::vector<Eigen::MatrixXd>>& q_nu,const std::vector<double>& strengthvec);
void printDiagonals(const std::vector<Eigen::MatrixXd>& V_it);

std::vector<MatrixIndex> findNonZeroDiagonal(const std::vector<Eigen::MatrixXd>& V_it);

Eigen::MatrixXd transqjtoqm(Eigen::MatrixXd qj,int t, int mu);

double computeElementwiseProductSum(const Eigen::MatrixXd& mat1,
                                     const Eigen::MatrixXd& mat2);

Eigen::MatrixXd buildResultMatrix(const Eigen::MatrixXd& qmitmu,
                                 const std::vector<basism>& basm,
                                 const std::vector<std::vector<Eigen::MatrixXd>>& ystrm);

Eigen::MatrixXd buildResultMatrix_1(const Eigen::MatrixXd& qmitmu,
                                 const std::vector<basism>& basm,
                                 const std::vector<std::vector<Eigen::MatrixXd>>& ystrm,
                                 const std::vector<basism>& basm_1,
                                 const std::vector<std::vector<Eigen::MatrixXd>>& ystrm_1);

std::map<int, std::map<int, Eigen::MatrixXd>> qmatjcal(
    const std::vector<std::vector<Eigen::MatrixXd>>& q_pi,
    const Eigen::MatrixXd& transformation_matrix,
    const std::vector<Eigen::MatrixXd>& V_it,
    const std::vector<basis>& bas,
    const std::vector<basism>& basm,
    const std::vector<std::vector<Eigen::MatrixXd>>&ystrm);

std::map<int, std::map<int, Eigen::MatrixXd>> qmatjcal_1(
    const std::vector<std::vector<Eigen::MatrixXd>>& q_pi,
    const Eigen::MatrixXd& transformation_matrix,
    const Eigen::MatrixXd& transformation_matrix_1,
    const std::vector<Eigen::MatrixXd>& V_it,
    const std::vector<basis>& bas,
    const std::vector<basism>& basm,
    const std::vector<basism>& basm_1,
    const std::vector<std::vector<Eigen::MatrixXd>>&ystrm,
    const std::vector<std::vector<Eigen::MatrixXd>>&ystrm_1);

Eigen::MatrixXd bemecal(
    const Eigen::MatrixXd& qti,
    const Eigen::MatrixXd& transformation_matrix,
    const std::vector<basis>& bas,
    const std::vector<basism>& basm,
    const std::vector<std::vector<Eigen::MatrixXd>>&ystrm,
    const int& t);

Eigen::MatrixXd bemecal_1(
    const Eigen::MatrixXd& qti,
    const Eigen::MatrixXd& transformation_matrix,
    const Eigen::MatrixXd& transformation_matrix_1,
    const std::vector<basis>& bas,
    const std::vector<basism>& basm,
    const std::vector<basism>& basm_1,
    const std::vector<std::vector<Eigen::MatrixXd>>&ystrm,
    const std::vector<std::vector<Eigen::MatrixXd>>&ystrm_1,
    const int& t);

bool isZero(double value, double epsilon );

void updateMatrices(std::vector<basis>& allbasisp, std::vector<std::vector<double>>& overlapp,std::vector<std::vector<std::vector<std::vector<double>>>>& ystrallp);

std::vector<Eigen::MatrixXd> qmatcalall(
    const std::vector<Eigen::MatrixXd>& qti,
    const Eigen::MatrixXd& transformation_matrix,
    const std::vector<basis>& bas,
    const std::vector<basism>& basm,
    const std::vector<std::vector<Eigen::MatrixXd>>&ystrm,
    const std::vector<int>& t);

std::vector<Eigen::MatrixXd> qmatcalall_1(
    const std::vector<Eigen::MatrixXd>& qti,
    const Eigen::MatrixXd& transformation_matrix,
    const Eigen::MatrixXd& transformation_matrix_1,
    const std::vector<basis>& bas,
    const std::vector<basism>& basm,
    const std::vector<basism>& basm_1,
    const std::vector<std::vector<Eigen::MatrixXd>>&ystrm,
    const std::vector<std::vector<Eigen::MatrixXd>>&ystrm_1,
    const std::vector<int>& t);

void calham(const int& nuculnum,
    const int& m1cal,
    const double& alpha1,
    const std::vector<double>& gvec,
    const std::vector<int>&qorderp1,
    const std::vector<double>& energyp,
    const std::vector<double>& strength,
    const std::vector<std::vector<double> >& ystrgetp,
    const std::vector<std::vector<int>> &rorderp,
    const std::vector<Eigen::MatrixXd>& V_it,
    const std::vector<std::map<int, Matrix4D>>& buildVValue2,
    const std::vector<std::vector<Eigen::MatrixXd>>& q_pi,
    std::vector<basism>& mtest,
    std::vector<basism>& mtest_1,
    Eigen::MatrixXd& schmitmat1,
    std::vector<basis>& allbasisp,
    std::vector<std::vector<double>>& hampcc,
    std::vector<Eigen::MatrixXd>& bemematcal,
    std::vector<int>& singleindex,
    std::vector<int>& tvec,
    std::map<int,std::map<int,Eigen::MatrixXd>>& qmatpicha,
    std::vector<int>& jvecp,
    std::vector<int>& parityvecp,
    std::vector<std::vector<Eigen::MatrixXd>>& ystrm1,
    std::vector<std::vector<Eigen::MatrixXd>>& ystrm1_1,
    Eigen::MatrixXd & m1,
    Eigen::MatrixXd & m1_1
    );








#endif //MAINCAL_H
