//
// Created by wang- on 25-7-11.
//
#define EIGEN_NO_DEBUG
#define EIGEN_UNROLLING_LIMIT 0
// #define EIGEN_USE_MKL_ALL
// #define EIGEN_USE_BLAS
// #define EIGEN_USE_LAPACKE
#include <iostream>
#include <fstream>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include "mtcal.h"
#include "qcal.h"
#include <cstdlib>
#include <lambda_lanczos/lambda_lanczos.hpp>
#include <Eigen/Dense>
#include <vector>
#include"writeout.h"
#include"bcscal.h"
#include <sstream>
#include <string>
#include <map>
#include <utility>
#include <maincal.h>
#include "inputmod.h"
int run(const std::string& inputfile, std::vector<basis> allbasisp1, ) {
    // 这里写你原来 main 的逻辑
    // 比如：处理输入、计算、保存文件等
    std::cout << "Program called with " << argc << " arguments." << std::endl;

    for (int i = 0; i < argc; i++) {
        std::cout << "argv[" << i << "] = " << argv[i] << std::endl;
    }

    // 假设你原本要处理的逻辑在这里……
    std::cout << "Running main logic...\n";

    return 0;
}
int main() {
    Data data = readMultipleArraysFromFile("../inputfile");
    std::map<int,Matrix4D>vpnmat= getvpnval(data.pnData,data.efcstrength3);
    std::vector<std::vector<Eigen::MatrixXd>> q_pi;
    std::vector<Eigen::MatrixXd> V_it;
    std::vector<std::vector<Eigen::MatrixXd>> q_nu;
    getsvdresult( vpnmat,q_pi,  // 输出: q_pi[i](j_alpha, j_beta)
        V_it,               // 输出: V_it(i, t)
        q_nu,data.efcstrength3);
    printDiagonals(V_it);


    //cal ham1
    rvecall=data.rorder1[0];
    rparityvec=data.rorder1.back();
    std::vector<std::vector<int>>rorderp=data.rorder1;
    rorderp.pop_back();
    Eigen::MatrixXd schmitmat1;
    std::vector<basis> allbasisp1;
    std::vector<std::vector<double>> hampcc1;
    std::vector<Eigen::MatrixXd> bemematcal1;
    std::vector<int> singleindex1;
    std::vector<int> tvec1;
    std::map<int,std::map<int,Eigen::MatrixXd>> qmatpicha1;
    std::vector<int> jvecp1;
    std::vector<int> parityvecp1;
    std::vector<double> gvec={data.eg1[1],data.eg1[2]};
    std::vector<basism> allbasism01;
    std::vector<basism> allbasism11;
    calham(data.num1,data.mcal,data.alpha1,
        gvec,data.bej,data.energyp,data.efcstrength1,data.ystrget1,data.rorder1,
        V_it,data.efc1,q_pi,allbasism01,allbasism11,schmitmat1,allbasisp1,hampcc1,bemematcal1,singleindex1,
        tvec1,qmatpicha1,jvecp1,parityvecp1);
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile << "allbasisp1" << std::endl;
    printBasisVectorOneLine(allbasisp1);
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile << "allbasism01" << std::endl;
    writeBasismvecToFile(allbasism01);
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile << "allbasism11" << std::endl;
    writeBasismvecToFile(allbasism11);
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile << "hampcc1" << std::endl;
    writeMatrixToFile(hampcc1,outfile);


    //cal ham2
    auto nucleustemp=nucleus;
    nucleus=nucleus2;

    rvecall=data.rorder2[0];
    rparityvec=data.rorder2.back();
    std::vector<std::vector<int>>rorderp2=data.rorder2;
    rorderp2.pop_back();
    Eigen::MatrixXd schmitmat2;
    std::vector<basis> allbasisp2;
    std::vector<std::vector<double>> hampcc2;
    std::vector<Eigen::MatrixXd> bemematcal2;
    std::vector<int> singleindex2;
    std::vector<int> tvec2;
    std::map<int,std::map<int,Eigen::MatrixXd>> qmatpicha2;
    std::vector<int> jvecp2;
    std::vector<int> parityvecp2;
    std::vector<double> gvec2={data.eg2[1],data.eg2[2]};
    std::vector<basism> allbasism02={};
    std::vector<basism> allbasism12={};
    calham(data.num2,data.mcal,data.alpha2,
        gvec2,data.bej,data.energyn,data.efcstrength2,data.ystrget2,data.rorder2,
        V_it,data.efc2,q_nu,allbasism02,allbasism12,schmitmat2,allbasisp2,hampcc2,bemematcal2,singleindex2,
        tvec2,qmatpicha2,jvecp2,parityvecp2);
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile << "allbasisp2" << std::endl;
    printBasisVectorOneLine(allbasisp2);
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile << "allbasism02" << std::endl;
    writeBasismvecToFile(allbasism02);
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile << "allbasism12" << std::endl;
    writeBasismvecToFile(allbasism12);
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile << "hampcc2" << std::endl;
    writeMatrixToFile(hampcc2,outfile);

    //get ham matrix and get result
    std::map<int, std::vector<CoupledBasis>> coupleBasesall2;
    std::map<int, std::vector<std::vector<double>>> eigenre;
    std::map<int,std::vector<double>> eigenvalue;
    std::vector<Eigen::MatrixXd> bematre;
    calcouple(allbasisp1,allbasisp2,hampcc1,hampcc2,
        jvecp1,jvecp2,parityvecp1,parityvecp2,
        bemematcal1,bemematcal2,qmatpicha1,qmatpicha2,
        data.eg1,data.eg2,V_it,data.bej,coupleBasesall2,
        eigenre,eigenvalue,bematre);

    return 0;
}