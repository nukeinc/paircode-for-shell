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
#include "maincal.h"
#include "inputmod.h"
int main() {
    Data data = readMultipleArraysFromFile("../inputfilesingle");
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
    std::vector<std::vector<Eigen::MatrixXd>>ystrm1;
    std::vector<std::vector<Eigen::MatrixXd>>ystrm1_1;
    Eigen::MatrixXd m01;
    Eigen::MatrixXd m11;
    std::vector<int> blockform1;
    std::vector<int> jjv;
    std::map<std::pair<int, int>, int> myMap1;
    std::map<int, std::vector<double>> eigenre1;
    std::map<std::pair<int, int>, std::vector<double>> engre;
    calhamsingle(data.num1,data.mcal,data.alpha1,
        gvec,data.bej,data.energyp,data.efcstrength1,
        data.ystrget1,data.rorder1,data.efc1,
        allbasism01,allbasism11,schmitmat1,allbasisp1,hampcc1,bemematcal1,singleindex1,
        tvec1,jvecp1,parityvecp1,ystrm1,ystrm1_1,
        m01,m11,blockform1,jjv,myMap1,eigenre1,engre);
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
    std::vector<std::vector<Eigen::MatrixXd>>ystrm2;
    std::vector<std::vector<Eigen::MatrixXd>>ystrm2_1;
    Eigen::MatrixXd m1;
    Eigen::MatrixXd m1_1;
    std::vector<int> blockform2;
    std::vector<int> jjv2;
    std::map<std::pair<int, int>, int> myMap2;
    std::map<int, std::vector<double>> eigenre2;
    std::map<std::pair<int, int>, std::vector<double>> engre2;
    calhamsingle(data.num2,data.mcal,data.alpha2,
        gvec2,data.bej,data.energyn,data.efcstrength2,
        data.ystrget2,data.rorder2,
        data.efc2,allbasism02,allbasism12,schmitmat2,
        allbasisp2,hampcc2,bemematcal2,singleindex2,
        tvec2,jvecp2,parityvecp2,ystrm2,
        ystrm2_1,m1,m1_1,blockform2,jjv2,myMap2,eigenre2,engre2);
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
    outfile << "hampcc2" << std::endl;
    writeMatrixToFile(hampcc2,outfile);

    std::vector<int> trannumvec = {0, 1};
    std::vector<int> tranvec;
    for (int i = 0; i < trannumvec.size(); i++) {
        tranvec.push_back(rvecall[trannumvec[i]]);
    }

    std::vector<std::vector<std::vector<double>>> transtrvec;
    transtrvec = transtr(trannumvec, data.ystrget2);
    std::vector<Eigen::MatrixXd> doublemat=caldoublemat(allbasisp1,allbasism01,allbasism11,
        ystrm1,ystrm1_1,m01,m11,schmitmat1,allbasisp2,allbasism02,allbasism12,ystrm2,ystrm2_1,
        m1,m1_1,schmitmat2,tranvec,transtrvec);
    for (int y=0;y<doublemat.size();++y)
    {
        std::vector<std::vector<double> > matrixa2ch=eigenToNestedVector(doublemat[y]);
        for (auto& entry : engre) {
            const std::pair<int, int>& key = entry.first;
            std::vector<double>& values = entry.second;
            for (auto& entry1:engre2)
            {
                const std::pair<int, int>& key1 = entry1.first;
                std::vector<double>& values1 = entry1.second;
                double result=0;
                for (int i=0;i<values.size();++i)
                {
                    int mnum=myMap1[{key.first,i}];
                    for (int j=0;j<values1.size();++j)
                    {
                        int mnum1=myMap2[{key1.first,j}];
                        result += values[i]*values1[j]*matrixa2ch[mnum][mnum1];
                        if ( isnan(result))
                        {
                            std::cout<< matrixa2ch[mnum][mnum1]<<" "<<values[i]<<" "<<values1[j]<<std::endl;
                        }
                    }
                }

                if (!outfile.is_open())
                {
                    // 如果文件未打开，则尝试打开文件
                    outfile.open("basis_output.txt",std::ios::app);
                }
                if (outfile.is_open()) {
                    outfile <<"核子角动量"<<y
                            << "初态角动量 " << key1.first
                            << " 初态位置 " << key1.second
                            << " 末态角动量 " << key.first
                            << " 末态位置 " << key.second
                            << " 谱因子 S " << std::pow(result, 2)
                            << std::endl;
                    outfile.close(); // 关闭文件
                } else {
                    std::cerr << "无法打开文件！" << std::endl;
                }




            }
            // 操作键（std::pair<int, int>）
        }
    }



    return 0;
}