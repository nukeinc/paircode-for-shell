//
// Created by 王恒毅 on 2025/9/29.
//
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
int main()
{
    std::string inputFile = "../inputfile_6_6";
    Data data = readMultipleArraysFromFile(inputFile);
    HamResult result1 = computeHamiltonians(data,outfile);
    int userInput = 0;
    std::cout << "请输入 1 以继续计算，其他键退出: ";
    userInput = 1;
    // std::cin >> userInput;
    if (userInput == 1)
    {
        std::string inputFile2 = "../inputfile_6_8";
        Data data2 = readMultipleArraysFromFile(inputFile2);
        HamResult result2 = computeHamiltonians(data2,outfile);
        int transinput=0;
        std::cout<<"计算转移反应输入1质子，2中子，输入其他退出：";
        transinput = 2;
        // std::cin>>transinput;
        if (transinput==1)
        {
            rvecall=result1.rvecall1;
            nucleus=result1.nucleus1;
            std::vector<int> trannumvec={0,1};
            std::vector<int>tranvec={};
            for (int i=0;i<trannumvec.size();i++)
            {
                tranvec.push_back(rvecall[trannumvec[i]]);
            }
            std::vector<std::vector<std::vector<double> > > transtrvec;
            transtrvec=transtr(trannumvec,data.ystrget1);
            std::vector<Eigen::MatrixXd> doublemat=
                caldoublemat(result1.allbasisp1,result1.allbasism01,result1.allbasism11,
                result1.ystrm01,result1.ystrm11,result1.changemat01,result1.changemat11,
                result1.schmitmat1,result2.allbasisp1,result2.allbasism01,result2.allbasism11,
                result2.ystrm01,result2.ystrm11,result2.changemat01,result2.changemat11,
                result2.schmitmat1,tranvec,transtrvec);
        }
        else if (transinput == 2) {
            // ---- 中子转移 ----
            rvecall = result1.rvecall2;
            nucleus = result1.nucleus2;

            std::vector<int> trannumvec = {0, 1};
            std::vector<int> tranvec;
            for (int i = 0; i < trannumvec.size(); i++) {
                tranvec.push_back(rvecall[trannumvec[i]]);
            }

            // std::vector<std::vector<std::vector<double>>> transtrvec;
            // transtrvec = transtr(trannumvec, data2.ystrget2);
            std::vector<int>rget={0,4};
            std::vector<Eigen::MatrixXd>  qstrm0={};
            std::vector<Eigen::MatrixXd> qstrm1={};
            std::vector<int> rvec={};
            std::vector<std::vector<std::vector<double>>> transtrvec={};
            outfile << "allbasisp2" << std::endl;
            transtrnocouple(rget, qstrm0, qstrm1, transtrvec,rvec);
            tranvec=rvec;

            std::vector<Eigen::MatrixXd> doublemat=
                caldoublemat(
                result1.allbasisp2, result1.allbasism02, result1.allbasism12,
                result1.ystrm02, result1.ystrm12, result1.changemat02, result1.changemat12,
                result1.schmitmat2,
                result2.allbasisp2, result2.allbasism02, result2.allbasism12,
                result2.ystrm02, result2.ystrm12, result2.changemat02, result2.changemat12,
                result2.schmitmat2,
                rvec, transtrvec
            );
            std::cout<< doublemat[0]<<std::endl;
            std::vector<Eigen::MatrixXd> resultsf=
                getdoublematresult(result1.coupleBases,result2.coupleBases,
                result1.eigenRe,result2.eigenRe,
                result1.jvecp1,result2.jvecp1,
                result1.jvecp2,result2.jvecp2,
                tranvec,doublemat,transinput,
                result1.allbasisp1,result1.allbasism01,result1.allbasism11,
                result1.ystrm01,result1.ystrm11,result1.changemat01,result1.changemat11,
                result1.schmitmat1,result2.allbasisp1,result2.allbasism01,result2.allbasism11,
                result2.ystrm01,result2.ystrm11,result2.changemat01,result2.changemat11,
                result2.schmitmat1);
            writeResultMatricesToFile("doublesfmat.txt",resultsf,
                result1.eigenRe,result2.eigenRe,rvec,transtrvec);

        }

    }

}