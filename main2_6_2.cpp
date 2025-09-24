//
// Created by wang- on 2025/3/27.
//
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


// 定义一个函数来计算矩阵的特征值和特征向量
/*void computeEigenvaluesAndEigenvectors(const std::vector<std::vector<double>>& matrix_data) {
    // 将 std::vector<std::vector<double>> 转换为 Eigen::MatrixXd
    Eigen::MatrixXd A(matrix_data.size(), matrix_data[0].size());
    for (size_t i = 0; i < matrix_data.size(); ++i) {
        for (size_t j = 0; j < matrix_data[i].size(); ++j) {
            A(i, j) = matrix_data[i][j];
        }
    }

    // 使用 EigenSolver 计算特征值
    Eigen::EigenSolver<Eigen::MatrixXd> solver(A);

    // 获取特征值（复数）
    Eigen::VectorXd eigenValues = solver.eigenvalues().real();  // 取实部
    std::cout << "Eigenvalues: \n" << eigenValues << std::endl;

    // 获取特征向量（复数矩阵）
    Eigen::MatrixXd eigenVectors = solver.eigenvectors().real();  // 取实部
    std::cout << "Eigenvectors: \n" << eigenVectors << std::endl;
}*/

using std::setprecision;
using lambda_lanczos::LambdaLanczos;

template<typename T>
using vector = std::vector<T>;
// TIP To <b>Run</b> code, press <shortcut actionId="Run"/> or
// click the <icon src="AllIcons.Actions.Execute"/> icon in the gutter.
void printMatrix(const std::vector<std::vector<double>>& matrix)
{
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            std::cout << matrix[i][j] << "  ";
        }
        std::cout << std::endl;
    }
}



// 打印单个 basis 到一行
void printBasisOneLine(const basis& b)
{
    // 一行输出，字段之间用空格区分
    std::cout
        << "j=" << b.j
        << " jn=" << b.jn
        << " r=" << vecToStr(b.r)
        << " rn=" << vecToStr(b.rn)
        << " sj=" << vecToStr(b.sj)
        << " rparity=" << vecToStr(b.rparity)
        << std::endl;
}

// 打印 std::vector<basis>，每个 basis 一行
void printBasisVectorOneLine(const std::vector<basis>& basisVec)
{
    for (size_t i = 0; i < basisVec.size(); ++i) {
        printBasisOneLine(basisVec[i]);
    }
}
bool isZero(double value, double epsilon = 1e-6) {
    return std::fabs(value) < epsilon;
}

void updateMatrices(std::vector<basis>& allbasisp, std::vector<std::vector<double>>& overlapp,std::vector<std::vector<std::vector<std::vector<double>>>>& ystrallp) {
    // 遍历 overlapp 矩阵的对角线元素
    for (size_t a = 0; a < overlapp.size(); ++a) {
        // 由于我们可能会在删除后改变 overlapp 的大小，所以每次循环前需要检查大小
        if (a >= overlapp.size()) {
            break; // 如果索引 a 超出了矩阵大小，则终止循环
        }

        if (isZero(overlapp[a][a])) {
            // 如果 overlapp[a][a] 接近零，删除 allbasisp 中的对应元素
            allbasisp.erase(allbasisp.begin() + a);  // 删除 allbasisp 中第 a 个元素
            ystrallp.erase(ystrallp.begin() + a);

            // 同时更新 overlapp 矩阵，删除第 a 行和第 a 列
            for (size_t i = 0; i < overlapp.size(); ++i) {
                overlapp[i].erase(overlapp[i].begin() + a);  // 删除第 i 行的第 a 列
            }
            overlapp.erase(overlapp.begin() + a);  // 删除第 a 行

            // 由于删除了一个元素，a 需要递减，避免跳过下一个元素
            --a; // 递减 a，保持索引不跳过元素
        }
    }
    return ;
}

int main()
{
    std::cout << "precal`wigner_6j`..." << std::endl;
    // util::wigner.initialize_wigner_6j() ;  // 预计算所有 Wigner 6j 系数
    std::cout << "Precompute done.\n";
    // util::wigner.precompute_wigner_6j(); // **一次性计算所有 `wigner_6j` 值**
    std::cout << "pre ok" << std::endl;

    if (!outfile.is_open()) {
        std::cerr << "无法打开文件！" << std::endl;
        return 1;
    }
    std::vector<std::vector<double>> Zmatrixp = {
        {-0.3138929, -0.2683482, 0.0, 0.0, 0.0},
        {-0.3747607, 0.2247637, 0.0, 0.0, 0.0},
        {0.0, 0.0, -0.0738756, -0.5707341, 0.0},
        {0.0, 0.0, 0.3942333, -0.1069502, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.4006788}
    };
    nucleus={
        {0, 4, 3},
        // {0, 0, 1}
    };
    int size1=nucleus.size();
    allnucleus={
        {0, 4, 3},
        // {0, 0, 1}
    };
    std::vector<std::vector<double> > ystrgetp;
    // ystrgetp = {
    //     {0, 0, 0,  0.6},
    //     {1, 1, 0,  -0.8},
    //     {0, 0, 1, -0.7},
    //     {0, 1, 1,  0.2},
    //     {1, 0, 1, 0.2},
    //     {1, 1, 1, -0.66}
    // };
    ystrgetp = {
        {0, 0, 0,  1},
        {0, 0, 1, 1},
    };
    normalizeBySecondColumn(ystrgetp);
    std::vector<int> rorderp={0,4};
    std::vector<std::vector<int>>rorderrp={rorderp,{10,10}};
    std::vector<std::pair<std::vector<int>, std::vector<int>>> catchpair1;
    catchpair1=generateValidPairs(6,rorderrp,51);
    std::vector<basis> allbasisp;
    allbasisp=calculateBasisWithRange(6, catchpair1);
    // for (size_t i = 0; i < allbasisp.size(); ++i)
    // {
    //     if (!outfile.is_open())
    //     {
    //         // 如果文件未打开，则尝试打开文件
    //         outfile.open("basis_output.txt",std::ios::app);
    //     }
    //     outfile<<"____________________begin_____________________"<<std::endl;
    //     outfile<<"rjl1"<<std::endl;
    //     outfile<<"n1"<<i<<std::endl;
    //     printBasisOneLinefile(outfile,allbasisp[i]);
    // }
    std::vector<std::vector<double> > qstrgetp;
    std::vector<int> qorderp={4};
    qstrgetp=qstrall(qorderp,1);

    int sizeq=qorderp.size();
    std::vector<std::vector<std::vector<double>>> qstrallp(sizeq, std::vector<std::vector<double>>(size1,
    std::vector<double>(size1, 0)));

    getystrbasis(qstrallp, qstrgetp, {0});
    std::vector<std::vector<double> > astrgetp;
    std::vector<int> aorderp={0};
    astrgetp=qstrall(aorderp,1);

    int sizea=aorderp.size();
    std::vector<std::vector<std::vector<double>>> astrallp(sizea, std::vector<std::vector<double>>(size1,
    std::vector<double>(size1, 0)));
    getystrbasis(astrallp, astrgetp, {0});
    int sizep=allbasisp.size();
    int sizepr=allbasisp[0].r.size();
    std::vector<std::vector<std::vector<std::vector<double>>>> ystrallp(
        sizep,  // 第一维：大小为 sizep
        std::vector<std::vector<std::vector<double>>>(
            sizepr,  // 第二维：大小为 sizepr
            std::vector<std::vector<double>>(
                size1,  // 第三维：大小为 size1
                std::vector<double>(size1, 0.0)  // 第四维：大小为 size1，元素初始化为 0.0
            )
        )
    );
    calculateystr(ystrallp,ystrgetp,allbasisp);
    std::vector<double> energyp={-0.02825000 ,0.94409507};
    std::vector<std::vector<double>> matslip;
    std::vector<double> ktp={-0.02};
    std::vector<double> gsp={-0.13};
    printBasisVectorOneLine(allbasisp);
    std::vector<std::vector<double> > overlapp;
    overlapp=overlapmatcal(allbasisp,ystrallp);
    updateMatrices(allbasisp, overlapp,ystrallp);
    printBasisVectorOneLine(allbasisp);
    printMatrix(overlapp);
    // std::ofstream outfile("basis_output.txt");
    std::vector<std::vector<double>>schmitmatp;
    schmitmatp=gramSchmidtInOverlap(overlapp,allbasisp,ystrallp);
    printMatrix(schmitmatp);
    printMatrix(overlapp);
    printBasisVectorOneLine(allbasisp);
    std::cout << "fffff1 " << std::endl;
    // schmitmatp=Zmatrixp;
    std::vector<std::vector<double> > schmitmatpinv;
    schmitmatpinv=transposeMatrix(schmitmatp);
    printMatrix(schmitmatpinv);
    std::vector<std::vector<double> > schmitmatpinvre=multiplyMatrices(schmitmatp,schmitmatpinv);
    printMatrix(schmitmatpinvre);
    // printMatrix(overlapp);
    std::vector<std::vector<double> > overlappp;
    overlappp=multiplyMatrices(schmitmatp,overlapp);
    printMatrix(overlappp);
    overlappp=multiplyMatrices(overlappp,schmitmatpinv);
    printMatrix(overlappp);
    // sqrtOverlap(overlappp);

    auto start1 = std::chrono::high_resolution_clock::now();
    matslip=matcalsim(allbasisp,ystrallp,qorderp,aorderp,energyp,qstrallp,ktp,astrallp,gsp);
    std::vector<std::vector<double> > matslipch=multiplyMatrices(schmitmatp,matslip);
    matslipch=multiplyMatrices(matslipch,schmitmatpinv);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start1);
    std::cout << "Time taken by function: " << duration.count() << " seconds" << std::endl;
    printMatrix(matslip);
    printMatrix(matslipch);
    printBasisVectorOneLine(allbasisp);
    std::vector<std::vector<std::vector<double>>> matsliqp;
    matsliqp=matspl(allbasisp,ystrallp,{4},qstrallp);
    printMatrix(matsliqp[0]);
    std::vector<std::vector<std::vector<double>>>matsliqpch;
    for(int i=0;i<matsliqp.size();++i)
    {
        std::vector<std::vector<double>> ch1=multiplyMatrices(schmitmatp,matsliqp[i]);
        ch1=multiplyMatrices(ch1,schmitmatpinv);
        matsliqpch.push_back(ch1);
    }
    printMatrix(matsliqpch[0]);
    std::vector<int> jvecp;
    for (int i=0;i<allbasisp.size();++i)
    {
        jvecp.push_back(allbasisp[i].sj.back());
    }




    std::vector<std::vector<double>> Zmatrixn = {
        {-0.4266297, -0.2739617, 0.0, 0.0, 0.0},
        {-0.5720773, 0.2043084, 0.0, 0.0, 0.0},
        {0.0, 0.0, -0.0357049, -0.6387638, 0.0},
        {0.0, 0.0, 0.4036302, -0.0565047, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.4197558}
    };

    nucleus=
    {
        {0, 10,  9},
     {2,  6,  7},
     {2,  6,  5},
     {4,  2,  3},
     {4,  2,  1},
     {0, 12, 13}
    };
    int size2=nucleus.size();
    allnucleus={
        {0, 10,  9},
    {2,  6,  7},
    {2,  6,  5},
    {4,  2,  3},
    {4,  2,  1},
    {0, 12, 13}
    };
    std::vector<std::vector<double> > ystrgetn;
    ystrgetn = {
        {5,  5,  0, -0.00674229},
    {0,  0,  0,  0.11867957},
    {1,  1,  0,  0.97301617},
    {2,  2,  0,  0.08850505},
    {3,  3,  0,  0.16530680},
    {4,  4,  0,  0.06284695},
    {5,  5,  1, -0.09111491},
    {0,  0,  1, -0.10107346},
    {0,  1,  1,  0.03785131},
    {0,  2,  1, -0.07718849},
    {1,  0,  1, -0.03785131},
    {1,  1,  1, -0.83192200},
    {1,  2,  1, -0.05598418},
    {1,  3,  1, -0.32844468},
    {2,  0,  1, -0.07718849},
    {2,  1,  1,  0.05598418},
    {2,  2,  1, -0.07230329},
    {2,  3,  1,   0.04714734},
    {2,  4,  1, -0.07054737},
    {3,  1,  1, -0.32844468},
    {3,  2,  1, -0.04714734},
    {3,  3,  1, -0.12946939},
    {3,  4,  1, -0.08999021},
    {4,  2,  1, -0.07054737},
    {4,  3,  1,  0.08999021}
    };
    normalizeBySecondColumn(ystrgetn);
    std::vector<int> rordern={0,4};
    std::vector<std::vector<int>>rorderrn={rordern,{10,10}};
    std::vector<std::pair<std::vector<int>, std::vector<int>>> catchpair2;
    catchpair2=generateValidPairs(4,rorderrn,50);
    std::vector<basis> allbasisn;
    allbasisn=calculateBasisWithRange(  4, catchpair2);
    std::vector<int> qordern={4};
    std::vector<std::vector<double> > qstrgetn;
    qstrgetn=qstrall(qordern,1);
    int sizeqn=qordern.size();
    std::vector<std::vector<std::vector<double>>> qstralln(sizeqn, std::vector<std::vector<double>>(size2,
    std::vector<double>(size2, 0)));
    getystrbasis(qstralln, qstrgetn, {0});
    std::vector<int> aordern={0};
    std::vector<std::vector<double> > astrgetn;
    astrgetn=qstrall(aordern,1);

    int sizean=aordern.size();
    std::vector<std::vector<std::vector<double>>> astralln(sizean, std::vector<std::vector<double>>(size2,
    std::vector<double>(size2, 0)));
    getystrbasis(astralln, astrgetn, {0});
    int sizen=allbasisn.size();
    int sizenr=allbasisn[0].r.size();
    std::vector<std::vector<std::vector<std::vector<double>>>> ystralln(
        sizen,  // 第一维：大小为 sizep
        std::vector<std::vector<std::vector<double>>>(
            sizenr,  // 第二维：大小为 sizepr
            std::vector<std::vector<double>>(
                size2,  // 第三维：大小为 size1
                std::vector<double>(size2, 0.0)  // 第四维：大小为 size1，元素初始化为 0.0
            )
        )
    );
    calculateystr(ystralln,ystrgetn,allbasisn);
    std::vector<double> energyn={1.42644173    ,-0.16711269   ,1.75830984 ,   0.54762673,   1.05662673,  2.57981581};
    std::vector<std::vector<double>> matslin;
    std::vector<double> ktn={-0.02};
    std::vector<double> gsn={-0.13};
    std::vector<std::vector<double> > overlapn;
    overlapn=overlapmatcal(allbasisn,ystralln);
    printMatrix(overlapn);
    updateMatrices(allbasisn, overlapn,ystralln);
    std::vector<std::vector<double>>schmitmatn;
    schmitmatn=gramSchmidtInOverlap(overlapn,allbasisn,ystralln);
    // schmitmatn=Zmatrixn;
    std::vector<std::vector<double> > schmitmatninv;
    schmitmatninv=transposeMatrix(schmitmatn);
    matslin=matcalsim(allbasisn,ystralln,qordern,aordern,energyn,qstralln,ktn,astralln,gsn);
    printMatrix(matslin);
    std::vector<std::vector<std::vector<double>>> matsliqn;
    matsliqn=matspl(allbasisn,ystralln,qordern,qstralln);
    printMatrix(matsliqn[0]);
    std::vector<std::vector<double> > overlapnn=multiplyMatrices(schmitmatn,overlapn);
    overlapnn=multiplyMatrices(overlapnn,schmitmatninv);
    printMatrix(overlapnn);

    std::vector<std::vector<double> > matslinch=multiplyMatrices(schmitmatn,matslin);
    matslinch=multiplyMatrices(matslinch,schmitmatninv);
    printMatrix(matslinch);

    std::vector<std::vector<std::vector<double>>>matsliqnch;
    for(int i=0;i<matsliqn.size();++i)
    {
        std::vector<std::vector<double>> ch1=multiplyMatrices(schmitmatn,matsliqn[i]);
        ch1=multiplyMatrices(ch1,schmitmatninv);
        matsliqnch.push_back(ch1);
    }
    printMatrix(matsliqnch[0]);
    std::vector<int> jvecn;
    for (int i=0;i<allbasisn.size();++i)
    {
        jvecn.push_back(allbasisn[i].sj.back());
    }
    std::map<int, std::vector<CoupledBasis>> coupleBasesall2;
    coupleBasesall2= coupleBases2(jvecn,jvecp);
    double gnnn=-0.06;
    std::map<int, std::vector<double>> eigenre;
    for (const auto& [key, value] : coupleBasesall2)
    {

        int siz1=coupleBasesall2[key].size();
        for (int tnum:qorderp)
        {
            std::vector<std::vector<double>> ham(siz1, std::vector<double>(siz1, 0));
            std::vector<std::vector<double>> ham2(siz1, std::vector<double>(siz1, 0));
            std::vector<std::vector<double>> ham3(siz1, std::vector<double>(siz1, 0));
            for (int i=0;i< siz1;i++)
            {
                for (int j=0;j< siz1;j++)
                {
                    int n1=coupleBasesall2[key][i].ni;
                    int n2=coupleBasesall2[key][j].ni;
                    int p1=coupleBasesall2[key][i].pi;
                    int p2=coupleBasesall2[key][j].pi;
                    int jn=jvecn[n1];
                    int jnp=jvecn[n2];
                    int jp=jvecp[p1];
                    int jpp=jvecp[p2];
                    if (n1==n2)
                    {
                        ham[i][j]= ham[i][j]+matslipch[p1][p2];
                    }
                    if (p1==p2)
                    {
                        ham[i][j]=ham[i][j]+matslinch[n1][n2];
                    }
                    int deltann=deltatwo(jn,jnp);
                    int deltapp=deltatwo(jp,jpp);
                    ham[i][j]=ham[i][j]*deltann*deltapp;
                    std::vector<int> sigl={jpp,jn,key,tnum};
                    int sig=sign_func(sigl);
                    double c2=util::wigner_6j(jnp,jpp,key,jp,jn,tnum)*mysqrt(jnp+1)*mysqrt(jpp+1)*gnnn;
                    double c3=matsliqnch[0][n1][n2]*matsliqpch[0][p1][p2];
                    ham2[i][j]=c2*c3*sig;
                    ham3[i][j]=ham[i][j]+ham2[i][j];
                }
            }
            if (key<5)
            {
                printMatrix(ham3);
            }
            // printMatrix(ham3);
            int nmat=ham3.size();
            auto mv_mul1 = [&](const std::vector<double>& in, std::vector<double>& out) {
                for(int i = 0; i < nmat; ++i) {
                    for(int j = 0; j < nmat; ++j) {
                        out[i] +=ham3[i][j]*in[j];
                    }
                }
            };

            LambdaLanczos<double> engine1(mv_mul1, nmat, false, 2); // true means to calculate the largest eigenvalue.
            std::vector<double> eigenvalues1;
            std::vector<std::vector<double>> eigenvectors1;
            engine1.run(eigenvalues1, eigenvectors1);
            eigenre[key]=eigenvectors1[0];
            std::cout << "J= " <<key<< std::endl;
            std::cout << "Eigenvalues: " << std::endl;
            for (size_t i = 0; i < eigenvalues1.size(); ++i) {
                std::cout << "Eigenvalue " << i + 1 << ": " << eigenvalues1[i] << std::endl;
            }

            // 输出特征向量
            std::cout << "Eigenvectors: " << std::endl;
            for (size_t i = 0; i < eigenvectors1.size(); ++i) {
                std::cout << "Eigenvector " << i + 1 << ": [ ";
                for (size_t j = 0; j < eigenvectors1[i].size(); ++j) {
                    std::cout << eigenvectors1[i][j] << " ";
                }
                std::cout << "]" << std::endl;
            }
            //computeEigenvaluesAndEigenvectors(matchange1);

            double mmm=1;
        }
    //
    //
    }



    // double alpha1=4.3;
    // std::vector<std::vector<double> > qstrgetp1;
    // std::vector<int> qorderp1={4};
    // qstrgetp1=qstrall(qorderp1,alpha1);
    //
    // int sizeq1=qorderp1.size();
    // std::vector<std::vector<std::vector<double>>> qstrallp1(sizeq, std::vector<std::vector<double>>(size1,
    // std::vector<double>(size1, 0)));
    //
    // getystrbasis(qstrallp1, qstrgetp1, {0});
    //
    // std::vector<int> qordern1={4};
    // std::vector<std::vector<double> > qstrgetn1;
    // qstrgetn1=qstrall(qordern1,1);
    // int sizeqn1=qordern1.size();
    // std::vector<std::vector<std::vector<double>>> qstralln1(sizeqn1, std::vector<std::vector<double>>(size2,
    // std::vector<double>(size2, 0)));
    // getystrbasis(qstralln1, qstrgetn1, {0});
    //
    // std::vector<std::vector<std::vector<double>>> matsliqp1;
    // matsliqp1=matspl(allbasisp,ystrallp,{4},qstrallp1);
    // std::vector<std::vector<std::vector<double>>>matsliqpch1;
    // for(int i=0;i<matsliqp1.size();++i)
    // {
    //     std::vector<std::vector<double>> ch1=multiplyMatrices(schmitmatp,matsliqp1[i]);
    //     ch1=multiplyMatrices(ch1,schmitmatpinv);
    //     matsliqpch1.push_back(ch1);
    // }
    //
    // std::vector<std::vector<std::vector<double>>> matsliqn1;
    // matsliqn1=matspl(allbasisn,ystralln,qordern,qstralln1);
    // std::vector<std::vector<std::vector<double>>>matsliqnch1;
    // for(int i=0;i<matsliqn1.size();++i)
    // {
    //     std::vector<std::vector<double>> ch1=multiplyMatrices(schmitmatn,matsliqn1[i]);
    //     ch1=multiplyMatrices(ch1,schmitmatninv);
    //     matsliqnch1.push_back(ch1);
    // }
    //
    // int key1=7;
    // int key2=3;
    // double ep=1.5;
    // double ev=0.5;
    // double tfi=0;
    // int siz1=coupleBasesall2[key1].size();
    // int siz2=coupleBasesall2[key2].size();
    // int siz11=eigenre[key1].size();
    // int siz22=eigenre[key2].size();
    //
    // // 第一个循环，遍历 eigenre[key1]
    // for (int i = 0; i < siz11; ++i) {
    //     double value1 = eigenre[key1][i];
    //     int n1=coupleBasesall2[key1][i].ni;
    //     int p1=coupleBasesall2[key1][i].pi;
    //     // 第二个循环，遍历 eigenre[key2]
    //     for (int j = 0; j < siz22; ++j) {
    //         double value2 = eigenre[key2][j];
    //         int n2=coupleBasesall2[key2][j].ni;
    //         int p2=coupleBasesall2[key2][j].pi;
    //         double qp1=matsliqpch1[0][p1][p2]*value1*value2*ep;
    //         double qn1=matsliqnch1[0][n1][n2]*value1*value2*ev;
    //         tfi+=qp1+qn1;
    //
    //     }
    // }











    // for (const auto& [key, value] : coupleBasesall2)
    // {
    //
    //     int siz1=coupleBasesall2[key].size();
    //     for (int tnum:qorderp)
    //     {
    //         std::vector<std::vector<double>> ham(siz1, std::vector<double>(siz1, 0));
    //         std::vector<std::vector<double>> ham2(siz1, std::vector<double>(siz1, 0));
    //         std::vector<std::vector<double>> ham3(siz1, std::vector<double>(siz1, 0));
    //         for (int i=0;i< siz1;i++)
    //         {
    //             for (int j=0;j< siz1;j++)
    //             {
    //                 int n1=coupleBasesall2[key][i].ni;
    //                 int n2=coupleBasesall2[key][j].ni;
    //                 int p1=coupleBasesall2[key][i].pi;
    //                 int p2=coupleBasesall2[key][j].pi;
    //                 int jn=jvecn[n1];
    //                 int jnp=jvecn[n2];
    //                 int jp=jvecp[p1];
    //                 int jpp=jvecp[p2];
    //                 if (n1==n2)
    //                 {
    //                     ham[i][j]= ham[i][j]+matslipch[p1][p2];
    //                 }
    //                 if (p1==p2)
    //                 {
    //                     ham[i][j]=ham[i][j]+matslinch[n1][n2];
    //                 }
    //                 std::vector<int> sigl={jpp,jn,key,tnum};
    //                 int sig=sign_func(sigl);
    //                 double c2=util::wigner_6j(jnp,jpp,key,jp,jn,tnum)*mysqrt(jnp+1)*mysqrt(jpp+1)*gnnn;
    //                 double c3=matsliqnch[0][n1][n2]*matsliqpch[0][p1][p2];
    //                 ham2[i][j]=c2*c3;
    //                 ham3[i][j]=ham[i][j]+ham2[i][j];
    //             }
    //         }
    //         printMatrix(ham3);
    //         int nmat=ham3.size();
    //         auto mv_mul1 = [&](const vector<double>& in, vector<double>& out) {
    //             for(int i = 0; i < nmat; ++i) {
    //                 for(int j = 0; j < nmat; ++j) {
    //                     out[i] +=ham3[i][j]*in[j];
    //                 }
    //             }
    //         };
    //
    //         LambdaLanczos<double> engine1(mv_mul1, nmat, false, 2); // true means to calculate the largest eigenvalue.
    //         vector<double> eigenvalues1;
    //         vector<vector<double>> eigenvectors1;
    //         engine1.run(eigenvalues1, eigenvectors1);
    //         //computeEigenvaluesAndEigenvectors(matchange1);
    //
    //         double mmm=1;
    //     }
    //
    //
    // }

    // std::map<int, std::vector<CoupledBasis>> coupleBasesall;
    // coupleBasesall=coupleBases(allbasisn,allbasisp,overlapp,overlapn);
    // int Nall=coupleBasesall.size();
    // std::map<int,std::vector<std::vector<double>>> schmitall;
    // for (const auto& [key, value] : coupleBasesall) {
    //     std::vector<std::vector<double>>schmitmat;
    //     std::vector<std::vector<double>>copoverget;
    //     copoverget=cpomatcal(overlapn,overlapp,coupleBasesall[key]);
    //     schmitall[key]=copoverget;
    //     schmitall[key]=gramSchmidtInOverlap(copoverget);
    //     std::vector<std::vector<double>>l1=transposeMatrix(schmitall[key]);
    //     std::vector<std::vector<double>> l2=multiplyMatrices(schmitall[key],copoverget);
    //     std::vector<std::vector<double>> l3=multiplyMatrices(l2,l1);
    // }
    // std::map<int,std::vector<std::vector<double>>> cpmat;
    // std::vector<double> ktt={-0.06};
    // cpmat=cpmatcal(matslip,matslin,coupleBasesall,allbasisn,allbasisp,qorderp,qordern,matsliqn,matsliqp,ktt);
    // for (const auto& [key, value] : coupleBasesall)
    // {
    //     std::vector<std::vector<double>> hmat=cpmat[key];
    //     printMatrix( hmat);
    //     std::vector<std::vector<double>> schmitmat;
    //     std::vector<std::vector<double>> schmitmatinv;
    //     schmitmat=schmitall[key];
    //     printMatrix(schmitmat);
    //     schmitmatinv=transposeMatrix(schmitmat);
    //     std::vector<std::vector<double>> matchange;
    //     matchange=multiplyMatrices(schmitmat,hmat);
    //     std::vector<std::vector<double>> matchange1;
    //     matchange1=multiplyMatrices(matchange,schmitmatinv);
    //     printMatrix(matchange1);
    //
    //     int nmat=matchange1.size();
    //     auto mv_mul = [&](const vector<double>& in, vector<double>& out) {
    //         for(int i = 0; i < nmat; ++i) {
    //             for(int j = 0; j < nmat; ++j) {
    //                 out[i] +=matchange1[i][j]*in[j];
    //             }
    //         }
    //     };
    //
    //     LambdaLanczos<double> engine(mv_mul, nmat, false, 2); // true means to calculate the largest eigenvalue.
    //     vector<double> eigenvalues;
    //     vector<vector<double>> eigenvectors;
    //     engine.run(eigenvalues, eigenvectors);
    //
    //     auto mv_mul1 = [&](const vector<double>& in, vector<double>& out) {
    //         for(int i = 0; i < nmat; ++i) {
    //             for(int j = 0; j < nmat; ++j) {
    //                 out[i] +=hmat[i][j]*in[j];
    //             }
    //         }
    //     };
    //
    //     LambdaLanczos<double> engine1(mv_mul1, nmat, false, 2); // true means to calculate the largest eigenvalue.
    //     vector<double> eigenvalues1;
    //     vector<vector<double>> eigenvectors1;
    //     engine1.run(eigenvalues1, eigenvectors1);
    //     //computeEigenvaluesAndEigenvectors(matchange1);
    //
    //     double mmm=1;
    // }




    if (outfile.is_open()) {
        outfile.close();
    }

    return 0;
}


// TIP See CLion help at <a
// href="https://www.jetbrains.com/help/clion/">jetbrains.com/help/clion/</a>.
//  Also, you can try interactive lessons for CLion by selecting
//  'Help | Learn IDE Features' from the main menu.