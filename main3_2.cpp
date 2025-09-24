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
#include"bcscal.h"


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
std::vector<std::vector<std::vector<double>>> splitBlockDiagonalMatrix(const std::vector<std::vector<double>>& A, const std::vector<int>& block_sizes) {
    int N = A.size();
    int start_row = 0;
    std::vector<std::vector<std::vector<double>>> blocks; // 存储所有对角块

    for (size_t i = 0; i < block_sizes.size(); ++i) {
        int block_size = block_sizes[i];
        std::vector<std::vector<double>> diagonalBlock(block_size, std::vector<double>(block_size));

        // 提取当前对角块
        for (int j = 0; j < block_size; ++j) {
            for (int k = 0; k < block_size; ++k) {
                diagonalBlock[j][k] = A[start_row + j][start_row + k];
            }
        }

        // 将当前对角块添加到 blocks 中
        blocks.push_back(diagonalBlock);
        start_row += block_size;
    }

    return blocks; // 返回所有对角块
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
        {1, 2, 3},
        {0, 6, 5},
        {1, 2, 1},
        {0, 8, 9}
    };

    allnucleus={
        {1, 2, 3},
        {0, 6, 5},
        {1, 2, 1},
        {0, 8, 9}
    };
    std::vector<double> energyp={0.0 ,0.78  ,1.56,4.52};
    int size1=nucleus.size();
    std::vector<std::vector<double>> j_array(2, std::vector<double>(energyp.size()));
    std::vector<double> Omega;

    // 填充二维数组
    for (size_t i = 0; i < energyp.size(); ++i) {
        Omega.push_back(nucleus[i].j);
        j_array[0][i] = Omega[i];           // 第一行为能量值
        j_array[1][i] =  energyp[i];    // 第二行为Omega值

    }
    double G=0.331;
    double N=4;
    auto yystrall1 = calculateYStrAll(j_array, G, N, Omega);

    std::vector<std::vector<double> > ystrgetp=yystrall1;

    normalizeBySecondColumn(ystrgetp);
    std::vector<int> rorderp={0,4};
    std::vector<std::vector<int>>rorderrp={rorderp,{10,10}};
    std::vector<std::pair<std::vector<int>, std::vector<int>>> catchpair1;
    catchpair1=generateValidPairs(4,rorderrp,51);
    std::vector<basis> allbasisp;
    allbasisp=calculateBasisWithRange(4, catchpair1);
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
    std::vector<std::vector<double>> matslip;
    std::vector<double> ktp={-0.226};
    std::vector<double> gsp={-0.3331};
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
    std::vector<int> jvecp;
    std::vector<int> blockform1;
    std::vector<int> jjv;
    for (int i=0;i<allbasisp.size();++i)
    {
        int ii;
        if (i==0)
        {
            ii=0;
        }
        jvecp.push_back(allbasisp[i].sj.back());
        if (i>0)
        {
            if (jvecp[i]!=jvecp[i-1])
            {
                blockform1.push_back(i-ii);
                ii=i;
                jjv.push_back(allbasisp[i-1].sj.back());
            }
        }
        if (i==allbasisp.size()-1)
        {
            blockform1.push_back(i+1-ii);
            jjv.push_back(allbasisp[i].sj.back());

        }

    }
    std::vector<std::vector<std::vector<double>> > blockh;
    blockh=splitBlockDiagonalMatrix(matslipch,blockform1);
    std::map<int, std::vector<double>> eigenre1;
    std::map<std::pair<int, int>, std::vector<double>> engre;
    for (size_t num=0;num<blockh.size();++num)
    {
        int key=jjv[num];
        std::vector<std::vector<double> > ham3=blockh[num];
        int nmat=ham3.size();
        auto mv_mul1 = [&](const vector<double>& in, vector<double>& out) {
            for(int i = 0; i < nmat; ++i) {
                for(int j = 0; j < nmat; ++j) {
                    out[i] +=ham3[i][j]*in[j];
                }
            }
        };

        LambdaLanczos<double> engine1(mv_mul1, nmat, false, 2); // true means to calculate the largest eigenvalue.
        vector<double> eigenvalues1;
        vector<vector<double>> eigenvectors1;
        engine1.run(eigenvalues1, eigenvectors1);
        eigenre1[key]=eigenvectors1[0];
        for (size_t lever=0;lever<eigenvectors1.size();++lever)
        {
            engre[{key,lever}]=eigenvectors1[lever];
        }
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
    }
    std::map<std::pair<int, int>, int> myMap1;
    for (size_t num=0;num<blockh.size();++num)
    {
        int mapn;
        for(int i=0;i<blockh[num].size();++i)
        {
            if (num==0 && i==0)
            {
                mapn=0;
            }
            myMap1[{jjv[num],i}]=mapn;
            mapn=mapn+1;
        }
    }



    // 填充二维数组
    for (size_t i = 0; i < energyp.size(); ++i) {
        Omega.push_back(nucleus[i].j);
        j_array[0][i] = Omega[i];           // 第一行为能量值
        j_array[1][i] =  energyp[i];    // 第二行为Omega值

    }
    double G1=0.331;
    double N1=6;
    auto yystrall11 = calculateYStrAll(j_array, G1, N1, Omega);

    std::vector<std::vector<double> > ystrgetp1=yystrall11;

    normalizeBySecondColumn(ystrgetp1);
    std::vector<int> rorderp1={0,4};
    std::vector<std::vector<int>>rorderrp1={rorderp1,{10,10}};
    std::vector<std::pair<std::vector<int>, std::vector<int>>> catchpair11;
    catchpair11=generateValidPairs(6,rorderrp1,51);
    std::vector<basis> allbasisp1;
    allbasisp1=calculateBasisWithRange(6, catchpair11);
    std::vector<std::vector<double> > qstrgetp1;
    std::vector<int> qorderp1={4};
    qstrgetp1=qstrall(qorderp1,1);

    int sizeq1=qorderp1.size();
    std::vector<std::vector<std::vector<double>>> qstrallp1(sizeq1, std::vector<std::vector<double>>(size1,
    std::vector<double>(size1, 0)));

    getystrbasis(qstrallp1, qstrgetp1, {0});
    std::vector<std::vector<double> > astrgetp1;
    std::vector<int> aorderp1={0};
    astrgetp=qstrall(aorderp1,1);

    int sizea1=aorderp1.size();
    std::vector<std::vector<std::vector<double>>> astrallp1(sizea1, std::vector<std::vector<double>>(size1,
    std::vector<double>(size1, 0)));
    getystrbasis(astrallp1, astrgetp1, {0});
    int sizep1=allbasisp1.size();
    int sizepr1=allbasisp1[0].r.size();
    std::vector<std::vector<std::vector<std::vector<double>>>> ystrallp1(
        sizep1,  // 第一维：大小为 sizep
        std::vector<std::vector<std::vector<double>>>(
            sizepr1,  // 第二维：大小为 sizepr
            std::vector<std::vector<double>>(
                size1,  // 第三维：大小为 size1
                std::vector<double>(size1, 0.0)  // 第四维：大小为 size1，元素初始化为 0.0
            )
        )
    );
    calculateystr(ystrallp1,ystrgetp1,allbasisp1);
    std::vector<std::vector<double>> matslip1;
    std::vector<double> ktp1={-0.226};
    std::vector<double> gsp1={-0.3331};
    printBasisVectorOneLine(allbasisp1);
    std::vector<std::vector<double> > overlapp1;
    overlapp1=overlapmatcal(allbasisp1,ystrallp1);
    updateMatrices(allbasisp1, overlapp1,ystrallp1);
    printBasisVectorOneLine(allbasisp1);
    printMatrix(overlapp1);
    // std::ofstream outfile("basis_output.txt");
    std::vector<std::vector<double>>schmitmatp1;
    schmitmatp1=gramSchmidtInOverlap(overlapp1,allbasisp1,ystrallp1);
    printMatrix(schmitmatp1);
    printMatrix(overlapp1);
    printBasisVectorOneLine(allbasisp1);
    std::cout << "fffff1 " << std::endl;
    // schmitmatp=Zmatrixp;
    std::vector<std::vector<double> > schmitmatpinv1;
    schmitmatpinv1=transposeMatrix(schmitmatp1);
    printMatrix(schmitmatpinv1);
    std::vector<std::vector<double> > schmitmatpinvre1=multiplyMatrices(schmitmatp1,schmitmatpinv1);
    printMatrix(schmitmatpinvre1);
    // printMatrix(overlapp);
    std::vector<std::vector<double> > overlappp1;
    overlappp1=multiplyMatrices(schmitmatp1,overlapp1);
    printMatrix(overlappp1);
    overlappp1=multiplyMatrices(overlappp1,schmitmatpinv1);
    printMatrix(overlappp1);
    // sqrtOverlap(overlappp);

    auto start11 = std::chrono::high_resolution_clock::now();
    matslip1=matcalsim(allbasisp1,ystrallp1,qorderp1,aorderp1,energyp,qstrallp1,ktp1,astrallp1,gsp1);
    std::vector<std::vector<double> > matslipch1=multiplyMatrices(schmitmatp1,matslip1);
    matslipch1=multiplyMatrices(matslipch1,schmitmatpinv1);
    auto end1 = std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::seconds>(end1 - start11);
    std::cout << "Time taken by function: " << duration1.count() << " seconds" << std::endl;
    printMatrix(matslip1);
    printMatrix(matslipch1);
    printBasisVectorOneLine(allbasisp1);
    std::vector<int> jvecp1;
    std::vector<int> blockform11;
    std::vector<int> jjv1;
    for (int i=0;i<allbasisp1.size();++i)
    {
        int ii;
        if (i==0)
        {
            ii=0;
        }
        jvecp1.push_back(allbasisp1[i].sj.back());
        if (i>0)
        {
            if (jvecp1[i]!=jvecp1[i-1])
            {
                blockform11.push_back(i-ii);
                ii=i;
                jjv1.push_back(allbasisp1[i-1].sj.back());
            }
        }
        if (i==allbasisp1.size()-1)
        {
            blockform11.push_back(i+1-ii);
            jjv1.push_back(allbasisp1[i].sj.back());

        }

    }
    std::vector<std::vector<std::vector<double>> > blockh1;
    blockh1=splitBlockDiagonalMatrix(matslipch1,blockform11);
    std::map<int, std::vector<double>> eigenre11;
    std::map<std::pair<int, int>, std::vector<double>> engre1;
    for (size_t num=0;num<blockh1.size();++num)
    {
        int key=jjv1[num];
        std::vector<std::vector<double> > ham3=blockh1[num];
        int nmat=ham3.size();
        auto mv_mul1 = [&](const vector<double>& in, vector<double>& out) {
            for(int i = 0; i < nmat; ++i) {
                for(int j = 0; j < nmat; ++j) {
                    out[i] +=ham3[i][j]*in[j];
                }
            }
        };

        LambdaLanczos<double> engine11(mv_mul1, nmat, false, 2); // true means to calculate the largest eigenvalue.
        vector<double> eigenvalues1;
        vector<vector<double>> eigenvectors1;
        engine11.run(eigenvalues1, eigenvectors1);
        eigenre11[key]=eigenvectors1[0];
        for (size_t lever=0;lever<eigenvectors1.size();++lever)
        {
            engre1[{key,lever}]=eigenvectors1[lever];
        }

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
    }

    std::map<std::pair<int, int>, int> myMap2;
    for (size_t num=0;num<blockh1.size();++num)
    {
        int mapn;
        for(int i=0;i<blockh1[num].size();++i)
        {
            if (num==0 && i==0)
            {
                mapn=0;
            }
            myMap2[{jjv1[num],i}]=mapn;
            mapn=mapn+1;
        }
    }

    // 打开输出文件
    std::ofstream outFile("output_table.csv");

    // 定义 n1 和 n2 的符号映射
    std::map<int, std::string> nSymbolMap = {
        {0, "2p3/2"},
        {1, "1f5/2"},
        {2, "2p1/2"},
        {3, "1g9/2"}
    };

    // 写入表头：所有 (n1符号, n2符号)r_i 组合
    outFile << "初态角动量,初态位置,末态角动量,末态位置";
    for (int n1 = 0; n1 <= 3; ++n1) {
        for (int n2 = 0; n2 <= 3; ++n2) {
            for (int r_i = 0; r_i <= 4; r_i += 4) {
                outFile << ",(" << nSymbolMap[n1] << "*" << nSymbolMap[n2] << ")" << r_i;
            }
        }
    }
    outFile << "\n";

    // 遍历所有 (n1, n2, r_i) 组合
    std::map<std::tuple<int, int, int, int>, std::map<std::tuple<int, int, int>, double>> resultMap;
    for (int n1 = 3; n1 <= 3; ++n1) {
        for (int n2 = 0; n2 <= 3; ++n2) {
            for (int r_i = 0; r_i <= 4; r_i += 4) {
                // 初始化矩阵
                int basis1size = allbasisp.size();
                int basis2size = allbasisp1.size();

                std::vector<std::vector<double>> ystadd1(size1, std::vector<double>(size1, 0));
                if (n1==n2)
                {
                    ystadd1[n1][n2] = 1;
                }else
                {
                    ystadd1[n1][n2] = 1/mysqrt(2);
                    int sigys=-std::pow(-1,(nucleus[n1].j+nucleus[n2].j+r_i)/2);
                    ystadd1[n2][n1]=sigys*1/mysqrt(2);

                }


                std::vector<std::vector<double>> matrixa2(basis1size, std::vector<double>(basis2size, 0));

                // 填充 matrixa2
                for (size_t i = 0; i < allbasisp.size(); ++i) {
                    std::vector<int> li1 = genmvec(jvecp[i], r_i);
                    int j1v = allbasisp[i].sj.back();
                    for (size_t j = 0; j < allbasisp1.size(); ++j) {
                        int j2v = allbasisp1[j].sj.back();
                        if (is_in_array(li1, j2v)) {
                            basis basisin = allbasisp[i];
                            basisin.sj.push_back(j2v);
                            basisin.r.push_back(r_i);
                            std::vector<std::vector<std::vector<double>>> ystrin(sizeq, std::vector<std::vector<double>>(size1, std::vector<double>(size1, 0)));
                            ystrin = ystrallp[i];
                            ystrin.push_back(ystadd1);
                            int termd = deltatwo(n1, n2);
                            double termd1 = 1 / mysqrt(termd + 1);
                            double cgterm = std::pow(-1, (j1v + r_i - j2v)/2);
                            double m_value = overlap(basisin, allbasisp1[j], ystrin, ystrallp1[j]);
                            matrixa2[i][j] = m_value * termd1;
                        }
                    }
                }

                // 计算 matrixa2ch
                std::vector<std::vector<double>> matrixa2ch(basis1size, std::vector<double>(basis2size, 0));
                matrixa2ch = multiplyMatrices(schmitmatp, matrixa2);
                matrixa2ch = multiplyMatrices(matrixa2ch, schmitmatpinv1);

                // 计算谱因子并存储到 resultMap
                for (auto& entry : engre) {
                    const std::pair<int, int>& key = entry.first;
                    std::vector<double>& values = entry.second;
                    for (auto& entry1 : engre1) {
                        const std::pair<int, int>& key1 = entry1.first;
                        std::vector<double>& values1 = entry1.second;
                        double result = 0;
                        for (int i = 0; i < values.size(); ++i) {
                            int mnum = myMap1[{key.first, i}];
                            for (int j = 0; j < values1.size(); ++j) {
                                int mnum1 = myMap2[{key1.first, j}];
                                result += values[i] * values1[j] * matrixa2ch[mnum][mnum1];
                            }
                        }
                        // 存储结果
                        resultMap[{key.first, key.second, key1.first, key1.second}][{n1, n2, r_i}] = result;
                    }
                }
            }
        }
    }

    // 将结果写入文件
    for (const auto& row : resultMap) {
        const auto& state = row.first;
        outFile << std::get<0>(state) << "," << std::get<1>(state) << "," << std::get<2>(state) << "," << std::get<3>(state);
        for (int n1 = 0; n1 <= 3; ++n1) {
            for (int n2 = 0; n2 <= 3; ++n2) {
                for (int r_i = 0; r_i <= 4; r_i += 4) {
                    outFile << "," << row.second.at({n1, n2, r_i});
                }
            }
        }
        outFile << "\n";
    }

    // 关闭文件
    outFile.close();
    std::cout << "结果已写入 output_table.csv 文件。" << std::endl;







    if (outfile.is_open()) {
        outfile.close();
    }

    return 0;
}


// TIP See CLion help at <a
// href="https://www.jetbrains.com/help/clion/">jetbrains.com/help/clion/</a>.
//  Also, you can try interactive lessons for CLion by selecting
//  'Help | Learn IDE Features' from the main menu.