//
// Created by wang- on 2025/3/25.
//
#include <iostream>
#include <fstream>
#include <array>
#include <cmath>
#include <iomanip>
#include "mtcal.h"
#include "qcal.h"

#include <lambda_lanczos/lambda_lanczos.hpp>
using std::cout;
using std::endl;
using std::setprecision;
using lambda_lanczos::LambdaLanczos;
std::vector<std::vector<double>> matrixcc(std::vector<std::vector<double>>mat1,std::vector<std::vector<double>>mat2,
std::vector<std::vector<double>>mat3)
{
    std::vector<std::vector<double>> result;
    result=multiplyMatrices(mat1,mat2);
    result=multiplyMatrices(result,mat3);
    return result;
};

std::vector<std::vector<double>> matrixAdd(const std::vector<std::vector<double>>& A, const
std::vector<std::vector<double>>&
B) {
    // 获取矩阵 A 和 B 的行数和列数
    int rows = A.size();
    int cols = A[0].size();

    // 创建结果矩阵 C，大小与 A 和 B 相同，初始值为 0
    std::vector<std::vector<double>> C(rows, std::vector<double>(cols, 0));

    // 矩阵加法操作
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            C[i][j] = A[i][j] + B[i][j];  // 对应元素相加
        }
    }

    return C;  // 返回结果矩阵 C
}

// 函数：矩阵数乘（浮点数版本）
std::vector<std::vector<double>> matrixScalarMultiply(const std::vector<std::vector<double>>& matrix, double scalar) {
    int rows = matrix.size();
    int cols = matrix[0].size();

    // 创建结果矩阵，大小与原矩阵相同
    std::vector<std::vector<double>> result(rows, std::vector<double>(cols, 0.0));

    // 数乘操作：将矩阵中的每个元素乘以标量
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[i][j] = matrix[i][j] * scalar;
        }
    }

    return result;
}



int main()

{

    std::vector<std::vector<double> > ystrget;
    nucleus={
        {0, 8, 7},
        {0, 4, 5},
    };
    allnucleus={
        {0, 8, 7},
        {0, 4, 5},
    };
    int size1=nucleus.size();
    ystrget = {
        {0, 0, 0,  0.6},
        {1, 1, 0,  0.8},
        {0, 0, 1, -0.7},
        {0, 1, 1,  0.2},
        {1, 0, 1, -0.2},
        {1, 1, 1, 0.66},
    };
    basis rjl;
    //rjl={0,0,{4,4,4},{1,1,1},{2,6,10}};
    rjl={7,0,{0,4},{0,1},{7,5},{1,1},1};
    // rjl={5,1,{4},{1},{3}};
    std::vector<int> ro1=rjl.rn;
    basis rjlinv;
    // rjlinv={7,0,{4},{1},{3}};

    //rjlinv={0,0,{4,4,4},{1,1,1},{2,6,10}};
    rjlinv={5,1,{0,0},{0,0},{5,5},{1,1},1};
    std::vector<int> ro2=rjlinv.rn;
    std::vector<int> rorder={0,4};
    int size2=2;

    std::vector<std::vector<std::vector<double>>> ystr1(size2, std::vector<std::vector<double>>(size1, std::vector<double>(size1, 0)));
    std::vector<std::vector<std::vector<double>>> ystr2(size2, std::vector<std::vector<double>>(size1, std::vector<double>(size1, 0)));

    getystrbasis(ystr1, ystrget, ro1);
    getystrbasis(ystr2, ystrget, ro2);
    double m111=overlap(rjl,rjlinv,ystr1,ystr2);
    std::vector<std::vector<double> > qstrgetp;
    std::vector<int> qorderp={4};
    qstrgetp=qstrall(qorderp);
    int sizeq=qorderp.size();
    std::vector<std::vector<std::vector<double>>> qstrallp(sizeq, std::vector<std::vector<double>>(size1,
    std::vector<double>(size1, 0)));
    getystrbasis(qstrallp, qstrgetp, {0});
    int tnum=0;
    std::vector<double> energy={1,1};
    // double mmm=sigqt(rjl,rjlinv,ystr1,ystr2,qstrallp,0,qorderp);
    // double mmm2=sigqt(rjlinv,rjl,ystr2,ystr1,qstrallp,0,qorderp);
    double mmm3=h0cal(rjl,rjlinv,ystr1,ystr2,energy);
    double n=qtqttest(rjl,rjlinv,ystr1,ystr2,qstrallp,tnum,qorderp);
    double n2=qtqttest(rjlinv,rjl,ystr2,ystr1,qstrallp,tnum,qorderp);
    /*
    std::vector<std::vector<double>>y1r2=
        {

            {-0.783923,  0.163586,  -0.132294,  0,  0},
            {-0.163586,  -0.49317,  -0.0687668,  -0.122677,  0},
            {-0.132294,  0.0687668,  -0.0596549,  0.0532782,  0},
            {0,  -0.122677,  -0.0532782,  0,  0},
            {0,  0,  0,  0,  -0.0658248}

    };
    std::vector<std::vector<double>>y1r1=
        {
        { {3.16119, 0, 0, 0, 0},
        {0, 2.28718, 0, 0, 0},
        {0, 0, 0.99166, 0, 0},
        {0, 0, 0, 0.859289295, 0},
        {0, 0, 0, 0, 0.267042}
        }};
    std::vector<std::vector<double>>y2r1={
        {0.53584976, 0, 0, 0, 0},
        {0, 0.837711, 0, 0, 0},
        {0, 0, 0.0683803, 0, 0},
        {0, 0, 0, 0.0541206, 0},
        {0, 0, 0, 0, -0.0591689}
    };
    std::vector<std::vector<double>>y2r2={
        {-0.783923,  0.163586,  -0.132294,  0,  0},
            {-0.163586,  -0.49317,  -0.0687668,  -0.122677,  0},
            {-0.132294,  0.0687668,  -0.0596549,  0.0532782,  0},
            {0,  -0.122677,  -0.0532782,  0,  0},
            {0,  0,  0,  0,  -0.0658248}
    };
    ystr1[0]=y1r1;
    ystr1[1]=y1r2;
    ystr2[0]=y2r1;
    ystr2[1]=y2r2;
    */
    /*double m12=pair2cal(ystr2,ystr1,rjlinv,rjl);

    std::vector<double> energyp={-0.09824028 ,0.84320220  ,2.25299294};

    double m1m=h0cal(rjl,rjlinv,ystr1,ystr2,energyp);
    double m2m=h0cal(rjlinv,rjl,ystr2,ystr1,energyp);
    double m3m=h0cal(rjlinv,rjlinv,ystr2,ystr2,energyp);

    std::vector<std::vector<double> > qstrgetp;
    std::vector<int> qorderp={4};
    qstrgetp=qstrall(qorderp);
    int sizeq=qorderp.size();
    std::vector<std::vector<std::vector<double>>> qstrallp(sizeq, std::vector<std::vector<double>>(size1,
    std::vector<double>(size1, 0)));
    getystrbasis(qstrallp, qstrgetp, {0});
    int tnum=0;
    double m1=sigqt(rjl,rjlinv,ystr1,ystr2,qstrallp,0,qorderp);
    double m2=sigqt(rjlinv,rjl,ystr2,ystr1,qstrallp,0,qorderp);
    */
    std::vector<std::vector<double> > astrgetp;
    std::vector<int> aorderp={0};
    astrgetp=qstrall(aorderp);
    int sizea=aorderp.size();
    std::vector<std::vector<std::vector<double>>> astrallp(sizea, std::vector<std::vector<double>>(size1,
    std::vector<double>(size1, 0)));
    getystrbasis(astrallp, astrgetp, {0});

    double na=atateven(rjl,rjlinv,ystr1,ystr2,astrallp,0,aorderp);
    /*
    double n=qtqt(rjl,rjlinv,ystr1,ystr2,qstrallp,tnum,qorderp);
    double n2=qtqt(rjlinv,rjl,ystr2,ystr1,qstrallp,tnum,qorderp);

    //std::cout<<"n"<<n<<std::endl;
    //std::cout<<"n2"<<n2<<std::endl;
    */
    return 0;


}

