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

using namespace std;

// Function to calculate the first equation (N)
std::tuple<double, double> calculate_bcs_parameters(
    const std::vector<std::vector<double>>& j_array,
    double G,
    int N,
    double lowlimit = 1e-10,
    int max_iterations = 1000)
{
    // 初始化参数
    double lambda1 = 1.0;
    double Delta1 = 1.0;
    double lambda2, Delta2;
    int iteration = 0;

    do {
        // 保存前次迭代结果
        lambda2 = lambda1;
        Delta2 = Delta1;

        // 计算新的lambda1
        double lambda_sum = 0.0;
        double summation1 = 0.0;

        for (size_t x1 = 0; x1 < j_array[0].size(); ++x1) {
            const double j3 = j_array[0][x1];  // 第三行对应索引2
            const double j5 = j_array[1][x1];  // 第五行对应索引4

            // 第一项求和
            const double term1 = (j5 / std::sqrt(std::pow(j5 - lambda2, 2) + std::pow(Delta2, 2)) - 1);
            lambda_sum += (j3 + 1.0) / 2.0 * term1;

            // 第二项求和
            summation1 += (j3 + 1.0) / 2.0 / std::sqrt(std::pow(j5 - lambda2, 2) + std::pow(Delta2, 2));
        }

        // 更新lambda
        lambda1 = (lambda_sum + N) / summation1;

        // 计算新的Delta1
        double delta_sum = 0.0;
        for (size_t x1 = 0; x1 < j_array[0].size(); ++x1) {
            const double j3 = j_array[0][x1];
            const double j5 = j_array[1][x1];

            delta_sum += (j3 + 1.0) / 2.0 / std::sqrt(1.0 + std::pow((j5 - lambda2)/Delta2, 2));
        }
        Delta1 = delta_sum * G / 2.0;

        // 防止无限循环
        if (++iteration > max_iterations) {
            throw std::runtime_error("BCS未收敛，达到最大迭代次数");
        }

    } while (std::max(std::abs(lambda1 - lambda2),
                     std::abs(Delta1 - Delta2)) > lowlimit);

    return {lambda1, Delta1};
}

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
calculateUVArrays(double lambda, double Delta, const std::vector<std::vector<double>>& j_array) {
    // 检查输入合法性



    // 分配 u_array 和 v_array
    std::vector<std::vector<double>> u_array(2, std::vector<double>(j_array[0].size()));
    std::vector<std::vector<double>> v_array(2, std::vector<double>(j_array[0].size()));

    // 计算 u_array 和 v_array
    for (size_t i = 0; i < j_array[0].size(); ++i) {
        v_array[0][i] = j_array[0][i];  // j_array(3,i)
        double term = (j_array[1][i] - lambda) / std::sqrt(std::pow(j_array[1][i] - lambda, 2) + std::pow(Delta, 2));
        v_array[1][i] = std::sqrt(0.5 * (1 - term));

        u_array[0][i] = j_array[0][i];  // j_array(3,i)
        u_array[1][i] = std::sqrt(1 - std::pow(v_array[1][i], 2));
    }

    // 返回结果
    return {u_array, v_array};
}
std::vector<std::vector<double>> calculateYStrAll(const std::vector<std::vector<double>>& j_array, double G, int N, const std::vector<double>& Omega) {
    // Step 1: 计算 BCS 参数
    int size1=nucleus.size();
    auto [lambda2, delta2] = calculate_bcs_parameters(j_array, G, N);

    // Step 2: 计算 u_array 和 v_array
    auto [u_array, v_array] = calculateUVArrays(lambda2, delta2, j_array);

    // Step 3: 初始化结果容器
    std::vector<std::vector<double>> yystrall;

    // Step 4: 计算 yystrgetp
    std::vector<int> yorderp = {0, 4};
    auto yystrgetp = qstrall(yorderp, 1);

    // Step 5: 初始化 yystrallp
    std::vector<std::vector<std::vector<double>>> yystrallp(yorderp.size(), std::vector<std::vector<double>>(size1, std::vector<double>(size1, 0)));

    // Step 6: 填充 yystrallp 和 yystrall
    for (size_t i = 0; i < u_array[1].size(); ++i) {
        yystrallp[0][i][i] = std::sqrt(Omega[i] + 1) * v_array[1][i] / u_array[1][i];
        yystrall.push_back({static_cast<double>(i), static_cast<double>(i), 0.0, yystrallp[0][i][i]});
    }

    // Step 7: 计算 qstrgetp
    std::vector<int> qorderp = {4};
    auto qstrgetp = qstrall(qorderp, 1);

    // Step 8: 初始化 qstrallp
    std::vector<std::vector<std::vector<double>>> qstrallp(qorderp.size(), std::vector<std::vector<double>>(size1, std::vector<double>(size1, 0)));

    // Step 9: 计算 qstrallp
    getystrbasis(qstrallp, qstrgetp, {0});

    // Step 10: 填充 yystrallp 和 yystrall
    for (size_t i = 0; i < u_array[1].size(); ++i) {
        for (size_t j = 0; j < u_array[1].size(); ++j) {
            yystrallp[1][i][j] = 0.5 * qstrallp[0][i][j] * ((yystrallp[0][i][i] / std::sqrt(Omega[i] + 1)) + (yystrallp[0][j][j] / std::sqrt(Omega[j] + 1)));
            yystrall.push_back({static_cast<double>(i), static_cast<double>(j), 1.0, yystrallp[1][i][j]});
        }
    }

    // Step 11: 返回 yystrall
    return yystrall;
}

int main() {
    nucleus={
        {1, 2, 3},
        {0, 6, 5},
        {1, 2, 1}
    };

    allnucleus={
        {1, 2, 3},
        {0, 6, 5},
        {1, 2, 1}
    };
    int size1=nucleus.size();

    // Example data (you should input your own data)
    vector<double> e = {0, 0.7, 1.08};    // Energy levels
    vector<double> Omega = {3, 5, 1}; // Omega_i values
    // 创建二维数组（两行）

    std::vector<std::vector<double>> j_array(2, std::vector<double>(e.size()));

    // 填充二维数组
    for (size_t i = 0; i < e.size(); ++i) {
        j_array[0][i] = Omega[i];           // 第一行为能量值
        j_array[1][i] =  e[i];    // 第二行为Omega值
    }
    double G = -0.33;                        // Constant G
    double N = 4.0;                        // Known value of N


    // Call the function to solve the equations using Newton's method
    std::tuple<double, double>m=calculate_bcs_parameters( j_array, G, N);
    double lambda2=std::get<0>(m);
    double delta2 = std::get<1>(m);


    auto [u_array, v_array]=calculateUVArrays(lambda2, delta2, j_array);
    std::vector<std::vector<double> > yystrgetp;
    std::vector<int> yorderp={0,4};
    yystrgetp=qstrall(yorderp,1);

    int sizeyy=yorderp.size();
    std::vector<std::vector<std::vector<double>>> yystrallp(sizeyy, std::vector<std::vector<double>>(size1,
    std::vector<double>(size1, 0)));
    std::vector<std::vector<double>> yystrall;

    for (size_t i = 0; i < u_array[1].size(); ++i) {
        yystrallp[0][i][i]=sqrt(Omega[i]+1)*v_array[1][i]/u_array[1][i];
        double ii=i;
        std::vector<double> x={ii,ii,0,yystrallp[0][i][i]};
        yystrall.push_back(x);

    }






    std::vector<std::vector<double> > qstrgetp;
    std::vector<int> qorderp={4};
    qstrgetp=qstrall(qorderp,1);

    int sizeq=qorderp.size();
    std::vector<std::vector<std::vector<double>>> qstrallp(sizeq, std::vector<std::vector<double>>(size1,
    std::vector<double>(size1, 0)));

    getystrbasis(qstrallp, qstrgetp, {0});
    for (size_t i = 0; i < u_array[1].size(); ++i) {
        for (size_t j = 0; j < u_array[1].size(); ++j)
        {
            yystrallp[1][i][j]=0.5*qstrallp[0][i][j]*((yystrallp[0][i][i]/sqrt(Omega[i]+1)) + (yystrallp[0][j][j]/sqrt(Omega[j]+1)));
            double ii=i;
            double jj=j;
            std::vector<double> x={ii,jj,1,yystrallp[1][i][j]};
            yystrall.push_back(x);
        }

    }
    auto yystrall1 = calculateYStrAll(j_array, G, N, Omega);
    // 输出结果
    for (const auto& row : yystrall) {
        for (double val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
    std::cout <<"----------"<< std::endl;
    // 输出结果
    for (const auto& row : yystrall1) {
        for (double val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }



    return 0;
}
