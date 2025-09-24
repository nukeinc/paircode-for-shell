//
// Created by wang- on 2025/1/13.
//
// 阶乘函数
#include "qcal.h"
#include <cmath>
#include <algorithm> // std::min

#include <cmath>
#include <algorithm> // std::min
#include "WignerSymbol.hpp"
#include "moe.h"
#include "moe.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
// 阶乘函数 (支持分数)
double factorial(double x) {
    if (x < 0) return 0; // 阶乘参数必须非负
    return std::tgamma(x + 1); // 使用伽玛函数计算阶乘
}

// 双阶乘函数 (支持分数)
double double_factorial(double x) {
    if (x <= 0) return 1; // 双阶乘停止条件
    double result = 1.0;
    for (double i = x; i > 0; i -= 2) {
        result *= i;
    }
    return result;
}

// 求和部分的计算
double calculate_sum(int n, int n_prime, int l, int l_prime, int lambda) {
    double sum = 0.0;

    // 计算常量部分
    int l_l_prime_lambda = l + l_prime + lambda;

    // 遍历 q 的所有可能值
    for (int q = 0; q <= std::min(n, n_prime); ++q) {
        // 计算每一项中的分母因子
        int term1 = n - q;
        int term2 = n_prime - q;
        double term3 = q + (l_prime - l + lambda) / 2.0 - n;
        double term4 = q + (l - l_prime + lambda) / 2.0 - n_prime;

        // 确保所有参数非负，否则跳过该项
        if (term1 < 0 || term2 < 0 || term3 < 0 || term4 < 0) {
            continue;
        }

        // 计算分子部分 (l + l' + λ + 2q + 1)!!
        double numerator = double_factorial(l_l_prime_lambda + 2 * q + 1);

        // 计算分母部分
        double denominator = std::pow(2, q) * factorial(q) *
                             factorial(term1) * factorial(term2) *
                             factorial(term3) * factorial(term4);

        // 累加到总和
        sum += numerator / denominator;
    }

    return sum;
}

// 主公式的计算
double calculate_integral(int n, int n_prime, int l, int l_prime, int lambda, double alpha) {
    // (-1)^(n - n')
    int sign_factor = (n - n_prime) % 2 == 0 ? 1 : -1;

    // 开平方项
    double sqrt_factor = std::sqrt(
        (factorial(n) * factorial(n_prime) * std::pow(2, n + n_prime - lambda)) /
        (double_factorial(2 * n + 2 * l + 1) * double_factorial(2 * n_prime + 2 * l_prime + 1))
    );

    // 对称性项 (1 + (-1)^(l + l' + λ)) / 2
    int symmetry_factor = (l + l_prime + lambda) % 2 == 0 ? 1 : 0;

    // 阶乘部分
    double factorial_part = factorial((l_prime - l + lambda) / 2.0) *
                            factorial((l - l_prime + lambda) / 2.0);

    // 求和部分
    double summation = calculate_sum(n, n_prime, l, l_prime, lambda);

    // 最终结果
    return sign_factor * sqrt_factor * symmetry_factor * factorial_part * summation / std::pow(alpha, lambda);
}

double computeQ(int j11, int j22, double tt, int n11, int l11, int n22, int l22,double alpha) {
    double j1=j11/2.0;
    double j2=j22/2.0;
    double t=tt/2.0;
    int n1=n11/2;
    int n2=n22/2;
    int l1=l11/2;
    int l2=l22/2;
    // 计算前面的符号部分
    double prefactor = std::pow(-1, j1 + t - 0.5);

    // 计算根号部分
    double sqrtPart = std::sqrt((2 * j1 + 1) * (2 * j2 + 1) / (4 * M_PI * (2 * t + 1)));

    // 计算 Clebsch-Gordan 系数
    double clebsch = util::CG(j11, j22, tt, 1, -1,0);
    // double alphaa=4.33158;

    // 计算径向积分
    double radialIntegral = calculate_integral(n1, n2, l1, l2, t,alpha);

    // 计算最终结果
    double result = prefactor * sqrtPart * clebsch * (1 + std::pow(-1, l1 + l2 + t)) / 2.0 * radialIntegral;
    if (t==0)
    {
        if (j11==j22)
        {
            result=sqrt(2.0*j1+1)/2;
        }
        else
        {
            return 0;
        }
    }
    return result;
}

std::vector<std::vector<double>> qstrall(std::vector<int> qorder,double alpha)
{
    std::vector<std::vector<double>> result;
    int size1 = nucleus.size();
    for (int q = 0; q < qorder.size(); ++q)
    {

        for (int i = 0; i < size1; ++i)
        {

            for (int j = 0; j < size1; ++j)
            {
                std::vector<double> q1={};
                q1.push_back(i);
                q1.push_back(j);
                q1.push_back(q);
                double q1num=computeQ(nucleus[i].j,nucleus[j].j,qorder[q],nucleus[i].n,nucleus[i].l,nucleus[j].n,nucleus[j].l, alpha);
                q1.push_back(q1num);
                result.push_back(q1);
            }
        }

    }
    return result;
}

double computem1q(int l,int cj,int dj)
{
    double lhalf = l / 2.0;
    double cjhalf = cj / 2.0;
    double djhalf = dj / 2.0;
    double sigfac = std::pow(-1, lhalf  + 0.5+djhalf);
    double qrtterm=std::sqrt(lhalf*(lhalf+1)/3.0);
    double jiaoterm=std::sqrt((cj+1)*(dj+1)*(l+1));
    double jterm=util::wigner_6j(cj,dj,2,l,l,1);
    double result= sigfac * qrtterm * jiaoterm * jterm;
    return result;
}

double computem1qp(int l,int cj,int dj)
{
    double lhalf = l / 2.0;
    double cjhalf = cj / 2.0;
    double djhalf = dj / 2.0;
    double sigfac = std::pow(-1, lhalf  + 0.5+cjhalf);
    double qrtterm=1/std::sqrt(2);
    double jiaoterm=std::sqrt((cj+1)*(dj+1));
    double jterm=util::wigner_6j(cj,dj,2,1,1,l);
    double result= sigfac * qrtterm * jiaoterm * jterm;
    return result;
}

std::vector<std::vector<double>> m1qstrall(double gl,double gs)
{
    std::vector<std::vector<double>> result;
    int size1 = nucleus.size();
    int q=0;


    for (int i = 0; i < size1; ++i)
    {

        for (int j = 0; j < size1; ++j)
        {
            if (nucleus[i].l != nucleus[j].l)
            {
                continue; // 只考虑 l 相同的情况
            }
            std::vector<double> q1={};
            q1.push_back(i);
            q1.push_back(j);
            q1.push_back(q);
            int cj=nucleus[i].j;
            int dj=nucleus[j].j;
            int l=nucleus[i].l;
            double re1=computem1q(l,cj,dj);
            double re2=computem1qp(l,cj,dj);
            double q1num=re1*gl+re2*gs;
            q1.push_back(q1num);
            result.push_back(q1);
        }
    }


    return result;
}
