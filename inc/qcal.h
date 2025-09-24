//
// Created by wang- on 2025/1/13.
//

#ifndef QCAL_H
#define QCAL_H
#include "../src/qcal.cpp"
#include "WignerSymbol.hpp"
double factorial(double x);
double double_factorial(double x);
double calculate_sum(int n, int n_prime, int l, int l_prime, int lambda);
double calculate_integral(int n, int n_prime, int l, int l_prime, int lambda, double alpha);
double computeQ(int j11, int j22, double tt, int n11, int l11, int n22, int l22,double alpha);
std::vector<std::vector<double>> m1qstrall(double gl,double gs);
#endif //QCAL_H
