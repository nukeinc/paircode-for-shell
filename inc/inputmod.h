//
// Created by wang- on 25-7-31.
//

#ifndef INPUTMOD_H
#define INPUTMOD_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cctype>
#include <algorithm>
#include <map>
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
#include"writeout.h"
#include"bcscal.h"
#include <sstream>
#include <utility>
#include "maincal.h"

struct Data {
    std::vector<Nucleus> nucleus1;
    std::vector<std::vector<int> >allnucleus1;
    std::vector<Nucleus> nucleus2;
    std::vector<std::vector<int> >allnucleus2;
    std::vector<double> energyp;
    std::vector<double> energyn;
    std::vector<std::vector<int>> rorder1;
    std::vector<std::vector<int>> rorder2;
    std::vector<std::vector<double>> ystrget1; // 新增yst1数据
    std::vector<std::vector<double>> ystrget2; // 新增yst2数据
    std::vector<double> efcstrength1;
    std::vector<double> efcstrength2;
    std::vector<double> efcstrength3;
    std::vector<std::map<int, Matrix4D>>  efc1; // 存储efc1数据
    std::vector<std::map<int, Matrix4D>> efc2; // 存储efc2数据
    std::vector<std::map<int, Matrix4D>> pnData; // 存储V_value数据
    std::vector<int> bej;
    std::vector<double> eg1;
    std::vector<double>eg2;
    int num1;
    int num2;
    int mcal;
    double alpha1;
    double alpha2;

};
#include "inputmod.cpp"
void processystrgetRow(const std::vector<double>& row, Data& data,int f);
std::string toLower(const std::string& str);









#endif //INPUTMOD_H
