//
// Created by wang- on 25-7-31.
//

#include "inputmod.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cctype>
#include <algorithm>
// 辅助函数：转换为小写
#include "maincal.h"
std::string toLower(const std::string& str) {
    std::string lowerStr = str;
    std::transform(lowerStr.begin(), lowerStr.end(), lowerStr.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return lowerStr;
}

std::vector<std::map<int, Matrix4D>> buildVValue1(const std::vector<DataRow>& data, const int num=1) {
    // 1. 找出最大的t值
    int max_t = 0;
    for (const auto& row : data) {
        if (row.t > max_t) {
            max_t = row.t;
        }
    }

    // 2. 初始化V_value结构
    std::vector<std::map<int, Matrix4D>> V_value(max_t + 1);

    // 3. 遍历数据填充矩阵
    for (const auto& row : data) {
        int sizeab=nucleus.size();
        if (num==2)
        {
            sizeab=nucleus2.size();
        }

        // 检查是否需要初始化该J对应的Matrix4D
        if (V_value[row.t][row.J].empty()) {
            // 获取所有数据的最大c和d值来确定Matrix4D的维度
            int max_a = 0;
            int max_b = 0;
            for (const auto& r : data) {
                if (r.t == row.t && r.J == row.J) {
                    max_a = std::max(max_b, r.a);
                    max_b = std::max(max_b, r.b);
                }
            }

            // 初始化Matrix4D：c×d的矩阵，每个矩阵初始为0×0
            V_value[row.t][row.J] = Matrix4D(
                sizeab,
                std::vector<Eigen::MatrixXd>(
                    sizeab,
                    Eigen::MatrixXd::Zero(sizeab, sizeab)
                )
            );
        }
        int j1=0;
        int j2=0;
        int j3=0;
        int j4=0;
        if (num==1)
        {
            j1=nucleus[row.a].j;
            j2=nucleus[row.b].j;
            j3=nucleus[row.c].j;
            j4=nucleus[row.d].j;
        }else if (num==2)
        {
            j1=nucleus2[row.a].j;
            j2=nucleus2[row.b].j;
            j3=nucleus2[row.c].j;
            j4=nucleus2[row.d].j;

        }

        double f2134=-1*std::pow(-1,(j1+j2-row.J)/2);
        double f1243=-1*std::pow(-1,(j3+j4-row.J)/2);
        double f2143=std::pow(-1,(j1+j2+j3+j4)/2);

        // 存入value
        V_value[row.t][row.J][row.a][row.b](row.c,row.d) = row.value;
        V_value[row.t][row.J][row.c][row.d](row.a,row.b) = row.value;
        V_value[row.t][row.J][row.b][row.a](row.c,row.d) = f2134 * row.value;
        V_value[row.t][row.J][row.c][row.d](row.b,row.a) = f2134 * row.value;
        V_value[row.t][row.J][row.a][row.b](row.d,row.c) = f1243 * row.value;
        V_value[row.t][row.J][row.d][row.c](row.a,row.b) = f1243 * row.value;
        V_value[row.t][row.J][row.b][row.a](row.d,row.c) = f2143 * row.value;
        V_value[row.t][row.J][row.d][row.c](row.b,row.a) = f2143 * row.value;
    }

    return V_value;
}


std::vector<std::map<int, Matrix4D>> buildVValuepn1(const std::vector<DataRow>& data) {
    // 1. 找出最大的t值
    int max_t = 0;
    for (const auto& row : data) {
        if (row.t > max_t) {
            max_t = row.t;
        }
    }

    // 2. 初始化V_value结构
    std::vector<std::map<int, Matrix4D>> V_value(max_t + 1);

    // 3. 遍历数据填充矩阵
    for (const auto& row : data) {
        int sizeab=nucleus.size();
        int sizecd=nucleus2.size();
        // 检查是否需要初始化该J对应的Matrix4D
        if (V_value[row.t][row.J].empty()) {
            // 获取所有数据的最大c和d值来确定Matrix4D的维度
            int max_a = 0;
            int max_b = 0;
            int max_c = 0;
            int max_d = 0;
            for (const auto& r : data) {
                if (r.t == row.t && r.J == row.J) {
                    max_a = std::max(max_a, r.a);
                    max_b = std::max(max_b, r.b);
                    max_c = std::max(max_c, r.c);
                    max_d = std::max(max_d, r.d);
                }
            }
            // sizeab=std::max(max_a, max_b);
            // sizecd=std::max(max_c, max_d);

            // 初始化Matrix4D：c×d的矩阵，每个矩阵初始为0×0
            V_value[row.t][row.J] = Matrix4D(
                sizeab,
                std::vector<Eigen::MatrixXd>(
                    sizeab,
                    Eigen::MatrixXd::Zero(sizecd, sizecd)
                )
            );
        }
        int j1=nucleus[row.a].j;
        int j2=nucleus[row.b].j;
        int j3=nucleus2[row.c].j;
        int j4=nucleus2[row.d].j;


        // 存入value
        V_value[row.t][row.J][row.a][row.b](row.c,row.d) = row.value;
        V_value[row.t][row.J][row.b][row.a](row.d,row.c) = row.value;
    }

    return V_value;
}

// 处理efc1数据行
void processystrgetRow(const std::vector<double>& row, Data& data,int f) {
    // 创建行的副本（避免修改原始数据）
    std::vector<double> processedRow = row;

    // 将前三个值转换为0-based索引
    for (int i = 0; i < 3; i++) {
        processedRow[i] -= 1;
    }

    // 添加到efc1数据
    if (f==1)
    {
        data.ystrget1.push_back(processedRow);
    }
    if (f==2)
    {
        data.ystrget2.push_back(processedRow);
    }
}



// 主文件读取函数
Data readMultipleArraysFromFile(const std::string& filename) {
    std::ifstream file(filename);  // 打开文件
    Data result;                   // 用于存储不同类型数据的结构体
    std::string currentSection;    // 当前正在处理的数据部分
    std::string line;              // 用来存储每一行
    // 临时存储 pn 数据行
    std::vector<DataRow> pnRows;
    std::vector<DataRow> efc1;
    std::vector<DataRow> efc2;

    if (!file.is_open()) {
        // 检查文件是否成功打开
        throw std::runtime_error("无法打开文件: " + filename);
    }

    // 逐行读取文件
    while (std::getline(file, line))
    {
        // 移除行尾的回车符（Windows兼容）
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        // 跳过空行
        if (line.empty()) {
            currentSection.clear();
            continue;
        }

        // 跳过注释行（以#开头）
        if (line[0] == '#') {
            continue;
        }

        // 检查是否是关键字行
        std::istringstream iss(line);
        std::string firstWord;
        if (iss >> firstWord) {
            // 转换为小写进行关键字匹配
            std::string lowerWord = toLower(firstWord);
            if (lowerWord =="nucleus1")
            {
                currentSection = "nucleus1";
                int n_val, l_val, j_val;
                if (iss >> n_val >> l_val >> j_val) {
                    result.nucleus1.push_back(Nucleus(n_val, l_val, j_val));
                    result.allnucleus1.push_back({n_val, l_val, j_val});
                }
                continue;
            }
            else if (lowerWord =="nucleus2")
            {
                currentSection = "nucleus2";
                int n_val, l_val, j_val;
                if (iss >> n_val >> l_val >> j_val) {
                    result.nucleus2.push_back(Nucleus(n_val, l_val, j_val));
                    result.allnucleus2.push_back({n_val, l_val, j_val});
                }
                continue;
            }
            else if (lowerWord == "energy1") {
                currentSection = "energyp";
                // 检查同一行是否有数据
                double value;
                while (iss >> value) {
                    result.energyp.push_back(value);
                }
                continue;
            }
            else if(lowerWord == "energy2") {
                currentSection = "energyn";
                // 检查同一行是否有数据
                double value;
                while (iss >> value) {
                    result.energyn.push_back(value);
                }
                continue;
            }

            else if(lowerWord == "efcstrength1") {
                currentSection = "efcstrength1";
                // 检查同一行是否有数据
                double value;
                while (iss >> value) {
                    result.efcstrength1.push_back(value);
                }
                continue;
            }
            else if(lowerWord == "efcstrength2") {
                currentSection = "efcstrength2";
                // 检查同一行是否有数据
                double value;
                while (iss >> value) {
                    result.efcstrength2.push_back(value);
                }
                continue;
            }
            else if(lowerWord == "efcstrength3") {
                currentSection = "efcstrength3";
                // 检查同一行是否有数据
                double value;
                while (iss >> value) {
                    result.efcstrength3.push_back(value);
                }
                continue;
            }
            else if(lowerWord == "eg1") {
                currentSection = "eg1";
                // 检查同一行是否有数据
                double value;
                while (iss >> value) {
                    result.eg1.push_back(value);
                }
                continue;
            }
            else if(lowerWord == "eg2") {
                currentSection = "eg2";
                // 检查同一行是否有数据
                double value;
                while (iss >> value) {
                    result.eg2.push_back(value);
                }
                continue;
            }
            else if(lowerWord == "bej") {
                currentSection = "bej";
                // 检查同一行是否有数据
                int value;
                while (iss >> value) {
                    result.bej.push_back(value);
                }
                continue;
            }
            else if(lowerWord == "mcal") {
                currentSection = "mcal";
                // 检查同一行是否有数据
                int value;
                while (iss >> value) {
                    result.mcal=value;
                }
                continue;
            }
            else if(lowerWord == "num1") {
                currentSection = "num1";
                // 检查同一行是否有数据
                int value;
                while (iss >> value) {
                    result.num1=value;
                }
                continue;
            }
            else if(lowerWord == "num2") {
                currentSection = "num2";
                // 检查同一行是否有数据
                int value;
                while (iss >> value) {
                    result.num2=value;
                }
                continue;
            }
            else if(lowerWord == "alpha1") {
                currentSection = "alpha1";
                // 检查同一行是否有数据
                double value;
                while (iss >> value) {
                    result.alpha1=value;
                }
                continue;
            }
            else if(lowerWord == "alpha2") {
                currentSection = "alpha2";
                // 检查同一行是否有数据
                double value;
                while (iss >> value) {
                    result.alpha2=value;
                }
                continue;
            }
            else if (lowerWord =="rorder1")
            {
                currentSection = "rorder1";
                int value;
                std::vector<int> row;
                while (iss >> value)
                {
                    row.push_back(value);
                }
                if (!row.empty()) {
                    result.rorder1.push_back(row);
                }
                continue;
            }
            else if (lowerWord =="rorder2")
            {
                currentSection = "rorder2";
                int value;
                std::vector<int> row;
                while (iss >> value)
                {
                    row.push_back(value);
                }
                if (!row.empty()) {
                    result.rorder2.push_back(row);
                }
                continue;
            }
            else if (lowerWord == "ystrget1") {
                currentSection = "ystrget1";
                // 尝试读取同一行的数据
                std::vector<double> row;
                double value;
                while (iss >> value) {
                    row.push_back(value);
                }
                // 如果一行有4个值，处理并添加
                if (row.size() == 4) {
                    processystrgetRow(row, result,1);
                }
                continue;
            }
            else if (lowerWord == "ystrget2") {
                currentSection = "ystrget2";
                // 尝试读取同一行的数据
                std::vector<double> row;
                double value;
                while (iss >> value) {
                    row.push_back(value);
                }
                // 如果一行有4个值，处理并添加
                if (row.size() == 4) {
                    processystrgetRow(row, result,2);
                }
                continue;
            }
            else if (lowerWord == "efc1")
            {
                currentSection = "efc1";
                continue;
            }
            else if (lowerWord == "efc2")
            {
                currentSection = "efc2";
                continue;
            }
            else if (lowerWord == "pnefc") {
                currentSection = "pnefc";
                continue;
            }
        }

        // 处理数据行
        if (!currentSection.empty()) {
            std::istringstream dataIss(line);
            if (currentSection == "nucleus1")
            {
                int n_val, l_val, j_val;
                if (dataIss >> n_val >> l_val >> j_val) {
                    result.nucleus1.push_back(Nucleus(n_val, l_val, j_val));
                    result.allnucleus1.push_back({n_val, l_val, j_val});
                } else {
                    std::cerr << "警告: nucleus数据行需要3个整数 - " << line << std::endl;
                }
            }
            else if (currentSection == "nucleus2")
            {
                int n_val, l_val, j_val;
                if (dataIss >> n_val >> l_val >> j_val) {
                    result.nucleus2.push_back(Nucleus(n_val, l_val, j_val));
                    result.allnucleus2.push_back({n_val, l_val, j_val});
                } else {
                    std::cerr << "警告: nucleus数据行需要3个整数 - " << line << std::endl;
                }
            }
            else if (currentSection == "energyp") {
                double value;
                while (dataIss >> value) {
                    result.energyp.push_back(value);
                }

                // 检查是否有解析错误
                if (dataIss.fail() && !dataIss.eof()) {
                    std::cerr << "警告: 粒子能量数据解析错误 - " << line << std::endl;
                    dataIss.clear();
                }
            }
            else if (currentSection == "energyn") {
                double value;
                while (dataIss >> value) {
                    result.energyn.push_back(value);
                }

                // 检查是否有解析错误
                if (dataIss.fail() && !dataIss.eof()) {
                    std::cerr << "警告: 结构系数数据解析错误 - " << line << std::endl;
                    dataIss.clear();
                }
            }
            else if (currentSection == "bej") {
                int value;
                while (dataIss >> value) {
                    result.bej.push_back(value);
                }

                // 检查是否有解析错误
                if (dataIss.fail() && !dataIss.eof()) {
                    std::cerr << "警告: 结构系数数据解析错误 - " << line << std::endl;
                    dataIss.clear();
                }
            }
            else if (currentSection == "eg1") {
                double value;
                while (dataIss >> value) {
                    result.eg1.push_back(value);
                }

                // 检查是否有解析错误
                if (dataIss.fail() && !dataIss.eof()) {
                    std::cerr << "警告: 粒子能量数据解析错误 - " << line << std::endl;
                    dataIss.clear();
                }
            }
            else if (currentSection == "eg2") {
                double value;
                while (dataIss >> value) {
                    result.eg2.push_back(value);
                }

                // 检查是否有解析错误
                if (dataIss.fail() && !dataIss.eof()) {
                    std::cerr << "警告: 粒子能量数据解析错误 - " << line << std::endl;
                    dataIss.clear();
                }
            }
            else if (currentSection == "efcstrength1") {
                double value;
                while (dataIss >> value) {
                    result.efcstrength1.push_back(value);
                }

                // 检查是否有解析错误
                if (dataIss.fail() && !dataIss.eof()) {
                    std::cerr << "警告: 结构系数数据解析错误 - " << line << std::endl;
                    dataIss.clear();
                }
            }
            else if (currentSection == "efcstrength2") {
                double value;
                while (dataIss >> value) {
                    result.efcstrength2.push_back(value);
                }

                // 检查是否有解析错误
                if (dataIss.fail() && !dataIss.eof()) {
                    std::cerr << "警告: 结构系数数据解析错误 - " << line << std::endl;
                    dataIss.clear();
                }
            }
            else if (currentSection == "efcstrength3") {
                double value;
                while (dataIss >> value) {
                    result.efcstrength3.push_back(value);
                }

                // 检查是否有解析错误
                if (dataIss.fail() && !dataIss.eof()) {
                    std::cerr << "警告: 结构系数数据解析错误 - " << line << std::endl;
                    dataIss.clear();
                }
            }
            else if (currentSection == "num1") {
                double value;
                while (dataIss >> value) {
                    result.num1= value;
                }

                // 检查是否有解析错误
                if (dataIss.fail() && !dataIss.eof()) {
                    std::cerr << "警告: 结构系数数据解析错误 - " << line << std::endl;
                    dataIss.clear();
                }
            }
            else if (currentSection == "num2") {
                double value;
                while (dataIss >> value) {
                    result.num2= value;
                }

                // 检查是否有解析错误
                if (dataIss.fail() && !dataIss.eof()) {
                    std::cerr << "警告: 结构系数数据解析错误 - " << line << std::endl;
                    dataIss.clear();
                }
            }
            else if (currentSection == "alpha1") {
                double value;
                while (dataIss >> value) {
                    result.alpha1= value;
                }

                // 检查是否有解析错误
                if (dataIss.fail() && !dataIss.eof()) {
                    std::cerr << "警告: 结构系数数据解析错误 - " << line << std::endl;
                    dataIss.clear();
                }
            }
            else if (currentSection == "alpha2") {
                double value;
                while (dataIss >> value) {
                    result.alpha2= value;
                }

                // 检查是否有解析错误
                if (dataIss.fail() && !dataIss.eof()) {
                    std::cerr << "警告: 结构系数数据解析错误 - " << line << std::endl;
                    dataIss.clear();
                }
            }
            else if (currentSection == "mcal") {
                int value;
                while (dataIss >> value) {
                    result.mcal= value;
                }

                // 检查是否有解析错误
                if (dataIss.fail() && !dataIss.eof()) {
                    std::cerr << "警告: 结构系数数据解析错误 - " << line << std::endl;
                    dataIss.clear();
                }
            }
            else if (currentSection == "rorder1") {
                int value;
                std::vector<int> row;
                while (dataIss >> value) {
                    row.push_back(value);
                }
                result.rorder1.push_back(row);
                // 检查是否有解析错误
                if (dataIss.fail() && !dataIss.eof()) {
                    std::cerr << "警告: qorder1解析错误 - " << line << std::endl;
                    dataIss.clear();
                }
            }
            else if (currentSection == "rorder2") {
                int value;
                std::vector<int> row;
                while (dataIss >> value) {
                    row.push_back(value);
                }
                result.rorder2.push_back(row);
                // 检查是否有解析错误
                if (dataIss.fail() && !dataIss.eof()) {
                    std::cerr << "警告: qorder2解析错误 - " << line << std::endl;
                    dataIss.clear();
                }
            }
            else if (currentSection == "ystrget1") {
                std::vector<double> row;
                double value;
                while (dataIss >> value) {
                    row.push_back(value);
                }

                if (row.size() == 4) {
                    processystrgetRow(row, result,1);
                } else if (!row.empty()) {
                    std::cerr << "警告: ystrget数据行应有4个数值，实际找到 "
                              << row.size() << " 个 - " << line << std::endl;
                }
            }
            else if (currentSection == "ystrget2") {
                std::vector<double> row;
                double value;
                while (dataIss >> value) {
                    row.push_back(value);
                }

                if (row.size() == 4) {
                    processystrgetRow(row, result,2);
                } else if (!row.empty()) {
                    std::cerr << "警告: ystrget数据行应有4个数值，实际找到 "
                              << row.size() << " 个 - " << line << std::endl;
                }
            }
            else if (currentSection == "efc1") {
                DataRow row;
                if (dataIss >> row.t >> row.a >> row.b >> row.c >> row.d >> row.J >> row.value) {
                    // 转换为0-based索引
                    row.t -= 1;
                    row.a -= 1;
                    row.c -= 1;
                    row.b -= 1;
                    row.d -= 1;
                    row.value=row.value*4; // 乘以4
                    efc1.push_back(row);
                } else {
                    std::cerr << "警告: efc1数据行解析失败 - " << line << std::endl;
                }

            }
            else if (currentSection == "efc2") {
                DataRow row;
                if (dataIss >> row.t >> row.a >> row.b >> row.c >> row.d >> row.J >> row.value) {
                    // 转换为0-based索引
                    row.t -= 1;
                    row.a -= 1;
                    row.c -= 1;
                    row.b -= 1;
                    row.d -= 1;
                    row.value=row.value*4; // 乘以4
                    efc2.push_back(row);
                } else {
                    std::cerr << "警告: efc2数据行解析失败 - " << line << std::endl;
                }

            }
            else if (currentSection == "pnefc") {
                DataRow row;
                if (dataIss >> row.t >> row.a >> row.c >> row.b >> row.d >> row.J >> row.value) {
                    // 转换为0-based索引
                    row.t -= 1;
                    row.a -= 1;
                    row.c -= 1;
                    row.b -= 1;
                    row.d -= 1;
                    pnRows.push_back(row);
                } else {
                    std::cerr << "警告: pn数据行解析失败 - " << line << std::endl;
                }

            }
        }
    }
    nucleus=result.nucleus1;
    nucleus2=result.nucleus2;
    allnucleus=result.allnucleus1;
    allnucleus2=result.allnucleus2;
    if (!pnRows.empty()) {
        result.pnData = buildVValuepn1(pnRows);
    }
    if (!efc1.empty()) {
        result.efc1 = buildVValue1(efc1,1);
    }
    if (!efc2.empty()) {
        result.efc2 = buildVValue1(efc2,2);
    }
    file.close();  // 关闭文件
    return result; // 返回读取到的数据
}
