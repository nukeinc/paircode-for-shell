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



// int main() {
//     // 开始计时
//     auto start = std::chrono::high_resolution_clock::now();
//
//     double sixj = util::wigner_6j(20, 20, 24, 13, 37, 29);
//
//     // 结束计时
//     auto end = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
//
//     // 输出结果，控制不同的精度
//     std::cout << "科学计数法: " << std::scientific << std::setprecision(16) << sixj << std::endl;
//
//     // 恢复默认格式
//     std::cout << std::defaultfloat;
//     std::cout << "默认格式: " << std::setprecision(6) << sixj << std::endl;
//
//     std::cout << "计算时间: " << duration.count() / 1000000.0 << " 毫秒" << std::endl;
//
//     return 0;
// }

int main() {
    auto start = std::chrono::high_resolution_clock::now();

    // 打开输出文件
    std::ofstream outfile("sixj_table.txt");
    if (!outfile.is_open()) {
        std::cerr << "无法打开输出文件!" << std::endl;
        return 1;
    }

    outfile << std::scientific << std::setprecision(15);
    std::cout << std::scientific << std::setprecision(15);

    int count = 0;
    const int max_j = 30;  // 最大角动量值

    // 遍历所有可能的组合
    for (int j1 = 0; j1 <= max_j; ++j1) {
        for (int j2 = 0; j2 <= max_j; ++j2) {
            for (int j3 = 0; j3 <= max_j; ++j3) {
                for (int j4 = 0; j4 <= max_j; ++j4) {
                    for (int j5 = 0; j5 <= max_j; ++j5) {
                        for (int j6 = 0; j6 <= max_j; ++j6) {
                            // 计算6j符号
                            double sixj = util::wigner_6j(j1, j2, j3, j4, j5, j6);

                            // 检查是否非零（使用适当的阈值）
                            if (std::abs(sixj) > 1e-15) {
                                outfile << "6j(" << j1 << ", " << j2 << ", " << j3
                                        << ", " << j4 << ", " << j5 << ", " << j6
                                        << ") = " << sixj << std::endl;
                                count++;

                                // 每1000个非零结果在控制台输出一次进度
                                if (count % 1000 == 0) {
                                    std::cout << "已找到 " << count << " 个非零结果..." << std::endl;
                                }
                            }
                        }
                    }
                }
            }
        }
        // 显示进度
        std::cout << "已完成 j1 = " << j1 << "/" << max_j << std::endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    // 输出统计信息
    std::cout << "\n=== 计算完成 ===" << std::endl;
    std::cout << "总非零结果数: " << count << std::endl;
    std::cout << "计算时间: " << duration.count() << " 秒" << std::endl;
    std::cout << "结果已保存到 sixj_table.txt" << std::endl;

    outfile.close();
    return 0;
}

