//
// Created by wang- on 2025/1/14.
//

#ifndef WRITEOUT_H
#define WRITEOUT_H
#include <iostream>
#include <fstream>
#include"moe.h"


inline void writeMatrixToFile(const std::vector<std::vector<double>>& matrix, std::ofstream& outfile) {


    for (const auto& row : matrix) {
        for (const auto& element : row) {
            outfile << element << " ";  // 使用空格分隔元素
        }
        outfile << "\n";  // 每行矩阵数据换行
    }

    outfile.close();

}

inline std::string vecToStr(const std::vector<int>& v) {
    std::string s = "{";
    for (size_t i = 0; i < v.size(); ++i) {
        s += std::to_string(v[i]) + (i + 1 < v.size() ? "," : "");
    }
    s += "}";
    return s;
}

// 简化的 printBasisOneLinefile 函数
inline void printBasisOneLinefile(std::ofstream& outfile, const basis& b) {
    if (outfile.is_open()) {
        outfile << "j=" << b.j
                << " jn=" << b.jn
                << " r=" << vecToStr(b.r)
                << " rn=" << vecToStr(b.rn)
                << " sj=" << vecToStr(b.sj)<<std::endl
                ;
    } else {
        std::cerr << "文件未打开，无法写入数据！" << std::endl;
    }
}





#endif //WRITEOUT_H
