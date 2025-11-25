//
// Created by wang- on 25-1-7.
//
#define EIGEN_NO_DEBUG
#define EIGEN_UNROLLING_LIMIT 0
#include "../inc/moe.h"
#include <unordered_set>
#include <functional>
#include <set>
#include <tuple>
#include <map>
#include <numeric>
#include <algorithm>
#include <vector>
#include"../inc/getbasis.h"
#include <Eigen/Dense>
#include <chrono>
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Sparse>
#include <array>
#include <mutex>
#include <memory>
#include <functional>

// #include "maincal.h"
#include "../inc/moe.h"
// 生成所有合法的粒子对组合，并跟踪每个 r 的序号

bool isrpaircan(const std::vector<std::vector<int>>& r,std::vector<int> count)
{
    int m1=1;
    for (int i=1;i<r.size();++i)
    {
        std::vector<int> discuss=r[i];
        int min=discuss[0];
        int max=discuss[1];
        int rnumsum=0;
        for (int j=2;j<discuss.size();++j)
        {
            rnumsum+=count[discuss[j]];
        }
        if (rnumsum < min || rnumsum > max) {
            return false; // 任一约束不满足即返回 false
        }
    }
    return true; // 合法组合
}


std::vector<std::pair<std::vector<int>, std::vector<int>>> generateValidPairs(int m,
    const std::vector<std::vector<int>>& r, int maxSum) {
    std::vector<std::pair<std::vector<int>, std::vector<int>>> allResults; // 存储所有合法配对及其对应的序号
    std::vector<int> currentR;                // 当前配对的耦合结果
    std::vector<int> currentRN;               // 当前配对的序号
    std::vector<int> count(r[0].size(), 0);   // 每种 r 值的计数

    // 回溯生成配对方案
    std::function<void(int, int)> backtrack = [&](int index, int sum) {
        if (currentR.size() == m / 2) {
            if (isrpaircan(r,count) && sum <= maxSum) {
                allResults.emplace_back(currentR, currentRN); // 存储当前合法配对及其序号
            }
            return;
        }

        // 遍历 r[0] 中的所有可能值
        for (int i = 0; i < r[0].size(); ++i) {
            int last=0;
            if (currentR.size() != 0)
            {
                last=currentRN.back();
            }
            if (isrpaircan(r,count) && i>=last)
            { // 检查限制条件
                currentR.push_back(r[0][i]);  // 插入当前的 r 值
                currentRN.push_back(i);       // 插入当前的 r 序号
                count[i]++;
                backtrack(index + 1, sum + r[0][i]); // 递归处理下一步
                currentR.pop_back();                // 回溯
                currentRN.pop_back();               // 回溯
                count[i]--;
            }
        }
    };

    backtrack(0, 0);
    return allResults;
}
std::vector<basis> calculateBasisWithRange( int m,
    const std::vector<std::pair<std::vector<int>, std::vector<int>>>& validPairs) {
    std::set<basis> allBases;
    // 使用集合存储，避免重复
    int nucnum = nucleus.size();
    int dienum=0;

    // 奇数核子情况：单核子角动量
    if (m % 2 == 1) {
        for (int i = 0; i < nucnum; i++) { // 遍历所有可能的单核子角动量
            for (const auto& [r, rn] : validPairs) { // 遍历所有合法的 r 组合及其序号
                std::function<void(std::vector<int>&, int)> backtrackSj;
                int parity;
                int j = nucleus[i].j; // 提取核子角动量
                parity=nucleus[i].parity;
                basis b;       // 当前 basis 结构
                b.j = j;       // 单粒子角动量
                b.jn = i;      // 单粒子的序号
                b.r = r;       // 存储粒子对角动量
                for (int element : rn)
                {
                    parity=parity*rparityvec[element]; // 计算 parity
                }
                b.rn = rn;     // 存储对应的 r 序号
                std::vector<int> sj; // 存储当前耦合角动量序列
                b.parity=parity;

                // 回溯生成所有 sj 的可能组合
                backtrackSj = [&](std::vector<int>& currentSj, int step) {
                    if (step == r.size()) { // 如果所有 r 已经处理完，存储结果
                        b.sj = currentSj;
                        allBases.insert(b); // 存储到集合中，自动去重
                        // dienum++;
                        // std::cout<<"dienum="<<dienum<<std::endl;
                        // printBasisOneLine(b);
                        // std::cout<<"now allbases"<<std::endl;
                        //
                        // int size=allBases.size();
                        // std::cout<<"size="<<size<<std::endl;
                        return;
                    }

                    // 第一个 sj 的值直接等于第一个 r
                    if (step == 0) {
                        int presj = b.j;
                        int rValue = r[step];
                        int minSj=std::abs(presj-rValue);
                        int maxSj = presj + rValue;
                        for (int newSj = minSj; newSj <= maxSj; newSj=newSj+2)
                        {
                            currentSj.push_back(newSj);     // 添加当前值
                            backtrackSj(currentSj, step + 1); // 递归处理下一步
                            currentSj.pop_back();          // 回溯
                        }
                    } else {
                        // 计算 sj[step] 的可能范围
                        int prevSj = currentSj.back(); // 上一个 sj
                        int rValue = r[step];

                        int minSj = std::abs(prevSj - rValue);
                        int maxSj = prevSj + rValue;

                        // 遍历可能的 sj[step] 值
                        for (int newSj = minSj; newSj <= maxSj; newSj=newSj+2) {
                            currentSj.push_back(newSj);     // 添加当前值
                            backtrackSj(currentSj, step + 1); // 递归处理下一步
                            currentSj.pop_back();          // 回溯
                        }
                    }
                };

                std::vector<int> currentSj;
                backtrackSj(currentSj, 0);
            }
        }
    } else { // 偶数核子情况：不考虑单核子
        for (const auto& [r, rn] : validPairs) { // 遍历所有合法的 r 组合及其序号
            std::function<void(std::vector<int>&, int)> backtrackSj;

            basis b;       // 当前 basis 结构
            b.j = 0;       // 没有单核子时 j 设为 0
            b.jn = 0;      // 没有单核子时序号设为 0
            b.r = r;       // 存储粒子对角动量
            b.rn = rn;     // 存储对应的 r 序号
            std::vector<int> sj; // 存储当前耦合角动量序列
            int parcal=1;
            for (int element : rn)
            {
                parcal=parcal*rparityvec[element];
            }
            b.parity=parcal; // 计算 parity

            // 回溯生成所有 sj 的可能组合
            backtrackSj = [&](std::vector<int>& currentSj, int step) {
                if (step == r.size()) { // 如果所有 r 已经处理完，存储结果
                    b.sj = currentSj;
                    // if ((b.sj.back()/2)%2==0)
                    // {
                    //     allBases.insert(b);
                    // } // 存储到集合中，自动去重
                    allBases.insert(b);
                    return;
                }

                // 第一个 sj 的值直接等于第一个 r
                if (step == 0) {
                    currentSj.push_back(r[0]);
                    backtrackSj(currentSj, step + 1);
                    currentSj.pop_back(); // 回溯
                } else {

                    // 计算 sj[step] 的可能范围
                    int prevSj = currentSj.back(); // 上一个 sj
                    int rValue = r[step];

                    int minSj = std::abs(prevSj - rValue);
                    int maxSj = prevSj + rValue;

                    // 遍历可能的 sj[step] 值
                    for (int newSj = minSj; newSj <= maxSj; newSj=newSj+2) {
                        if (step==1)
                        {
                            // if ((newSj/2)%2==1)
                            // {
                            //     continue;
                            // }
                        }
                        currentSj.push_back(newSj);     // 添加当前值
                        backtrackSj(currentSj, step + 1); // 递归处理下一步
                        currentSj.pop_back();          // 回溯
                    }
                }
            };

            std::vector<int> currentSj;
            backtrackSj(currentSj, 0);
        }
    }
    std::vector<basis>allbasisc(allBases.begin(), allBases.end());
    // std::cout <<"dienum="<< dienum << std::endl;
    return allbasisc;
}
// 计算所有可能的 basis
// std::vector<basis> calculateBasisWithRange( int m,
//     const std::vector<std::pair<std::vector<int>, std::vector<int>>>& validPairs) {
//     std::set<basis> allBases;
//     // 使用集合存储，避免重复
//     int nucnum = nucleus.size();
//
//     // 奇数核子情况：单核子角动量
//     if (m % 2 == 1) {
//         for (int i = 0; i < nucnum; i++) { // 遍历所有可能的单核子角动量
//             for (const auto& [r, rn] : validPairs) { // 遍历所有合法的 r 组合及其序号
//                 std::function<void(std::vector<int>&, int)> backtrackSj;
//                 int j = nucleus[i].j; // 提取核子角动量
//
//                 basis b;       // 当前 basis 结构
//                 b.j = j;       // 单粒子角动量
//                 b.jn = i;      // 单粒子的序号
//                 b.r = r;       // 存储粒子对角动量
//                 b.rn = rn;     // 存储对应的 r 序号
//                 std::vector<int> sj; // 存储当前耦合角动量序列
//
//                 // 回溯生成所有 sj 的可能组合
//                 backtrackSj = [&](std::vector<int>& currentSj, int step) {
//                     if (step == r.size()) { // 如果所有 r 已经处理完，存储结果
//                         b.sj = currentSj;
//                         allBases.insert(b); // 存储到集合中，自动去重
//                         return;
//                     }
//
//                     // 第一个 sj 的值直接等于第一个 r
//                     if (step == 0) {
//                         int presj = b.j;
//                         int rValue = r[step];
//                         int minSj=std::abs(presj-rValue);
//                         int maxSj = presj + rValue;
//                         for (int newSj = minSj; newSj <= maxSj; newSj=newSj+2)
//                         {
//                             currentSj.push_back(newSj);     // 添加当前值
//                             backtrackSj(currentSj, step + 1); // 递归处理下一步
//                             currentSj.pop_back();          // 回溯
//                         }
//                     } else {
//                         // 计算 sj[step] 的可能范围
//                         int prevSj = currentSj.back(); // 上一个 sj
//                         int rValue = r[step];
//
//                         int minSj = std::abs(prevSj - rValue);
//                         int maxSj = prevSj + rValue;
//
//                         // 遍历可能的 sj[step] 值
//                         for (int newSj = minSj; newSj <= maxSj; newSj=newSj+2) {
//                             currentSj.push_back(newSj);     // 添加当前值
//                             backtrackSj(currentSj, step + 1); // 递归处理下一步
//                             currentSj.pop_back();          // 回溯
//                         }
//                     }
//                 };
//
//                 std::vector<int> currentSj;
//                 backtrackSj(currentSj, 0);
//             }
//         }
//     } else { // 偶数核子情况：不考虑单核子
//         for (const auto& [r, rn] : validPairs) { // 遍历所有合法的 r 组合及其序号
//             std::function<void(std::vector<int>&, int)> backtrackSj;
//
//             basis b;       // 当前 basis 结构
//             b.j = 0;       // 没有单核子时 j 设为 0
//             b.jn = 0;      // 没有单核子时序号设为 0
//             b.r = r;       // 存储粒子对角动量
//             b.rn = rn;     // 存储对应的 r 序号
//             std::vector<int> sj; // 存储当前耦合角动量序列
//
//             // 回溯生成所有 sj 的可能组合
//             backtrackSj = [&](std::vector<int>& currentSj, int step) {
//                 if (step == r.size()) { // 如果所有 r 已经处理完，存储结果
//                     b.sj = currentSj;
//                     // if ((b.sj.back()/2)%2==0)
//                     // {
//                     //     allBases.insert(b);
//                     // } // 存储到集合中，自动去重
//                     allBases.insert(b);
//                     return;
//                 }
//
//                 // 第一个 sj 的值直接等于第一个 r
//                 if (step == 0) {
//                     currentSj.push_back(r[0]);
//                     backtrackSj(currentSj, step + 1);
//                     currentSj.pop_back(); // 回溯
//                 } else {
//
//                     // 计算 sj[step] 的可能范围
//                     int prevSj = currentSj.back(); // 上一个 sj
//                     int rValue = r[step];
//
//                     int minSj = std::abs(prevSj - rValue);
//                     int maxSj = prevSj + rValue;
//
//                     // 遍历可能的 sj[step] 值
//                     for (int newSj = minSj; newSj <= maxSj; newSj=newSj+2) {
//                         if (step==1)
//                         {
//                             if ((newSj/2)%2==1)
//                             {
//                                 continue;
//                             }
//                         }
//                         currentSj.push_back(newSj);     // 添加当前值
//                         backtrackSj(currentSj, step + 1); // 递归处理下一步
//                         currentSj.pop_back();          // 回溯
//                     }
//                 }
//             };
//
//             std::vector<int> currentSj;
//             backtrackSj(currentSj, 0);
//         }
//     }
//     std::vector<basis>allbasisc(allBases.begin(), allBases.end());
//     return allbasisc;
// }

void getystrbasis(
    std::vector<std::vector<std::vector<double>>>& ystr1,
    const std::vector<std::vector<double>>& ystrget,
    const std::vector<int>& ro1) {

    // 遍历 ro1
    for (size_t i = 0; i < ro1.size(); ++i) {
        int firstDimValue = ro1[i]; // 取出 ro1 的值

        // 遍历 ystrget 找到第三列等于 firstDimValue 的所有行
        for (const auto& row : ystrget) {
            if (static_cast<int>(row[2]) == firstDimValue) {
                // 根据第一列和第二列确定三维数组的第二维和第三维
                int secondDim = static_cast<int>(row[0]);
                int thirdDim = static_cast<int>(row[1]);

                // 将第四列的值赋值到 ystr1
                ystr1[i][secondDim][thirdDim] = row[3]; // 直接赋值浮点数
            }
        }
    }
}



 void calculateystr(
    std::vector<std::vector<std::vector<std::vector<double>>>>&ystrall,
    const std::vector<std::vector<double>>& ystrget,
    std::vector<basis> allbasis)
{
    int nb=allbasis.size();
    for (int i=0;i<nb ;i++)
    {
        std::vector<int> ro1=allbasis[i].rn;
        getystrbasis(ystrall[i],ystrget,ro1);
    }
}

double coupleoverlap(const std::vector<std::vector<double>> overlapmatn,
    const std::vector<std::vector<double>> overlapmatp,
    CoupledBasis cpball1,CoupledBasis cpball2)
{
    double cpoverlap = 0.0;
    if (cpball1.J!=cpball2.J)
    {
        return 0.0;
    }
    cpoverlap=overlapmatn[cpball1.ni][cpball2.ni]*overlapmatp[cpball1.pi][cpball2.pi];
    return cpoverlap;
}

std::vector<std::vector<double>> cpomatcal(const std::vector<std::vector<double>> overlapmatn,
    const std::vector<std::vector<double>> overlapmatp,std::vector<CoupledBasis>& inputBases)
{
    int n = inputBases.size();
    std::vector<std::vector<double>> cpomatcal(n,std::vector<double>(n,0));
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<n;j++)
        {
            double cpoverlap=0.0;
            cpoverlap=coupleoverlap(overlapmatn,overlapmatp,inputBases[i],inputBases[j]);
            cpomatcal[i][j]=cpoverlap;
        }
    }
    return cpomatcal;
}



// std::vector<std::vector<double>> gramSchmidtInOverlap(
//     const std::vector<std::vector<double>>& cpomat,
//     std::vector<basis>& allbasisp,
//     std::vector<std::vector<std::vector<std::vector<double>>>>& ystrallp)
//
// {
//     int n = (int)cpomat.size();
//     if(n == 0) return {};
//
//     // 检查 cpomat 是否是方阵
//     for (int i = 0; i < n; i++){
//         if ((int)cpomat[i].size() != n) {
//             throw std::runtime_error("gramSchmidtInOverlap: cpomat not square.");
//         }
//     }
//
//     // 初始化 C 为单位阵 (每个新基先等于对应的旧基)
//     std::vector<std::vector<double>> C(n, std::vector<double>(n, 0.0));
//     for(int i = 0; i < n; i++){
//         C[i][i] = 1.0;  // 第 i 行表示: |\psi_i> 初始就是 |\phi_i>
//     }
//
//     // 开始做类似于 "v_k -= sum_{j<k} (v_k·u_j) * u_j" 的过程
//     // 只是现在内积是由 cpomat 给出的,
//     // 并且 v_k, u_j 通过系数阵 C 的行表示.
//     for(int k = 0; k < n; k++){
//         // (1) 对 k行, 从 j=0..k-1 做投影减法
//         for(int j = 0; j < k; j++){
//             // overlap_kj = <psi_k | psi_j>
//             //            = sum_{alpha,beta} C[k][alpha] * cpomat[alpha][beta] * C[j][beta]
//             double overlap_kj = 0.0;
//             for(int alpha = 0; alpha < n; alpha++){
//                 for(int beta = 0; beta < n; beta++){
//                     overlap_kj += C[k][alpha] * cpomat[alpha][beta] * C[j][beta];
//                 }
//             }
//             // 从 k行向量中减去 overlap_kj * (j行向量)
//             for(int alpha = 0; alpha < n; alpha++){
//                 C[k][alpha] -= overlap_kj * C[j][alpha];
//             }
//         }
//
//         // (2) 计算 k行向量的范数(norm_k) = sqrt( <psi_k | psi_k> )
//         double norm_k_sq = 0.0;
//         for(int alpha = 0; alpha < n; alpha++){
//             for(int beta = 0; beta < n; beta++){
//                 norm_k_sq += C[k][alpha] * cpomat[alpha][beta] * C[k][beta];
//             }
//         }
//         if(norm_k_sq < 1e-14){
//             // 说明该向量几乎为0, 无法继续正交化
//             // 你可以选择跳过或报错
//             // throw std::runtime_error("Gram-Schmidt: found near-zero vector => linear dependence or singular overlap.");
//             std::cout<<"=0dek="<<k<<std::endl;
//         }
//         double norm_k = std::sqrt(norm_k_sq);
//
//         // (3) 归一化:  C[k][alpha] /= norm_k
//         for(int alpha = 0; alpha < n; alpha++){
//             C[k][alpha] /= norm_k;
//         }
//     }
//     return C;
// }

std::vector<std::vector<double>> gramSchmidtInOverlap(
    std::vector<std::vector<double>>& cpomat,  // 允许修改 cpomat
    std::vector<basis>& allbasisp,
    std::vector<std::vector<std::vector<std::vector<double>>>>& ystrallp)
{
    int n = (int)cpomat.size();
    if (n == 0) return {};

    // 检查 cpomat 是否是方阵
    for (int i = 0; i < n; i++) {
        if ((int)cpomat[i].size() != n) {
            throw std::runtime_error("gramSchmidtInOverlap: cpomat not square.");
        }
    }

    // 初始化 C 为单位阵 (每个新基先等于对应的旧基)
    std::vector<std::vector<double>> C(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        C[i][i] = 1.0;  // 第 i 行表示: |\psi_i> 初始就是 |\phi_i>
    }

    // 存储有效向量的索引以及无效的基矢索引
    std::vector<int> validIndices;
    std::vector<int> invalidIndices;

    // 开始做类似于 "v_k -= sum_{j<k} (v_k·u_j) * u_j" 的过程
    for (int k = 0; k < n; k++) {
        // (1) 对 k行, 从 j=0..k-1 做投影减法
        for (int j = 0; j < k; j++) {
            // 如果 j 行已经被视为线性相关，跳过该行的投影计算
            if (std::find(validIndices.begin(), validIndices.end(), j) == validIndices.end()) {
                continue;
            }

            double overlap_kj = 0.0;
            for (int alpha = 0; alpha < n; alpha++) {
                for (int beta = 0; beta < n; beta++) {
                    overlap_kj += C[k][alpha] * cpomat[alpha][beta] * C[j][beta];
                }
            }
            // 从 k行向量中减去 overlap_kj * (j行向量)
            for (int alpha = 0; alpha < n; alpha++) {
                C[k][alpha] -= overlap_kj * C[j][alpha];
            }
        }

        // (2) 计算 k行向量的范数(norm_k)
        double norm_k_sq = 0.0;
        for (int alpha = 0; alpha < n; alpha++) {
            for (int beta = 0; beta < n; beta++) {
                norm_k_sq += C[k][alpha] * cpomat[alpha][beta] * C[k][beta];
            }
        }

        // 判断是否为零向量（范数接近零）
        if (norm_k_sq < 1e-14) {  // 如果接近零
            std::cerr << "Warning: Found near-zero vector at k=" << k << ", possibly due to linear dependence.\n";

            // 记录无效基矢的索引
            invalidIndices.push_back(k);

            continue;  // 跳过该向量
        }

        // 计算范数并归一化
        double norm_k = std::sqrt(norm_k_sq);
        for (int alpha = 0; alpha < n; alpha++) {
            C[k][alpha] /= norm_k;
        }

        // 记录有效的基向量索引
        validIndices.push_back(k);
    }

    // 删除无效向量，构建最终的正交基
    std::vector<std::vector<double>> orthogonalBasis;
    for (int idx : validIndices) {
        orthogonalBasis.push_back(C[idx]);  // 只保存有效的向量
    }
    std::cout << "fffff3 " << std::endl;
    // 统一删除无效的基矢及对应的数据
    for (int i = invalidIndices.size() - 1; i >= 0; i--) {
        std::cout << "fffff3i "<< i << std::endl;
        int k = invalidIndices[i];

        // 删除 allbasisp 中的基矢
        allbasisp.erase(allbasisp.begin() + k);

        ystrallp.erase(ystrallp.begin() + k);
        // 删除 ystrallp 中对应的第 0 维部分


        // 删除 cpomat 中的行列
        cpomat.erase(cpomat.begin() + k);
        int xxx=cpomat.size();

        for (int j = 0; j < cpomat.size(); j++) {

            cpomat[j].erase(cpomat[j].begin() + k);
            // orthogonalBasis[j].erase(orthogonalBasis[j].begin() + k);
        }
        for (int jj = 0; jj < orthogonalBasis.size(); jj++) {


            orthogonalBasis[jj].erase(orthogonalBasis[jj].begin() + k);
        }

    }


    // 删除无效向量，构建最终的正交基


    return orthogonalBasis;
}

std::vector<std::vector<double>> gramSchmidtInOverlap(
    std::vector<std::vector<double>>& cpomat,
    std::vector<basis>& allbasisp,
    std::vector<std::vector<std::vector<std::vector<double>>>>& ystrallp,
    const std::vector<std::pair<int, int>>& blocks)  // 块定义为 [start, end)
{
    int n = (int)cpomat.size();
    if (n == 0) return {};

    // 检查 cpomat 是否是方阵
    for (int i = 0; i < n; i++) {
        if ((int)cpomat[i].size() != n) {
            throw std::runtime_error("gramSchmidtInOverlap: cpomat not square.");
        }
    }

    std::vector<std::vector<double>> C(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        C[i][i] = 1.0;
    }

    std::vector<int> validIndices;
    std::vector<int> invalidIndices;

    // 预计算每个向量属于哪个块
    std::vector<int> block_index(n, -1);
    for (int b = 0; b < blocks.size(); b++) {
        int start = blocks[b].first;
        int end = blocks[b].second;
        for (int i = start; i < end; i++) {  // 注意：i < end
            block_index[i] = b;
        }
    }

    // 对每个块单独处理
    for (const auto& block : blocks) {
        int start = block.first;
        int end = block.second;

        for (int k = start; k < end; k++) {  // 注意：k < end
            // 投影减法
            for (int j : validIndices) {
                // 只处理同一块内的向量
                if (block_index[j] != block_index[k]) continue;

                double overlap_kj = 0.0;
                // 只计算块内的重叠积分
                for (int alpha = start; alpha < end; alpha++) {  // alpha < end
                    for (int beta = start; beta < end; beta++) {  // beta < end
                        overlap_kj += C[k][alpha] * cpomat[alpha][beta] * C[j][beta];
                    }
                }

                // 向量减法
                for (int alpha = 0; alpha < n; alpha++) {
                    C[k][alpha] -= overlap_kj * C[j][alpha];
                }
            }

            // 计算范数（只在块内计算）
            double norm_k_sq = 0.0;
            for (int alpha = start; alpha < end; alpha++) {  // alpha < end
                for (int beta = start; beta < end; beta++) {  // beta < end
                    norm_k_sq += C[k][alpha] * cpomat[alpha][beta] * C[k][beta];
                }
            }

            if (norm_k_sq < 1e-10) {
                std::cerr << "Warning: Found near-zero vector at k=" << k << ", possibly due to linear dependence.\n";
                invalidIndices.push_back(k);
                continue;
            }

            double norm_k = std::sqrt(norm_k_sq);
            for (int alpha = 0; alpha < n; alpha++) {
                C[k][alpha] /= norm_k;
            }
            validIndices.push_back(k);
        }
    }

    // 构建正交基
    std::vector<std::vector<double>> orthogonalBasis;
    for (int idx : validIndices) {
        orthogonalBasis.push_back(C[idx]);
    }

    // 删除无效基矢（从后往前删除以避免索引问题）
    std::sort(invalidIndices.begin(), invalidIndices.end(), std::greater<int>());
    for (int k : invalidIndices) {
        allbasisp.erase(allbasisp.begin() + k);
        ystrallp.erase(ystrallp.begin() + k);
        cpomat.erase(cpomat.begin() + k);
        for (auto& row : cpomat) {
            row.erase(row.begin() + k);
        }
        for (auto& vec : orthogonalBasis) {
            vec.erase(vec.begin() + k);
        }
    }

    return orthogonalBasis;
}

std::vector<std::pair<int, int>> findblocks(std::vector<basis>& allbasisp)
{
    // 找到所有的块
    std::vector<std::pair<int, int>> blocks;
    if (allbasisp.empty()) return blocks;  // 如果没有基矢，返回空块

    int start = 0;
    int currentJ = allbasisp[0].sj.back();  // 初始 J 值

    for (int i = 0; i < allbasisp.size(); i++) {
        if (allbasisp[i].sj.back() != currentJ) {
            blocks.emplace_back(start, i);  // 添加当前块
            start = i;                      // 更新起始位置
            currentJ = allbasisp[i].sj.back();      // 更新 J 值
        }
    }
    blocks.emplace_back(start, allbasisp.size());  // 添加最后一个块

    return blocks;
}


// std::vector<std::vector<double>> gramSchmidtInOverlap(
//     std::vector<std::vector<double>>& cpomat,
//     std::vector<basis>& allbasisp,
//     std::vector<std::vector<std::vector<std::vector<double>>>>& ystrallp)
// {
//     int n = cpomat.size();
//     if (n == 0) return {};
//     std::vector<std::vector<double>> C1(n, std::vector<double>(n, 0.0));
//     for (int i = 0; i < n; ++i) {
//         C1[i][i] = 1.0;
//     }
//
//     std::vector<int> invalidIndices;
//
//     for (int k = 0; k < n; ++k) {
//         for (int j = 0; j < k; ++j) {
//             double overlap = 0.0;
//             for (int a = 0; a < n; ++a) {
//                 for (int b = 0; b < n; ++b) {
//                     overlap += C1[k][a] * cpomat[a][b] * C1[j][b];
//                 }
//             }
//             for (int a = 0; a < n; ++a) {
//                 C1[k][a] -= overlap * C1[j][a];
//             }
//         }
//
//         double norm_sq = 0.0;
//         for (int a = 0; a < n; ++a) {
//             for (int b = 0; b < n; ++b) {
//                 norm_sq += C1[k][a] * cpomat[a][b] * C1[k][b];
//             }
//         }
//
//         if (norm_sq < 1e-15) {
//             invalidIndices.push_back(k);
//             continue;
//         }
//         double norm = std::sqrt(norm_sq);
//         for (int a = 0; a < n; ++a) {
//             C1[k][a] /= norm;
//         }
//
//     }
//
//     // 1. 检测并标记无效基矢（线性相关）
//     // for (int k = 0; k < n; ++k) {
//     //     double norm_sq = 0.0;
//     //     for (int i = 0; i < n; ++i) {
//     //         for (int j = 0; j < n; ++j) {
//     //             norm_sq += (i == k ? 1.0 : 0.0) * cpomat[i][j] * (j == k ? 1.0 : 0.0);
//     //         }
//     //     }
//     //     if (norm_sq < 1e-9) {
//     //         invalidIndices.push_back(k);
//     //     }
//     // }
//
//     // 2. 逆序删除无效基矢
//     for (auto it = invalidIndices.rbegin(); it != invalidIndices.rend(); ++it) {
//         int k = *it;
//         allbasisp.erase(allbasisp.begin() + k);
//         ystrallp.erase(ystrallp.begin() + k);
//         cpomat.erase(cpomat.begin() + k);
//         for (auto& row : cpomat) {
//             row.erase(row.begin() + k);
//         }
//     }
//
//     // 更新矩阵维度
//     int new_n = cpomat.size();
//     if (new_n == 0) return {};
//
//     // 3. 在缩减后的矩阵上执行正交化
//     std::vector<std::vector<double>> C(new_n, std::vector<double>(new_n, 0.0));
//     for (int i = 0; i < new_n; ++i) {
//         C[i][i] = 1.0;
//     }
//
//     std::vector<int> validIndices;
//     for (int k = 0; k < new_n; ++k) {
//         for (int j = 0; j < k; ++j) {
//             double overlap = 0.0;
//             for (int a = 0; a < new_n; ++a) {
//                 for (int b = 0; b < new_n; ++b) {
//                     overlap += C[k][a] * cpomat[a][b] * C[j][b];
//                 }
//             }
//             for (int a = 0; a < new_n; ++a) {
//                 C[k][a] -= overlap * C[j][a];
//             }
//         }
//
//         double norm_sq = 0.0;
//         for (int a = 0; a < new_n; ++a) {
//             for (int b = 0; b < new_n; ++b) {
//                 norm_sq += C[k][a] * cpomat[a][b] * C[k][b];
//             }
//         }
//
//         if (norm_sq < 1e-10) {
//             throw std::runtime_error("Unexpected linear dependence after pruning.");
//         }
//
//         double norm = std::sqrt(norm_sq);
//         for (int a = 0; a < new_n; ++a) {
//             C[k][a] /= norm;
//         }
//         validIndices.push_back(k);
//     }
//
//     // 4. 构建正交基
//     std::vector<std::vector<double>> orthogonalBasis;
//     for (int idx : validIndices) {
//         orthogonalBasis.push_back(C[idx]);
//     }
//
//     return orthogonalBasis;
// }


bool approxEqualRows(const std::vector<double>& rowA,
                     const std::vector<double>& rowB,
                     double tol = 1e-14)
{
    if (rowA.size() != rowB.size()) return false;
    for (size_t i = 0; i < rowA.size(); ++i) {
        if (std::fabs(rowA[i] - rowB[i]) > tol) {
            return false;
        }
    }
    return true;
}
std::map<int, std::vector<CoupledBasis>> coupleBases(
    const std::vector<basis>& neutronBases,
    const std::vector<basis>& protonBases,std::vector<std::vector<double> > overlapp,std::vector<std::vector<double> > overlapn
) {
    // 定义一个 map，以 J 值为键，每个键对应一个 vector
    std::map<int, std::vector<CoupledBasis>> coupledBasesMap;

    // 遍历中子和质子的基组
    for (size_t neutronIndex = 0; neutronIndex < neutronBases.size(); ++neutronIndex) {
        bool duplicatedNeutron = false;
        for (size_t i = 0; i < neutronIndex; ++i) {
            if (approxEqualRows(overlapn[neutronIndex] , overlapn[i])) {
                // 说明第 neutronIndex 行与之前的第 i 行完全相同
                duplicatedNeutron = true;
                break;
            }
        }
        // 如果发现重复，就跳过本次 neutronIndex
        if (duplicatedNeutron) {
            continue;
        }
        for (size_t protonIndex = 0; protonIndex < protonBases.size(); ++protonIndex) {
            bool duplicatedProton = false;
            for (size_t j = 0; j < protonIndex; ++j) {
                if (approxEqualRows(overlapp[protonIndex] , overlapp[j])) {
                    duplicatedProton = true;
                    break;
                }
            }
            if (duplicatedProton) {
                continue;
            }
            // 获取 sj 的最后一个值
            int sj_neutron = neutronBases[neutronIndex].sj.back();
            int sj_proton = protonBases[protonIndex].sj.back();

            // 计算耦合后的 sj 最后一个值的范围
            int minSj = std::abs(sj_neutron - sj_proton);
            int maxSj = sj_neutron + sj_proton;

            // 遍历可能的耦合结果，步长为 2
            for (int sj_coupled = minSj; sj_coupled <= maxSj; sj_coupled += 2) {
                // 创建新的 CoupledBasis
                CoupledBasis newCoupledBasis;
                newCoupledBasis.ni = neutronIndex; // 存储中子索引
                newCoupledBasis.pi = protonIndex;   // 存储质子索引
                newCoupledBasis.J = sj_coupled;             // 设置耦合后的总角动量

                if (overlapn[neutronIndex][neutronIndex]*overlapp[protonIndex][protonIndex]<0.000000001)
                {
                    continue;
                }

                // 将新的 CoupledBasis 添加到对应 J 的 vector 中
                coupledBasesMap[sj_coupled].push_back(newCoupledBasis);
            }
        }
    }

    return coupledBasesMap; // 返回 map，其中每个键是不同的 J 值
}

std::vector<std::vector<double>> invertMatrix(const std::vector<std::vector<double>>& matrix) {
    int n = matrix.size();

    // 检查矩阵是否是方阵
    for (const auto& row : matrix) {
        if (row.size() != n) {
            throw std::invalid_argument("Matrix must be square");
        }
    }

    // 创建增广矩阵
    std::vector<std::vector<double>> augMatrix(n, std::vector<double>(2 * n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            augMatrix[i][j] = matrix[i][j]; // 原始矩阵
        }
        augMatrix[i][n + i] = 1.0; // 单位矩阵
    }

    // 高斯-乔丹消元
    for (int i = 0; i < n; ++i) {
        // 找主元
        double pivot = augMatrix[i][i];
        if (pivot == 0.0) {
            throw std::runtime_error("Matrix is singular and cannot be inverted");
        }

        // 将当前行标准化
        for (int j = 0; j < 2 * n; ++j) {
            augMatrix[i][j] /= pivot;
        }

        // 消去其他行的当前列
        for (int k = 0; k < n; ++k) {
            if (k != i) {
                double factor = augMatrix[k][i];
                for (int j = 0; j < 2 * n; ++j) {
                    augMatrix[k][j] -= factor * augMatrix[i][j];
                }
            }
        }
    }

    // 提取逆矩阵
    std::vector<std::vector<double>> inverseMatrix(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            inverseMatrix[i][j] = augMatrix[i][n + j];
        }
    }

    return inverseMatrix;
}
std::vector<std::vector<double>> multiplyMatrices(
    const std::vector<std::vector<double>>& A,
    const std::vector<std::vector<double>>& B) {

    // 获取矩阵维度
    size_t m = A.size();        // A 的行数
    size_t n = A[0].size();     // A 的列数
    size_t p = B[0].size();     // B 的列数

    // 检查矩阵维度是否符合乘法规则
    if (n != B.size()) {
        throw std::invalid_argument("Matrix dimensions do not allow multiplication");
    }

    // 初始化结果矩阵
    std::vector<std::vector<double>> C(m, std::vector<double>(p, 0.0));

    // 矩阵乘法核心计算
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < p; ++j) {
            for (size_t k = 0; k < n; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return C;
}

// std::vector<std::vector<double>> multiplyMatrices(
//     const std::vector<std::vector<double>>& A,
//     const std::vector<std::vector<double>>& B)
// {
//     // 1. 增强维度校验
//     if (A.empty() || B.empty() || A[0].empty() || B[0].empty()) {
//         throw std::invalid_argument("Input matrices cannot be empty");
//     }
//
//     const size_t m = A.size();
//     const size_t n = A[0].size();
//     const size_t p = B[0].size();
//     const size_t b_rows = B.size();
//
//     // 校验矩阵维度一致性
//     if (n != b_rows) {
//         throw std::invalid_argument(
//             "A.columns(" + std::to_string(n) +
//             ") != B.rows(" + std::to_string(b_rows) + ")"
//         );
//     }
//
//     // 校验每行元素数量一致性
//     for (const auto& row : A) {
//         if (row.size() != n) {
//             throw std::invalid_argument("A contains inconsistent row lengths");
//         }
//     }
//     for (const auto& row : B) {
//         if (row.size() != p) {
//             throw std::invalid_argument("B contains inconsistent row lengths");
//         }
//     }
//
//     // 2. 优化内存布局（转置B矩阵）
//     std::vector<std::vector<double>> B_T(p, std::vector<double>(b_rows));
//     for (size_t i = 0; i < b_rows; ++i) {
//         for (size_t j = 0; j < p; ++j) {
//             B_T[j][i] = B[i][j];
//         }
//     }
//
//     // 3. 优化循环顺序（i-k-j）并启用编译器优化
//     std::vector<std::vector<double>> C(m, std::vector<double>(p, 0.0));
//     for (size_t i = 0; i < m; ++i) {
//         const auto& a_row = A[i];  // 缓存行引用
//         auto& c_row = C[i];
//         for (size_t k = 0; k < n; ++k) {
//             const double a_ik = a_row[k];  // 缓存标量值
//             const auto& b_T_row = B_T[k];  // 缓存转置行
//             for (size_t j = 0; j < p; ++j) {
//                 c_row[j] += a_ik * b_T_row[j];
//             }
//         }
//     }
//
//     return C;
// }


/*

std::vector<std::vector<double>> matall(std::vector<CoupledBasis> couplebasis,
    std::vector<std::vector<double>> slimatp,std::vector<std::vector<double>> slimatn,
    std::vector<std::vector<double>> pqmat,std::vector<std::vector<double>> nqmat,
    std::vector<basis> rjlp,std::vector<basis>rjln)
{
    int n = couplebasis.size();
    std::vector<std::vector<double>> matresult(n, std::vector<double>(n, 0.0));
    for (int k = 0; k < n; ++k)
    {
        for (int j = 0; j < n; ++j)
        {
            CoupledBasis ba1=couplebasis[k];
            CoupledBasis ba2=couplebasis[j];
            int n1=ba1.ni;
            int n2=ba2.ni;
            int p1=ba1.pi;
            int p2=ba2.pi;
            double c1=slimatn[n1][n2];
            double c2=slimatp[p1][p2];

        }
    }
}*/
std::vector<std::vector<double>> transposeMatrix(const std::vector<std::vector<double>>& matrix) {
    size_t rows = matrix.size();
    size_t cols = matrix[0].size();

    // 创建一个转置后的矩阵
    std::vector<std::vector<double>> transposed(cols, std::vector<double>(rows));

    // 填充转置矩阵
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            transposed[j][i] = matrix[i][j];
        }
    }

    return transposed;
}

std::map<int, std::vector<CoupledBasis>> coupleBases2(
    const std::vector<int>& jvecn,
    const std::vector<int>& jvecp,
    const std::vector<int>& parityvecp,
    const std::vector<int>& parityvecn
) {
    // 定义一个 map，以 J 值为键，每个键对应一个 vector
    std::map<int, std::vector<CoupledBasis>> coupledBasesMap;

    // 遍历中子和质子的基组
    for (size_t neutronIndex = 0; neutronIndex < jvecn.size(); ++neutronIndex) {
        for (size_t protonIndex = 0; protonIndex < jvecp.size(); ++protonIndex) {
            bool duplicatedProton = false;
            // 获取 sj 的最后一个值
            int sj_neutron = jvecn[neutronIndex];
            int sj_proton = jvecp[protonIndex];
            int parn=parityvecn[neutronIndex];
            int parp = parityvecp[protonIndex];
            int par= parn * parp; // 耦合后的宇称

            // 计算耦合后的 sj 最后一个值的范围
            int minSj = std::abs(sj_neutron - sj_proton);
            int maxSj = sj_neutron + sj_proton;
            int mx = std::max(minSj,-4); // 最大值
            int mn = std::min(maxSj, 4); // 最小值
            // 遍历可能的耦合结果，步长为 2
            for (int sj_coupled = mx; sj_coupled <= mn; sj_coupled += 2) {
                // 创建新的 CoupledBasis
                CoupledBasis newCoupledBasis;
                newCoupledBasis.ni = neutronIndex; // 存储中子索引
                newCoupledBasis.pi = protonIndex;   // 存储质子索引
                newCoupledBasis.J = sj_coupled;             // 设置耦合后的总角动量


                // 将新的 CoupledBasis 添加到对应 J 的 vector 中
                coupledBasesMap[sj_coupled*par].push_back(newCoupledBasis);
            }
        }
    }

    return coupledBasesMap; // 返回 map，其中每个键是不同的 J 值
}

// std::vector<basis> mktransfer( int n1,int n2,int getj,std::vector<basis>basisori)
// {
//     std::vector<basis> newBasis;
//     for (size_t i = 0; i < basisori.size(); ++i)
//     {
//         int jj_i=basisori[i].sj.back();
//         std::vector<int> li1 = genmvec(jj_i, getj);
//
//     }
// }
void writeMatrixToFile(const Eigen::MatrixXd& matrix) {
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    if (outfile.is_open()) {
        outfile << matrix << std::endl;  // Eigen 原生支持的流输出
        outfile.close();
    } else {
        std::cerr << "Error: Unable to open file "  << std::endl;
    }
}

// void writeTensor3ToFile(const Eigen::Tensor<double, 3>& tensor) {
//     if (!outfile.is_open())
//     {
//         // 如果文件未打开，则尝试打开文件
//         outfile.open("basis_output.txt",std::ios::app);
//     }
//
//
//     outfile << "=== Tensor3 Output ===\n";
//     for (int i = 0; i < tensor.dimension(2); ++i) {
//         outfile << "Slice " << i << ":\n";
//         for (int j = 0; j < tensor.dimension(0); ++j) {
//             for (int k = 0; k < tensor.dimension(1); ++k) {
//                 outfile << tensor(j,k,i) << " ";
//             }
//             outfile << "\n";
//         }
//         outfile << "\n"<<std::endl;
//     }
// }
void writeTensor3ToFile(const Eigen::Tensor<double, 3>& tensor) {
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }


    outfile << "=== Tensor3 Output ===\n";
    for (int i = 0; i < tensor.dimension(0); ++i) {
        for (int j = 0; j < tensor.dimension(1); ++j) {
            for (int k = 0; k < tensor.dimension(2); ++k) {
                if (tensor(i,j,k)!=0)
                {
                    outfile <<"i,j,k="<<i<<","<<j<<","<<k<<","<< tensor(i,j,k) <<"\n"<<std::endl;
                }
            }
            // outfile << "\n";
        }
        // outfile << "\n"<<std::endl;
    }
}

void writeSimpleTensorToFile(const Eigen::Tensor<double, 2>& tensor) {
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }

    // 遍历张量的每一行
    for (int i = 0; i < tensor.dimension(0); ++i) {
        // 遍历当前行的每一列
        for (int j = 0; j < tensor.dimension(1); ++j) {
            outfile << tensor(i, j) << " ";  // 用空格分隔元素
        }
        outfile << "\n";  // 每行末尾换行
    }
    outfile.close();
}

void writeBasismToFile(const basism& bm, const std::string& filename = "basis_output.txt", bool append = true) {
    std::ofstream outfile;
    std::ios_base::openmode mode = append ? std::ios::app : std::ios::out;
    outfile.open(filename, mode);

    if (!outfile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }

    // 辅助函数：输出vector
    auto writeVector = [&outfile](const std::vector<int>& vec, const std::string& name) {
        outfile << name << ": [";
        for (size_t i = 0; i < vec.size(); ++i) {
            if (i != 0) outfile << ", ";
            outfile << vec[i];
        }
        outfile << "]\n";
    };

    // 辅助函数：输出2D vector
    auto write2DVector = [&outfile](const std::vector<std::vector<int>>& vec, const std::string& name) {
        outfile << name << ": [\n";
        for (const auto& inner : vec) {
            outfile << "  [";
            for (size_t i = 0; i < inner.size(); ++i) {
                if (i != 0) outfile << ", ";
                outfile << inner[i];
            }
            outfile << "]\n";
        }
        outfile << "]\n";
    };

    // 输出各个成员
    outfile << "--- basism entry ---\n";
    writeVector(bm.r, "r");
    writeVector(bm.rn, "rn");
    writeVector(bm.rparity, "rparity");
    outfile << "parity: " << bm.parity << "\n";
    write2DVector(bm.mvec, "mvec");
    writeVector(bm.rvecnum, "rvecnum");
    writeVector(bm.m_rvec, "m_rvec");
    writeVector(bm.bigmvec, "bigmvec");
    outfile << "\n"<<std::endl;  // 添加空行分隔不同条目

    outfile.close();
}

void writeBasismvecToFile(const std::vector<basism>& basismvec, const std::string& filename = "basis_output.txt", bool append = true) {
    std::ofstream outfile;
    std::ios_base::openmode mode = append ? std::ios::app : std::ios::out;
    outfile.open(filename, mode);

    if (!outfile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }
    int i=0;

    for (const auto& bm : basismvec) {
        outfile<< "--- basism num ="<<i<<"-----"<<std::endl;
        writeBasismToFile(bm, filename, append);
        i=i+1;
    }

    outfile.close();
}

using Vec3D = std::vector<std::vector<std::vector<int>>>;

bool is_abs_b_less_than_abs_a(const std::vector<int>& a, const std::vector<int>& b) {
    if (a.size() != b.size()) return false;  // 确保长度一致

    return std::equal(b.begin(), b.end(), a.begin(),
        [](int bi, int ai) {
            return std::abs(bi) <= std::abs(ai);  // 比较绝对值
        });
}
std::vector<int> generateArray(int m) {
    std::vector<int> result;
    for (int i = -m; i <= m; i += 2) {
        result.push_back(i);
    }
    return result;
}

std::vector<int> countFrequencies(const std::vector<int>& m, const std::vector<int>& n) {
    // 创建一个哈希表统计m中所有数字的出现频率
    std::unordered_map<int, int> freq_map;
    for (int num : m) {
        freq_map[n[num]]++;
    }

    // 为n中的每个数字查询频率，存入结果数组
    std::vector<int> result;
    result.reserve(n.size());
    for (int target : n) {
        result.push_back(freq_map[target]);
    }

    return result;
}

void backtrack(const std::vector<int>& nums, std::vector<std::vector<int>>& res,
               std::vector<int>& path, int start, int n) {
    if (path.size() == n) {
        res.push_back(path);
        return;
    }
    for (int i = start; i < nums.size(); ++i) {
        path.push_back(nums[i]);
        backtrack(nums, res, path, i, n); // 允许重复，所以传入i而非i+1
        path.pop_back();
    }
}


std::vector<std::vector<int>> uniqueCombinations(const std::vector<int>& nums, int n) {
    std::vector<int> sorted_nums = nums;
    std::sort(sorted_nums.begin(), sorted_nums.end()); // 先排序
    std::vector<std::vector<int>> res;
    std::vector<int> path;
    backtrack(sorted_nums, res, path, 0, n);
    return res;
}

// 递归构建n的三维结构
void buildNestedArrays(const Vec3D& m, Vec3D& n, size_t currentDim = 0) {
    if (currentDim >= m.size()) return;

    if (currentDim == 0) {
        // 第一维：直接复制m的第一层到n的每一组
        for (const auto& secondDim : m[currentDim]) {
            n.push_back({secondDim});  // 初始化n的第一层
        }
    } else {
        // 后续维度：将当前层数据追加到n的已有结构中
        Vec3D newN;
        for (size_t i = 0; i < n.size(); ++i) {
            for (const auto& secondDim : m[currentDim]) {
                auto temp = n[i];      // 复制现有结构
                temp.push_back(secondDim);  // 追加新层
                newN.push_back(temp);  // 保存扩展后的结构
            }
        }
        n = std::move(newN);  // 更新n
    }

    buildNestedArrays(m, n, currentDim + 1);  // 处理下一维
}


// std::vector<basism> calculateBasisform(std::vector<basis> allBases)
// {
//     std::set<basism> allbasesm;
//     int x=0;
//     for (const auto& b : allBases)
//     {
//         basism onebasism;
//         std::vector<int> rv;
//         rv.push_back(b.j);
//         rv.insert(rv.end(), b.r.begin(), b.r.end());
//         std::vector<int> rnv;
//         rnv.push_back(b.jn);
//         rnv.insert(rnv.end(), b.rn.begin(), b.rn.end());
//         std::vector<int> sj;
//         sj.push_back(b.j);
//         sj.insert(sj.end(), b.sj.begin(), b.sj.end());
//         onebasism.rvecnum=countFrequencies(b.rn, rvecall);
//         std::vector<std::vector<std::vector<int>>> mgeteveryr;
//         onebasism.parity=b.parity;
//         onebasism.r=rv;
//         onebasism.rn=b.rn;
//         onebasism.sj=sj;
//         onebasism.rn=rnv;
//         onebasism.basisn=x;
//         for (size_t i = 0; i < rvecall.size(); ++i)
//         {
//             std::vector<int> mvec=generateArray(rvecall[i]);
//             std::vector<std::vector<int>> currentm;
//             currentm=uniqueCombinations(mvec,onebasism.rvecnum[i]);
//             mgeteveryr.push_back(currentm);
//         }
//         Vec3D mforbasis;
//         buildNestedArrays(mgeteveryr, mforbasis);
//         for (size_t i = 0; i < mforbasis.size(); ++i)
//         {
//             for (int l = -onebasism.r[0]; l <=onebasism.r[0] ; ++l)
//             {
//                 basism bdie;
//                 bdie=onebasism;
//                 bdie.bigmvec.push_back(l);
//                 bdie.mvec=mforbasis[i];
//                 std::vector<int>msum;
//                 for (size_t j = 0; j < bdie.mvec.size(); ++j)
//                 {
//                     msum.push_back(std::accumulate(bdie.mvec[j].begin(), bdie.mvec[j].end(), 0));
//                 }
//                 bdie.m_rvec=msum;
//                 int msumm=0;
//
//                 for (size_t j = 0; j < bdie.mvec.size(); ++j)
//                 {
//                     for (size_t k = 0; k < bdie.mvec[j].size(); ++k)
//                     {
//                         msumm+=bdie.mvec[j][k];
//                         bdie.bigmvec.push_back(msumm);
//                     }
//                 }
//                 int sum_this = std::accumulate(bdie.m_rvec.begin(), bdie.m_rvec.end(), 0);
//                 // 根据 r[0] 的值进行筛选
//                 if (!bdie.r.empty()) {  // 确保 r 不为空
//                     if (is_abs_b_less_than_abs_a(bdie.sj,bdie.bigmvec))
//                     {
//                         if (bdie.r[0] == 0) {
//                             // 保留 m_rvec 和为 0 的项
//                             if (sum_this == 0 ) {
//                                 allbasesm.insert(bdie);  // 和为 0 的排在前面
//                             }
//                         } else {
//                             // 保留 m_rvec 和为 1 的项
//                             if (sum_this == 1 ) {
//                                 allbasesm.insert(bdie);  // 和为 1 的排在前面
//                             }
//                         }
//                     }
//                 }
//             }
//         }
//         x++;
//     }
//     std::vector<basism>allbasisc(allbasesm.begin(), allbasesm.end());
//     return allbasisc;
// }

std::vector<basism> calculateBasisform(int nucnum,std::vector<std::vector<int>> rorderrp)
{
    std::set<basism> allbasesm;
    std::vector<std::pair<std::vector<int>, std::vector<int>>> catchpair1;
    catchpair1=generateValidPairs(nucnum,rorderrp,51);
    for (const auto& b : catchpair1)
    {
        std::vector<int> nuclearvec={0};
        if (nucnum%2 != 0 ) {
            nuclearvec={};
            for (auto & i:nucleus)
            {
                nuclearvec.push_back(i.j);
            }
        }
        int m=0;
        int partyorder=0;
        for (const auto& jval : nuclearvec)
        {

            int parity=1;
            if (nucnum%2 != 0 )
            {
                parity=nucleus[partyorder].parity;
            }
            basism onebasism;
            std::vector<int> rv;
            rv.push_back(jval);
            std::vector<int>b1=b.first;
            std::vector<int>b2=b.second;
            for (int i : b2)
            {
                parity=rparityvec[i]*parity;
            }
            onebasism.parity=parity;
            rv.insert(rv.end(), b1.begin(), b1.end());
            std::vector<int> rnv;
            rnv.push_back(m);
            rnv.insert(rnv.end(), b2.begin(), b2.end());
            int n=rvecall.size();
            std::vector<int> arr;
            arr.reserve(n); // 预分配空间
            for (int i = 0; i < n; ++i) {
                arr.push_back(i);
            }
            onebasism.rvecnum=countFrequencies(b2, arr);
            std::vector<std::vector<std::vector<int>>> mgeteveryr;
            onebasism.r=rv;
            onebasism.rn=b2;
            onebasism.rn=rnv;
            for (size_t i = 0; i < rvecall.size(); ++i)
            {
                std::vector<int> mvec=generateArray(rvecall[i]);
                std::vector<std::vector<int>> currentm;
                currentm=uniqueCombinations(mvec,onebasism.rvecnum[i]);
                mgeteveryr.push_back(currentm);
            }
            Vec3D mforbasis;
            buildNestedArrays(mgeteveryr, mforbasis);
            for (size_t i = 0; i < mforbasis.size(); ++i)
            {
                for (int l = -onebasism.r[0]; l <=onebasism.r[0] ; ++l)
                {
                    basism bdie;
                    bdie=onebasism;
                    bdie.bigmvec.push_back(l);
                    bdie.mvec=mforbasis[i];
                    std::vector<int>msum;
                    msum.push_back(l);
                    for (size_t j = 0; j < bdie.mvec.size(); ++j)
                    {
                        msum.push_back(std::accumulate(bdie.mvec[j].begin(), bdie.mvec[j].end(), 0));
                    }
                    bdie.m_rvec=msum;
                    int msumm=l;

                    for (size_t j = 0; j < bdie.mvec.size(); ++j)
                    {
                        for (size_t k = 0; k < bdie.mvec[j].size(); ++k)
                        {
                            msumm+=bdie.mvec[j][k];
                            bdie.bigmvec.push_back(msumm);
                        }
                    }
                    int sum_this = std::accumulate(bdie.m_rvec.begin(), bdie.m_rvec.end(), 0);
                    // 根据 r[0] 的值进行筛选
                    if (!bdie.r.empty()) {  // 确保 r 不为空


                        if (bdie.r[0] == 0) {
                            // 保留 m_rvec 和为 0 的项
                            if (sum_this == 0 ) {
                                allbasesm.insert(bdie);  // 和为 0 的排在前面
                            }
                        } else {
                            // 保留 m_rvec 和为 1 的项
                            if (sum_this == 1 ) {
                                allbasesm.insert(bdie);  // 和为 1 的排在前面
                            }
                        }

                    }
                }
            }
            m=m+1;
            partyorder=partyorder+1;

        }
    }
    std::vector<basism>allbasisc(allbasesm.begin(), allbasesm.end());
    return allbasisc;
}


std::vector<basism> calculateBasisform_1(int nucnum,std::vector<std::vector<int>>rorderrp)
{
    std::set<basism> allbasesm;
    std::vector<std::pair<std::vector<int>, std::vector<int>>> catchpair1;
    catchpair1=generateValidPairs(nucnum,rorderrp,51);
    for (const auto& b : catchpair1)
    {
        std::vector<int> nuclearvec={0};
        if (nucnum%2 != 0 ) {
            nuclearvec={};
            for (auto & i:nucleus)
            {
                nuclearvec.push_back(i.j);
            }
        }
        int m=0;
        int partyorder=0;
        for (const auto& jval : nuclearvec)
        {

            int parity=1;
            if (nucnum%2 != 0 )
            {
                parity=nucleus[partyorder].parity;
            }
            basism onebasism;
            std::vector<int> rv;
            rv.push_back(jval);
            std::vector<int>b1=b.first;
            std::vector<int>b2=b.second;
            for (int i : b2)
            {
                parity=rparityvec[i]*parity;
            }
            onebasism.parity=parity;
            rv.insert(rv.end(), b1.begin(), b1.end());
            std::vector<int> rnv;
            rnv.push_back(m);
            rnv.insert(rnv.end(), b2.begin(), b2.end());
            int n=rvecall.size();
            std::vector<int> arr;
            arr.reserve(n); // 预分配空间
            for (int i = 0; i < n; ++i) {
                arr.push_back(i);
            }
            onebasism.rvecnum=countFrequencies(b2, arr);
            std::vector<std::vector<std::vector<int>>> mgeteveryr;
            onebasism.r=rv;
            onebasism.rn=b2;
            onebasism.rn=rnv;
            for (size_t i = 0; i < rvecall.size(); ++i)
            {
                std::vector<int> mvec=generateArray(rvecall[i]);
                std::vector<std::vector<int>> currentm;
                currentm=uniqueCombinations(mvec,onebasism.rvecnum[i]);
                mgeteveryr.push_back(currentm);
            }
            Vec3D mforbasis;
            buildNestedArrays(mgeteveryr, mforbasis);
            for (size_t i = 0; i < mforbasis.size(); ++i)
            {
                for (int l = -onebasism.r[0]; l <=onebasism.r[0] ; ++l)
                {
                    basism bdie;
                    bdie=onebasism;
                    bdie.bigmvec.push_back(l);
                    bdie.mvec=mforbasis[i];
                    std::vector<int>msum;
                    msum.push_back(l);
                    for (size_t j = 0; j < bdie.mvec.size(); ++j)
                    {
                        msum.push_back(std::accumulate(bdie.mvec[j].begin(), bdie.mvec[j].end(), 0));
                    }
                    bdie.m_rvec=msum;
                    int msumm=l;

                    for (size_t j = 0; j < bdie.mvec.size(); ++j)
                    {
                        for (size_t k = 0; k < bdie.mvec[j].size(); ++k)
                        {
                            msumm+=bdie.mvec[j][k];
                            bdie.bigmvec.push_back(msumm);
                        }
                    }
                    int sum_this = std::accumulate(bdie.m_rvec.begin(), bdie.m_rvec.end(), 0);
                    // 根据 r[0] 的值进行筛选
                    if (!bdie.r.empty()) {  // 确保 r 不为空


                        if (bdie.r[0] == 0) {
                            // 保留 m_rvec 和为 0 的项
                            if (sum_this == 2 ) {
                                allbasesm.insert(bdie);  // 和为 0 的排在前面
                            }
                        }

                    }
                }
            }
            m=m+1;
            partyorder=partyorder+1;

        }
    }
    std::vector<basism>allbasisc(allbasesm.begin(), allbasesm.end());
    return allbasisc;
}

using namespace std;


// void backtrack(vector<vector<int>>& res, vector<int>& nums, int start) {
//     if (start == nums.size()) {
//         res.push_back(nums);  // 找到一个排列
//         return;
//     }
//     for (int i = start; i < nums.size(); i++) {
//         swap(nums[start], nums[i]);      // 交换当前元素（允许重复数字交换）
//         backtrack(res, nums, start + 1); // 递归处理后续位置
//         swap(nums[start], nums[i]);      // 恢复交换（回溯）
//     }
// }
//
// vector<vector<int>> permuteUnique(vector<int>& nums) {
//     vector<vector<int>> res;
//     backtrack(res, nums, 0);
//     return res;
// }

// void backtrack(std::vector<std::vector<int>>& res, std::vector<int>& nums, int start) {
//     if (start == nums.size()) {
//         res.push_back(nums);  // 找到一个排列
//         return;
//     }
//     for (int i = start; i < nums.size(); i++) {
//         // 如果当前数字已经在 start 位置交换过，跳过
//         if (i > start && nums[i] == nums[start]) {
//             continue;
//         }
//         std::swap(nums[start], nums[i]);     // 交换
//         backtrack(res, nums, start + 1);     // 递归处理后续位置
//         std::swap(nums[start], nums[i]);     // 恢复交换（回溯）
//     }
// }
void backtrack(std::vector<std::vector<int>>& res, std::vector<int>& nums, int start) {
    if (start == nums.size()) {
        res.push_back(nums);
        return;
    }
    std::unordered_set<int> seen;  // 记录当前层已使用的值
    for (int i = start; i < nums.size(); i++) {
        if (seen.count(nums[i])) continue;  // 跳过重复值
        seen.insert(nums[i]);
        std::swap(nums[start], nums[i]);
        backtrack(res, nums, start + 1);
        std::swap(nums[start], nums[i]);  // 回溯
    }
}

std::vector<std::vector<int>> permuteUnique(std::vector<int>& nums) {
    std::vector<std::vector<int>> res;
    std::sort(nums.begin(), nums.end());  // 先排序，让相同数字相邻
    backtrack(res, nums, 0);
    return res;
}

// Eigen::MatrixXd ComputeTransformationMatrix(
//     const std::vector<basism>& allbasisc ,const std::vector<basis>& allbasis  // 核子对类型
// )
// {
//     int n = allbasisc.size();
//     int m = allbasis.size();
//     Eigen::MatrixXd cjmat= Eigen::MatrixXd::Zero(n, m);
//     for (size_t i = 0; i < allbasisc.size(); ++i)
//     {
//         for (size_t j = 0; j < allbasis.size(); ++j)
//         {
//             double cjv=0;
//             std::vector<int> sub_vec;
//
//             std::copy(allbasisc[i].r.begin() + 1, allbasisc[i].r.end(), std::back_inserter(sub_vec));
//             if (allbasis[j].r == sub_vec && allbasisc[i].r[0] == allbasis[j].j)
//             {
//                 int sjnum=0;
//                 basism mb=allbasisc[i];
//                 basis jb=allbasis[j];
//                 int bigmvalue=mb.bigmvec[0];
//                 int sj0=jb.j;
//                 int m0=mb.bigmvec[0];
//                 double sumall=1;
//                 for (size_t k=0;k<mb.mvec.size(); ++k)
//                 {
//                     int sjfirst=sj0;
//                     if (k!=0 && !jb.sj.empty())
//                     {
//                         sjnum=sjnum+mb.rvecnum[k-1];
//                     }
//                     if (sjnum!=0)
//                     {
//                         sjfirst=jb.sj[sjnum-1];
//                     }
//                     vector<vector<int>> mvecall=permuteUnique(mb.mvec[k]);
//                     double sumcj=0;
//                     for (const auto& mvec : mvecall)
//                     {
//                         int m00=m0;
//                         int sjnum1=sjnum;
//                         double cj=1;
//                         for (size_t l = 0; l < mvec.size(); ++l)
//                         {
//                             int m01=m00+mvec[l];
//                             int sj1=jb.sj[sjnum1];
//                             if (l==0)
//                             {
//                                 cj=util::CG(sjfirst,rvecall[k],sj1,m00,mvec[l],m01);
//                             }else
//                             {
//                                 cj=cj*util::CG(jb.sj[sjnum1-1],rvecall[k],sj1,m00,mvec[l],m01);
//                             }
//                             sjnum1++;
//                             m00=m01;
//                         }
//                         sumcj+=cj;
//                     }
//                     m0=mb.m_rvec[k];
//                     sumall=sumall*sumcj;
//                 }
//                 cjmat(i,j)= sumall;
//
//             }
//         }
//     }
//     return cjmat;
// }

Eigen::MatrixXd ComputeTransformationMatrix(
    const std::vector<basism>& allbasisc ,const std::vector<basis>& allbasis  // 核子对类型
)
{
    int n = allbasisc.size();
    int m = allbasis.size();
    Eigen::MatrixXd cjmat= Eigen::MatrixXd::Zero(n, m);
    for (size_t i = 0; i < allbasisc.size(); ++i)
    {
        for (size_t j = 0; j < allbasis.size(); ++j)
        {
            double cjv=0;
            std::vector<int> sub_vec;

            std::copy(allbasisc[i].r.begin() + 1, allbasisc[i].r.end(), std::back_inserter(sub_vec));
            if (allbasis[j].r == sub_vec && allbasisc[i].r[0] == allbasis[j].j)
            {
                int sjnum=0;
                basism mb=allbasisc[i];
                basis jb=allbasis[j];
                int bigmvalue=mb.bigmvec[0];
                int sj0=jb.j;
                int m0=0;
                double sumall=1;
                for (size_t k=0;k<mb.mvec.size(); ++k)
                {
                    m0+=mb.m_rvec[k];
                    int sjfirst=sj0;
                    if (k!=0 && !jb.sj.empty())
                    {
                        sjnum=sjnum+mb.rvecnum[k-1];
                    }
                    if (sjnum!=0)
                    {
                        sjfirst=jb.sj[sjnum-1];
                    }
                    vector<vector<int>> mvecall=permuteUnique(mb.mvec[k]);
                    double sumcj=0;
                    for (const auto& mvec : mvecall)
                    {
                        int m00=m0;
                        int sjnum1=sjnum;
                        double cj=1;
                        for (size_t l = 0; l < mvec.size(); ++l)
                        {
                            int m01=m00+mvec[l];
                            int sj1=jb.sj[sjnum1];
                            if (l==0)
                            {
                                cj=util::CG(sjfirst,rvecall[k],sj1,m00,mvec[l],m01);
                            }else
                            {
                                cj=cj*util::CG(jb.sj[sjnum1-1],rvecall[k],sj1,m00,mvec[l],m01);
                            }
                            sjnum1++;
                            m00=m01;
                        }
                        sumcj+=cj;
                    }

                    sumall=sumall*sumcj;
                }
                cjmat(i,j)= sumall;

            }
        }
    }
    return cjmat;
}


Eigen::MatrixXd ComputeTransformationMatrixeven(
    const std::vector<basism>& allbasisc ,const std::vector<basis>& allbasis  // 核子对类型
)
{
    int n = allbasisc.size();
    int m = allbasis.size();
    Eigen::MatrixXd cjmat= Eigen::MatrixXd::Zero(n, m);
    for (size_t i = 0; i < allbasisc.size(); ++i)
    {
        for (size_t j = 0; j < allbasis.size(); ++j)
        {
            double cjv=0;
            std::vector<int> sub_vec;

            std::copy(allbasisc[i].rn.begin() + 1, allbasisc[i].rn.end(), std::back_inserter(sub_vec));
            if (allbasis[j].rn == sub_vec && allbasisc[i].r[0] == allbasis[j].j)
            {
                int sjnum=0;
                basism mb=allbasisc[i];
                basis jb=allbasis[j];
                int bigmvalue=mb.bigmvec[0];
                int sj0=jb.j;
                int m0=0;
                double sumall=1;
                for (size_t k=0;k<mb.mvec.size(); ++k)
                {
                    m0+=mb.m_rvec[k];
                    int sjfirst=sj0;
                    if (k!=0 && !jb.sj.empty())
                    {
                        sjnum=sjnum+mb.rvecnum[k-1];
                    }
                    if (sjnum!=0)
                    {
                        sjfirst=jb.sj[sjnum-1];
                    }
                    vector<vector<int>> mvecall=permuteUnique(mb.mvec[k]);
                    double sumcj=0;
                    for (const auto& mvec : mvecall)
                    {
                        int m00=m0;
                        int sjnum1=sjnum;
                        double cj=1;
                        for (size_t l = 0; l < mvec.size(); ++l)
                        {
                            int m01=m00+mvec[l];
                            int sj1=jb.sj[sjnum1];
                            if (l==0)
                            {
                                cj=util::CG(sjfirst,rvecall[k],sj1,m00,mvec[l],m01);
                            }else
                            {
                                cj=cj*util::CG(jb.sj[sjnum1-1],rvecall[k],sj1,m00,mvec[l],m01);
                            }
                            sjnum1++;
                            m00=m01;
                        }
                        sumcj+=cj;
                    }

                    sumall=sumall*sumcj;
                }
                cjmat(i,j)= sumall;

            }
        }
    }
    return cjmat;
}


void removeZeroRows(
    const Eigen::MatrixXd& m1,
    std::vector<basism>& mtest,
    std::vector<std::vector<Eigen::MatrixXd>>& ystrm1,
    double tolerance = 1e-10)
{
    // 步骤1: 查找全零行的行号
    std::vector<int> zeroRows;
    for (int i = 0; i < m1.rows(); ++i) {
        if (m1.row(i).isZero(tolerance)) {
            zeroRows.push_back(i);
        }
    }

    // 如果没有全零行，直接返回
    if (zeroRows.empty()) return;

    // 步骤2: 按降序排序以便安全删除
    std::sort(zeroRows.begin(), zeroRows.end(), std::greater<int>());

    // 步骤3: 删除对应元素
    for (int row : zeroRows) {
        // 确保索引在有效范围内
        if (row < mtest.size()) {
            mtest.erase(mtest.begin() + row);
            std::cout<<"remove row "<<row<<std::endl;
        }else
        {
            std::cout<<"error in removeZeroRows"<<std::endl;
        }
        if (row < ystrm1.size()) {
            ystrm1.erase(ystrm1.begin() + row);
        }else
        {
            std::cout<<"error in removeZeroRows"<<std::endl;
        }
    }
}

Eigen::MatrixXd convertToMatrixFast(const std::vector<std::vector<double>>& data) {
    if (data.empty()) return Eigen::MatrixXd();

    size_t rows = data.size();
    size_t cols = data[0].size();

    // 检查所有行的列数是否一致
    for (const auto& row : data) {
        if (row.size() != cols) {
            throw std::invalid_argument("All rows must have the same number of columns.");
        }
    }

    // 将数据复制到连续内存
    std::vector<double> flat_data;
    for (const auto& row : data) {
        flat_data.insert(flat_data.end(), row.begin(), row.end());
    }

    // 使用 Eigen::Map 直接映射内存（避免拷贝）
    return Eigen::Map<Eigen::MatrixXd>(flat_data.data(), cols, rows).transpose();
}

std::vector<int> generateSequence(int n) {
    std::vector<int> seq(n ); // 创建大小为 n+1 的 vector
    std::iota(seq.begin(), seq.end(), 0); // 从 0 开始填充
    return seq;
}

void printMatrix(const std::vector<std::vector<double>>& matrix)
{
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            std::cout << matrix[i][j] << "  ";
        }
        std::cout << std::endl;
    }
}

std::vector<std::vector<Eigen::MatrixXd>> computeystrm(const std::vector<basism>& allbasisc,
    const std::vector<std::vector<double> >& ystrgetp)
{
    int ii=0;
    nucleusm.clear();
    for (auto & nucleu : nucleus)
    {
        int jv=nucleu.j;
        int n=nucleu.n;
        int l=nucleu.l;
        for (int m=-nucleu.j;m<=nucleu.j;m+=2)
        {
            nucleusm.emplace_back(n,l,jv,m,ii);
        }
        ii++;
    }
    basis basisr;
    int n = rvecall.size();
    std::vector<int> vec(n );
    std::iota(vec.begin(), vec.end(), 0);
    basisr.rn=vec;
    std::vector<basis> basis1;
    int size1 = nucleusm.size();

    std::vector<std::vector<std::vector<std::vector<double>>>> ystrallp(
        1,  // 第一维：大小为 sizep
        std::vector<std::vector<std::vector<double>>>(
            n,  // 第二维：大小为 sizepr
            std::vector<std::vector<double>>(
                nucleus.size(),  // 第三维：大小为 size1
                std::vector<double>(nucleus.size(), 0.0)  // 第四维：大小为 size1，元素初始化为 0.0
            )
        )
    );
    basis1.push_back(basisr);

    calculateystr(ystrallp,ystrgetp,basis1);

    std::vector<std::vector<std::vector<double>>> ystr=ystrallp[0];
    // for (std::vector<std::vector<double>>yymat:ystr)
    // {
    //     std::cout<<"ymatget"<<std::endl;
    //     printMatrix(yymat);
    // }
    std::vector<Eigen::MatrixXd> ystrm;
    for (size_t m=0;m<ystr.size();m++)
    {
        ystrm.emplace_back(convertToMatrixFast(ystr[m])); // 将每个二维向量转换为 Eigen::MatrixXd
    }
    std::vector<std::vector<Eigen::MatrixXd>> ystrbasismresult;
    for (size_t i=0; i<allbasisc.size(); ++i)
    {
        std::vector<Eigen::MatrixXd> ystrbasism= std::vector<Eigen::MatrixXd>();
        int r0=allbasisc[i].r[0];
        if (r0==0)
        {
            Eigen::MatrixXd I_dynamic = Eigen::MatrixXd::Identity(size1, size1);
            Eigen::MatrixXd imat;
            ystrbasism.push_back(imat);
        }else
        {
            Eigen::VectorXd column_vector(size1);
            column_vector.setZero();
            for (int j = 0; j < size1; ++j)
            {
                int num=nucleusm[j].num;
                int j0=nucleusm[j].j;
                int m0=nucleusm[j].m;
                if (j0==allbasisc[i].r[0] && m0==allbasisc[i].bigmvec[0])
                {
                    column_vector(j)=1;
                }
            }
            ystrbasism.emplace_back(column_vector);

        }
        for (size_t j=1; j<allbasisc[i].r.size(); ++j)
        {
            Eigen::MatrixXd ystrm1(size1,size1);
            int j3=allbasisc[i].r[j];
            int m3=allbasisc[i].bigmvec[j]- allbasisc[i].bigmvec[j-1];
            int jn=allbasisc[i].rn[j];
            for (int k=0; k<size1; ++k)
            {
                int num1=nucleusm[k].num;
                int j1=nucleusm[k].j;
                int m1=nucleusm[k].m;
                for (int l=0; l<size1; ++l)
                {
                    int num2=nucleusm[l].num;
                    int j2=nucleusm[l].j;
                    int m2=nucleusm[l].m;
                    double cg_value = util::CG(j1, j2, j3, m1, m2, m3);
                    double yy=ystr[jn][num1][num2];
                    double matv=cg_value*yy;
                    ystrm1(k,l)=matv;

                }
            }


            ystrbasism.push_back(ystrm1);
        }
        ystrbasismresult.push_back(ystrbasism);
    }
    return ystrbasismresult;
}


std::vector<std::vector<std::vector<int>>> computeystrmvec(const std::vector<basism>& allbasisc,
    const std::vector<std::vector<double> >& ystrgetp,std::vector<std::vector<Eigen::MatrixXd>>& ystrmvec)
{
    ystrmvec.clear();
    int ii=0;
    nucleusm.clear();
    for (auto & nucleu : nucleus)
    {
        int jv=nucleu.j;
        int n=nucleu.n;
        int l=nucleu.l;
        for (int m=-nucleu.j;m<=nucleu.j;m+=2)
        {
            nucleusm.emplace_back(n,l,jv,m,ii);
        }
        ii++;
    }
    basis basisr;
    int n = rvecall.size();
    std::vector<int> vec(n );
    std::iota(vec.begin(), vec.end(), 0);
    basisr.rn=vec;
    std::vector<basis> basis1;
    int size1 = nucleusm.size();

    std::vector<std::vector<std::vector<std::vector<double>>>> ystrallp(
        1,  // 第一维：大小为 sizep
        std::vector<std::vector<std::vector<double>>>(
            n,  // 第二维：大小为 sizepr
            std::vector<std::vector<double>>(
                nucleus.size(),  // 第三维：大小为 size1
                std::vector<double>(nucleus.size(), 0.0)  // 第四维：大小为 size1，元素初始化为 0.0
            )
        )
    );
    basis1.push_back(basisr);

    calculateystr(ystrallp,ystrgetp,basis1);

    std::vector<std::vector<std::vector<double>>> ystr=ystrallp[0];
    std::vector<Eigen::MatrixXd> ystrm;
    for (size_t m=0;m<ystr.size();m++)
    {
        ystrm.emplace_back(convertToMatrixFast(ystr[m])); // 将每个二维向量转换为 Eigen::MatrixXd
    }
    for (int i=0;i < ystr.size();i++)
    {
        int rnn=rvecall[i];
        std::vector<Eigen::MatrixXd> rnnystrvec= {};
        for (int mi=0; mi<(2*rnn); mi=mi+2)
        {
            Eigen::MatrixXd ystrm1(size1,size1);
            int j3=rnn;
            int m3=mi-rnn;
            int jn=i;
            for (int k=0; k<size1; ++k)
            {
                int num1=nucleusm[k].num;
                int j1=nucleusm[k].j;
                int m1=nucleusm[k].m;
                for (int l=0; l<size1; ++l)
                {
                    int num2=nucleusm[l].num;
                    int j2=nucleusm[l].j;
                    int m2=nucleusm[l].m;
                    double cg_value = util::CG(j1, j2, j3, m1, m2, m3);
                    double yy=ystr[jn][num1][num2];
                    double matv=cg_value*yy;
                    ystrm1(k,l)=matv;

                }
            }
            rnnystrvec.emplace_back(ystrm1);

        }
        ystrmvec.emplace_back(rnnystrvec);
    }

    std::vector<std::vector<std::vector<int>>> ystrbasismresult;
    for (size_t i=0; i<allbasisc.size(); ++i)
    {
        std::vector<std::vector<int>> ystrbasism={};
        int r0=allbasisc[i].r[0];
        if (r0==0)
        {
            int iii=-1;
            ystrbasism.push_back({iii});
        }else
        {
            std::vector<int> vecc= {allbasisc[i].r[0],allbasisc[i].bigmvec[0]};
            ystrbasism.emplace_back(vecc);
        }
        for (size_t j=1; j<allbasisc[i].r.size(); ++j)
        {
            std::vector<int> ystrm1={};
            int j3=allbasisc[i].r[j];
            int m3=allbasisc[i].bigmvec[j]- allbasisc[i].bigmvec[j-1];
            int jn=allbasisc[i].rn[j];
            ystrm1={jn,m3};
            ystrbasism.push_back(ystrm1);
        }
        ystrbasismresult.push_back(ystrbasism);
    }
    return ystrbasismresult;
}

// 定义一个 4D 矩阵：维度为 (dim1 × dim2 × dim3 × dim4)
using Matrix4D = std::vector<std::vector<Eigen::MatrixXd>>;
using Matrix4D_sp = std::vector<std::vector<Eigen::SparseMatrix<double>>>;


Matrix4D_sp createZeroMatrix4D_sp(int dim1, int dim2, int rows, int cols) {
    Matrix4D_sp result;
    result.reserve(dim1); // 预分配第一维
    for (int i = 0; i < dim1; ++i) {
        result.emplace_back(); // 添加第二维向量
        result[i].reserve(dim2); // 预分配第二维
        for (int j = 0; j < dim2; ++j) {
            // 创建全零稀疏矩阵
            result[i].emplace_back(rows, cols);
        }
    }
    return result;
}

Eigen::Tensor<double, 4> matrix4d_to_tensor(const Matrix4D& mat4d) {
    // 获取维度大小
    const int dim0 = mat4d.size();          // 第一维大小
    const int dim1 = dim0 > 0 ? mat4d[0].size() : 0;  // 第二维大小

    // 检查内部矩阵是否一致
    if (dim0 == 0 || dim1 == 0) {
        return Eigen::Tensor<double, 4>(0,0,0,0);
    }

    const int dim2 = mat4d[0][0].rows();    // 第三维大小（矩阵行数）
    const int dim3 = mat4d[0][0].cols();    // 第四维大小（矩阵列数）

    // 创建四维张量
    Eigen::Tensor<double, 4> tensor(dim0, dim1, dim2, dim3);

    // 复制数据
    for (int i = 0; i < dim0; ++i) {
        for (int j = 0; j < dim1; ++j) {
            // 检查每个矩阵维度是否一致
            //if (mat4d[i][j].rows() != dim2 || mat4d[i][j].cols() != dim3) {
            //throw std::runtime_error("Inconsistent matrix dimensions in Matrix4D");
            //}

            for (int k = 0; k < dim2; ++k) {
                for (int l = 0; l < dim3; ++l) {
                    tensor(i, j, k, l) = mat4d[i][j](k, l);
                }
            }
        }
    }
    return tensor;
}


Matrix4D tensor_to_matrix4d(const Eigen::Tensor<double, 4>& tensor_expr) {
    // 强制表达式求值为实际张量
    Eigen::Tensor<double, 4> tensor = tensor_expr.eval();

    const int dim0 = tensor.dimension(0);
    const int dim1 = tensor.dimension(1);
    const int dim2 = tensor.dimension(2);
    const int dim3 = tensor.dimension(3);

    Matrix4D result(
        dim0,
        std::vector<Eigen::MatrixXd>(
            dim1,
            Eigen::MatrixXd::Zero(dim2, dim3)
        )
    );

    // 手动计算步长（兼容所有 Eigen 版本）
    const Eigen::Index stride0 = dim1 * dim2 * dim3;  // 第一维度步长
    const Eigen::Index stride1 = dim2 * dim3;         // 第二维度步长

    // 获取数据指针
    const double* base_ptr = tensor.data();

    for (int i = 0; i < dim0; ++i) {
        for (int j = 0; j < dim1; ++j) {
            // 计算内存偏移量
            const Eigen::Index offset =
                i * stride0 +
                j * stride1;

            // 创建内存映射
            Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
                matrix_map(base_ptr + offset, dim2, dim3);

            // 直接赋值
            result[i][j] = matrix_map;
        }
    }

    return result;
}


void printNonZeroElements(const Matrix4D& tensor) {
    const std::vector<std::string> labels = {"a", "b", "c", "d"};

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Non-zero elements report:\n"<<std::endl;
    std::cout << "-------------------------\n"<<std::endl;

    for (size_t a = 0; a < tensor.size(); ++a) {
        for (size_t b = 0; b < tensor[a].size(); ++b) {
            const auto& matrix = tensor[a][b];
            for (int c = 0; c < matrix.rows(); ++c) {
                for (int d = 0; d < matrix.cols(); ++d) {
                    double val = matrix(c, d);
                    if (val != 0.0) {
                        std::cout <<a<<',' << b<< ','<<c<<','<<d<<','<<val<<"\n"<<std::endl;
                    }
                }
            }
        }
    }
}
// 假设 galltest 是已定义和初始化的四维数组
void checkInequality(const std::vector<std::vector<Matrix4D>>& galltest) {
    // 获取 galltest 的尺寸
    int size_l = galltest.size();

    // 遍历所有索引 l, m, a, b, c, d
    for (int l = 0; l < size_l; ++l) {
        for (int m = 0; m < galltest[l].size(); ++m) {
            for (int a = 0; a < galltest[l][m].size(); ++a) {
                for (int b = 0; b < galltest[l][m][a].size(); ++b) {
                    for (int c = 0; c < galltest[l][m][a].size(); ++c) {
                        for (int d = 0; d < galltest[l][m][a].size(); ++d) {
                            // 获取 galltest[l][m][a][b](c, d)
                            double element1 = galltest[l][m][a][b](c, d);
                            // 获取 galltest[m][l][c][d](a, b)
                            double element2 = galltest[m][l][c][d](a, b);

                            // 判断两个元素是否不等
                            if (std::abs(element1-element2)> 1e-8) { // 使用一个小的阈值来比较浮点数
                                std::cout << "Elements are not equal!" << std::endl;
                                std::cout << "galltest[" << l << "][" << m << "][" << a << "][" << b << "]("
                                          << c << "," << d << ")"<< element1<<" != "
                                          << "galltest[" << m << "][" << l << "][" << c << "][" << d << "]("
                                          << a << "," << b << ")"<< element2 << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }
}

// 假设 galltest 是已定义和初始化的四维数组
void checkInequalityt(const Matrix4D& fall) {
    // 获取 galltest 的尺寸
    int size_l = fall.size();

    // 遍历所有索引 l, m, a, b, c, d
    for (int l = 0; l < size_l; ++l) {
        for (int m = 0; m < fall[l].size(); ++m) {
            for (int a = 0; a < fall[l][m].cols(); ++a) {
                for (int b = 0; b < fall[l][m].cols(); ++b) {
                    // 获取 galltest[l][m][a][b](c, d)
                    double element1 = fall[l][m](a, b);
                    // 获取 galltest[m][l][c][d](a, b)
                    double element2 = fall[m][l](b, a);

                    // 判断两个元素是否不等
                    if (std::abs(element1-element2)> 1e-8) { // 使用一个小的阈值来比较浮点数
                        std::cout << "Elements are not equal!" << std::endl;
                        std::cout << "galltest[" << l << "][" << m << "][" << a << "][" << b << "]"<< element1<<" != "
                                  << "galltest[" << m << "][" << l << "][" << b << "][" << a << "]"<< element2 << std::endl;


                    }
                }
            }
        }
    }
}
Matrix4D multiplyScalar(const Matrix4D& mat4d, double scalar) {
    Matrix4D result = mat4d;  // 先复制原数据
    for (auto& row : result) {
        for (auto& matrix : row) {
            matrix *= scalar;  // Eigen 的 MatrixXd 支持 *= 运算符
        }
    }
    return result;
}

using SparseMatrix = Eigen::SparseMatrix<double>;
using SparseMatrix4D = std::vector<std::vector<SparseMatrix>>;

SparseMatrix4D convertToSparse(const Matrix4D& denseMat) {
    SparseMatrix4D sparseMat;
    sparseMat.resize(denseMat.size());

    for (size_t a = 0; a < denseMat.size(); ++a) {
        sparseMat[a].resize(denseMat[a].size());
        for (size_t b = 0; b < denseMat[a].size(); ++b) {
            // 仅存储非零元素
            sparseMat[a][b] = denseMat[a][b].sparseView();
        }
    }
    return sparseMat;
}

Matrix4D convertToDense(const SparseMatrix4D& sparseMat) {
    Matrix4D denseMat;
    denseMat.resize(sparseMat.size());
    for (size_t a = 0; a < sparseMat.size(); ++a) {
        denseMat[a].resize(sparseMat[a].size());
        for (size_t b = 0; b < sparseMat[a].size(); ++b) {
            // 核心转换：稀疏矩阵转稠密矩阵
            denseMat[a][b] = Eigen::MatrixXd(sparseMat[a][b]);
        }
    }
    return denseMat;
}

// Matrix4D addMatrix4D(const Matrix4D& mat1, const Matrix4D& mat2) {
//     // 首先检查维度是否匹配
//     // if (mat1.size() != mat2.size()) {
//     //     throw std::invalid_argument("Outer dimension mismatch!");
//     // }
//     // for (size_t i = 0; i < mat1.size(); ++i) {
//     //     if (mat1[i].size() != mat2[i].size()) {
//     //         throw std::invalid_argument("Inner dimension mismatch!");
//     //     }
//     //     for (size_t j = 0; j < mat1[i].size(); ++j) {
//     //         if (mat1[i][j].rows() != mat2[i][j].rows() ||
//     //             mat1[i][j].cols() != mat2[i][j].cols()) {
//     //             throw std::invalid_argument("Matrix size mismatch at (" + std::to_string(i) + "," + std::to_string(j) + ")");
//     //             }
//     //     }
//     // }
//
//     // 执行加法
//     Matrix4D result = mat1;  // 拷贝 mat1 的结构
//     for (size_t i = 0; i < mat1.size(); ++i) {
//         for (size_t j = 0; j < mat1[i].size(); ++j) {
//             result[i][j] = mat1[i][j] + mat2[i][j];  // Eigen 矩阵加法
//         }
//     }
//     return result;
// }


Eigen::MatrixXd sparseAdd(const Eigen::MatrixXd& m1, const Eigen::MatrixXd& m2) {
    Eigen::MatrixXd result = m1; // 复制m1

    // 只更新m2中的非零元素
    for (int j = 0; j < m2.cols(); ++j) {
        for (int i = 0; i < m2.rows(); ++i) {
            if (m2(i, j) != 0.0) {
                result(i, j) += m2(i, j);
            }
        }
    }
    return result;
}

SparseMatrix4D addMatrix4D(const Matrix4D& mat1, const Matrix4D& mat2) {
    // 首先检查维度是否匹配
    // if (mat1.size() != mat2.size()) {
    //     throw std::invalid_argument("Outer dimension mismatch!");
    // }
    // for (size_t i = 0; i < mat1.size(); ++i) {
    //     if (mat1[i].size() != mat2[i].size()) {
    //         throw std::invalid_argument("Inner dimension mismatch!");
    //     }
    //     for (size_t j = 0; j < mat1[i].size(); ++j) {
    //         if (mat1[i][j].rows() != mat2[i][j].rows() ||
    //             mat1[i][j].cols() != mat2[i][j].cols()) {
    //             throw std::invalid_argument("Matrix size mismatch at (" + std::to_string(i) + "," + std::to_string(j) + ")");
    //             }
    //     }
    // }

    // 执行加法
    SparseMatrix4D result = convertToSparse(mat1);  // 拷贝 mat1 的结构
    SparseMatrix4D sparseMat2=convertToSparse(mat2);
    for (size_t i = 0; i < mat1.size(); ++i) {
        for (size_t j = 0; j < mat1[i].size(); ++j)
        {
            Eigen::SparseMatrix<double> sparseSum = result[i][j] + sparseMat2[i][j];
            result[i][j] = sparseSum;
        }
    }
    return result;
}

SparseMatrix4D addMatrix4D(const SparseMatrix4D& mat1, const Matrix4D& mat2) {
    SparseMatrix4D result = mat1;  // 拷贝 mat1 的结构
    SparseMatrix4D sparseMat2=convertToSparse(mat2);
    for (size_t i = 0; i < mat1.size(); ++i) {
        for (size_t j = 0; j < mat1[i].size(); ++j)
        {
            Eigen::SparseMatrix<double> sparseSum = result[i][j] + sparseMat2[i][j];
            result[i][j] = sparseSum;
        }
    }
    return result;
}

SparseMatrix4D addMatrix4D(const SparseMatrix4D& mat1, const SparseMatrix4D& mat2) {
    SparseMatrix4D result = mat1;  // 拷贝 mat1 的结构
    SparseMatrix4D sparseMat2=mat2;
    for (size_t i = 0; i < mat1.size(); ++i) {
        for (size_t j = 0; j < mat1[i].size(); ++j)
        {
            Eigen::SparseMatrix<double> sparseSum = result[i][j] + sparseMat2[i][j];
            result[i][j] = sparseSum;
        }
    }
    return result;
}

Eigen::MatrixXd currentq(Eigen::MatrixXd p1p,Eigen::MatrixXd p1,Eigen::MatrixXd p2)
{
    Eigen::MatrixXd resultp;
    Eigen::MatrixXd t1=p1*p1p*p2;
    Eigen::MatrixXd t2=p2*p1p*p1;
    Eigen::MatrixXd t31=p1*p1p;
    double t311=t31.trace()*0.5;
    Eigen::MatrixXd t3=t311*p2;
    resultp=t1+t2-t3;
    return resultp;
}

Matrix4D currentq1(Eigen::MatrixXd p1p,Eigen::MatrixXd p1,Eigen::MatrixXd p2,Eigen::MatrixXd p3)
{
    Eigen::MatrixXd re1=currentq(p1p,p1,p2);
    int nucnum=nucleusm.size();
    Matrix4D resultq(nucnum,
    std::vector<Eigen::MatrixXd>(nucnum,
    Eigen::MatrixXd::Zero(nucnum, nucnum))
    );
    for (int i=0;i<nucnum;++i)
    {
        for (int j=0;j<nucnum;++j)
        {
            double t1=re1(i,j);
            if (t1==0)
            {
                continue;
            }
            for (int k=0;k<nucnum;++k)
            {
                for (int l=0;l<nucnum;++l)
                {
                    double t2=p3(k,l);
                    double t3=t1*t2;
                    resultq[i][j](k,l)=4*t3;

                }
            }
        }
    }
    return resultq;
}
Eigen::MatrixXd currentq2_2(Eigen::MatrixXd p1p,Eigen::MatrixXd p1,Eigen::MatrixXd p3)
{
    Eigen::MatrixXd resultp;
    Eigen::MatrixXd t1=p1*p1p*p3;
    Eigen::MatrixXd t2=p3*p1p*p1;
    Eigen::MatrixXd t31=p3*p1p;
    double t311=t31.trace()*0.5;
    Eigen::MatrixXd t3=t311*p1;
    resultp=t1+t2-t3;
    return resultp;
}

Matrix4D currentq2(Eigen::MatrixXd p1p,Eigen::MatrixXd p1,Eigen::MatrixXd p2,Eigen::MatrixXd p3)
{
    Eigen::MatrixXd re1=currentq2_2(p1p,p1,p3);
    int nucnum=nucleusm.size();
    Matrix4D resultq(nucnum,
    std::vector<Eigen::MatrixXd>(nucnum,
    Eigen::MatrixXd::Zero(nucnum, nucnum))
    );
    for (int i=0;i<nucnum;++i)
    {
        for (int j=0;j<nucnum;++j)
        {
            double t1=re1(i,j);
            if (t1==0)
            {
                continue;
            }
            for (int k=0;k<nucnum;++k)
            {
                for (int l=0;l<nucnum;++l)
                {
                    double t2=p2(k,l);
                    double t3=t1*t2;
                    resultq[i][j](k,l)=4*t3;

                }
            }
        }
    }
    return resultq;
}

Eigen::MatrixXd currentq3_3(Eigen::MatrixXd p1p,Eigen::MatrixXd p2,Eigen::MatrixXd p3)
{
    Eigen::MatrixXd resultp;
    Eigen::MatrixXd t1=p2*p1p*p3;
    Eigen::MatrixXd t2=p3*p1p*p2;
    Eigen::MatrixXd t31=p2*p1p;
    double t311=t31.trace()*0.5;
    Eigen::MatrixXd t3=t311*p3;
    resultp=t1+t2-t3;
    return resultp;
}
Matrix4D currentq3(Eigen::MatrixXd p1p,Eigen::MatrixXd p1,Eigen::MatrixXd p2,Eigen::MatrixXd p3)
{
    Eigen::MatrixXd re1=currentq3_3(p1p,p2,p3);
    int nucnum=nucleusm.size();
    Matrix4D resultq(nucnum,
    std::vector<Eigen::MatrixXd>(nucnum,
    Eigen::MatrixXd::Zero(nucnum, nucnum))
    );
    for (int i=0;i<nucnum;++i)
    {
        for (int j=0;j<nucnum;++j)
        {
            double t1=re1(i,j);
            if (t1==0)
            {
                continue;
            }
            for (int k=0;k<nucnum;++k)
            {
                for (int l=0;l<nucnum;++l)
                {
                    double t2=p1(k,l);
                    double t3=t1*t2;
                    resultq[k][l](i,j)=4*t3;

                }
            }
        }
    }
    return resultq;
}


Eigen::MatrixXd calsigp(Eigen::MatrixXd pk,Eigen::MatrixXd pn,Eigen::MatrixXd pnp)
{
    Eigen::MatrixXd t1=pk*pnp;
    double t2=-t1.trace()*2;
    // if (t2==0)
    // {
    //     std::cout<<"t2 is zero"<<std::endl;
    // }
    Eigen::MatrixXd t3=t2*pn;
    return t3;
}

Eigen::MatrixXd cal2p(Eigen::MatrixXd pk,Eigen::MatrixXd pnp,Eigen::MatrixXd pi)
{
    Eigen::MatrixXd t1=pk*pnp*pi;
    Eigen::MatrixXd t2=pi*pnp*pk;
    Eigen::MatrixXd t3=t1+t2;
    Eigen::MatrixXd t4=4*t3;
    return t4;
}



SparseMatrix4D calq(const basism& bas1,const basism& bas2,const std::vector<Eigen::MatrixXd>&ystrm1,const std::vector<Eigen::MatrixXd>&ystrm2)
{
    // auto start = std::chrono::high_resolution_clock::now();
    int nucnum=nucleusm.size();
    int nnp=ystrm2.size();
    int nn=ystrm1.size();
    Matrix4D resultmat(nucnum,
    std::vector<Eigen::MatrixXd>(nucnum,
    Eigen::MatrixXd::Zero(nucnum, nucnum))  // 最内层是 5x5 矩阵
    );
    int rsize=ystrm1.size();
    // 检查输入是否有效
    if (ystrm1.size() - ystrm2.size()!=2) {
        throw std::invalid_argument("ystrm1 and ystrm2 must have the same size!");
    }
    if (bas1.r[0]==0)
    {
        if (bas1.r.size()==4)
        {
            Eigen::MatrixXd p1p=ystrm2[1];
            Eigen::MatrixXd p1=ystrm1[1];
            Eigen::MatrixXd p2=ystrm1[2];
            Eigen::MatrixXd p3=ystrm1[3];

            SparseMatrix4D resultq1=addMatrix4D(currentq1(p1p,p1,p2,p3),currentq2(p1p,p1,p2,p3));
            SparseMatrix4D resultq2=addMatrix4D(resultq1,currentq3(p1p,p1,p2,p3));
            return resultq2;
        }else if (bas1.r.size()==3)
        {
            Eigen::MatrixXd p1=ystrm1[1];
            Eigen::MatrixXd p2=ystrm1[2];
            for (int i=0;i<nucnum;++i)
            {
                for (int j=0;j<nucnum;++j)
                {
                    double t1=p1(i,j);
                    for (int k=0;k<nucnum;++k)
                    {
                        for (int l=0;l<nucnum;++l)
                        {
                            double t2=p2(k,l);
                            double t3=t1*t2;
                            resultmat[i][j](k,l)=t3;
                        }
                    }
                }
            }
            SparseMatrix4D resultmatsp=convertToSparse(resultmat);
            return resultmatsp;
        }else
        {
            Matrix4D re1(nucnum,
            std::vector<Eigen::MatrixXd>(nucnum,Eigen::MatrixXd::Zero(nucnum, nucnum))  // 最内层是 5x5 矩阵
            );
            SparseMatrix4D term1=convertToSparse(re1);
            SparseMatrix4D term2=term1;
            for (int k=1;k<rsize;++k)
            {
                basism bas11=bas1;
                basism bas22=bas2;
                std::vector<Eigen::MatrixXd>ystrm11=ystrm1;
                std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
                Eigen::MatrixXd pk=ystrm1[k];
                ystrm11.erase(ystrm11.begin() + k);
                Eigen::MatrixXd pn=ystrm11.back();
                Eigen::MatrixXd pnp=ystrm2[rsize-3];
                Eigen::MatrixXd pkre=calsigp(pk,pn,pnp);
                // ystrm11[nn-1]=pkre;
                ystrm22.pop_back();
                bas11.r.pop_back();
                bas22.r.pop_back();
                ystrm11.back()=pkre;
                SparseMatrix4D tt1=calq(bas11,bas22,ystrm11,ystrm22);
                term1=addMatrix4D(term1,tt1);
                // std::cout<<"term1"<<std::endl;
                // std::cout<<"k:"<<k<<std::endl;
                // printNonZeroElements(term1);
            }
            // std::cout<<"term1"<<std::endl;
            // printNonZeroElements(term1);
            for (int k=2;k<rsize;++k)
            {
                for (int i=1;i<=k-1;++i)
                {
                    basism bas11=bas1;
                    basism bas22=bas2;
                    std::vector<Eigen::MatrixXd>ystrm11=ystrm1;
                    std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
                    Eigen::MatrixXd pk=ystrm1[k];
                    Eigen::MatrixXd pnp=ystrm2[rsize-3];
                    Eigen::MatrixXd pi=ystrm1[i];
                    Eigen::MatrixXd pkre=cal2p(pk,pnp,pi);
                    ystrm11.erase(ystrm11.begin() + k);
                    ystrm11.erase(ystrm11.begin() + i);

                    ystrm11.push_back(pkre);
                    ystrm22.pop_back();
                    bas11.r.pop_back();
                    bas22.r.pop_back();
                    SparseMatrix4D tt1=calq(bas11,bas22,ystrm11,ystrm22);
                    term2=addMatrix4D(term2,tt1);
                    // std::cout<<"term2"<<std::endl;
                    // std::cout<<"k:"<<k<<"i"<<i<<std::endl;
                    // printNonZeroElements(term2);

                }
            }
            // std::cout<<"term2"<<std::endl;
            // printNonZeroElements(term2);
            term1= addMatrix4D(term1, term2);
            // std::cout<<"resultmat"<<std::endl;
            // printNonZeroElements(resultmat);
            return term1;
        }
    }else
    {
        // if (bas1.r.size()==3)
        // {
        //     Eigen::VectorXd s=ystrm1[0];
        //     Eigen::VectorXd sp=ystrm2[0];
        //     double t1=sp.transpose()*s;
        //     Eigen::MatrixXd p1= ystrm1[1];
        //     Eigen::MatrixXd p2= ystrm1[2];
        //     Eigen::MatrixXd t2_1= s*sp.transpose()*p1;
        //     Eigen::MatrixXd t2_2= p1*sp*s.transpose();
        //     Eigen::MatrixXd t2=t2_1+t2_2;
        //     Eigen::MatrixXd t3_1= s*sp.transpose()*p2;
        //     Eigen::MatrixXd t3_2= p2*sp*s.transpose();
        //     Eigen::MatrixXd t3=t3_1+t3_2;
        //     for (int a=0;a<nucnum;++a)
        //     {
        //         for (int b=0;b<nucnum;++b)
        //         {
        //             for (int c=0;c<nucnum;++c)
        //             {
        //                 for (int d=0;d<nucnum;++d)
        //                 {
        //                     double term1=t1*p1(a,b)*p2(c,d);
        //                     double term2=t2(a,b)*p2(c,d);
        //                     double term3=t3(a,b)*p1(c,d);
        //                     resultmat[a][b](c,d)=term1-term2-term3;
        //                 }
        //             }
        //         }
        //     }
        //     return resultmat;
        // }
        if (bas1.r.size() == 3)
        {
            Eigen::VectorXd s = ystrm1[0];
            Eigen::VectorXd sp = ystrm2[0];
            const double t1 = sp.dot(s);  // 更高效的向量点积

            Eigen::MatrixXd p1 = ystrm1[1];
            Eigen::MatrixXd p2 = ystrm1[2];

            // 优化矩阵乘法
            Eigen::MatrixXd t2 = s * sp.transpose() * p1 + p1 * sp * s.transpose();
            Eigen::MatrixXd t3 = s * sp.transpose() * p2 + p2 * sp * s.transpose();


            // 使用矩阵运算替代循环
            #pragma omp parallel for collapse(2)
            for (int a = 0; a < nucnum; ++a) {
                for (int b = 0; b < nucnum; ++b) {
                    // 计算标量因子
                    const double scalar1 = t1 * p1(a, b) - t2(a, b);
                    const double scalar2 = t3(a, b);

                    // 直接计算整个矩阵
                    resultmat[a][b] = scalar1 * p2 - scalar2 * p1;
                }
            }
            // auto end = std::chrono::high_resolution_clock::now();
            // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            // std::cout << "calgtime: " <<"\n"<< duration.count() << "μs\n"<<std::endl;
            SparseMatrix4D term1=convertToSparse(resultmat);


            return term1;
        }
        else
        {
            Matrix4D re1(nucnum,
            std::vector<Eigen::MatrixXd>(nucnum,Eigen::MatrixXd::Zero(nucnum, nucnum))  // 最内层是 5x5 矩阵
            );
            SparseMatrix4D term1=convertToSparse(re1);
            SparseMatrix4D term2=term1;
            SparseMatrix4D term3=term2;
            for (int k=1;k<rsize;++k)
            {
                basism bas11=bas1;
                basism bas22=bas2;
                std::vector<Eigen::MatrixXd>ystrm11=ystrm1;
                std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
                Eigen::MatrixXd pk=ystrm1[k];
                ystrm11.erase(ystrm11.begin() + k);
                Eigen::MatrixXd pn=ystrm11.back();
                Eigen::MatrixXd pnp=ystrm2[rsize-3];
                Eigen::MatrixXd pkre=calsigp(pk,pn,pnp);
                ystrm11.back()=pkre;
                ystrm22.pop_back();
                bas11.r.pop_back();
                bas22.r.pop_back();
                SparseMatrix4D tt1=calq(bas11,bas22,ystrm11,ystrm22);
                term1=addMatrix4D(term1,tt1);
            }
            for (int k=2;k<rsize;++k)
            {
                for (int i=1;i<=k-1;++i)
                {
                    basism bas11=bas1;
                    basism bas22=bas2;
                    std::vector<Eigen::MatrixXd>ystrm11=ystrm1;
                    std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
                    Eigen::MatrixXd pk=ystrm1[k];
                    Eigen::MatrixXd pnp=ystrm2[rsize-3];
                    Eigen::MatrixXd pi=ystrm1[i];
                    Eigen::MatrixXd pkre=cal2p(pk,pnp,pi);
                    ystrm11.erase(ystrm11.begin() + k);
                    ystrm11.erase(ystrm11.begin() + i);

                    ystrm11.push_back(pkre);
                    ystrm22.pop_back();
                    bas11.r.pop_back();
                    bas22.r.pop_back();
                    SparseMatrix4D tt1=calq(bas11,bas22,ystrm11,ystrm22);
                    term2=addMatrix4D(term2,tt1);

                }
            }
            for (int k=1;k<rsize;++k)
            {
                basism bas11=bas1;
                basism bas22=bas2;
                std::vector<Eigen::MatrixXd>ystrm11=ystrm1;
                std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
                Eigen::MatrixXd pk=ystrm1[k];
                ystrm11.erase(ystrm11.begin() + k);
                Eigen::MatrixXd pnp=ystrm2[rsize-3];
                Eigen::VectorXd s= ystrm1[0];
                Eigen::MatrixXd pkre=4*(pk*pnp*s);
                ystrm11[0]=pkre;
                ystrm22.pop_back();
                bas11.r.pop_back();
                bas22.r.pop_back();
                SparseMatrix4D tt1=calq(bas11,bas22,ystrm11,ystrm22);
                term1=addMatrix4D(term1,tt1);
            }
            term1= addMatrix4D(term1, term2);
            return term1;

        }
    }

}

// Matrix4D gchange(Matrix4D gor)
// {
//
//     double six=1.0/6.0;
//     int nucnum=nucleusm.size();
//     Matrix4D result(nucnum,
//            std::vector<Eigen::MatrixXd>(nucnum,Eigen::MatrixXd::Zero(nucnum, nucnum))  // 最内层是 5x5 矩阵
//            );
//     for (int a=0;a<nucnum;++a)
//     {
//         for (int b=0;b<nucnum;++b)
//         {
//             for (int c=0;c<nucnum;++c)
//             {
//                 for (int d=0;d<nucnum;++d)
//                 {
//                     result[a][b](c,d)=six*(gor[a][b](c,d)-gor[a][c](b,d)+gor[a][d](b,c)
//                         +gor[b][c](a,d)-gor[b][d](a,c)+gor[c][d](a,b));
//                 }
//             }
//         }
//     }
//     return result;
// }

// Matrix4D gchange(const Matrix4D& gor)//_optimized
// {
//     const double six = 1.0 / 6.0;
//     const int nucnum = nucleusm.size();
//
//     // 预分配结果张量
//     Matrix4D result(nucnum,
//            std::vector<Eigen::MatrixXd>(nucnum, Eigen::MatrixXd::Zero(nucnum, nucnum)));
//
//     // 并行化外层循环
//     #pragma omp parallel for collapse(2) schedule(dynamic)
//     for (int a = 0; a < nucnum; ++a) {
//         for (int b = 0; b < nucnum; ++b) {
//             // 获取当前结果矩阵的引用
//             Eigen::MatrixXd& res_ab = result[a][b];
//
//             // 第一项：gor[a][b](c,d)
//             res_ab = gor[a][b];
//
//             // 添加其他项
//             for (int c = 0; c < nucnum; ++c) {
//                 for (int d = 0; d < nucnum; ++d) {
//                     // 第二项：-gor[a][c](b,d)
//                     res_ab(c, d) -= gor[a][c](b, d);
//
//                     // 第三项：+gor[a][d](b,c)
//                     res_ab(c, d) += gor[a][d](b, c);
//
//                     // 第四项：+gor[b][c](a,d)
//                     res_ab(c, d) += gor[b][c](a, d);
//
//                     // 第五项：-gor[b][d](a,c)
//                     res_ab(c, d) -= gor[b][d](a, c);
//
//                     // 第六项：+gor[c][d](a,b)
//                     res_ab(c, d) += gor[c][d](a, b);
//                 }
//             }
//
//             // 乘以系数1/6
//             res_ab *= six;
//         }
//     }
//
//     return result;
// }

Matrix4D gchange(const Matrix4D& gor)//_optimized
{
    const double six = 1.0 / 6.0;
    const int nucnum = nucleusm.size();

    // 预分配结果张量
    Matrix4D result(nucnum,
           std::vector<Eigen::MatrixXd>(nucnum, Eigen::MatrixXd::Zero(nucnum, nucnum)));

    // 并行化外层循环
// #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int a = 0; a < nucnum; ++a) {
        for (int b = a; b < nucnum; ++b) {
         // 添加其他项
            for (int c = b; c < nucnum; ++c) {
                for (int d = c; d < nucnum; ++d) {
                    double resum=gor[a][b](c,d);
                    // 第二项：-gor[a][c](b,d)

                    resum -= gor[a][c](b, d);

                    // 第三项：+gor[a][d](b,c)
                    resum += gor[a][d](b, c);

                    // 第四项：+gor[b][c](a,d)
                    resum += gor[b][c](a, d);

                    // 第五项：-gor[b][d](a,c)
                    resum -= gor[b][d](a, c);

                    // 第六项：+gor[c][d](a,b)
                    resum += gor[c][d](a, b);
                    resum*= six;
                    // 反对称化处理
                    if (resum==0)
                    {
                        continue;
                    }
                    result[a][b](c,d)=resum;
                    result[a][b](d,c)=-resum;
                    result[a][c](b,d)=-resum;
                    result[a][c](d,b)=resum;
                    result[a][d](b,c)=resum;
                    result[a][d](c,b)=-resum;
                    result[b][a](c,d)=-resum;
                    result[b][a](d,c)=resum;
                    result[b][c](a,d)=resum;
                    result[b][c](d,a)=-resum;
                    result[b][d](a,c)=-resum;
                    result[b][d](c,a)=resum;
                    result[c][a](b,d)=resum;
                    result[c][a](d,b)=-resum;
                    result[c][b](a,d)=-resum;
                    result[c][b](d,a)=resum;
                    result[c][d](a,b)=resum;
                    result[c][d](b,a)=-resum;
                    result[d][a](b,c)=-resum;
                    result[d][a](c,b)=resum;
                    result[d][b](a,c)=resum;
                    result[d][b](c,a)=-resum;
                    result[d][c](a,b)=-resum;
                    result[d][c](b,a)=resum;

                }
            }
        }
    }

    return result;
}

using Tensor4D = Eigen::Tensor<double, 4>;

Tensor4D gchange_tensor(const Matrix4D& gor)
{
    const double six = 1.0 / 6.0;
    const int nucnum = nucleusm.size();

    // 创建并初始化四维张量 [nucnum, nucnum, nucnum, nucnum]
    Tensor4D result(nucnum, nucnum, nucnum, nucnum);
    result.setZero();

    // 只处理四个下标互不相同的情况
    for (int a = 0; a < nucnum; ++a) {
        for (int b = a ; b < nucnum; ++b) {
            for (int c = b ; c < nucnum; ++c) {
                for (int d = c ; d < nucnum; ++d) {
                    double resum = gor[a][b](c, d);
                    resum -= gor[a][c](b, d);  // -gor[a][c](b,d)
                    resum += gor[a][d](b, c);   // +gor[a][d](b,c)
                    resum += gor[b][c](a, d);   // +gor[b][c](a,d)
                    resum -= gor[b][d](a, c);   // -gor[b][d](a,c)
                    resum += gor[c][d](a, b);   // +gor[c][d](a,b)
                    resum *= six;

                    if (resum == 0) continue;

                    // 直接对张量的24个排列位置赋值
                    result(a, b, c, d) = resum;
                    result(a, b, d, c) = -resum;
                    result(a, c, b, d) = -resum;
                    result(a, c, d, b) = resum;
                    result(a, d, b, c) = resum;
                    result(a, d, c, b) = -resum;

                    result(b, a, c, d) = -resum;
                    result(b, a, d, c) = resum;
                    result(b, c, a, d) = resum;
                    result(b, c, d, a) = -resum;
                    result(b, d, a, c) = -resum;
                    result(b, d, c, a) = resum;

                    result(c, a, b, d) = resum;
                    result(c, a, d, b) = -resum;
                    result(c, b, a, d) = -resum;
                    result(c, b, d, a) = resum;
                    result(c, d, a, b) = resum;
                    result(c, d, b, a) = -resum;

                    result(d, a, b, c) = -resum;
                    result(d, a, c, b) = resum;
                    result(d, b, a, c) = resum;
                    result(d, b, c, a) = -resum;
                    result(d, c, a, b) = -resum;
                    result(d, c, b, a) = resum;
                }
            }
        }
    }

    return result;
}

// Matrix4D gchange(const Matrix4D& gor)//_optimized
// {
//     constexpr double six = 1.0 / 6.0;
//     const int nucnum = nucleusm.size();
//     Matrix4D result(nucnum, std::vector<Eigen::MatrixXd>(nucnum,
//                  Eigen::MatrixXd::Zero(nucnum, nucnum)));
//
//     // 预定义索引映射表 (避免重复计算排列)
//     constexpr int perm[24][4] = {
//         {0,1,2,3}, {0,1,3,2}, {0,2,1,3}, {0,2,3,1}, {0,3,1,2}, {0,3,2,1},
//         {1,0,2,3}, {1,0,3,2}, {1,2,0,3}, {1,2,3,0}, {1,3,0,2}, {1,3,2,0},
//         {2,0,1,3}, {2,0,3,1}, {2,1,0,3}, {2,1,3,0}, {2,3,0,1}, {2,3,1,0},
//         {3,0,1,2}, {3,0,2,1}, {3,1,0,2}, {3,1,2,0}, {3,2,0,1}, {3,2,1,0}
//     };
//     constexpr int sign[24] = {
//         1, -1, -1, 1, 1, -1,
//         -1, 1, 1, -1, -1, 1,
//         1, -1, -1, 1, 1, -1,
//         -1, 1, 1, -1, -1, 1
//     };
//
//     for (int a = 0; a < nucnum; ++a) {
//         auto& res_a = result[a];
//         auto& gor_a = gor[a];
//
//         for (int b = a; b < nucnum; ++b) {
//             auto& res_ab = res_a[b];
//             auto& gor_ab = gor_a[b];
//
//             for (int c = b; c < nucnum; ++c) {
//                 auto& res_ac = res_a[c];
//                 auto& gor_ac = gor_a[c];
//
//                 for (int d = c; d < nucnum; ++d) {
//                     double resum = gor_ab(c, d);
//                     resum -= gor_ac(b, d);
//                     resum += gor[a][d](b, c);
//                     resum += gor[b][c](a, d);
//                     resum -= gor[b][d](a, c);
//                     resum += gor[c][d](a, b);
//                     resum *= six;
//                     if (resum==0) {
//                         continue;  // 如果结果为0，跳过
//                     }
//
//                     int idx[4] = {a, b, c, d};
//
//                     // 修复：使用 const int* 而不是 int*
//                     for (int p = 0; p < 24; ++p) {
//                         const int* pm = perm[p];  // 关键修复：添加 const
//                         result[idx[pm[0]]][idx[pm[1]]](idx[pm[2]], idx[pm[3]]) = sign[p] * resum;
//                     }
//                 }
//             }
//         }
//     }
//
//     return result;
// }

// // 定义排列组合类型
// using Permutation = std::array<int, 5>; // [i, j, k, l, sign]
// using IndexMap = std::array<Permutation, 6>;
//
// // 全局索引映射表智能指针
// static std::shared_ptr<IndexMap[]> global_index_map;
// static int cached_nucnum = 0;
// static std::mutex init_mutex;
//
// // 初始化全局索引映射表
// void initialize_global_index_map(int nucnum) {
//     std::lock_guard<std::mutex> lock(init_mutex);
//     if (nucnum == cached_nucnum && global_index_map)
//         return;
//
//     const size_t total_elements = static_cast<size_t>(nucnum) * nucnum * nucnum * nucnum;
//     global_index_map = std::make_shared<IndexMap[]>(total_elements);
//     cached_nucnum = nucnum;
//
//     // 预定义排列组合
//     const std::array<Permutation, 6> base_permutations = {{
//         {{0, 1, 2, 3, +1}},  // q_{αβγδ}
//         {{0, 2, 1, 3, -1}},  // -q_{αγβδ}
//         {{0, 3, 1, 2, +1}},  // +q_{αδβγ}
//         {{1, 2, 0, 3, +1}},  // +q_{βγαδ}
//         {{1, 3, 0, 2, -1}},  // -q_{βδαγ}
//         {{2, 3, 0, 1, +1}}   // +q_{γδαβ}
//     }};
//
//     // 填充索引映射表
//     for (int a = 0; a < nucnum; ++a) {
//         for (int b = 0; b < nucnum; ++b) {
//             for (int c = 0; c < nucnum; ++c) {
//                 for (int d = 0; d < nucnum; ++d) {
//                     const size_t idx = static_cast<size_t>(a) * nucnum*nucnum*nucnum +
//                                       static_cast<size_t>(b) * nucnum*nucnum +
//                                       static_cast<size_t>(c) * nucnum + d;
//
//                     for (int p = 0; p < 6; ++p) {
//                         const auto& base = base_permutations[p];
//                         auto& dest = global_index_map[idx][p];
//
//                         // 应用排列
//                         dest[0] = (base[0] == 0) ? a : (base[0] == 1) ? b : (base[0] == 2) ? c : d;
//                         dest[1] = (base[1] == 0) ? a : (base[1] == 1) ? b : (base[1] == 2) ? c : d;
//                         dest[2] = (base[2] == 0) ? a : (base[2] == 1) ? b : (base[2] == 2) ? c : d;
//                         dest[3] = (base[3] == 0) ? a : (base[3] == 1) ? b : (base[3] == 2) ? c : d;
//                         dest[4] = base[4];  // 符号因子
//                     }
//                 }
//             }
//         }
//     }
// }
//
// // 针对 nucnum=32 的专用优化版本
// Matrix4D gchange(const Matrix4D& gor)//_optimized
// {
//     constexpr int nucnum = 32;  // 固定为32
//     constexpr double six = 1.0 / 6.0;
//     constexpr size_t total = static_cast<size_t>(nucnum) * nucnum * nucnum * nucnum;
//
//     // 初始化全局索引映射表（如果尚未初始化）
//     if (!global_index_map || cached_nucnum != nucnum) {
//         initialize_global_index_map(nucnum);
//     }
//
//     // 预分配结果张量
//     Matrix4D result(nucnum,
//            std::vector<Eigen::MatrixXd>(nucnum, Eigen::MatrixXd::Zero(nucnum, nucnum)));
//
//     // 创建一维视图函数
//     auto get_element = [nucnum](const Matrix4D& tensor, int a, int b, int c, int d) {
//         return tensor[a][b](c, d);
//     };
//
//     auto set_element = [nucnum](Matrix4D& tensor, int a, int b, int c, int d, double value) {
//         tensor[a][b](c, d) = value;
//     };
//
//     // 主计算循环
//     for (int a = 0; a < nucnum; ++a) {
//         for (int b = 0; b < nucnum; ++b) {
//             for (int c = 0; c < nucnum; ++c) {
//                 for (int d = 0; d < nucnum; ++d) {
//                     const size_t idx = static_cast<size_t>(a) * nucnum*nucnum*nucnum +
//                                       static_cast<size_t>(b) * nucnum*nucnum +
//                                       static_cast<size_t>(c) * nucnum + d;
//
//                     double sum = 0.0;
//                     for (int p = 0; p < 6; ++p) {
//                         const auto& perm = global_index_map[idx][p];
//                         sum += perm[4] * get_element(gor, perm[0], perm[1], perm[2], perm[3]);
//                     }
//
//                     set_element(result, a, b, c, d, sum * six);
//                 }
//             }
//         }
//     }
//
//     return result;
// }

std::vector<Eigen::MatrixXd> tensor3DToMatrixVector(
    const Eigen::Tensor<double, 3>& tensor)
{
    // 检查张量维度
    static_assert(Eigen::Tensor<double, 3>::NumDimensions == 3,
                 "Input must be a 3D tensor");

    // 获取各维度大小
    const int dim0 = tensor.dimension(0);
    const int dim1 = tensor.dimension(1);
    const int dim2 = tensor.dimension(2);

    // 预分配结果向量
    std::vector<Eigen::MatrixXd> result(dim0,Eigen::MatrixXd::Zero(dim1, dim2));


    // 逐层转换
    for (int i = 0; i < dim0; ++i) {
        for (int j = 0; j < dim1; ++j)
        {
            for (int k = 0; k < dim2; ++k)
            {
                result[i](j,k)= tensor(i, j, k);
            }
        }
    }

    return result;
}
Eigen::Tensor<double, 3> calt(const basism& bas1,const basism& bas2,const std::vector<Eigen::MatrixXd>&ystrm1,const std::vector<Eigen::MatrixXd>&ystrm2)
{
    // auto start = std::chrono::high_resolution_clock::now();

    int nucnum=nucleusm.size();
    int rsize=ystrm1.size();
    Eigen::Tensor<double, 3> resultmat(nucnum,nucnum,nucnum);
    resultmat.setZero();
    if (bas1.r.size()==3)
    {
        Eigen::VectorXd s=ystrm1[0];
        Eigen::MatrixXd p1=ystrm1[1];
        Eigen::MatrixXd p2=ystrm1[2];
        Eigen::MatrixXd p1p=ystrm2[1];

        Eigen::MatrixXd t1_1=p1*p1p*p2;
        Eigen::MatrixXd t1_2_1=p1*p1p;
        Eigen::MatrixXd t1_2=0.5*t1_2_1.trace()*p2;
        Eigen::MatrixXd t1=t1_1-t1_2;

        Eigen::MatrixXd t2_1=p2*p1p*p1;
        Eigen::MatrixXd t2_2_1=p2*p1p;
        Eigen::MatrixXd t2_2=0.5*t2_2_1.trace()*p1;
        Eigen::MatrixXd t2=t2_1-t2_2;
        Eigen::MatrixXd t12=t1+t2;

        Eigen::MatrixXd t3_1=p1*p1p*s;

        Eigen::MatrixXd t4_1=p2*p1p*s;

        for (int a=0;a<nucnum;++a)
        {
            for (int b=0;b<nucnum;++b)
            {
                for (int c=0;c<nucnum;++c)
                {
                    double result;
                    double term12=4*s[a]*t12(b,c);
                    double term3=4*t3_1(a)*p2(b,c);
                    double term4=4*t4_1(a)*p1(b,c);
                    result=term12+term3+term4;
                    resultmat(a,b,c)=result;
                }
            }
        }
        // auto end = std::chrono::high_resolution_clock::now();
        // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        // std::cout << "calttime: " <<"\n"<< duration.count() << "μs\n"<<std::endl;
        return resultmat;

    }else if (bas1.r.size()==2)
    {
        Eigen::VectorXd s=ystrm1[0];
        Eigen::MatrixXd p1=ystrm1[1];
        for (int a=0;a<nucnum;++a)
        {
            for (int b=0;b<nucnum;++b)
            {
                for (int c=0;c<nucnum;++c)
                {
                    resultmat(a,b,c)=s[a]*p1(b,c);
                }
            }
        }
        return resultmat;
    }else
    {
        int nn=ystrm1.size();
        Eigen::Tensor<double, 3> term1(nucnum,nucnum,nucnum);
        term1.setZero(); // 显式置零
        Eigen::Tensor<double, 3> term2(nucnum,nucnum,nucnum);
        term2.setZero(); // 显式置零
        Eigen::Tensor<double, 3> term3(nucnum,nucnum,nucnum);
        term3.setZero(); // 显式置零
        for (int k=1;k<rsize;++k)
        {
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<< "calt,term1, k: " << k << "\n"<<std::endl;
            basism bas11=bas1;
            basism bas22=bas2;
            std::vector<Eigen::MatrixXd>ystrm11=ystrm1;
            std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
            Eigen::MatrixXd pk=ystrm1[k];
            ystrm11.erase(ystrm11.begin() + k);
            Eigen::MatrixXd pn=ystrm11.back();
            Eigen::MatrixXd pnp=ystrm2[rsize-2];
            Eigen::MatrixXd pkre=calsigp(pk,pn,pnp);
            ystrm11.back()=pkre;
            ystrm22.pop_back();
            bas11.r.pop_back();
            bas22.r.pop_back();
            Eigen::Tensor<double, 3> tt1=calt(bas11,bas22,ystrm11,ystrm22);

            term1= term1 + tt1; // 使用 += 操作符进行加法
            // writeTensor3ToFile(tt1);
        }
        for (int k=2;k<rsize;++k)
        {
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<< "calt,term2, k: " << k << "\n"<<std::endl;
            for (int i=1;i<=k-1;++i)
            {
                // if (!outfile.is_open())
                // {
                //     // 如果文件未打开，则尝试打开文件
                //     outfile.open("basis_output.txt",std::ios::app);
                // }
                // outfile<< "calt,term2, i: " << i << "\n"<<std::endl;
                basism bas11=bas1;
                basism bas22=bas2;
                std::vector<Eigen::MatrixXd>ystrm11=ystrm1;
                std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
                Eigen::MatrixXd pk=ystrm1[k];
                Eigen::MatrixXd pnp=ystrm2[rsize-2];
                Eigen::MatrixXd pi=ystrm1[i];
                Eigen::MatrixXd pkre=cal2p(pk,pnp,pi);
                ystrm11.erase(ystrm11.begin() + k);
                ystrm11.erase(ystrm11.begin() + i);

                ystrm11.push_back(pkre);
                ystrm22.pop_back();
                bas11.r.pop_back();
                bas22.r.pop_back();
                Eigen::Tensor<double, 3> tt1=calt(bas11,bas22,ystrm11,ystrm22);
                term2= term2 + tt1; // 使用 += 操作符进行加法
                // writeTensor3ToFile(tt1);

            }
        }
        for (int k=1;k<rsize;++k)
        {
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<< "calt,term3, k: " << k << "\n"<<std::endl;
            basism bas11=bas1;
            basism bas22=bas2;
            std::vector<Eigen::MatrixXd>ystrm11=ystrm1;
            std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
            Eigen::MatrixXd pk=ystrm1[k];
            ystrm11.erase(ystrm11.begin() + k);
            Eigen::MatrixXd pnp=ystrm2[rsize-2];
            Eigen::VectorXd s= ystrm1[0];
            Eigen::MatrixXd pkre=4*(pk*pnp*s);
            ystrm11[0]=pkre;
            ystrm22.pop_back();
            bas11.r.pop_back();
            bas22.r.pop_back();
            Eigen::Tensor<double, 3> tt1=calt(bas11,bas22,ystrm11,ystrm22);
            term3= term3 + tt1; // 使用 += 操作符进行加法
            // writeTensor3ToFile(tt1);
        }
        resultmat= term1+term2+term3;
        return resultmat;

    }
}

Eigen::Tensor<double,3> tchange(Eigen::Tensor<double,3> tor)
{
    double six=1.0/3.0;
    int nucnum=nucleusm.size();
    Eigen::Tensor<double, 3> result(nucnum,nucnum,nucnum);
    result.setZero();
    result.setZero(); // 显式置零
    double resum=0;
    for (int a=0;a<nucnum;++a)
    {
        for (int b=a;b<nucnum;++b)
        {
            for (int c=b;c<nucnum;++c)
            {
                double t1=tor(a,b,c);
                double t2=tor(b,a,c);
                double t3=tor(c,a,b);
                resum=six*(t1-t2+t3);
                if (resum==0)
                {
                    continue;
                }
                result(a,b,c)=resum;
                result(a,c,b)=-resum;
                result(b,c,a)=resum;
                result(c,b,a)=-resum;
                result(c,a,b)=resum;
                result(b,a,c)=-resum;


            }
        }
    }
    return result;
}

Eigen::MatrixXd tensor_to_matrix(const Eigen::Tensor<double, 2>& tensor)
{
    const int rows = tensor.dimension(0);
    const int cols = tensor.dimension(1);

    // 创建内存映射（零拷贝）
    return Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
        tensor.data(),  // 原始数据指针
        rows,           // 行数
        cols            // 列数
    );
}




Eigen::MatrixXd calb(const basism& bas1,const basism& bas2,const std::vector<Eigen::MatrixXd>&ystrm1,const std::vector<Eigen::MatrixXd>&ystrm2)
{
    // auto start = std::chrono::high_resolution_clock::now();

    int nucnum=nucleusm.size();
    int nnp=ystrm2.size();
    int nn=ystrm1.size();
    Eigen::MatrixXd resultmat= Eigen::MatrixXd::Zero(nucnum, nucnum);
    int rsize=ystrm1.size();
    // 检查输入是否有效
    if (ystrm1.size() - ystrm2.size()!=1) {
        throw std::invalid_argument("ystrm1 and ystrm2 must have the same size!");
    }
    if (bas1.r[0]==0)
    {
        if (bas1.r.size()==3)
        {
            Eigen::MatrixXd p1p=ystrm2[1];
            Eigen::MatrixXd p1=ystrm1[1];
            Eigen::MatrixXd p2=ystrm1[2];
            Eigen::MatrixXd t1=p1*p1p*p2;
            Eigen::MatrixXd t2=p2*p1p*p1;
            Eigen::MatrixXd t3=4*(t1+t2);
            Eigen::MatrixXd t11=p1*p1p;
            double t311=t11.trace()*2;
            Eigen::MatrixXd t22=t311*p2;
            Eigen::MatrixXd t4=p2*p1p;
            double t44=t4.trace()*2;
            Eigen::MatrixXd t55=t44*p1;
            Eigen::MatrixXd t66=t55+t22;
            Eigen::MatrixXd resultq1=t3-t66;
            return resultq1;
        }else if (bas1.r.size()==2)
        {
            resultmat=ystrm1[1];
            return resultmat;
        }else
        {
            // auto start = std::chrono::high_resolution_clock::now();
            basism bas22=bas2;
            Eigen::MatrixXd pn_1=ystrm2[rsize-2];
            std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
            ystrm22.pop_back();
            bas22.r.pop_back();
            Matrix4D qca=convertToDense(  calq (bas1,bas22,ystrm1,ystrm22));
            // std::cout<<"qca"<<std::endl;
            // printNonZeroElements(qca);
            Matrix4D qbar= gchange(qca);
            // auto end = std::chrono::high_resolution_clock::now();
            // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            // std::cout <<"time1:"<< duration.count() << std::endl;
            // Matrix4D qbar= qca;
            // std::cout<<"qbar"<<std::endl;
            // printNonZeroElements(qbar);
            // start = std::chrono::high_resolution_clock::now();
            for (int a=0;a<nucnum;++a)
            {
                for (int b=0;b<nucnum;++b)
                {
                    if (qbar[a][b].isZero(1e-12)) continue;
                    double matsum = qbar[a][b].cwiseProduct(pn_1).sum();
                    resultmat(a,b)=12*matsum;
                }
            }
            // end = std::chrono::high_resolution_clock::now();
            // duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            // std::cout <<"time2:"<< duration.count() << std::endl;
            return resultmat;

        }
    }else
    {
        if (bas1.r.size()==2)
        {
            Eigen::VectorXd s=ystrm1[0];
            Eigen::VectorXd sp=ystrm2[0];
            Eigen::MatrixXd p1= ystrm1[1];
            double t1=sp.transpose()*s;
            Eigen::MatrixXd t2=s*sp.transpose();
            Eigen::MatrixXd t3=sp*s.transpose();
            resultmat=(t1*p1)-(t2*p1)-(p1*t3);
            return resultmat;

        }else
        {
            Eigen::VectorXd sp=ystrm2[0];
            std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
            ystrm22[0]=Eigen::MatrixXd::Identity(nucnum,nucnum);
            basism bas22=bas2;
            bas22.r[0]=0;
            Eigen::Tensor<double, 3> t1=calt(bas1,bas22,ystrm1,ystrm22);
            Eigen::Tensor<double, 3> t2=tchange(t1);
            Eigen::Tensor<double, 2> tensor2dsum(nucnum, nucnum);
            tensor2dsum.setZero();

            // writeTensor3ToFile(t2);
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<<"sp"<<sp<<std::endl;
            for (int c=0;c<nucnum;++c)
            {
                if (sp(c)==0)
                {
                    continue; // 如果系数为0，跳过该切片
                }
                auto tslice=t2.chip(c,2); // 获取第 c 个切片
                tensor2dsum+=3*sp(c)*tslice;
                // writeSimpleTensorToFile(tensor2dsum);
            }
            // resultmat=tensor_to_matrix(tensor2dsum);
            for (int a=0;a<nucnum;++a)
            {
                for (int b=0;b<nucnum;++b)
                {

                    resultmat(a,b)=tensor2dsum(a,b);

                }
            }
            // auto end = std::chrono::high_resolution_clock::now();
            // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            // std::cout << "canbtime " <<"\n"<< duration.count() << "μs\n"<<std::endl;
            return resultmat;


        }
    }
}

// Matrix4D m2tom4(Eigen::MatrixXd pip, Eigen::MatrixXd bmat)
// {
//     int nucnum=nucleusm.size();
//     Matrix4D result(nucnum,
//            std::vector<Eigen::MatrixXd>(nucnum,Eigen::MatrixXd::Zero(nucnum, nucnum))  // 最内层是 5x5 矩阵
//            );
//     for (int a=0;a<nucnum;++a)
//     {
//         for (int b=0;b<nucnum;++b)
//         {
//             for (int c=0;c<nucnum;++c)
//             {
//                 for (int d=0;d<nucnum;++d)
//                 {
//                     result[a][b](c,d)=pip(a,b)*bmat(c,d);
//
//                 }
//             }
//         }
//     }
//     return result;
// }

Matrix4D m2tom4(const Eigen::MatrixXd& pip, const Eigen::MatrixXd& bmat)
{
    // auto start = std::chrono::high_resolution_clock::now();

    const int nucnum = nucleusm.size();

    // 预分配结果矩阵
    Matrix4D result(nucnum,
           std::vector<Eigen::MatrixXd>(nucnum, Eigen::MatrixXd::Zero(nucnum, nucnum)));

    // 优化1：缓存友好循环顺序
    // #pragma omp parallel for collapse(2)
    int i=0;
    for (int a = 0; a < nucnum; ++a) {
        for (int b = 0; b < nucnum; ++b) {
            if (abs(pip(a, b)) < 1e-16)
            {
                i+=1;
                continue; // 跳过零元素以减少计算量

            }
            // 直接计算整个矩阵，避免内层循环
            result[a][b].noalias() = pip(a, b) * bmat;
        }
    }
    std::cout <<"i="<< i << std::endl;
    // auto end = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    // std::cout << "m2tom4cal elapsed time: " << duration.count() << "μs\n";

    return result;
}

void m2tom4sp(const Eigen::MatrixXd& pip, const Eigen::MatrixXd& bmat, Matrix4D_sp& resultmat)
{
    auto start = std::chrono::high_resolution_clock::now();

    const int nucnum = nucleusm.size();

    // 预分配结果矩阵
    Eigen::SparseMatrix<double> bmat_sparse = bmat.sparseView();

    // 优化1：缓存友好循环顺序
    // #pragma omp parallel for collapse(2)
    for (int a = 0; a < nucnum; ++a) {
        for (int b = 0; b < nucnum; ++b) {
            if (abs(pip(a, b)) < 1e-16)
            {
                continue; // 跳过零元素以减少计算量

            }
            // 直接计算整个矩阵，避免内层循环
            resultmat[a][b] += pip(a, b)* 4 * bmat_sparse;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    // std::cout << "m2tom4cal elapsed time: " << duration.count() << "μs\n";
}

// Matrix4D m2tom4(const Eigen::MatrixXd& pip, const Eigen::MatrixXd& bmat)//_sparse_optimized
// {
//     const int nucnum = pip.rows();
//     const int bm_rows = bmat.rows();
//     const int bm_cols = bmat.cols();
//
//     // 1. 将稠密矩阵pip转换为稀疏矩阵
//     Eigen::SparseMatrix<double> sparse_pip = pip.sparseView();
//     sparse_pip.makeCompressed();  // 压缩格式以提高遍历效率
//
//
//     // 2. 初始化结果矩阵（全零矩阵）
//     Matrix4D result(nucnum, std::vector<Eigen::MatrixXd>(nucnum,
//                      Eigen::MatrixXd::Zero(bm_rows, bm_cols)));
//
//     // 3. 仅遍历稀疏矩阵中的非零元素
//     for (int k = 0; k < sparse_pip.outerSize(); ++k) {
//         for (Eigen::SparseMatrix<double>::InnerIterator it(sparse_pip, k); it; ++it) {
//             const int a = it.row();
//             const int b = it.col();
//
//             // 4. 直接计算矩阵乘法（避免临时对象）
//             result[a][b].noalias() = it.value() * bmat;
//         }
//     }
//
//     return result;
// }

// Matrix4D m2tom4(const Eigen::MatrixXd& pip, const Eigen::MatrixXd& bmat_dense)//_sparse_optimized
// {
//     const int nucnum = nucleusm.size();
//
//     // 将bmat转换为稀疏矩阵
//     Eigen::SparseMatrix<double> bmat = bmat_dense.sparseView();
//     bmat.makeCompressed();  // 压缩格式以提高效率
//
//     // 预分配结果矩阵
//     Matrix4D result(nucnum,
//            std::vector<Eigen::MatrixXd>(nucnum, Eigen::MatrixXd::Zero(nucnum, nucnum)));
//
//     // 预处理：收集所有非零元素位置
//     std::vector<std::tuple<int, int, double>> nonZeroEntries;
//     for (int k = 0; k < bmat.outerSize(); ++k) {
//         for (Eigen::SparseMatrix<double>::InnerIterator it(bmat, k); it; ++it) {
//             nonZeroEntries.emplace_back(it.row(), it.col(), it.value());
//         }
//     }
//
//     #pragma omp parallel for collapse(2) schedule(dynamic)
//     for (int a = 0; a < nucnum; ++a) {
//         for (int b = 0; b < nucnum; ++b) {
//             const double scalar = pip(a, b);
//             Eigen::MatrixXd& target = result[a][b];
//
//             // 只处理非零元素
//             for (const auto& entry : nonZeroEntries) {
//                 int c = std::get<0>(entry);
//                 int d = std::get<1>(entry);
//                 target(c, d) = scalar * std::get<2>(entry);
//             }
//         }
//     }
//
//     return result;
// }


// Matrix4D m4tom4(Eigen::MatrixXd pipa, Eigen::MatrixXd pipb,Matrix4D qmat)
// {
//     auto start = std::chrono::high_resolution_clock::now();
//     int nucnum=nucleusm.size();
//     Matrix4D result(nucnum,
//            std::vector<Eigen::MatrixXd>(nucnum,Eigen::MatrixXd::Zero(nucnum, nucnum))  // 最内层是 5x5 矩阵
//            );
//     for (int a=0;a<nucnum;++a)
//     {
//         for (int b=0;b<nucnum;++b)
//         {
//             for (int c=0;c<nucnum;++c)
//             {
//                 for (int d=0;d<nucnum;++d)
//                 {
//
//                     double sumre=0;
//                     for (int e=0;e<nucnum;++e)
//                     {
//                         for (int f=0;f<nucnum;++f)
//                         {
//                             sumre= sumre +qmat[e][f](c,d)*( pipa(a,e)*pipb(f,b)+pipa(f,b)*pipb(a,e));
//                         }
//                     }
//                     result[a][b](c,d)=sumre;
//
//                 }
//             }
//         }
//     }
//     auto end = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//     std::cout << "Elapsed time: " << duration.count() << "μs\n";
//     return result;
// }

// Matrix4D m4tom4(Eigen::MatrixXd pipa, Eigen::MatrixXd pipb,Matrix4D qmat)
// {
//     auto start = std::chrono::high_resolution_clock::now();
//     int nucnum=nucleusm.size();
//     Matrix4D result(nucnum,
//            std::vector<Eigen::MatrixXd>(nucnum,Eigen::MatrixXd::Zero(nucnum, nucnum))  // 最内层是 5x5 矩阵
//            );
//     #pragma omp parallel for collapse(2)
//     for (int c=0;c<nucnum;++c)
//     {
//         for (int d=0;d<nucnum;++d)
//         {
//             Eigen::MatrixXd term1=pipa*qmat[c][d]*pipb;
//             Eigen::MatrixXd term2=pipb*qmat[c][d]*pipa;
//             Eigen::MatrixXd term3(nucnum, nucnum);
//             term3.noalias() = term1 + term2;
//             for (int a=0;a<nucnum;++a)
//             {
//                 for (int b=0;b<nucnum;++b)
//                 {
//                     result[a][b](c,d)=term3(a,b);
//                 }
//             }
//
//
//         }
//     }
//     auto end = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//     std::cout << "Elapsed time: " << duration.count() << "μs\n";
//     return result;
// }

// Matrix4D m4tom4(const Eigen::MatrixXd& pipa_dense,
//                                  const Eigen::MatrixXd& pipb_dense,
//                                  const Matrix4D& qmat)
// {
//     auto start = std::chrono::high_resolution_clock::now();
//     const int nucnum = nucleusm.size();
//
//     // 将稠密矩阵转换为稀疏矩阵格式
//     Eigen::SparseMatrix<double> pipa = pipa_dense.sparseView();
//     Eigen::SparseMatrix<double> pipb = pipb_dense.sparseView();
//
//     // 预分配结果矩阵
//     Matrix4D result(nucnum,
//            std::vector<Eigen::MatrixXd>(nucnum, Eigen::MatrixXd::Zero(nucnum, nucnum)));
//
//     // 提前计算转置矩阵（稀疏矩阵转置效率高）
//     Eigen::SparseMatrix<double> pipa_transposed = pipa.transpose();
//     Eigen::SparseMatrix<double> pipb_transposed = pipb.transpose();
//
//     // #pragma omp parallel
//     // {
//         // 每个线程创建局部变量避免竞争
//         Eigen::SparseMatrix<double> term1, term2;
//
//         // #pragma omp for collapse(2) schedule(dynamic)
//         for (int c = 0; c < nucnum; ++c)
//         {
//             for (int d = 0; d < nucnum; ++d)
//             {
//                 // 将qmat[c][d]转换为稀疏矩阵（如果它还不是）
//                 Eigen::SparseMatrix<double> qmat_cd = qmat[c][d].sparseView();
//
//                 // 优化计算顺序：先计算小矩阵乘法
//                 // 使用稀疏矩阵乘法
//                 term1 = pipa * qmat_cd;
//                 term1 = term1 * pipb;
//
//                 term2 = pipb * qmat_cd;
//                 term2 = term2 * pipa;
//
//                 // 直接赋值到结果矩阵
//                 for (int a = 0; a < nucnum; ++a)
//                 {
//                     for (int b = 0; b < nucnum; ++b)
//                     {
//                         // 稀疏矩阵访问效率较低，但赋值次数少
//                         result[a][b](c, d) = term1.coeff(a, b) + term2.coeff(a, b);
//                     }
//                 }
//             }
//         }
//     // }
//
//     auto end = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//     std::cout << "Optimized sparse elapsed time: " << duration.count() << "μs\n";
//     return result;
// }

// Matrix4D m4tom4(const Eigen::MatrixXd& pipa_dense,
//                 const Eigen::MatrixXd& pipb_dense,
//                 const Matrix4D& qmat)
// {
//     auto start = std::chrono::high_resolution_clock::now();
//
//     const int nucnum = nucleusm.size();
//
//     // 将稠密矩阵转换为稀疏矩阵格式
//     Eigen::SparseMatrix<double> pipa = pipa_dense.sparseView();
//     Eigen::SparseMatrix<double> pipb = pipb_dense.sparseView();
//
//
//     // 预分配结果矩阵
//     Matrix4D result(nucnum, std::vector<Eigen::MatrixXd>(nucnum, Eigen::MatrixXd::Zero(nucnum, nucnum)));
//
//     // 避免转置操作，直接计算并缓存qmat与pipa、pipb的乘积
//     for (int c = 0; c < nucnum; ++c) {
//         for (int d = 0; d < nucnum; ++d) {
//             // 获取qmat的子矩阵
//             Eigen::SparseMatrix<double> qmat_cd = qmat[c][d].sparseView();
//
//             // 直接计算并缓存 term1 和 term2
//             Eigen::SparseMatrix<double> term1 = pipa * qmat_cd * pipb;
//             Eigen::SparseMatrix<double> term2 = pipb * qmat_cd * pipa;
//
//             // 将结果赋值到结果矩阵
//             for (int a = 0; a < nucnum; ++a) {
//                 for (int b = 0; b < nucnum; ++b) {
//                     // 稀疏矩阵访问效率较低，但赋值次数少
//                     result[a][b](c, d) = term1.coeff(a, b) + term2.coeff(a, b);
//                 }
//             }
//         }
//     }
//
//     auto end = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//     std::cout << "Optimized sparse elapsed time: " << duration.count() << "μs\n";
//     return result;
// }

Matrix4D m4tom4(const Eigen::MatrixXd& pipa_dense,
                          const Eigen::MatrixXd& pipb_dense,
                          const Matrix4D& qmat)//_optimized
{
    auto start = std::chrono::high_resolution_clock::now();
    const int nucnum = nucleusm.size();

    // 将稠密矩阵转换为稀疏矩阵格式（压缩存储）
    Eigen::SparseMatrix<double> pipa = pipa_dense.sparseView();
    Eigen::SparseMatrix<double> pipb = pipb_dense.sparseView();
    pipa.makeCompressed();
    pipb.makeCompressed();

    // 预分配结果矩阵
    Matrix4D result(nucnum,
                   std::vector<Eigen::MatrixXd>(nucnum,
                   Eigen::MatrixXd::Zero(nucnum, nucnum)));

    // 预计算 pipa 和 pipb 的非零元素位置
    std::vector<std::vector<int>> pipa_nnz(nucnum);
    std::vector<std::vector<int>> pipb_nnz(nucnum);

    #pragma omp parallel for
    for (int i = 0; i < nucnum; ++i) {
        for (Eigen::SparseMatrix<double>::InnerIterator it_pipa(pipa, i); it_pipa; ++it_pipa) {
            pipa_nnz[i].push_back(it_pipa.row());
        }
        for (Eigen::SparseMatrix<double>::InnerIterator it_pipb(pipb, i); it_pipb; ++it_pipb) {
            pipb_nnz[i].push_back(it_pipb.row());
        }
    }

    // 并行化外层循环
    // #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int c = 0; c < nucnum; ++c) {
        for (int d = 0; d < nucnum; ++d) {
            // 获取当前 qmat 矩阵（稀疏视图）
            const Eigen::MatrixXd& qmat_cd_dense = qmat[c][d];
            Eigen::SparseMatrix<double> qmat_cd = qmat_cd_dense.sparseView();
            qmat_cd.makeCompressed();

            // 只计算必要的部分
            Eigen::SparseMatrix<double> term1 = pipa * qmat_cd * pipb;
            Eigen::SparseMatrix<double> term2 = pipb * qmat_cd * pipa;

            // 只处理非零元素
            for (int k = 0; k < term1.outerSize(); ++k) {
                for (Eigen::SparseMatrix<double>::InnerIterator it(term1, k); it; ++it) {
                    const int a = it.row();
                    const int b = it.col();
                    result[a][b](c, d) += it.value();
                }
            }

            for (int k = 0; k < term2.outerSize(); ++k) {
                for (Eigen::SparseMatrix<double>::InnerIterator it(term2, k); it; ++it) {
                    const int a = it.row();
                    const int b = it.col();
                    result[a][b](c, d) += it.value();
                }
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    // std::cout << "Optimized sparse elapsed time: " << duration.count() << "μs\n";
    return result;
}




// void m4tom4sp(const Eigen::MatrixXd& pipa_dense,
//                           const Eigen::MatrixXd& pipb_dense,
//                           const Matrix4D& qmat, Matrix4D_sp& resultmat)//_optimized
// {
//     auto start = std::chrono::high_resolution_clock::now();
//     const int nucnum = nucleusm.size();
//
//     // 将稠密矩阵转换为稀疏矩阵格式（压缩存储）
//     Eigen::SparseMatrix<double> pipa = pipa_dense.sparseView();
//     Eigen::SparseMatrix<double> pipb = pipb_dense.sparseView();
//     pipa.makeCompressed();
//     pipb.makeCompressed();
//
//     // 预分配结果矩阵
//     Matrix4D result(nucnum,
//                    std::vector<Eigen::MatrixXd>(nucnum,
//                    Eigen::MatrixXd::Zero(nucnum, nucnum)));
//
//     // 预计算 pipa 和 pipb 的非零元素位置
//     std::vector<std::vector<int>> pipa_nnz(nucnum);
//     std::vector<std::vector<int>> pipb_nnz(nucnum);
//
//     #pragma omp parallel for
//     for (int i = 0; i < nucnum; ++i) {
//         for (Eigen::SparseMatrix<double>::InnerIterator it_pipa(pipa, i); it_pipa; ++it_pipa) {
//             pipa_nnz[i].push_back(it_pipa.row());
//         }
//         for (Eigen::SparseMatrix<double>::InnerIterator it_pipb(pipb, i); it_pipb; ++it_pipb) {
//             pipb_nnz[i].push_back(it_pipb.row());
//         }
//     }
//
//     // 并行化外层循环
//     #pragma omp parallel for collapse(2) schedule(dynamic)
//     for (int c = 0; c < nucnum; ++c) {
//         for (int d = 0; d < nucnum; ++d) {
//             // 获取当前 qmat 矩阵（稀疏视图）
//             const Eigen::MatrixXd& qmat_cd_dense = qmat[c][d];
//             Eigen::SparseMatrix<double> qmat_cd = qmat_cd_dense.sparseView();
//             qmat_cd.makeCompressed();
//
//             // 只计算必要的部分
//             Eigen::SparseMatrix<double> term1 = pipa * qmat_cd * pipb*96;
//             Eigen::SparseMatrix<double> term2 = pipb * qmat_cd * pipa*96;
//
//             // 只处理非零元素
//             for (int k = 0; k < term1.outerSize(); ++k) {
//                 for (Eigen::SparseMatrix<double>::InnerIterator it(term1, k); it; ++it) {
//                     const int a = it.row();
//                     const int b = it.col();
//                     resultmat [a][b].coeffRef(c, d) += it.value();
//                 }
//             }
//
//             for (int k = 0; k < term2.outerSize(); ++k) {
//                 for (Eigen::SparseMatrix<double>::InnerIterator it(term2, k); it; ++it) {
//                     const int a = it.row();
//                     const int b = it.col();
//                     resultmat[a][b].coeffRef(c, d) += it.value();
//                 }
//             }
//         }
//     }
//
//     auto end = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//     std::cout << "Optimized sparse elapsed time: " << duration.count() << "μs\n";
// }

void m4tom4sp(const Eigen::MatrixXd& pipa_dense,
                        const Eigen::MatrixXd& pipb_dense,
                        const Matrix4D& qmat, Matrix4D_sp& resultmat)//_optimized
{
    auto start = std::chrono::high_resolution_clock::now();
    const int nucnum = nucleusm.size();

    // 1. 预转换和压缩所有稀疏矩阵
    Eigen::SparseMatrix<double> pipa = pipa_dense.sparseView();
    Eigen::SparseMatrix<double> pipb = pipb_dense.sparseView();
    pipa.makeCompressed();
    pipb.makeCompressed();

    // 预转换qmat为稀疏矩阵并压缩
    std::vector<std::vector<Eigen::SparseMatrix<double>>> qmat_sp(nucnum);
    for (int c = 0; c < nucnum; ++c) {
        qmat_sp[c].resize(nucnum);
        for (int d = 0; d < nucnum; ++d) {
            qmat_sp[c][d] = qmat[c][d].sparseView();
            qmat_sp[c][d].makeCompressed();
        }
    }

    // 2. 优化非零元素位置预计算
    std::vector<std::unordered_set<int>> pipa_nnz_cols(nucnum);
    std::vector<std::unordered_set<int>> pipb_nnz_cols(nucnum);

    #pragma omp parallel for
    for (int col = 0; col < nucnum; ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(pipa, col); it; ++it) {
            pipa_nnz_cols[col].insert(it.row());
        }
        for (Eigen::SparseMatrix<double>::InnerIterator it(pipb, col); it; ++it) {
            pipb_nnz_cols[col].insert(it.row());
        }
    }

    // 3. 使用局部累加器减少锁竞争
    Matrix4D_sp local_result(nucnum);
    for (int a = 0; a < nucnum; ++a) {
        local_result[a].resize(nucnum);
        for (int b = 0; b < nucnum; ++b) {
            local_result[a][b] = Eigen::SparseMatrix<double>(nucnum, nucnum);
        }
    }

    // 4. 并行循环优化
    // #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int c = 0; c < nucnum; ++c) {
        for (int d = 0; d < nucnum; ++d) {
            const auto& qmat_cd = qmat_sp[c][d];

            // 5. 优化稀疏矩阵乘法 - 跳过全零矩阵
            if (qmat_cd.nonZeros() == 0) continue;

            // 6. 使用三元乘积避免中间存储
            Eigen::SparseMatrix<double> term1;
            Eigen::SparseMatrix<double> term2;

            // 直接计算 (pipa * qmat_cd) * pipb * 96
            term1 = Eigen::SparseMatrix<double>((pipa * qmat_cd).pruned() * pipb);
            term1 *= 96;

            // 直接计算 (pipb * qmat_cd) * pipa * 96
            term2 = Eigen::SparseMatrix<double>((pipb * qmat_cd).pruned() * pipa);
            term2 *= 96;

            // 7. 高效累加到局部结果
            for (int k = 0; k < term1.outerSize(); ++k) {
                for (Eigen::SparseMatrix<double>::InnerIterator it(term1, k); it; ++it) {
                    const int a = it.row();
                    const int b = it.col();
                    local_result[a][b].coeffRef(c, d) += it.value();
                }
            }

            for (int k = 0; k < term2.outerSize(); ++k) {
                for (Eigen::SparseMatrix<double>::InnerIterator it(term2, k); it; ++it) {
                    const int a = it.row();
                    const int b = it.col();
                    local_result[a][b].coeffRef(c, d) += it.value();
                }
            }
        }
    }

    // 8. 合并局部结果到最终输出
    #pragma omp parallel for collapse(2)
    for (int a = 0; a < nucnum; ++a) {
        for (int b = 0; b < nucnum; ++b) {
            resultmat[a][b] += local_result[a][b];
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    // std::cout << "Optimized sparse elapsed time: " << duration.count() << "μs\n";
}

Tensor4D m4tom4_tensor(const Eigen::MatrixXd& pipa_dense,
                       const Eigen::MatrixXd& pipb_dense,
                       const Tensor4D& qmat_tensor)
{
    auto start = std::chrono::high_resolution_clock::now();
    const int nucnum = pipa_dense.rows();  // 假设为方阵

    // 将稠密矩阵转换为稀疏矩阵
    Eigen::SparseMatrix<double> pipa = pipa_dense.sparseView();
    Eigen::SparseMatrix<double> pipb = pipb_dense.sparseView();
    pipa.makeCompressed();
    pipb.makeCompressed();

    // 创建结果张量并初始化为零
    Tensor4D result(nucnum, nucnum, nucnum, nucnum);
    result.setZero();

    // 预提取qmat张量的数据指针（列主序）
    const double* qmat_data = qmat_tensor.data();

    #pragma omp parallel
    {
        // 每个线程私有一个平面矩阵，避免内存分配竞争
        Eigen::MatrixXd plane(nucnum, nucnum);

        #pragma omp for collapse(2) schedule(dynamic)
        for (int c = 0; c < nucnum; ++c) {
            for (int d = 0; d < nucnum; ++d) {
                plane.setZero();

                // 直接映射qmat的(c,d)平面到Eigen矩阵（零拷贝）
                Eigen::Map<const Eigen::MatrixXd> qmat_cd_dense(
                    qmat_data +
                    c * nucnum * nucnum * nucnum +  // c维偏移
                    d * nucnum * nucnum,            // d维偏移
                    nucnum, nucnum
                );

                // 跳过全零平面
                if (qmat_cd_dense.isZero(0.0)) continue;

                // 转换为稀疏矩阵（仅当非零时才转换）
                Eigen::SparseMatrix<double> qmat_cd = qmat_cd_dense.sparseView();
                qmat_cd.makeCompressed();

                // 计算两项矩阵乘法
                Eigen::SparseMatrix<double> term1 = pipa * qmat_cd * pipb;
                Eigen::SparseMatrix<double> term2 = pipb * qmat_cd * pipa;

                // 将结果累加到平面矩阵
                for (int k = 0; k < term1.outerSize(); ++k) {
                    for (Eigen::SparseMatrix<double>::InnerIterator it(term1, k); it; ++it) {
                        plane(it.row(), it.col()) += it.value();
                    }
                }
                for (int k = 0; k < term2.outerSize(); ++k) {
                    for (Eigen::SparseMatrix<double>::InnerIterator it(term2, k); it; ++it) {
                        plane(it.row(), it.col()) += it.value();
                    }
                }

                // 将平面矩阵复制到结果张量的(c,d)平面
                Eigen::Map<Eigen::MatrixXd> result_plane(
                    &result(0, 0, c, d),
                    nucnum, nucnum
                );
                result_plane = plane;
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Tensor optimized elapsed time: " << duration.count() << "μs\n";
    return result;
}

Eigen::Tensor<double,4> m2tom4tensor(const Eigen::MatrixXd& pip, const Eigen::MatrixXd& bmat)
{
    const int nucnum = pip.rows();
    const int rows = bmat.rows();
    const int cols = bmat.cols();

    // 创建四维张量 [nucnum, nucnum, rows, cols]
    Eigen::Tensor<double, 4> result(nucnum, nucnum, rows, cols);

    // 将 bmat 广播到四维张量的后两个维度
    Eigen::TensorMap<const Eigen::Tensor<double, 2>> bmat_tensor(bmat.data(), rows, cols);
    auto bmat_4d = bmat_tensor.reshape(Eigen::array<Eigen::Index, 4>{1, 1, rows, cols})
                      .broadcast(Eigen::array<Eigen::Index, 4>{nucnum, nucnum, 1, 1});

    // 将 pip 重塑为四维张量 [nucnum, nucnum, 1, 1] 用于广播
    Eigen::TensorMap<const Eigen::Tensor<double, 2>> pip_tensor(pip.data(), nucnum, nucnum);
    auto pip_4d = pip_tensor.reshape(Eigen::array<Eigen::Index, 4>{nucnum, nucnum, 1, 1})
                     .broadcast(Eigen::array<Eigen::Index, 4>{1, 1, rows, cols});

    // 执行逐元素乘法：result = pip(a,b) * bmat
    result = pip_4d * bmat_4d;

    return result;
}

// Eigen::Tensor<double,4> m2tom4tensor(const Eigen::MatrixXd& pip, const Eigen::MatrixXd& bmat)//_optimized
// {
//     const int nucnum = pip.rows();
//     const int rows = bmat.rows();
//     const int cols = bmat.cols();
//
//     // 创建并初始化结果张量为全零
//     Eigen::Tensor<double, 4> result(nucnum, nucnum, rows, cols);
//     result.setZero();
//
//     // 提取bmat的TensorMap（避免拷贝，仅创建视图）
//     Eigen::TensorMap<const Eigen::Tensor<double, 2>> bmat_tensor(bmat.data(), rows, cols);
//
//     // 仅遍历pip中的非零元素
//     for (int a = 0; a < nucnum; ++a) {
//         for (int b = 0; b < nucnum; ++b) {
//             const double val = pip(a, b);
//             // 跳过零元素以利用pip的稀疏性
//             if (val == 0.0) continue;
//
//             // 获取当前(a,b)位置的二维切片（rows x cols平面）
//             auto slice = result.chip(a, 0).chip(b, 1);
//
//             // 直接计算：slice = val * bmat
//             slice = bmat_tensor * val;
//         }
//     }
//
//     return result;
// }

// Matrix4D m2tom4(const Eigen::MatrixXd& pip, const Eigen::MatrixXd& bmat)
// {
//     const int nucnum = pip.rows();
//     Matrix4D result(nucnum, std::vector<Eigen::MatrixXd>(nucnum,
//                  bmat * 1.0)); // 先拷贝bmat模板
//
//     for (int a = 0; a < nucnum; ++a) {
//         for (int b = 0; b < nucnum; ++b) {
//             result[a][b] *= pip(a, b); // 就地乘法
//         }
//     }
//     return result;
// }
//
// using Matrix4D = std::vector<std::vector<Eigen::MatrixXd>>;
//
// // 将Matrix4D转换为Eigen::Tensor<double,4>
// Eigen::Tensor<double, 4> matrix4d_to_tensor(const Matrix4D& mat4d) {
//     int n = mat4d.size();
//     Eigen::Tensor<double, 4> tensor(n, n, n, n);
//     for (int i = 0; i < n; ++i) {
//         for (int j = 0; j < n; ++j) {
//             for (int k = 0; k < n; ++k) {
//                 for (int l = 0; l < n; ++l) {
//                     tensor(i, j, k, l) = mat4d[i][j](k, l);
//                 }
//             }
//         }
//     }
//     return tensor;
// }
//
// // 将Eigen::Tensor<double,4>转换回Matrix4D
// Matrix4D tensor_to_matrix4d(const Eigen::Tensor<double, 4>& tensor) {
//     int n = tensor.dimension(0);
//     Matrix4D mat4d(n, std::vector<Eigen::MatrixXd>(n, Eigen::MatrixXd::Zero(n, n)));
//     for (int i = 0; i < n; ++i) {
//         for (int j = 0; j < n; ++j) {
//             for (int k = 0; k < n; ++k) {
//                 for (int l = 0; l < n; ++l) {
//                     mat4d[i][j](k, l) = tensor(i, j, k, l);
//                 }
//             }
//         }
//     }
//     return mat4d;
// }
//
// // 优化的张量运算版本
// Matrix4D m4tom4(const Eigen::MatrixXd& pipa,
//                          const Eigen::MatrixXd& pipb,
//                          const Matrix4D& qmat) {
//     auto start = std::chrono::high_resolution_clock::now();
//     int nucnum = pipa.rows();
//
//     // 将输入转换为Eigen Tensor
//     Eigen::Tensor<double, 2> pipa_tensor = Eigen::TensorMap<Eigen::Tensor<const double, 2>>(
//         pipa.data(), nucnum, nucnum);
//     Eigen::Tensor<double, 2> pipb_tensor = Eigen::TensorMap<Eigen::Tensor<const double, 2>>(
//         pipb.data(), nucnum, nucnum);
//     Eigen::Tensor<double, 4> qmat_tensor = matrix4d_to_tensor(qmat);
//
//     // 计算term1 = pipa(a,e) * pipb(f,b) -> shape (a,b,e,f)
//     Eigen::Tensor<double, 4> term1 =
//         pipa_tensor.reshape(Eigen::array<int,4>{nucnum, 1, nucnum, 1})
//                  .broadcast(Eigen::array<int,4>{1, nucnum, 1, nucnum}) *
//         pipb_tensor.reshape(Eigen::array<int,4>{1, nucnum, 1, nucnum})
//                  .broadcast(Eigen::array<int,4>{nucnum, 1, nucnum, 1});
//
//     // 计算term2 = pipa(f,b) * pipb(a,e) -> shape (a,b,e,f)
//     Eigen::Tensor<double, 4> term2 =
//         pipa_tensor.reshape(Eigen::array<int,4>{1, nucnum, 1, nucnum})
//                  .broadcast(Eigen::array<int,4>{nucnum, 1, nucnum, 1}) *
//         pipb_tensor.reshape(Eigen::array<int,4>{nucnum, 1, nucnum, 1})
//                  .broadcast(Eigen::array<int,4>{1, nucnum, 1, nucnum});
//
//     // 合并term1和term2
//     Eigen::Tensor<double, 4> combined = term1 + term2;
//
//     // 与qmat(e,f,c,d)进行双收缩：sum over e and f
//     Eigen::array<Eigen::IndexPair<int>, 2> contract_dims = {
//         Eigen::IndexPair<int>(2, 0),  // 收缩combined的e和qmat的e
//         Eigen::IndexPair<int>(3, 1)   // 收缩combined的f和qmat的f
//     };
//     Eigen::Tensor<double, 4> result_tensor = combined.contract(qmat_tensor, contract_dims);
//     auto end = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//     std::cout << "Elapsed time: " << duration.count() << "μs\n";
//
//     // 转换回Matrix4D格式
//     return tensor_to_matrix4d(result_tensor);
// }
// Matrix4D m4tom4(Eigen::MatrixXd pipa, Eigen::MatrixXd pipb,Matrix4D qmat)
// {
//     auto start = std::chrono::high_resolution_clock::now();
//     int nucnum=nucleusm.size();
//     Matrix4D result(nucnum,
//            std::vector<Eigen::MatrixXd>(nucnum,Eigen::MatrixXd::Zero(nucnum, nucnum))  // 最内层是 5x5 矩阵
//            );
//     for (int e=0;e<nucnum;++e)
//     {
//         for (int f=0;f<nucnum;++f)
//         {
//             Eigen::MatrixXd qmat_ef = qmat[e][f]; // 预取 qmat[e][f]
//             for (int a=0;a<nucnum;++a)
//             {
//                 double pipa_ae = pipa(a, e); // 预取 pipa(a, e)
//                 double pipb_ae = pipb(a, e); // 预取 pipb(f, b)
//                 for (int b=0;b<nucnum;++b)
//                 {
//                     double pipa_fb = qmat_ef(f,b);
//                     double pipb_fb = pipa(f, b);
//                     double sumre=0;
//                     for (int c=0;c<nucnum;++c)
//                     {
//                         for (int d=0;d<nucnum;++d)
//                         {
//                             result[a][b](c,d)= result[a][b](c,d) +qmat_ef(c,d)*( pipa_ae*pipb_fb+pipa_fb*pipb_ae);
//                         }
//                     }
//
//                 }
//             }
//         }
//     }
//     auto end = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//     std::cout << "Elapsed time: " << duration.count() << "μs\n";
//     return result;
// }
// Matrix4D m4tom4(const Eigen::MatrixXd& pipa,
//                          const Eigen::MatrixXd& pipb,
//                          const Matrix4D& qmat)
// {
//     const int nucnum = nucleusm.size();
//     Matrix4D result(nucnum,
//         std::vector<Eigen::MatrixXd>(nucnum,
//             Eigen::MatrixXd::Zero(nucnum, nucnum)));
//
//     const Eigen::MatrixXd pipbT = pipb.transpose(); // 预转置
//
//     for (int a = 0; a < nucnum; ++a) {
//         for (int b = 0; b < nucnum; ++b) {
//             auto& target = result[a][b];
//             const auto pipa_a = pipa.row(a); // 预取行
//             const auto pipb_b = pipbT.row(b); // 预取列转置行
//
//             for (int c = 0; c < nucnum; ++c) {
//                 for (int d = 0; d < nucnum; ++d) {
//                     double sumre = 0;
//
//                     // 循环展开(手动或编译器优化)
//                     for (int e = 0; e < nucnum; ++e) {
//                         const double pipa_ae = pipa_a[e];
//                         const double pipb_bf = pipb_b[e];
//
//                         for (int f = 0; f < nucnum; ++f) {
//                             sumre += qmat[e][f](c,d) *
//                                    (pipa_ae * pipbT(b,f) +
//                                     pipa(f,b) * pipbT(e,a));
//                         }
//                     }
//                     target(c,d) = sumre;
//                 }
//             }
//         }
//     }
//     return result;
// }

// Matrix4D m4tom4(const Eigen::MatrixXd& pipa,
//                          const Eigen::MatrixXd& pipb,
//                          const Matrix4D& qmat)
// {
//     const int nucnum = nucleusm.size();
//     Matrix4D result(nucnum,
//         std::vector<Eigen::MatrixXd>(nucnum,
//             Eigen::MatrixXd::Zero(nucnum, nucnum)));
//
//     // 1. 预计算所有转置和中间结果
//     const Eigen::MatrixXd pipbT = pipb.transpose().eval();  // 显式求值避免重复计算
//     std::vector<Eigen::VectorXd> pipa_rows(nucnum), pipbT_rows(nucnum);
//     for(int i=0; i<nucnum; ++i) {
//         pipa_rows[i] = pipa.row(i);
//         pipbT_rows[i] = pipbT.row(i);
//     }
//
//     // 2. 分块处理（提升缓存命中率）
//     const int block_size = std::min(16, nucnum);  // 根据L1缓存调整
//
//     // 3. 并行化+向量化
// #pragma omp parallel for collapse(2) schedule(dynamic)
//     for(int a=0; a<nucnum; ++a)
//         for(int b=0; b<nucnum; ++b) {
//             auto& target = result[a][b];
//             const auto& pipa_a = pipa_rows[a];
//             const auto& pipb_b = pipbT_rows[b];
//
//             // 4. 使用Eigen的向量化操作
//             for(int c=0; c<nucnum; ++c) {
//                 for(int d=0; d<nucnum; ++d) {
//                     double sumre = 0;
//                     // 5. 循环展开和SIMD优化
// #pragma omp simd reduction(+:sumre)
//                     for(int e=0; e<nucnum; ++e) {
//                         const double* qmat_row = &qmat[e][0](c,d);
//                         const double pipa_ae = pipa_a[e];
//                         const double pipb_be = pipb_b[e];
//
//                         for(int f=0; f<nucnum; ++f) {
//                             sumre += qmat_row[f*nucnum] *
//                                    (pipa_ae * pipbT_rows[f][b] +
//                                     pipa(f,b) * pipb_be);
//                         }
//                     }
//                     target(c,d) = sumre;
//                 }
//             }
//         }
//     return result;
// }
void printAsymmetricElements(const Matrix4D& tensor) {
    const std::vector<std::string> labels = {"a", "b", "c", "d"};
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Asymmetric elements report (where tensor[a][b](c,d) != tensor[c][d](a,b)):\n";
    std::cout << "-----------------------------------------------------------\n"<<std::endl;

    for (size_t a = 0; a < tensor.size(); ++a) {
        for (size_t b = 0; b < tensor[a].size(); ++b) {
            const auto& matrix_ab = tensor[a][b];
            for (int c = 0; c < matrix_ab.rows(); ++c) {
                for (int d = 0; d < matrix_ab.cols(); ++d) {
                    // 确保c/d作为外层索引不越界
                    if (c < tensor.size() && d < tensor[c].size()) {
                        const auto& matrix_cd = tensor[c][d];
                        // 确保a/b作为内层索引不越界
                        if (a < matrix_cd.rows() && b < matrix_cd.cols()) {
                            double val_ab = matrix_ab(c, d);
                            double val_cd = matrix_cd(a, b);

                            if (std::abs(val_ab - val_cd) > 1e-9) {  // 考虑浮点误差
                                std::cout << a <<','<<b << ","
                                          << c << "," << d << "," << val_ab
                                          << "," << c<<','<< d <<","
                                          << a << "," << b <<"," << val_cd << "\n"<< std::endl;
                            }
                        }
                    }
                }
            }
        }
    }
}

// Matrix4D calg(const basism& bas1,const basism& bas2,const std::vector<Eigen::MatrixXd>&ystrm1,const std::vector<Eigen::MatrixXd>&ystrm2)
// {
//     auto start = std::chrono::high_resolution_clock::now();
//
//     int nucnum=nucleusm.size();
//     int nnp=ystrm2.size();
//     Matrix4D term1(nucnum,
//            std::vector<Eigen::MatrixXd>(nucnum,Eigen::MatrixXd::Zero(nucnum, nucnum))
//            );
//     Matrix4D term2(nucnum,
//            std::vector<Eigen::MatrixXd>(nucnum,Eigen::MatrixXd::Zero(nucnum, nucnum))
//            );
//     Matrix4D term3(nucnum,
//            std::vector<Eigen::MatrixXd>(nucnum,Eigen::MatrixXd::Zero(nucnum, nucnum))
//            );
//     Eigen::Tensor<double, 4> result1(nucnum,nucnum,nucnum,nucnum);
//     result1.setZero();
//
//
//
//     // auto start = std::chrono::high_resolution_clock::now();
//     for (int k=1;k<nnp;++k)
//     {
//         basism bas22=bas2;
//         Eigen::MatrixXd pkp=ystrm2[k];
//         std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
//         ystrm22.erase(ystrm22.begin() + k);
//         bas22.r.erase( bas22.r.begin() + k);
//         Eigen::MatrixXd bcal=calb(bas1,bas22,ystrm1,ystrm22);
//         Matrix4D t1=m2tom4(pkp,bcal);
//         Eigen::Tensor<double,4>t11=m2tom4tensor(pkp,bcal);
//         // printAsymmetricElements(t1);
//         auto end = std::chrono::high_resolution_clock::now();
//         auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//         std::cout << "Elapsed time for term1_1: " <<"k"<<k<<"\n"<< duration.count() << "μs\n"<<std::endl;
//
//         term1=addMatrix4D(term1,t1);
//         result1=result1+t11;
//
//     }
//     auto end = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//     std::cout << "Elapsed time for term1: " << duration.count() << "μs\n"<<std::endl;
//     start = std::chrono::high_resolution_clock::now();
//     for (int k=2;k<nnp;++k)
//     {
//         for (int i=1;i<=k-1;++i)
//         {
//             basism bas22=bas2;
//             std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
//             Eigen::MatrixXd pkp=ystrm2[k];
//             Eigen::MatrixXd pip=ystrm2[i];
//             bas22.r.erase( bas22.r.begin() + k);
//             bas22.r.erase( bas22.r.begin() + i);
//             ystrm22.erase(ystrm22.begin() + k);
//             ystrm22.erase(ystrm22.begin() + i);
//             Matrix4D qcal=calq(bas1,bas22,ystrm1,ystrm22);
//
//
//
//             Matrix4D qbar= gchange(qcal);
//             end = std::chrono::high_resolution_clock::now();
//             duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//             std::cout << "Elapsed time for term2: " << duration.count() << "μs\n"<<std::endl;
//
//
//
//
//             // printAsymmetricElements(qbar);
//             // printNonZeroElements(qbar);
//             Matrix4D t2=m4tom4(pkp,pip,qbar);
//             // Matrix4D t3=m4tom4(pip,pkp,qbar);
//             // printAsymmetricElements(t2);
//             term2=addMatrix4D(term2,t2);
//             // term2=addMatrix4D(term2,t3);
//         }
//     }
//
//
//
//     term1=multiplyScalar(term1,4);
//     term2=multiplyScalar(term2,96);
//     // printAsymmetricElements(term1);
//     // printAsymmetricElements(term2);
//     Matrix4D resultmat=addMatrix4D(term1, term2);
//     start= std::chrono::high_resolution_clock::now();
//     if (bas1.r[0]!=0)
//     {
//         for (int k=1;k<nnp;++k)
//         {
//             basism bas22=bas2;
//             Eigen::MatrixXd pkp=ystrm2[k];
//             std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
//             ystrm22.erase(ystrm22.begin() + k);
//             bas22.r.erase( bas22.r.begin() + k);
//             Eigen::Tensor<double, 3> t1=calt(bas1,bas22,ystrm1,ystrm22);
//             Eigen::VectorXd sp=ystrm2[0];
//             Eigen::Tensor<double, 3> tmat=tchange(t1);
//             auto matrix_vec1 = tensor3DToMatrixVector(t1);
//             auto matrix_vec = tensor3DToMatrixVector(tmat);
//             // end = std::chrono::high_resolution_clock::now();
//             // duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//             // std::cout << "Elapsed time for term3_1: " << duration.count() << "μs\n";
//             // start= std::chrono::high_resolution_clock::now();
//             Eigen::MatrixXd tit(nucnum,nucnum);
//             // for (int a=0;a<nucnum;++a)
//             // {
//             //     double spa=sp(a);
//             //     Eigen::MatrixXd t3_11= Eigen::MatrixXd::Zero(nucnum, nucnum);
//             //     if (std::abs(spa)>1e-7)
//             //     {
//             //         t3_11=spa*pkp;
//             //     }
//             //
//             //     for (int b=0;b<nucnum;++b)
//             //     {
//             //
//             //         double spb=sp(b);
//             //         Eigen::ArrayXd t3_12 = Eigen::ArrayXd::Zero(nucnum);
//             //         if (std::abs(spb)>1e-7)
//             //         {
//             //             t3_12=spb*pkp.row(a).array();
//             //         }
//             //
//             //         tit.setZero();
//             //         for (int ap=0;ap<nucnum;++ap)
//             //         {
//             //             Eigen::MatrixXd tslice=matrix_vec[ap];
//             //             double spn=t3_11(b,ap)-t3_12(ap);
//             //             // double spn=5;
//             //             if (std::abs(spn)>1e-7)
//             //             {
//             //                 tit=tit+12*spn*matrix_vec[ap];
//             //             }
//             //
//             //         }
//             //         term3[a][b]=tit;
//             //     }
//             // }
//             Eigen::MatrixXd t3_11(nucnum, nucnum);
//             Eigen::ArrayXd t3_12(nucnum);
//             Eigen::MatrixXd temp_sum(nucnum, nucnum);
//
//             for (int a=0; a<nucnum; ++a) {
//                 const double spa = sp(a);
//                 t3_11.setZero();
//                 if (spa !=0)
//                 {
//                     t3_11.noalias() = spa * pkp;
//                 }
//
//
//                 for (int b=0; b<nucnum; ++b) {
//                     const double spb = sp(b);
//                     t3_12.setZero();
//
//                     if (spb != 0)
//                     {
//                         t3_12 = spb * pkp.row(a).array();
//                     }
//
//
//
//
//                     temp_sum.setZero();
//                     if (spb != 0 || spa != 0)
//                     {
//                         for (int ap=0; ap<nucnum; ++ap) {
//                             const double spn = t3_11(b,ap) - t3_12(ap);
//                             if (spn != 0)
//                             {
//                                 temp_sum.noalias() += 12 * spn * matrix_vec[ap];
//                             }
//                         }
//                     }
//
//                     term3[a][b] += temp_sum;
//                 }
//             }
//         }
//         end = std::chrono::high_resolution_clock::now();
//         duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//         std::cout << "Elapsed time for term3: " << duration.count() << "μs\n";
//         resultmat=addMatrix4D(resultmat,term3);
//     }
//     // auto end = std::chrono::high_resolution_clock::now();
//     // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//     // std::cout << "gcal time: " << duration.count() << "μs\n";
//
//     return resultmat;
//
// }


// Matrix4D calg(const basism& bas1,const basism& bas2,const std::vector<Eigen::MatrixXd>&ystrm1,const std::vector<Eigen::MatrixXd>&ystrm2)
// {
//     // auto start = std::chrono::high_resolution_clock::now();
//
//     int nucnum=nucleusm.size();
//     int nnp=ystrm2.size();
//     Matrix4D term1(nucnum,
//            std::vector<Eigen::MatrixXd>(nucnum,Eigen::MatrixXd::Zero(nucnum, nucnum))
//            );
//     Matrix4D term2(nucnum,
//            std::vector<Eigen::MatrixXd>(nucnum,Eigen::MatrixXd::Zero(nucnum, nucnum))
//            );
//     Matrix4D term3(nucnum,
//            std::vector<Eigen::MatrixXd>(nucnum,Eigen::MatrixXd::Zero(nucnum, nucnum))
//            );
//     auto start = std::chrono::high_resolution_clock::now();
//     for (int k=1;k<nnp;++k)
//     {
//         basism bas22=bas2;
//         Eigen::MatrixXd pkp=ystrm2[k];
//         std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
//         ystrm22.erase(ystrm22.begin() + k);
//         bas22.r.erase( bas22.r.begin() + k);
//         Eigen::MatrixXd bcal=calb(bas1,bas22,ystrm1,ystrm22);
//
//
//
//
//         Matrix4D t1=m2tom4(pkp,bcal);
//         // printAsymmetricElements(t1);
//         auto end = std::chrono::high_resolution_clock::now();
//         auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//         // std::cout << "Elapsed time for term1_1: " <<"k"<<k<<"\n"<< duration.count() << "μs\n";
//
//
//
//         term1=addMatrix4D(term1,t1);
//
//     }
//     auto end = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//     // std::cout << "Elapsed time for term1: " << duration.count() << "μs\n"<<std::endl;
//     start = std::chrono::high_resolution_clock::now();
//     for (int k=2;k<nnp;++k)
//     {
//         for (int i=1;i<=k-1;++i)
//         {
//             basism bas22=bas2;
//             std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
//             Eigen::MatrixXd pkp=ystrm2[k];
//             Eigen::MatrixXd pip=ystrm2[i];
//             bas22.r.erase( bas22.r.begin() + k);
//             bas22.r.erase( bas22.r.begin() + i);
//             ystrm22.erase(ystrm22.begin() + k);
//             ystrm22.erase(ystrm22.begin() + i);
//             Matrix4D qcal=calq(bas1,bas22,ystrm1,ystrm22);
//
//
//
//             Matrix4D qbar= gchange(qcal);
//             end = std::chrono::high_resolution_clock::now();
//             duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//             // std::cout << "Elapsed time for term2: " << duration.count() << "μs\n"<<std::endl;
//
//
//
//
//             // printAsymmetricElements(qbar);
//             // printNonZeroElements(qbar);
//             Matrix4D t2=m4tom4(pkp,pip,qbar);
//             // Matrix4D t3=m4tom4(pip,pkp,qbar);
//             // printAsymmetricElements(t2);
//             term2=addMatrix4D(term2,t2);
//             // term2=addMatrix4D(term2,t3);
//         }
//     }
//
//
//
//     term1=multiplyScalar(term1,4);
//     term2=multiplyScalar(term2,96);
//     // printAsymmetricElements(term1);
//     // printAsymmetricElements(term2);
//     Matrix4D resultmat=addMatrix4D(term1, term2);
//     start= std::chrono::high_resolution_clock::now();
//     if (bas1.r[0]!=0)
//     {
//         for (int k=1;k<nnp;++k)
//         {
//             basism bas22=bas2;
//             Eigen::MatrixXd pkp=ystrm2[k];
//             std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
//             ystrm22.erase(ystrm22.begin() + k);
//             bas22.r.erase( bas22.r.begin() + k);
//             Eigen::Tensor<double, 3> t1=calt(bas1,bas22,ystrm1,ystrm22);
//             Eigen::VectorXd sp=ystrm2[0];
//             Eigen::Tensor<double, 3> tmat=tchange(t1);
//             auto matrix_vec1 = tensor3DToMatrixVector(t1);
//             auto matrix_vec = tensor3DToMatrixVector(tmat);
//             // end = std::chrono::high_resolution_clock::now();
//             // duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//             // std::cout << "Elapsed time for term3_1: " << duration.count() << "μs\n";
//             // start= std::chrono::high_resolution_clock::now();
//             Eigen::MatrixXd tit(nucnum,nucnum);
//             // for (int a=0;a<nucnum;++a)
//             // {
//             //     double spa=sp(a);
//             //     Eigen::MatrixXd t3_11= Eigen::MatrixXd::Zero(nucnum, nucnum);
//             //     if (std::abs(spa)>1e-7)
//             //     {
//             //         t3_11=spa*pkp;
//             //     }
//             //
//             //     for (int b=0;b<nucnum;++b)
//             //     {
//             //
//             //         double spb=sp(b);
//             //         Eigen::ArrayXd t3_12 = Eigen::ArrayXd::Zero(nucnum);
//             //         if (std::abs(spb)>1e-7)
//             //         {
//             //             t3_12=spb*pkp.row(a).array();
//             //         }
//             //
//             //         tit.setZero();
//             //         for (int ap=0;ap<nucnum;++ap)
//             //         {
//             //             Eigen::MatrixXd tslice=matrix_vec[ap];
//             //             double spn=t3_11(b,ap)-t3_12(ap);
//             //             // double spn=5;
//             //             if (std::abs(spn)>1e-7)
//             //             {
//             //                 tit=tit+12*spn*matrix_vec[ap];
//             //             }
//             //
//             //         }
//             //         term3[a][b]=tit;
//             //     }
//             // }
//             Eigen::MatrixXd t3_11(nucnum, nucnum);
//             Eigen::ArrayXd t3_12(nucnum);
//             Eigen::MatrixXd temp_sum(nucnum, nucnum);
//
//             for (int a=0; a<nucnum; ++a) {
//                 const double spa = sp(a);
//                 t3_11.setZero();
//                 if (spa !=0)
//                 {
//                     t3_11.noalias() = spa * pkp;
//                 }
//
//
//                 for (int b=0; b<nucnum; ++b) {
//                     const double spb = sp(b);
//                     t3_12.setZero();
//
//                     if (spb != 0)
//                     {
//                         t3_12 = spb * pkp.row(a).array();
//                     }
//
//
//
//
//                     temp_sum.setZero();
//                     if (spb != 0 || spa != 0)
//                     {
//                         for (int ap=0; ap<nucnum; ++ap) {
//                             const double spn = t3_11(b,ap) - t3_12(ap);
//                             if (spn != 0)
//                             {
//                                 temp_sum.noalias() += 12 * spn * matrix_vec[ap];
//                             }
//                         }
//                     }
//
//                     term3[a][b] += temp_sum;
//                 }
//             }
//         }
//         end = std::chrono::high_resolution_clock::now();
//         duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//         // std::cout << "Elapsed time for term3: " << duration.count() << "μs\n";
//         resultmat=addMatrix4D(resultmat,term3);
//     }
//     // auto end = std::chrono::high_resolution_clock::now();
//     // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//     // std::cout << "gcal time: " << duration.count() << "μs\n";
//
//     return resultmat;
//
// }

Matrix4D_sp calg_sp(const basism& bas1,const basism& bas2,const std::vector<Eigen::MatrixXd>&ystrm1,const std::vector<Eigen::MatrixXd>&ystrm2)
{
    // auto start = std::chrono::high_resolution_clock::now();

    int nucnum=nucleusm.size();
    int nnp=ystrm2.size();
    Matrix4D_sp term1=createZeroMatrix4D_sp(nucnum,nucnum,nucnum,nucnum);
    Matrix4D_sp term2=term1;
    Matrix4D_sp term3=term1;
    Matrix4D_sp resultmatsp=term1;
    auto start = std::chrono::high_resolution_clock::now();
    for (int k=1;k<nnp;++k)
    {
        basism bas22=bas2;
        Eigen::MatrixXd pkp=ystrm2[k];
        std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
        ystrm22.erase(ystrm22.begin() + k);
        bas22.r.erase( bas22.r.begin() + k);
        Eigen::MatrixXd bcal=calb(bas1,bas22,ystrm1,ystrm22);
        m2tom4sp(pkp,bcal,resultmatsp);
        // printAsymmetricElements(t1);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        // std::cout << "Elapsed time for term1_1: " <<"k"<<k<<"\n"<< duration.count() << "μs\n";
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    // std::cout << "Elapsed time for term1: " << duration.count() << "μs\n"<<std::endl;
    start = std::chrono::high_resolution_clock::now();

    for (int k=2;k<nnp;++k)
    {
        for (int i=1;i<=k-1;++i)
        {
            basism bas22=bas2;
            std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
            Eigen::MatrixXd pkp=ystrm2[k];
            Eigen::MatrixXd pip=ystrm2[i];
            bas22.r.erase( bas22.r.begin() + k);
            bas22.r.erase( bas22.r.begin() + i);
            ystrm22.erase(ystrm22.begin() + k);
            ystrm22.erase(ystrm22.begin() + i);
            Matrix4D qcal=convertToDense(  calq(bas1,bas22,ystrm1,ystrm22));
            Matrix4D qbar= gchange(qcal);
            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            // std::cout << "Elapsed time for term2: " << duration.count() << "μs\n"<<std::endl;

            m4tom4sp(pkp,pip,qbar,resultmatsp);
        }
    }

    start= std::chrono::high_resolution_clock::now();
    if (bas1.r[0]!=0)
    {
        for (int k=1;k<nnp;++k)
        {
            basism bas22=bas2;
            Eigen::MatrixXd pkp=ystrm2[k];
            std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
            ystrm22.erase(ystrm22.begin() + k);
            bas22.r.erase( bas22.r.begin() + k);
            Eigen::Tensor<double, 3> t1=calt(bas1,bas22,ystrm1,ystrm22);
            Eigen::VectorXd sp=ystrm2[0];
            Eigen::Tensor<double, 3> tmat=tchange(t1);
            auto matrix_vec1 = tensor3DToMatrixVector(t1);
            auto matrix_vec = tensor3DToMatrixVector(tmat);
            Eigen::MatrixXd t3_11(nucnum, nucnum);
            Eigen::ArrayXd t3_12(nucnum);
            Eigen::MatrixXd temp_sum(nucnum, nucnum);



            for (int a=0; a<nucnum; ++a) {
                const double spa = sp(a);
                t3_11.setZero();
                if (spa !=0)
                {
                    t3_11.noalias() = spa * pkp;
                }
                for (int b=0; b<nucnum; ++b) {
                    const double spb = sp(b);
                    t3_12.setZero();

                    if (spb != 0)
                    {
                        t3_12 = spb * pkp.row(a).array();
                    }
                    temp_sum.setZero();
                    if (spb != 0 || spa != 0)
                    {
                        for (int ap=0; ap<nucnum; ++ap) {
                            const double spn = t3_11(b,ap) - t3_12(ap);
                            if (spn != 0)
                            {
                                Eigen::SparseMatrix<double> matrix_vec_ap = matrix_vec[ap].sparseView();
                                resultmatsp[a][b] += 12 * spn * matrix_vec_ap;
                            }
                        }
                    }
                }
            }
        }
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        // std::cout << "Elapsed time for term3: " << duration.count() << "μs\n";
    }
    // auto end = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    // std::cout << "gcal time: " << duration.count() << "μs\n";

    return resultmatsp;

}


// Matrix4D calg(const basism& bas1,const basism& bas2,const std::vector<Eigen::MatrixXd>&ystrm1,const std::vector<Eigen::MatrixXd>&ystrm2)
// {
//     auto start = std::chrono::high_resolution_clock::now();
//
//     int nucnum=nucleusm.size();
//     int nnp=ystrm2.size();
//     Eigen::Tensor<double, 4> result1(nucnum,nucnum,nucnum,nucnum);
//     result1.setZero();
//     Eigen::Tensor<double, 4> result2(nucnum,nucnum,nucnum,nucnum);
//     result2.setZero();
//     Eigen::Tensor<double, 4> result3(nucnum,nucnum,nucnum,nucnum);
//     result3.setZero();
//
//
//
//     // auto start = std::chrono::high_resolution_clock::now();
//     for (int k=1;k<nnp;++k)
//     {
//         basism bas22=bas2;
//         Eigen::MatrixXd pkp=ystrm2[k];
//         std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
//         ystrm22.erase(ystrm22.begin() + k);
//         bas22.r.erase( bas22.r.begin() + k);
//         Eigen::MatrixXd bcal=calb(bas1,bas22,ystrm1,ystrm22);
//         // Matrix4D t1=m2tom4(pkp,bcal);
//         Eigen::Tensor<double,4>t11=m2tom4tensor(pkp,bcal);
//         // printAsymmetricElements(t1);
//         auto end = std::chrono::high_resolution_clock::now();
//         auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//         std::cout << "Elapsed time for term1_1: " <<"k"<<k<<"\n"<< duration.count() << "μs\n"<<std::endl;
//
//         // term1=addMatrix4D(term1,t1);
//         result1=result1+t11;
//
//     }
//     auto end = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//     std::cout << "Elapsed time for term1: " << duration.count() << "μs\n"<<std::endl;
//     start = std::chrono::high_resolution_clock::now();
//     for (int k=2;k<nnp;++k)
//     {
//         for (int i=1;i<=k-1;++i)
//         {
//             basism bas22=bas2;
//             std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
//             Eigen::MatrixXd pkp=ystrm2[k];
//             Eigen::MatrixXd pip=ystrm2[i];
//             bas22.r.erase( bas22.r.begin() + k);
//             bas22.r.erase( bas22.r.begin() + i);
//             ystrm22.erase(ystrm22.begin() + k);
//             ystrm22.erase(ystrm22.begin() + i);
//             Matrix4D qcal=calq(bas1,bas22,ystrm1,ystrm22);
//             Matrix4D qbar= gchange(qcal);
//             Tensor4D qbar_tensor = matrix4d_to_tensor(qbar);
//
//             end = std::chrono::high_resolution_clock::now();
//             duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//             std::cout << "Elapsed time for term2: " << duration.count() << "μs\n"<<std::endl;
//
//             // printAsymmetricElements(qbar);
//             // printNonZeroElements(qbar);
//             // Matrix4D t2=m4tom4(pkp,pip,qbar);
//             Tensor4D t2_tensor = m4tom4_tensor(pkp,pip,qbar_tensor);
//             // Matrix4D t3=m4tom4(pip,pkp,qbar);
//             // printAsymmetricElements(t2);
//             // term2=addMatrix4D(term2,t2);
//             // term2=addMatrix4D(term2,t3);
//             result2=result2+t2_tensor;
//         }
//     }
//     result1=4.0*result1;
//     result2=96.0*result2;
//     Tensor4D result_tensor = result1 + result2;
//
//
//
//     // term1=multiplyScalar(term1,4);
//     // term2=multiplyScalar(term2,96);
//     // printAsymmetricElements(term1);
//     // printAsymmetricElements(term2);
//     start= std::chrono::high_resolution_clock::now();
//     if (bas1.r[0]!=0)
//     {
//         for (int k=1;k<nnp;++k)
//         {
//             basism bas22=bas2;
//             Eigen::MatrixXd pkp=ystrm2[k];
//             std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
//             ystrm22.erase(ystrm22.begin() + k);
//             bas22.r.erase( bas22.r.begin() + k);
//             Eigen::Tensor<double, 3> t1=calt(bas1,bas22,ystrm1,ystrm22);
//             Eigen::VectorXd sp=ystrm2[0];
//             Eigen::Tensor<double, 3> tmat=tchange(t1);
//             if (sp.isZero(1e-12) || pkp.isZero(1e-12)) continue;
//
//             // 3.2 使用张量广播进行向量化计算
//             // 重塑sp为矩阵便于广播
//             Eigen::MatrixXd sp_mat = sp * Eigen::RowVectorXd::Ones(nucnum);
//
//             // 计算系数张量 (a, b, ap)
//             Eigen::Tensor<double, 3> coeff_tensor(nucnum, nucnum, nucnum);
//             coeff_tensor.setZero();
//
//             // 正确方法：使用TensorMap创建视图
//             // 映射sp_mat为三维张量 (a, 1, ap)
//             Eigen::TensorMap<const Eigen::Tensor<double, 3>> sp_mat_map_a(
//                 sp_mat.data(), nucnum, 1, nucnum
//             );
//
//             // 映射pkp为三维张量 (1, b, ap) - 注意数据布局
//             Eigen::TensorMap<const Eigen::Tensor<double, 3>> pkp_map_b(
//                 pkp.data(), 1, nucnum, nucnum
//             );
//
//             // 第一部分: spa * pkp(b, ap)
//             coeff_tensor = sp_mat_map_a.broadcast(Eigen::array<long, 3>{1, nucnum, 1}) *
//                            pkp_map_b.broadcast(Eigen::array<long, 3>{nucnum, 1, 1});
//
//             // 第二部分: -spb * pkp(a, ap)
//             // 映射sp_mat为三维张量 (1, b, ap)
//             Eigen::TensorMap<const Eigen::Tensor<double, 3>> sp_mat_map_b(
//                 sp_mat.data(), 1, nucnum, nucnum
//             );
//
//             // 映射pkp为三维张量 (a, 1, ap)
//             Eigen::TensorMap<const Eigen::Tensor<double, 3>> pkp_map_a(
//                 pkp.data(), nucnum, 1, nucnum
//             );
//
//             coeff_tensor -= sp_mat_map_b.broadcast(Eigen::array<long, 3>{nucnum, 1, 1}) *
//                             pkp_map_a.broadcast(Eigen::array<long, 3>{1, nucnum, 1});
//
//             // 3.3 张量收缩: sum_ap [coeff(a,b,ap) * tmat(ap,c,d)]
//             Eigen::array<Eigen::IndexPair<int>, 1> contract_dims =
//                 {Eigen::IndexPair<int>(2, 0)};  // 收缩ap维度
//             Tensor4D increment = 12.0 * coeff_tensor.contract(tmat, contract_dims);
//
//             // 3.4 累加到结果
//             result3 += increment;
//         }
//         end = std::chrono::high_resolution_clock::now();
//         duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//         std::cout << "Elapsed time for term3: " << duration.count() << "μs\n";
//         result_tensor=result3 + result_tensor;
//
//     }
//     // auto end = std::chrono::high_resolution_clock::now();
//     // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//     // std::cout << "gcal time: " << duration.count() << "μs\n";
//     Matrix4D resultmat= tensor_to_matrix4d(result_tensor);
//
//     return resultmat;
//
// }

Eigen::MatrixXd calf(const basism& bas1,const basism& bas2,const std::vector<Eigen::MatrixXd>&ystrm1,const std::vector<Eigen::MatrixXd>&ystrm2)
{
    int nucnum=nucleusm.size();
    int sizer=ystrm1.size();
    Eigen::MatrixXd term1=Eigen::MatrixXd::Zero(nucnum, nucnum);
    Eigen::MatrixXd term2=Eigen::MatrixXd::Zero(nucnum, nucnum);
    for (int k=1;k<sizer;++k)
    {
        // if (!outfile.is_open())
        // {
        //     // 如果文件未打开，则尝试打开文件
        //     outfile.open("basis_output.txt",std::ios::app);
        // }
        // outfile<< "calf,term1, k: " << k << "\n"<<std::endl;
        basism bas22=bas2;
        Eigen::MatrixXd pkp=ystrm2[k];
        std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
        ystrm22.erase(ystrm22.begin() + k);
        bas22.r.erase( bas22.r.begin() + k);
        Eigen::MatrixXd bcal=calb(bas1,bas22,ystrm1,ystrm22);
        Eigen::MatrixXd bcalt= bcal.transpose();

        Eigen::MatrixXd t1;
        t1=pkp*bcalt;
        term1=t1+term1;
        // if (!outfile.is_open())
        // {
        //     // 如果文件未打开，则尝试打开文件
        //     outfile.open("basis_output.txt",std::ios::app);
        // }
        // outfile<<"term1matrix"<< ", k: " << k << "\n";
        // writeMatrixToFile(t1);

    }

    term1=term1*4;
    if (bas1.r[0]!=0)
    {
        // if (!outfile.is_open())
        // {
        //     // 如果文件未打开，则尝试打开文件
        //     outfile.open("basis_output.txt",std::ios::app);
        // }
        // outfile<< "calf, term2: "  << "\n"<<std::endl;
        basism bas22=bas2;
        Eigen::MatrixXd pnp=ystrm2.back();
        std::vector<Eigen::MatrixXd>ystrm22=ystrm2;
        ystrm22.pop_back();
        bas22.r.pop_back();
        Eigen::Tensor<double, 3> t1=calt(bas1,bas22,ystrm1,ystrm22);
        Eigen::Tensor<double,3>tchangecal=tchange(t1);
        Eigen::VectorXd sp=ystrm2[0];
        for (int a=0;a<nucnum;++a)
        {
            double spa=sp(a);
            for (int b=0;b<nucnum;++b)
            {
                double matsum=0;
                for (int c=0;c<nucnum;++c)
                {
                    for (int d=0;d<nucnum;++d)
                    {
                        matsum+=6*spa*tchangecal(b,c,d)*pnp(c,d);
                    }
                }
                term2(a,b)=matsum;
            }
        }
        // if (!outfile.is_open())
        // {
        //     // 如果文件未打开，则尝试打开文件
        //     outfile.open("basis_output.txt",std::ios::app);
        // }
        // outfile<< ", term2 size: " << term2.rows() << "x" << term2.cols() << "\n";
        // writeMatrixToFile(term2);
        term1=term1+term2;
    }
    return term1;

}

// std::vector<std::vector<Matrix4D>> calgall(const std::vector<basism>& allbasisc,const std::vector<std::vector<Eigen::MatrixXd>>& ystrm)
// {
//     // int n=allbasisc.size();
//     // std::vector<std::vector<Matrix4D>> result(n, std::vector<Matrix4D>(n));
//     // #pragma omp parallel for collapse(2)
//     // for (int i=0;i<n;++i)
//     // {
//     //     // std::cout << "Calculating g for basis " << i << " of " << n << std::endl;
//     //     for (int j=0;j<n;++j)
//     //     {
//     //         // std::cout << "Calculating g for basis pair (" << i << ", " << j << ")" << std::endl;
//     //         result[i][j]=calg(allbasisc[i],allbasisc[j],ystrm[i],ystrm[j]);
//     //         // printAsymmetricElements(result[i][j]);
//     //         // printNonZeroElements(result[i][j]);
//     //     }
//     // }
//     // return result;
//     int n = allbasisc.size();
//     // 使用一维连续存储优化缓存
//     std::vector<Matrix4D> result_flat(n * n);
//
// #pragma omp parallel for collapse(2) schedule(dynamic, 1)
//     for (int i = 0; i < n; ++i) {
//         for (int j = 0; j < n; ++j) {
//             // 使用局部变量避免竞态
//             Matrix4D tmp = calg(allbasisc[i], allbasisc[j], ystrm[i], ystrm[j]);
//             result_flat[i * n + j] = tmp;
//         }
//     }
//
//     // 转换为嵌套结构（如需）
//     std::vector<std::vector<Matrix4D>> result(n, std::vector<Matrix4D>(n));
//     for (int i = 0; i < n; ++i) {
//         for (int j = 0; j < n; ++j) {
//             result[i][j] = result_flat[i * n + j];
//         }
//     }
//     return result;
// }

std::vector<std::vector<Eigen::MatrixXd>> calfall(const std::vector<basism>& allbasisc,const std::vector<std::vector<Eigen::MatrixXd>>& ystrm)
{
    int n=allbasisc.size();
    std::vector<std::vector<Eigen::MatrixXd>> result(n, std::vector<Eigen::MatrixXd>(n));
    for (int i=0;i<n;++i)
    {
        for (int j=0;j<n;++j)
        {
            result[i][j]=calf(allbasisc[i],allbasisc[j],ystrm[i],ystrm[j]);
        }
    }
    return result;
}


Matrix4D calo(std::map<int, Matrix4D> V_value)
{
    int nucnum=nucleusm.size();
    Matrix4D result(nucnum,
           std::vector<Eigen::MatrixXd>(nucnum,Eigen::MatrixXd::Zero(nucnum, nucnum))  // 最内层是 5x5 矩阵
           );
    for (int a=0;a<nucnum;++a)
    {
        int ja=nucleusm[a].j;
        int numa=nucleusm[a].num;
        int ma= nucleusm[a].m;
        for (int b=0;b<nucnum;++b)
        {
            int jb=nucleusm[b].j;
            int numb=nucleusm[b].num;
            int mb= nucleusm[b].m;
            for (int c=0;c<nucnum;++c)
            {
                int jc=nucleusm[c].j;
                int numc=nucleusm[c].num;
                int mc= nucleusm[c].m;
                for (int d=0;d<nucnum;++d)
                {
                    int jd=nucleusm[d].j;
                    int numd=nucleusm[d].num;
                    int md= nucleusm[d].m;
                    double c1=mysqrt((1+deltatwo(numa,numb))*(1+deltatwo(numc,numd)))/4.0;
                    double matvalsum=0;
                    for (const auto& v_pair : V_value) {
                        int Jsum = v_pair.first;
                        const Matrix4D& v_matrix = v_pair.second;
                        int m1=ma+mb;
                        int m2=mc+md;
                        if (m1==m2)
                        {
                            double cj1=util::CG(ja,jb,Jsum,ma,mb,m1);
                            double cj2=util::CG(jc,jd,Jsum,mc,md,m2);
                            double matval=cj1*cj2*v_matrix[numa][numb](numc,numd)*c1;
                            matvalsum+=matval;

                        }
                    }
                    result[a][b](c,d)=matvalsum;
                    // result[c][d](a,b)=matvalsum;  // 确保对称性
                    // result[b][a](c,d)=std::pow(-1,ja+jb-Jsum)*matvalsum;  // 确保对称性
                }
            }
        }
    }
    return result;
}

Eigen::MatrixXd calqsige(const basism& bas1,const basism& bas2,std::vector<double>envec)
{
    int nucnum=nucleusm.size();
    Eigen::MatrixXd result=Eigen::MatrixXd::Zero(nucnum, nucnum);
    for (int a=0;a<nucnum;++a)
    {
        int ja=nucleusm[a].j;
        int numa=nucleusm[a].num;
        int ma= nucleusm[a].m;
        for (int b=0;b<nucnum;++b)
        {
            int jb=nucleusm[b].j;
            int numb=nucleusm[b].num;
            int mb= nucleusm[b].m;
            double c1=deltatwo(a,b);
            double matvalsum=envec[numa]*c1;
            result(a,b)=matvalsum;

        }
    }
    return result;
}

Matrix4D multip4dm(Matrix4D a,Matrix4D b)
{
    // auto start = std::chrono::high_resolution_clock::now();

    Matrix4D result= a;  // 先复制原数据
    for (int i=0;i<a.size();++i)
    {
        for (int j=0;j<a[i].size();++j)
        {
            result[i][j]=a[i][j].cwiseProduct(b[i][j]);
        }
    }
    // auto end = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    // std::cout << "calho time: " << duration.count() << "μs\n"<<std::endl;
    return result;

}

Matrix4D multip4dm_sp(const Matrix4D_sp& a, const Matrix4D& b)
{
    assert(a.size() == b.size() && "Outer dimension mismatch");

    Matrix4D result;
    result.resize(a.size());

    for (size_t i = 0; i < a.size(); ++i) {
        assert(a[i].size() == b[i].size() && "Middle dimension mismatch");
        result[i].resize(a[i].size());

        for (size_t j = 0; j < a[i].size(); ++j) {
            // 将稠密矩阵转换为稀疏矩阵
            Eigen::SparseMatrix<double> sparse_a = a[i][j];
            Eigen::SparseMatrix<double> sparse_b = b[i][j].sparseView();

            // 执行稀疏矩阵的逐元素乘法
            Eigen::SparseMatrix<double> sparse_result = sparse_a.cwiseProduct(sparse_b);

            // 将结果转换回稠密矩阵
            result[i][j] = Eigen::MatrixXd(sparse_result);
        }
    }
    return result;
}

Matrix4D multip4dm_sp(const Matrix4D& a, const Matrix4D& b)
{
    assert(a.size() == b.size() && "Outer dimension mismatch");

    Matrix4D result;
    result.resize(a.size());

    for (size_t i = 0; i < a.size(); ++i) {
        assert(a[i].size() == b[i].size() && "Middle dimension mismatch");
        result[i].resize(a[i].size());

        for (size_t j = 0; j < a[i].size(); ++j) {
            // 将稠密矩阵转换为稀疏矩阵
            Eigen::SparseMatrix<double> sparse_a = a[i][j].sparseView();
            Eigen::SparseMatrix<double> sparse_b = b[i][j].sparseView();

            // 执行稀疏矩阵的逐元素乘法
            Eigen::SparseMatrix<double> sparse_result = sparse_a.cwiseProduct(sparse_b);

            // 将结果转换回稠密矩阵
            result[i][j] = Eigen::MatrixXd(sparse_result);
        }
    }
    return result;
}

// Eigen::MatrixXd ham1(const std::vector<basism>& bas,const std::vector<std::vector<Eigen::MatrixXd>>&ystrm,
//     std::vector<std::vector<Matrix4D>> galltest,std::vector<std::map<int, Matrix4D>> buildVValue1,
//     std::vector<double>strength)
// {
//     int nucnum=nucleusm.size();
//     int sizebas=bas.size();
//     Eigen::MatrixXd result=Eigen::MatrixXd::Zero(sizebas, sizebas);
//     for (int t=0;t<buildVValue1.size();++t)
//     {
//         Matrix4D ocal=calo( buildVValue1[t]);
//         for (int l=0;l<sizebas;++l)
//         {
//             basism bas1=bas[l];
//             std::vector<Eigen::MatrixXd> ystrm1=ystrm[l];
//             for (int m=0;m<sizebas;++m)
//             {
//                 basism bas2=bas[m];
//                 std::vector<Eigen::MatrixXd> ystrm2=ystrm[m];
//
//                 // std::cout<<"t="<<t<<",l="<<l<<",m="<<m<<std::endl;
//                 // std::cout<<"ocal"<<std::endl;
//                 // printNonZeroElements(ocal);
//                 // std::cout<<"gs"<<std::endl;
//                 Matrix4D gs=galltest[l][m];
//                 // printNonZeroElements(gs);
//                 Matrix4D hcal=multip4dm(gs,ocal);
//                 double hsum=0;
//                 for (int b=0; b<nucnum;++b)
//                 {
//                     for (int a=0;a<=b;++a)
//                     {
//                         for (int d=0;d<nucnum;++d)
//                         {
//                             for (int c=0;c<=d;++c)
//                             {
//                                 hsum+=hcal[a][b](c,d);
//                             }
//                         }
//                     }
//                 }
//                 result(l,m)+=hsum*strength[t];
//             }
//
//         }
//     }
//     return result;
// }

// Eigen::MatrixXd ham1(const std::vector<basism>& bas,const std::vector<std::vector<Eigen::MatrixXd>>&ystrm,
//     std::vector<std::map<int, Matrix4D>> buildVValue1,
//     std::vector<double>strength)
// {
//     auto start = std::chrono::high_resolution_clock::now();
//     int nucnum=nucleusm.size();
//     int sizebas=bas.size();
//     const int t_size = buildVValue1.size();
//     Eigen::MatrixXd result=Eigen::MatrixXd::Zero(sizebas, sizebas);
//     std::vector<Matrix4D> ocal_cache(t_size);
//     #pragma omp parallel for schedule(dynamic)
//     for (int t = 0; t < t_size; ++t) {
//         ocal_cache[t] = calo(buildVValue1[t]);
//     }
//     auto end = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//     std::cout << "Elapsed time for ham1_1: " << duration.count() << "μs\n";
//
//
//     #pragma omp parallel for
//     for (int l=0;l<sizebas;++l)
//     {
//         basism bas1=bas[l];
//         std::vector<Eigen::MatrixXd> ystrm1=ystrm[l];
//         for (int m=0;m<sizebas;++m)
//         {
//             basism bas2=bas[m];
//             std::vector<Eigen::MatrixXd> ystrm2=ystrm[m];
//             if (bas1.parity != bas2.parity) {
//                 continue; // 如果基态的奇偶性不同，跳过计算
//             }
//
//             // std::cout<<"t="<<t<<",l="<<l<<",m="<<m<<std::endl;
//             // std::cout<<"ocal"<<std::endl;
//             // printNonZeroElements(ocal);
//             // std::cout<<"gs"<<std::endl;
//             Matrix4D gs=calg( bas[l],bas[m],ystrm[l],ystrm[m]);
//             // Matrix4D_sp gs_sp=calg_sp(bas[l],bas[m],ystrm[l],ystrm[m]);
//             for (int t=0;t<buildVValue1.size();++t)
//             {
//                 Matrix4D ocal=ocal_cache[t];
//
//
//                 // printNonZeroElements(gs);
//                 Matrix4D hcal=multip4dm(gs,ocal);
//                 double hsum=0;
//                 // auto start1 = std::chrono::high_resolution_clock::now();
//                 for (int b=0; b<nucnum;++b)
//                 {
//                     for (int a=0;a<=b;++a)
//                     {
//                         for (int d=0;d<nucnum;++d)
//                         {
//                             for (int c=0;c<=d;++c)
//                             {
//                                 hsum+=hcal[a][b](c,d);
//                             }
//                         }
//                     }
//                 }
//                 // auto end1 = std::chrono::high_resolution_clock::now();
//                 // auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);
//                 // std::cout << "Elapsed time for hamsum: " << duration1.count() << "μs\n";
//                 result(l,m)+=hsum*strength[t];
//             }
//             // std::cout<<"l="<<l<<std::endl;
//             // std::cout<<"result="<<result<<std::endl;
//
//         }
//     }
//     end = std::chrono::high_resolution_clock::now();
//     duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//     std::cout << "Elapsed time for ham1: " << duration.count() << "μs\n";
//     return result;
// }

// 在并行区域外声明进度计数器
std::string format_time(int seconds) {
    if (seconds < 0) {
        return "negative time";
    }
    if (seconds == 0) {
        return "0s";
    }

    int hours = seconds / 3600;
    int minutes = (seconds % 3600) / 60;
    int secs = seconds % 60;

    std::string result;

    if (hours > 0) {
        result += std::to_string(hours) + "h ";
    }
    if (minutes > 0 || hours > 0) {
        if (hours > 0) {
            result += std::to_string(hours) + "h ";
        }
        if (minutes > 0) {
            result += std::to_string(minutes) + "m ";
        }
        result += std::to_string(secs) + "s";
    } else {
        // 没有小时，分钟也没有（即分钟为0），那么只显示秒
        result += std::to_string(secs) + "s";
    }

    return result;
}
Eigen::MatrixXd ham1(const std::vector<basism>& bas,const std::vector<std::vector<Eigen::MatrixXd>>&ystrm,
    std::vector<std::map<int, Matrix4D>> buildVValue1,
    std::vector<double>strength)
{
    auto start = std::chrono::high_resolution_clock::now();

    int nucnum=nucleusm.size();
    int sizebas=bas.size();
    const int t_size = buildVValue1.size();
    Eigen::MatrixXd result=Eigen::MatrixXd::Zero(sizebas, sizebas);
    std::vector<Matrix4D> ocal_cache(t_size);
    #pragma omp parallel for schedule(dynamic)
    for (int t = 0; t < t_size; ++t) {
        ocal_cache[t] = calo(buildVValue1[t]);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Elapsed time for ham1_1: " << duration.count() << "μs\n";


    size_t total_pairs = sizebas * sizebas;
    std::atomic<size_t> completed_pairs(0);
    const size_t update_interval = std::max<size_t>(10, total_pairs / 100);
    auto start_time = std::chrono::steady_clock::now();

    std::cout << "Starting calculation of " << total_pairs << " pairs...\n";


    #pragma omp parallel for collapse(2)
    for (int l=0;l<sizebas;++l)
    {
        basism bas1=bas[l];
        std::vector<Eigen::MatrixXd> ystrm1=ystrm[l];
        for (int m=0;m<sizebas;++m)
        {
            basism bas2=bas[m];
            std::vector<Eigen::MatrixXd> ystrm2=ystrm[m];
            if (bas1.parity != bas2.parity) {
                continue; // 如果基态的奇偶性不同，跳过计算
            }

            // std::cout<<"t="<<t<<",l="<<l<<",m="<<m<<std::endl;
            // std::cout<<"ocal"<<std::endl;
            // printNonZeroElements(ocal);
            // std::cout<<"gs"<<std::endl;

            Matrix4D_sp gs_sp=calg_sp(bas[l],bas[m],ystrm[l],ystrm[m]);
            for (int t=0;t<buildVValue1.size();++t)
            {
                Matrix4D ocal=ocal_cache[t];


                // printNonZeroElements(gs);
                Matrix4D hcal=multip4dm_sp(gs_sp,ocal);
                double hsum=0;
                // auto start1 = std::chrono::high_resolution_clock::now();
                for (int b=0; b<nucnum;++b)
                {
                    for (int a=0;a<=b;++a)
                    {
                        for (int d=0;d<nucnum;++d)
                        {
                            for (int c=0;c<=d;++c)
                            {
                                hsum+=hcal[a][b](c,d);
                            }
                        }
                    }
                }
                // auto end1 = std::chrono::high_resolution_clock::now();
                // auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);
                // std::cout << "Elapsed time for hamsum: " << duration1.count() << "μs\n";
                result(l,m)+=hsum*strength[t];
            }
            // std::cout<<"l="<<l<<std::endl;
            // std::cout<<"result="<<result<<std::endl;

            // 更新进度计数器
            size_t current = ++completed_pairs;

            // 主线程显示进度
            if (omp_get_thread_num() == 0) {
                if (current % update_interval == 0 || current == total_pairs) {
                    auto now = std::chrono::steady_clock::now();
                    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count();

                    std::cout << "\rProcessed: " << current << "/" << total_pairs
                              << " (" << std::fixed << std::setprecision(1)
                              << 100.0 * current / total_pairs << "%)"
                              << " | Time: " << format_time(elapsed);
                    std::cout.flush();
                }
            }

        }
    }
    // 完成后
    auto end_time = std::chrono::steady_clock::now();
    auto total_time = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();

    std::cout << "\n\nCalculation completed in " << format_time(total_time) << "!\n";
    std::cout << "Total pairs: " << total_pairs << "\n";
    std::cout << "Pairs per second: "
              << static_cast<double>(total_pairs) / total_time << "\n";
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile << "Elapsed time for ham1: " << duration.count() << "μs\n"<<std::endl;
    std::cout << "Elapsed time for ham1: " << duration.count() << "μs\n";
    return result;
}

Eigen::MatrixXd hamsige(const std::vector<basism>& bas,const std::vector<std::vector<Eigen::MatrixXd>>&ystrm,std::vector<double> sigenvec)
{
    int nucnum=nucleusm.size();
    int sizebas=bas.size();
    Eigen::MatrixXd result=Eigen::MatrixXd::Zero(sizebas, sizebas);

    #pragma omp parallel for
    for (int l=0;l<sizebas;++l)
    {
        basism bas1=bas[l];
        std::vector<Eigen::MatrixXd> ystrm1=ystrm[l];
        for (int m=0;m<sizebas;++m)
        {
            basism bas2=bas[m];
            std::vector<Eigen::MatrixXd> ystrm2=ystrm[m];
            Eigen::MatrixXd qcal=calqsige(bas1,bas2,sigenvec);
            Eigen::MatrixXd gs=calf(bas1,bas2,ystrm1,ystrm2);
            Eigen::MatrixXd hcal=gs.cwiseProduct(qcal);
            for (int a=0; a<nucnum;++a)
            {
                result(l,m)+=hcal(a,a);
            }
        }

    }
    return result;

}

Eigen::MatrixXd hamchange(const std::vector<basism>& basm,const std::vector<basis>&bas,
    Eigen::MatrixXd transjm,int t,int mu,Eigen::MatrixXd ham1)
{
    int sizebas=bas.size();
    int nucnum=nucleusm.size();
    int sizebasm=basm.size();
    Eigen::MatrixXd result=Eigen::MatrixXd::Zero(sizebas, sizebas);
    Eigen::MatrixXd transjmt=transjm.transpose();
    Eigen::MatrixXd t1=transjmt*ham1*transjm;
    int m11=0;
    if (basm[0].r[0]!=0)
    {
        m11=1;
    }
    for (int l=0;l<sizebas;++l)
    {
        int j1=bas[l].sj.back();
        for (int m=0;m<sizebas;++m)
        {
            int j2=bas[m].sj.back();
            double cj=util::CG(j1,t,j2,m11,mu,m11);
            double rem=0;
            if (cj!=0)
            {
                rem=t1(l,m)/cj;
            }
            result(l,m)=rem;  // 计算结果
        }
    }
    return result;
}

Eigen::MatrixXd hamchange_1(const std::vector<basism>& basm,const std::vector<basis>&bas,
    Eigen::MatrixXd transjm,int t,int mu,Eigen::MatrixXd ham1,const std::vector<basism>& basm_1,Eigen::MatrixXd transjm_1)
{
    int sizebas=bas.size();
    int nucnum=nucleusm.size();
    int sizebasm=basm.size();
    Eigen::MatrixXd result=Eigen::MatrixXd::Zero(sizebas, sizebas);
    Eigen::MatrixXd transjmt=transjm_1.transpose();

    Eigen::MatrixXd t1_1=transjmt*ham1;
    Eigen::MatrixXd t1=t1_1*transjm;
    // std::cout<<"transjmt="<<std::endl;
    // std::cout<<transjmt<<std::endl;
    // std::cout<<"transjm_1="<<std::endl;
    // std::cout<<transjm_1<<std::endl;
    // std::cout<<"t1_1mat="<<std::endl;
    // std::cout<<t1_1<<std::endl;
    // std::cout<<"t1mat="<<std::endl;
    // std::cout<<t1<<std::endl;



    for (int l=0;l<sizebas;++l)
    {
        int j1=bas[l].sj.back();

        for (int m=0;m<sizebas;++m)
        {

            int j2=bas[m].sj.back();
            double cj=util::CG(j1,t,j2,2,-2,0);
            double rem=0;
            if (cj!=0)
            {
                rem=t1(l,m)/cj;
            }
            result(l,m)=rem;  // 计算结果
        }
    }
    return result;
}


double caloverlapm(basism bas1,basism bas2,std::vector<Eigen::MatrixXd> ystrm1,std::vector<Eigen::MatrixXd> ystrm2)
{
    double result=0;
    Eigen::MatrixXd pnp=ystrm2.back();
    std::vector<Eigen::MatrixXd> ystrm22=ystrm2;
    ystrm22.pop_back();
    basism bas22=bas2;
    bas22.r.pop_back();
    Eigen::MatrixXd bmat=calb(bas1,bas22,ystrm1,ystrm22);
    Eigen::MatrixXd pbmat =pnp*bmat;
    result=-2.0*pbmat.trace();
    // std::cout<<"bmat="<<result<<std::endl;
    // std::cout<<"pbmat="<<pbmat<<std::endl;
    // std::cout<<"result="<<result<<std::endl;
    return result;
}



Eigen::MatrixXd caloverlapmmat(const std::vector<basism>& bas1,const std::vector<basism>& bas2,
    const std::vector<std::vector<Eigen::MatrixXd>>& ystrm1,const std::vector<std::vector<Eigen::MatrixXd>>& ystrm2)
{
    int nucnum=bas1.size();
    std::cout<< "num=" << nucnum <<std::endl;
    // 进度跟踪
    std::atomic<int> completed(0);
    const int total_tasks = nucnum + (nucnum * (nucnum - 1)) / 2; // 对角 + 上三角

    auto update_progress = [&]() {
        int current = ++completed;
        #pragma omp critical
        {
            float percent = 100.f * current / total_tasks;
            int bars = static_cast<int>(percent / 2);
            std::cout << "\r["
                      << std::string(bars, '=') << ">"
                      << std::string(50 - bars, ' ') << "] "
                      << std::fixed << std::setprecision(1) << percent << "%"
                      << " (" << current << "/" << total_tasks << ")"
                      << std::flush;
        }
    };
    Eigen::MatrixXd result=Eigen::MatrixXd::Zero(nucnum, nucnum);
    // 预计算对角线并标记零行
    std::vector<bool> zero_row(nucnum, false);

    // 第一阶段：并行计算对角线
    #pragma omp parallel for
    for (int a = 0; a < nucnum; ++a) {
        basism bas11 = bas1[a];
        std::vector<Eigen::MatrixXd> ystrm11 = ystrm1[a];
        double diag_val = caloverlapm(bas11, bas11, ystrm11, ystrm11);

        if (std::abs(diag_val) < 1e-12) { // 使用适当的小量阈值
            zero_row[a] = true;
        }
        result(a,a) = diag_val;
        update_progress();
    }

    // 第二阶段：并行计算非零行的非对角元素
    #pragma omp parallel for schedule(dynamic)  // dynamic适用于负载不均衡的情况
    for (int a = 0; a < nucnum; ++a) {
        if (zero_row[a]) {
            // 整行置零（对角线已在第一阶段设置）
            for (int b = 0; b < nucnum; ++b) {
                if (b != a) result(a,b) = 0.0;
            }
            continue;
        }

        basism bas11 = bas1[a];
        std::vector<Eigen::MatrixXd> ystrm11 = ystrm1[a];

        // 只计算上三角部分（b > a），利用对称性
        for (int b = a + 1; b < nucnum; ++b) {
            if (zero_row[b]) {
                result(a,b) = 0.0;
                result(b,a) = 0.0;
                continue;
            }

            basism bas22 = bas2[b];
            std::vector<Eigen::MatrixXd> ystrm22 = ystrm2[b];

            if (bas11.parity != bas22.parity) {
                result(a,b) = 0.0;
                result(b,a) = 0.0;
            } else {
                double val = caloverlapm(bas11, bas22, ystrm11, ystrm22);
                result(a,b) = val;
                result(b,a) = val;  // 对称赋值
            }
            update_progress();
        }
    }
    return result;
}

Eigen::MatrixXd caloverlapmmatdif(const std::vector<basism>& bas1,const std::vector<basism>& bas2,
    const std::vector<std::vector<Eigen::MatrixXd>>& ystrm1,const std::vector<std::vector<Eigen::MatrixXd>>& ystrm2)
{
    int nucnum=bas1.size();
    int nucnum2=bas2.size();

    Eigen::MatrixXd result=Eigen::MatrixXd::Zero(nucnum, nucnum);



    // 第二阶段：并行计算非零行的非对角元素
    #pragma omp parallel for schedule(dynamic)  // dynamic适用于负载不均衡的情况
    for (int a = 0; a < nucnum; ++a)
    {

        basism bas11 = bas1[a];
        std::vector<Eigen::MatrixXd> ystrm11 = ystrm1[a];

        // 只计算上三角部分（b > a），利用对称性
        for (int b = 0; b < nucnum2; ++b)
        {
            basism bas22 = bas2[b];
            std::vector<Eigen::MatrixXd> ystrm22 = ystrm2[b];

            if (bas11.parity != bas22.parity) {
                result(a,b) = 0.0;
            } else {
                double val = caloverlapm(bas11, bas22, ystrm11, ystrm22);
                result(a,b) = val;
            }
        }
    }
    return result;
}

std::vector<std::vector<double>> eigenToNestedVector(const Eigen::MatrixXd& mat) {
    std::vector<std::vector<double>> result(mat.rows());
    for (int i = 0; i < mat.rows(); ++i) {
        result[i].resize(mat.cols());
        for (int j = 0; j < mat.cols(); ++j) {
            result[i][j] = mat(i, j); // 保证 n[i][j] == mat(i,j)
        }
    }
    return result;
}
Eigen::MatrixXd overlapdifchange(const std::vector<basism>& basm1,const std::vector<basis>&bas1,
    const std::vector<basism>& basm2,const std::vector<basis>&bas2,
    Eigen::MatrixXd transfermat1,Eigen::MatrixXd transfer2mat,Eigen::MatrixXd sfmat,int t)
{
    Eigen::MatrixXd transfer1t=transfermat1.transpose();
    Eigen::MatrixXd result1=transfer1t*sfmat*transfer2mat;
    return result1;
}

std::vector<std::vector<double> > overlapchange(Eigen::MatrixXd overm,Eigen::MatrixXd transjm)
{
    Eigen::MatrixXd transjmt=transjm.transpose();
    Eigen::MatrixXd t1=transjmt*overm*transjm;
    std::vector<std::vector<double>> result=eigenToNestedVector(t1);
    return result;
}

double calsf(const basism& bas1,const basism& bas2,std::vector<Eigen::MatrixXd> ystrm1,const std::vector<Eigen::MatrixXd>& ystrm2,int nuor)
{
    int nucnum=nucleusm.size();
    Eigen::VectorXd sp = Eigen::VectorXd::Zero(nucnum);
    sp(nuor) = 1.0;
    Eigen::VectorXd s=ystrm1[0];
    double result=0;
    int rsize=bas1.r.size();

    for (int l=0;l<rsize;++l)
    {
        basism bas11=bas1;
        std::vector<Eigen::MatrixXd> ystrm11=ystrm1;
        ystrm11[0]=Eigen::MatrixXd::Identity(nucnum, nucnum);
        bas11.r[0]=0;
        if (l==0)
        {
            double spa=sp.transpose()*s;
            ystrm11[1].array() *= spa;

        }else
        {
            Eigen::MatrixXd p1=ystrm11[l];
            Eigen::MatrixXd termp1=-((s*sp.transpose()*p1)+(p1*sp*s.transpose()));

            ystrm11[l]=termp1;
        }
        result+=caloverlapm(bas11,bas2,ystrm11,ystrm2);

    }
    return result;

}



Eigen::MatrixXd calsfmat(std::vector<basism> bas1,std::vector<basism> bas2,
    std::vector<std::vector<Eigen::MatrixXd>> ystrm1,std::vector<std::vector<Eigen::MatrixXd>> ystrm2,int nuor)
{
    int bas1num=bas1.size();
    int bas2num=bas2.size();
    int nucnum=nucleusm.size();
    Eigen::MatrixXd result=Eigen::MatrixXd::Zero(bas1num, bas2num);
    for (int l=0;l<bas1num;++l)
    {
        for (int m=0;m<bas2num;++m)
        {
            result(l,m)=calsf(bas1[l],bas2[m],ystrm1[l],ystrm2[m],nuor);
        }
    }

    return result;
}

Eigen::MatrixXd sfmatchange(const std::vector<basism>& basm1,const std::vector<basis>&bas1,
    const std::vector<basism>& basm2,const std::vector<basis>&bas2,
    Eigen::MatrixXd transfermat1,Eigen::MatrixXd transfer2mat,Eigen::MatrixXd sfmat,int nuor)
{
    double sigf=std::pow(-1, -nucleusm[nuor].j+1);
    Eigen::MatrixXd transfer2t=transfermat1.transpose();
    std::cout << "transfer2t: " << transfer2t.rows() << "x" << transfer2t.cols() << std::endl;
    std::cout << "sfmat: " << sfmat.rows() << "x" << sfmat.cols() << std::endl;
    std::cout << "transfermat1: " << transfer2mat.rows() << "x" << transfermat1.cols() << std::endl;
    Eigen::MatrixXd result1=transfer2t*sfmat*transfer2mat;
    int sizebas1=bas1.size();
    int sizebas2=bas2.size();
    int sizebasm1=basm1.size();
    int sizebasm2=basm2.size();
    int t=nucleusm[nuor].j;
    int mu=nucleusm[nuor].m;
    int m11=1;
    int m12=0;
    Eigen::MatrixXd result=Eigen::MatrixXd::Zero(sizebas1, sizebas2);
    for (int l=0;l<sizebas1;++l)
    {
        int j1=bas1[l].sj.back();
        for (int m=0;m<sizebas2;++m)
        {
            int j2=bas2[m].sj.back();
            double cj=util::CG(j2,t,j1,m12,-mu,m11);
            double rem=0;
            if (cj!=0)
            {
                rem=result1(l,m)/cj* sigf;
            }
            result(l,m)=rem;  // 计算结果
        }
    }

    return result;
}

Eigen::MatrixXd convertToEigenMatrix(const std::vector<std::vector<double>>& data) {
    // 检查输入是否为空
    if (data.empty() || data[0].empty()) {
        return Eigen::MatrixXd(0, 0);
    }

    size_t rows = data.size();
    size_t cols = data[0].size();

    // 检查所有行是否具有相同的列数
    for (const auto& row : data) {
        if (row.size() != cols) {
            throw std::invalid_argument("Input is not a rectangular matrix");
        }
    }

    // 创建Eigen矩阵并拷贝数据
    Eigen::MatrixXd mat(rows, cols);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            mat(i, j) = data[i][j];
        }
    }
    return mat;
}



std::vector<size_t> findMminus1PositionsSTL(const std::vector<Nucleusm>& nuclei,int i) {
    std::vector<size_t> positions;
    auto it = nuclei.begin();

    while ((it = std::find_if(it, nuclei.end(),
            [i](const Nucleusm& n) { return n.m == i; })) != nuclei.end()) {
        positions.push_back(std::distance(nuclei.begin(), it));
        ++it;
            }
    return positions;
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

struct NonZeroPosition {
    size_t matrix_index;  // 在vector中的索引
    int row;              // 行号
    int col;              // 列号
    double value;         // 元素值（可选）
};

using Positions = std::vector<NonZeroPosition>;

Positions findNonZeroPositions(const std::vector<Eigen::MatrixXd>& V_it) {
    Positions positions;

    for (size_t i = 0; i < V_it.size(); ++i) {
        const auto& matrix = V_it[i];

        for (int r = 0; r < matrix.rows(); ++r) {
            for (int c = 0; c < matrix.cols(); ++c) {
                if (matrix(r, c) != 0.0) {  // 判断是否非零
                    positions.push_back({i, r, c, matrix(r, c)});
                }
            }
        }
    }

    return positions;
}

std::vector<Eigen::MatrixXd> hamqcal(const std::vector<int>& qnpt,const std::vector<Eigen::MatrixXd>& q_nz_m)
{
    std::vector<Eigen::MatrixXd> hamqcalre={};

    for (int i=0;i<qnpt.size();++i)
    {
        int qn=qnpt[i];
        Eigen::MatrixXd hamqcalre1=Eigen::MatrixXd::Zero(qn, qn);
        for (int j=0;j<q_nz_m.size();++j)
        {
            Eigen::MatrixXd q_nz_m1=q_nz_m[j];
            hamqcalre1+=q_nz_m1;
        }
        hamqcalre.push_back(hamqcalre1);
    }
}
std::vector<Eigen::MatrixXd> hamqcal(const std::vector<int>& qnpt,const std::vector<Eigen::MatrixXd>& q_nz_m,
    const std::vector<basism>& bas,const std::vector<std::vector<Eigen::MatrixXd>>&ystrm)
{
    int nucnum=nucleusm.size();
    int sizebas=bas.size();
    std::vector<Eigen::MatrixXd> hamqcalre={};
    Eigen::MatrixXd term1=Eigen::MatrixXd::Zero(sizebas, sizebas);
    for (int i=0;i<qnpt.size();++i)
    {
        Eigen::MatrixXd hamqcalre1=Eigen::MatrixXd::Zero(sizebas, sizebas);
        hamqcalre.push_back(hamqcalre1);
    }
    for (int l=0;l<sizebas;++l)
    {
        basism bas1=bas[l];
        std::vector<Eigen::MatrixXd> ystrm1=ystrm[l];
        for (int m=0;m<sizebas;++m)
        {
            basism bas2=bas[m];
            std::vector<Eigen::MatrixXd> ystrm2=ystrm[m];
            Eigen::MatrixXd gs=calf(bas1,bas2,ystrm1,ystrm2);

            for (int i=0;i<qnpt.size();++i)
            {
                int qnt=qnpt[i];
                Eigen::MatrixXd q_nz_m1=q_nz_m[qnt];
                Eigen::MatrixXd hamqcalre1=gs.cwiseProduct(q_nz_m1);
                double relm=hamqcalre1.sum();
                hamqcalre[i](l,m)=relm;
            }
        }

    }
    return hamqcalre;

}

std::vector<std::vector<double>> readstrfile(const std::string& filename) {
    std::vector<std::vector<double>> data;

    // 打开文件
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("无法打开文件: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<double> row;
        double value;

        // 读取每行的四个数值
        while (iss >> value) {
            row.push_back(value);
        }

        // 确保每行有四个数值
        if (row.size() != 4) {
            throw std::runtime_error("文件格式错误：每行应有4个数值");
        }

        // 转换为0-based索引（前三个数值减1）
        row[0] -= 1;
        row[1] -= 1;
        row[2] -= 1;

        data.push_back(row);
    }

    file.close();
    return data;
}

// 打印单个 basis 到一行
void printBasisOneLine(const basis& b)
{
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    // 一行输出，字段之间用空格区分
    outfile
        << "j=" << b.j
        << " jn=" << b.jn
        << " r=" << vecToStr(b.r)
        << " rn=" << vecToStr(b.rn)
        << " sj=" << vecToStr(b.sj)
        << " rparity=" << vecToStr(b.rparity)
        << "parity=" << b.parity
        << std::endl;
}

// 打印 std::vector<basis>，每个 basis 一行
void printBasisVectorOneLine(const std::vector<basis>& basisVec)
{
    for (size_t i = 0; i < basisVec.size(); ++i) {
        printBasisOneLine(basisVec[i]);
    }
}

// 将3D张量的某个维度切片转换为列优先矩阵
Eigen::MatrixXd tensorSliceToMatrix(const Eigen::Tensor<double, 3>& tensor, int dim, int index) {
    // 获取切片后的张量维度
    Eigen::array<Eigen::Index, 2> dims;
    if (dim == 0) {
        dims = {tensor.dimension(1), tensor.dimension(2)}; // 切行 → 列×深度
    } else if (dim == 1) {
        dims = {tensor.dimension(0), tensor.dimension(2)}; // 切列 → 行×深度
    } else {
        dims = {tensor.dimension(0), tensor.dimension(1)}; // 切深度 → 行×列
    }

    // 提取切片并手动填充到列优先矩阵
    Eigen::MatrixXd matrix(dims[0], dims[1]);
    for (int i = 0; i < dims[0]; ++i) {
        for (int j = 0; j < dims[1]; ++j) {
            if (dim == 0) {
                matrix(i, j) = tensor(index, i, j); // 行切片
            } else if (dim == 1) {
                matrix(i, j) = tensor(i, index, j); // 列切片
            } else {
                matrix(i, j) = tensor(i, j, index); // 深度切片
            }
        }
    }
    return matrix;
}

Eigen::VectorXd calsvec
(const basism& bas1,const basism& bas2,
    std::vector<Eigen::MatrixXd> ystrm1,std::vector<Eigen::MatrixXd> ystrm2)
{
    int nucnum=nucleusm.size();
    Eigen::VectorXd result=Eigen::VectorXd::Zero(nucnum);
    Eigen::MatrixXd pnp=ystrm2.back();
    std::vector<Eigen::MatrixXd> ystrm22=ystrm2;
    ystrm22.pop_back();
    basism bas22=bas2;
    bas22.r.pop_back();
    Eigen::Tensor<double, 3> t1=calt(bas1,bas22,ystrm1,ystrm22);
    Eigen::Tensor<double, 3> t2=tchange(t1);
    // writeTensor3ToFile(t2);

    for (int i=0;i<nucnum;++i)
    {
        Eigen::MatrixXd stcal=tensorSliceToMatrix(t2,0,i);
        Eigen::MatrixXd mat1=pnp*stcal;

        double re1=-6* mat1.trace();
        result(i)=re1;
    }
    // std::cout<<"calsvec="<<result<<std::endl;
    return result;
}

double sigsfcal(const basism& bas1,const basism& bas2,
    std::vector<Eigen::MatrixXd> ystrm1,std::vector<Eigen::MatrixXd> ystrm2,int trannum)
{
    int sizen=ystrm2.size();
    if (bas1.r[0]==0)
    {
        double result=0;
        double term1=0;
        double term2=0;
        Eigen::VectorXd sp= ystrm2[0];
        for (int k=1;k<sizen;++k)
        {
            Eigen::MatrixXd pkp=ystrm2[k];
            std::vector<Eigen::MatrixXd> ystrm22=ystrm2;
            ystrm22[0]=Eigen::MatrixXd::Identity(nucleusm.size(), nucleusm.size());
            ystrm22.erase(ystrm22.begin()+k);
            basism bas22=bas2;
            bas22.r.erase(bas22.r.begin()+k);
            bas22.r[0]=0;
            Eigen::MatrixXd bmat=calb(bas1,bas22,ystrm1,ystrm22);

            // std::cout<<"pkp="<<pkp<<std::endl;
            // std::cout<<"bmat="<<bmat<<std::endl;
            // std::cout<<"sp="<<sp<<std::endl;
            // std::cout<<"re1"<<pkp*bmat<<std::endl;
            Eigen::MatrixXd term11=pkp*bmat*sp;
            term1+=term11(trannum);
            if (sp(trannum)!=0 && k==sizen-1)
            {
                Eigen::MatrixXd term12=pkp*bmat;
                term2+=term12.trace()*sp(trannum);
            }

        }
        result=4*term1-(2*term2);
        return result;
    }
    else
    {
        int nucnum=nucleusm.size();
        // std::cout<<"begincalsigsfcal"<<std::endl;
        double result=0;
        Eigen::VectorXd sp= ystrm1[0];
        for (int k=1;k<sizen;++k)
        {
            Eigen::MatrixXd pkp=ystrm2[k];
            for (int l=0;l<nucleusm.size();l++)
            {
                Eigen::VectorXd snew = Eigen::VectorXd::Zero(nucnum);
                snew(l) = 1.0;
                std::vector<Eigen::MatrixXd> ystrm22=ystrm2;
                ystrm22[0]=snew;
                ystrm22.erase(ystrm22.begin()+k);
                basism bas22=bas2;
                bas22.r.erase(bas22.r.begin()+k);
                bas22.r[0]=nucleusm[l].j;
                double overe=caloverlapm(bas1,bas22,ystrm1,ystrm22);
                double tempm=pkp(trannum,l)*overe;
                result+=tempm;
            }

            // std::vector<Eigen::MatrixXd> ystrm22=ystrm2;
            // ystrm22.erase(ystrm22.begin()+k);
            // basism bas22=bas2;
            // bas22.r.erase(bas22.r.begin()+k);
            // Eigen::VectorXd svec=calsvec(bas1,bas22,ystrm1,ystrm22);
            // Eigen::MatrixXd term11=pkp*svec;
            // result+=term11(trannum);
        }
        return 2*result;
    }
}

Eigen::MatrixXd sigsfmat(std::vector<basism> bas1,std::vector<basism> bas2,
    std::vector<std::vector<Eigen::MatrixXd>> ystrm1,std::vector<std::vector<Eigen::MatrixXd>> ystrm2,int trannum)
{
    int bas1num=bas1.size();
    int bas2num=bas2.size();
    Eigen::MatrixXd result=Eigen::MatrixXd::Zero(bas1num, bas2num);
    for (int l=0;l<bas1num;++l)
    {
        for (int m=0;m<bas2num;++m)
        {
            result(l,m)=sigsfcal(bas1[l],bas2[m],ystrm1[l],ystrm2[m],trannum);
        }
    }
    std::cout<<"sigsfcal"<<std::endl;

    std::cout<<result<<std::endl;

    return result;
}

Eigen::MatrixXd sigsfmatchange(const std::vector<basism>& basm1,const std::vector<basis>&bas1,
    const std::vector<basism>& basm2,const std::vector<basis>&bas2,
    Eigen::MatrixXd transfermat1,Eigen::MatrixXd transfer2mat,Eigen::MatrixXd sfmat,int trannum)
{
    Eigen::MatrixXd transfer1t=transfermat1.transpose();
    Eigen::MatrixXd result1=transfer1t*sfmat*transfer2mat;
    int sizebas1=bas1.size();
    int sizebas2=bas2.size();
    int sizebasm1=basm1.size();
    int sizebasm2=basm2.size();
    int t=nucleusm[trannum].j;
    int mu=nucleusm[trannum].m;
    int m11=basm1[0].bigmvec.back();
    int m12= basm2[0].bigmvec.back();
    Eigen::MatrixXd result=Eigen::MatrixXd::Zero(sizebas1, sizebas2);
#pragma omp parallel for
    for (int l=0;l<sizebas1;++l)
    {
        int j1=bas1[l].sj.back();
        for (int m=0;m<sizebas2;++m)
        {
            int j2=bas2[m].sj.back();
            double cj=util::CG(j1,t,j2,m11,mu,m12);
            double rem=0;
            if (cj!=0)
            {
                rem=result1(l,m)/cj;
            }
            result(l,m)=rem;  // 计算结果
        }
    }

    return result;
}

std::vector<Eigen::MatrixXd>calsfmatall(
    std::vector<basism> allbasism1,
    std::vector<std::vector<Eigen::MatrixXd>> ystrm1,
    Eigen::MatrixXd schmitmat1,
    std::vector<basis> allbasisp1,
    Eigen::MatrixXd changemat1,
    std::vector<basism> allbasism2,
    std::vector<std::vector<Eigen::MatrixXd>>ystrm2,
    Eigen::MatrixXd schmitmat2,
    std::vector<basis> allbasisp2,
    Eigen::MatrixXd changemat2)
{

    auto nucmmin1=findMminus1PositionsSTL(nucleusm,-1);
    std::vector<Eigen::MatrixXd> sfmatall;
    for (int i=0;i<nucmmin1.size();++i)
    {
        Eigen::MatrixXd sfmat=sigsfmat(allbasism1, allbasism2, ystrm1, ystrm2, nucmmin1[i]);
        Eigen::MatrixXd sfchange=sigsfmatchange(allbasism1,allbasisp1, allbasism2,allbasisp2,
            changemat1, changemat2, sfmat,nucmmin1[i]);
        Eigen::MatrixXd sfmatchsch= schmitmat1*sfchange* schmitmat2.transpose();
        if (!outfile.is_open())
        {
            // 如果文件未打开，则尝试打开文件
            outfile.open("basis_output.txt",std::ios::app);
        }
        outfile<< "sfmatchsch for nucleus " << i << ":\n";
        outfile<<"sfmat"<<"\n"
        <<sfmat << "\n"
        << "sfchange"<<"\n"
        << sfchange << "\n"
        << "sfmatchsch"<<"\n"
        << sfmatchsch << "\n"
        << "------------------------" << std::endl;
        sfmatall.push_back(sfmatchsch);
    }
    return sfmatall;
}

std::vector<Eigen::MatrixXd> getsfmatresult(
    std::map<int, std::vector<CoupledBasis>> couplebasis1,
    std::map<int,std::vector<CoupledBasis>> couplebasis2,
    std::map<int, std::vector<std::vector<double>>> eigenre1,
    std::map<int, std::vector<std::vector<double>>> eigenre2,
    std::vector<int> jvecp1,
    std::vector<int> jvecp2,
    std::vector<int> jvecn1,
    std::vector<int> jvecn2,
    std::vector<int> tvec,
    std::vector<Eigen::MatrixXd> transmat,
    int norp,
    std::vector<basis> pbas1,
    std::vector<basism>pbasm01,
    std::vector<basism>pbasm11,
    std::vector<std::vector<Eigen::MatrixXd>> pystr1,
    std::vector<std::vector<Eigen::MatrixXd>> pystr11,
    Eigen::MatrixXd pchangmat1,
    Eigen::MatrixXd pchangmat11,
    Eigen::MatrixXd pschammit1,
    std::vector<basis> pbas2,
    std::vector<basism>pbasm02,
    std::vector<basism>pbasm12,
    std::vector<std::vector<Eigen::MatrixXd>> pystr2,
    std::vector<std::vector<Eigen::MatrixXd>> pystr21,
    Eigen::MatrixXd pchangmat2,
    Eigen::MatrixXd pchangmat21,
    Eigen::MatrixXd pschammit2
    )
{
    Eigen::MatrixXd overlapp=caloverlapmmatdif(pbasm01,pbasm02,pystr1,pystr2);
    Eigen::MatrixXd overlappch=overlapdifchange(pbasm01,pbas1,pbasm02,pbas2,
        pchangmat1,pchangmat2,overlapp,0);
    overlappch=pschammit1*overlappch* pschammit2.transpose();
    int sizedft=transmat.size();
    int totalStates1 = 0;
    for (auto &[key, vec] : eigenre1) {
        totalStates1 += vec.size();  // 每个key的本征态数
    }
    int totalStates2 = 0;
    for (auto &[key, vec] : eigenre2) {
        totalStates2 += vec.size();
    }
    std::vector<Eigen::MatrixXd> resultmat;
    resultmat.resize(sizedft);
    for (int tmatnum = 0; tmatnum < sizedft; ++tmatnum) {
        resultmat[tmatnum] = Eigen::MatrixXd::Zero(totalStates1, totalStates2);
    }
    int offset1 = 0;
    if (norp==2)
    {

        for (const auto &[key, value]: eigenre1)
        {
            int jisum = std::abs(key);
            int siz1 = couplebasis1[key].size();
            int vitnum = 0;
            int sizeeigen1 = eigenre1[key].size();
            for (int u = 0; u < sizeeigen1; u++)
            {
                int offset2 = 0;
                int keynum2 = 0;
                for (const auto &[key2, value2]: eigenre2)
                {
                    int jfsum = std::abs(key2);
                    int siz2 = couplebasis2[key2].size();
                    std::vector<std::vector<std::vector<double> > > transmatch
                    (sizedft, std::vector<std::vector<double> >
                     (siz1, std::vector<double>(siz2, 0.0)));
                    for (int i = 0; i < siz1; i++)
                    {
                        for (int j = 0; j < siz2; j++)
                        {
                            int n1 = couplebasis1[key][i].ni;
                            int n2 = couplebasis2[key2][j].ni;
                            int p1 = couplebasis1[key][i].pi;
                            int p2 = couplebasis2[key2][j].pi;
                            int jn = jvecn1[n1];
                            int jnp = jvecn2[n2];
                            int jp = jvecp1[p1];
                            int jpp = jvecp2[p2];


                            int deltann = deltatwo(jn, jnp);
                            int deltapp = deltatwo(jp, jpp);
                            for (int tmatnum = 0; tmatnum < sizedft; ++tmatnum)
                            {
                                int t = tvec[tmatnum];
                                Eigen::MatrixXd transmat1 = transmat[tmatnum];
                                double calnu = deltapp * overlappch(p1, p2) * std::pow(-1, (jn - jnp + jfsum - jisum) / 2)
                                               * ufrom6_j(jpp, jn, jfsum, t, jisum, jnp) * transmat1(n1, n2);
                                transmatch[tmatnum][i][j] = calnu;
                            }
                        }
                    }
                    if (key==0)
                    {
                        if (key2==4)
                        {
                            std::cout<<"_______________________"<<std::endl;
                            printMatrix(transmatch[1]);
                        }
                    }


                    int sizeeigen2 = eigenre2[key2].size();

                    for (int v = 0; v < sizeeigen2; v++)
                    {
                        std::vector<double> eigencal1 = eigenre1[key][u];
                        std::vector<double> eigencal2 = eigenre2[key2][v];

                        for (int tmatnum = 0; tmatnum < sizedft; ++tmatnum)
                        {
                            int t = tvec[tmatnum] / 2;
                            int fact = 1;
                            double sum = 0.0;
                            Eigen::MatrixXd tranmatcal(eigencal1.size(), eigencal2.size());

                            for (int i = 0; i < eigencal1.size(); ++i)
                            {
                                for (int j = 0; j < eigencal2.size(); ++j)
                                {
                                    double cal1 = eigencal1[i];
                                    double cal2 = eigencal2[j];
                                    sum += cal2 * transmatch[tmatnum][i][j] * cal1;
                                }
                            }
                            sum = std::pow(sum, 2.0) ;
                            resultmat[tmatnum](offset1 + u, offset2 + v) = sum;
                        }
                    }
                    offset2 += sizeeigen2;
                }
            }
            offset1 += sizeeigen1;
        }
    }
    return resultmat;
}


Eigen::MatrixXd doublesfmatchange(const std::vector<basism>& basm1,const std::vector<basis>&bas1,
    const std::vector<basism>& basm2,const std::vector<basis>&bas2,
    Eigen::MatrixXd transfermat1,Eigen::MatrixXd transfer2mat,Eigen::MatrixXd sfmat,int t)
{
    Eigen::MatrixXd transfer1t=transfermat1.transpose();
    Eigen::MatrixXd result1=transfer1t*sfmat*transfer2mat;
    int sizebas1=bas1.size();
    int sizebas2=bas2.size();
    int sizebasm1=basm1.size();
    int sizebasm2=basm2.size();
    int m11=basm1[0].bigmvec.back();
    int m12= basm2[0].bigmvec.back();
    int mu=m12-m11;
    Eigen::MatrixXd result=Eigen::MatrixXd::Zero(sizebas1, sizebas2);
// #pragma omp parallel for
    for (int l=0;l<sizebas1;++l)
    {
        int j1=bas1[l].sj.back();
        for (int m=0;m<sizebas2;++m)
        {
            int j2=bas2[m].sj.back();
            double cj=util::CG(j1,t,j2,m11,mu,m12);
            double rem=0;
            if (cj!=0)
            {
                rem=result1(l,m)/cj;
            }
            result(l,m)=rem;  // 计算结果
        }
    }

    return result;
}


Eigen::MatrixXd doublesfmatele(const basism& bas1,const basism& bas2,
    std::vector<Eigen::MatrixXd> ystrm1,std::vector<Eigen::MatrixXd> ystrm2)
{
    std::cout << std::fixed << std::setprecision(4);
    int sizen=ystrm2.size();
    int nunum=nucleusm.size();
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(nunum,nunum);
    Eigen::MatrixXd term1=Eigen::MatrixXd::Zero(nunum,nunum);
    Eigen::MatrixXd term2 = Eigen::MatrixXd::Zero(nunum,nunum);
    for (int k = 1; k < sizen; ++k)
    {
        Eigen::MatrixXd pkp12 = ystrm2[k];
        std::vector<Eigen::MatrixXd> ystrm22 = ystrm2;
        ystrm22.erase(ystrm22.begin() + k);
        basism bas22 = bas2;
        bas22.r.erase(bas22.r.begin() + k);
        Eigen::MatrixXd pnp = ystrm22.back();
        std::vector<Eigen::MatrixXd> ystrm222 = ystrm22;
        ystrm222.pop_back();
        basism bas222 = bas22;
        bas222.r.pop_back();
        Eigen::MatrixXd bmat = calb(bas1, bas222, ystrm1, ystrm222);
        Eigen::MatrixXd pbmat = pnp.transpose() * bmat;
        term1 += pkp12 * pbmat.trace();
        // std::cout<<"term1"<<std::endl;
        // std::cout<<term1<<std::endl;
    }
    for (int k = 2; k < sizen; ++k)
    {
        Eigen::MatrixXd pkp = ystrm2[k];
        for (int i = 1; i < k; ++i)
        {
            Eigen::MatrixXd pip = ystrm2[i];
            std::vector<Eigen::MatrixXd> ystrm22 = ystrm2;
            basism bas22 = bas2;
            bas22.r.erase(bas22.r.begin() + k);
            bas22.r.erase(bas22.r.begin() + i);
            ystrm22.erase(ystrm22.begin() + k);
            ystrm22.erase(ystrm22.begin() + i);
            Eigen::MatrixXd bmat = calb(bas1, bas22, ystrm1, ystrm22);
            Eigen::MatrixXd midterm1 = pkp * bmat * pip;
            Eigen::MatrixXd midterm2 = pip * bmat * pkp;
            term2+=midterm1+midterm2;
            // std::cout<<term2<<std::endl;
        }
    }
    result=4*term1+8*term2;
    if (bas1.r[0]!=0)
    {
        Eigen::MatrixXd term3 = Eigen::MatrixXd::Zero(nunum,nunum);
        Eigen::VectorXd s0=ystrm2[0];
        for (int k = 1; k < sizen; ++k)
        {
            Eigen::MatrixXd pkp = ystrm2[k];
            std::vector<Eigen::MatrixXd> ystrm22 = ystrm2;
            ystrm22[0]=Eigen::MatrixXd::Identity(nunum, nunum);
            ystrm22.erase(ystrm22.begin()+k);
            basism bas22 = bas2;
            bas22.r.erase(bas22.r.begin()+k);
            bas22.r[0]=0;
            Eigen::VectorXd ss=calsvec(bas1, bas22, ystrm1, ystrm22);
            Eigen::RowVectorXd s0t=s0.transpose();
            Eigen::RowVectorXd sst=ss.transpose();
            Eigen::MatrixXd pkpt=pkp.transpose();
            term3+=s0*sst*pkpt-pkp*ss*s0t;
        }
        result+=term3*2;
    }
    return result;
}

Eigen::MatrixXd doubleftransqjtoqm(Eigen::MatrixXd qj,int t, int mu)
{
    int nucnum=nucleusm.size();
    int nucnu=nucleus.size();
    Eigen::MatrixXd qm=Eigen::MatrixXd::Zero(nucnum,nucnum);
    for (int a = 0; a < nucnum; ++a)
    {
        int ja=nucleusm[a].j;
        int numa=nucleusm[a].num;
        int ma=nucleusm[a].m;
        for (int b = 0; b < nucnum; ++b)
        {
            int jb=nucleusm[b].j;
            int numb=nucleusm[b].num;
            int mb=nucleusm[b].m;
            if (ma+mb==mu)
            {
                double cj=util::CG(ja, jb, t, ma, mb, mu);
                qm(a,b)=cj*qj(numa,numb);
            }
        }
    }
    return qm;
}

double computeElementwiseProductSum1(const Eigen::MatrixXd& mat1,
                                     const Eigen::MatrixXd& mat2) {
    // 假设 mat1 和 mat2 的大小相同
    double sum = 0.0;

    int rows = mat1.rows();
    int cols = mat1.cols();
    // std::cout<<"mat1 = "<<"/n"<<mat1<<std::endl;
    // std::cout<<"mat2 = "<<"/n"<<mat2<<std::endl;
    // Eigen::MatrixXd matre(rows,cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            sum += mat1(i, j) * mat2(i, j);
            // matre(i,j) = mat1(i,j)* mat2(i,j);
            // if (sum!=0)
            // {
            //     std::cout << "i = " << i<<","<<"j"<<j << std::endl;
            // }
        }
    }
    // std::cout<<"matre = "<<"/n"<<matre<<std::endl;

    return sum;
}

std::vector<Eigen::MatrixXd> doublesfmat(const std::vector<basism>& bas1,
    const std::vector<std::vector<Eigen::MatrixXd>>& ystrm1,
    const std::vector<basism>& bas2,
    const std::vector<std::vector<Eigen::MatrixXd>>& ystrm2,
    std::vector<Eigen::MatrixXd>& qstrm,
    const std::vector<int> & tvec)
{
    auto start = std::chrono::high_resolution_clock::now();

    int nucnum=nucleusm.size();
    int sizebas1=bas1.size();
    int sizebas2=bas2.size();
    int sizeqstr=qstrm.size();
    std::vector<Eigen::MatrixXd> resultmatvec;
    resultmatvec.resize(sizeqstr,
        Eigen::MatrixXd::Zero(sizebas1, sizebas2));
    // #pragma omp parallel for
    for (int l = 0; l < sizebas1; ++l) {
        const basism& basl = bas1[l];  // 用引用避免拷贝
        int m11=basl.bigmvec.back();
        const auto& ystrml = ystrm1[l]; // 用 const 引用
        for (int m = 0; m < sizebas2; ++m) {
            const basism& basmm = bas2[m];
            const auto& ystrmm = ystrm2[m];
            Eigen::MatrixXd doublesfmatone = doublesfmatele(basl, basmm, ystrml, ystrmm);
            for (int n = 0; n < sizeqstr; ++n)
            {
                Eigen::MatrixXd qtimu=qstrm[n];
                resultmatvec[n](l,m) =
                    computeElementwiseProductSum1(doublesfmatone, qtimu);
            }
        }
    }
    return resultmatvec;
}

std::vector<std::vector<std::vector<double>>> transtr(
    std::vector<int>rnnum,
    const std::vector<std::vector<double> >& ystrgetp)
{
    basis basisr;
    basisr.rn=rnnum;
    std::vector<basis> basis1;
    int size1 = nucleusm.size();

    std::vector<std::vector<std::vector<std::vector<double>>>> ystrallp(
        1,  // 第一维：大小为 sizep
        std::vector<std::vector<std::vector<double>>>(
            rnnum.size(),  // 第二维：大小为 sizepr
            std::vector<std::vector<double>>(
                nucleus.size(),  // 第三维：大小为 size1
                std::vector<double>(nucleus.size(), 0.0)  // 第四维：大小为 size1，元素初始化为 0.0
            )
        )
    );
    basis1.push_back(basisr);

    calculateystr(ystrallp,ystrgetp,basis1);

    std::vector<std::vector<std::vector<double>>> ystr=ystrallp[0];
    return ystr;
}

void transtrnocouple(
    std::vector<int>rget,std::vector<Eigen::MatrixXd> & qstrm0,
    std::vector<Eigen::MatrixXd>& qstrm1,std::vector<std::vector<std::vector<double>>>& transtrvec, std::vector<int>& rvec)
{
    std::cout << std::fixed << std::setprecision(4);
    int nucusize=nucleus.size();
    int nucnumsize=nucleusm.size();
    for (int i=0; i < rget.size(); ++i)
    {

        for (int j = 0; j < nucusize; ++j)
        {
            for (int k = 0; k < nucusize; ++k)
            {
                double fat=1;
                if (j==k)
                {
                    fat=1/sqrt(2.0);
                }
                Eigen::MatrixXd ystrtemp=Eigen::MatrixXd::Zero(nucusize, nucusize);
                ystrtemp(j,k) =fat*1;
                // ystrtemp(k,j)=fat*1;
                Eigen::MatrixXd ystrmtemp0=doubleftransqjtoqm(ystrtemp,rget[i],0);
                Eigen::MatrixXd ystrmtemp1=doubleftransqjtoqm(ystrtemp,rget[i],-2);
                qstrm0.push_back(ystrmtemp0);
                qstrm1.push_back(ystrmtemp1);
                std::cout<<"i,j,k"<<std::endl;
                std::cout<<i<<","<<j<<","<<k<<std::endl;
                std::cout<<"qstrm0"<<std::endl;
                std::cout<<ystrmtemp0<<std::endl;
                std::cout<<"qstrm1"<<std::endl;
                std::cout<<ystrmtemp1<<std::endl;
                std::vector<std::vector<double>>ystrr= eigenToNestedVector(ystrtemp);
                transtrvec.push_back(ystrr);

                rvec.push_back(rget[i]);
            }
        }
    }
    return;
}

std::vector<Eigen::MatrixXd> caldoublemat(
    std::vector<basis> bas1,
    std::vector<basism>basm01,
    std::vector<basism>basm11,
    std::vector<std::vector<Eigen::MatrixXd>> ystr1,
    std::vector<std::vector<Eigen::MatrixXd>> ystr11,
    Eigen::MatrixXd changmat1,
    Eigen::MatrixXd changmat11,
    Eigen::MatrixXd schammit1,
    std::vector<basis> bas2,
    std::vector<basism>basm02,
    std::vector<basism>basm12,
    std::vector<std::vector<Eigen::MatrixXd>> ystr2,
    std::vector<std::vector<Eigen::MatrixXd>> ystr21,
    Eigen::MatrixXd changmat2,
    Eigen::MatrixXd changmat21,
    Eigen::MatrixXd schammit2,
    std::vector<int> tranvec,
    std::vector<std::vector<std::vector<double>>> transtr

    )
{

    int sizedft=tranvec.size();
    std::vector<Eigen::MatrixXd> transtrmat={};
    std::vector<Eigen::MatrixXd> transtrmatm0={};
    std::vector<Eigen::MatrixXd> transtrmatm1={};
    for (int l = 0; l < sizedft; ++l)
    {
        transtrmat.push_back(convertToEigenMatrix(transtr[l]));
        Eigen::MatrixXd qstrm0=doubleftransqjtoqm(
            transtrmat.back(),tranvec[l],0);
        transtrmatm0.push_back(qstrm0);
        Eigen::MatrixXd qstrm1=doubleftransqjtoqm(
            transtrmat.back(),tranvec[l],-2);
        transtrmatm1.push_back(qstrm1);
    }

    std::vector<Eigen::MatrixXd> matj0={};
    std::vector<Eigen::MatrixXd> matj1={};
    std::vector<Eigen::MatrixXd> matj0sch={};
    std::vector<Eigen::MatrixXd> matj1sch={};


    std::vector<Eigen::MatrixXd> mat0=doublesfmat(basm01,ystr1,
        basm02,ystr2,transtrmatm0,tranvec);
    for (int l = 0; l < sizedft; ++l)
    {
        matj0.push_back(doublesfmatchange(basm01,bas1,
            basm02,bas2,
            changmat1,changmat2,
            mat0[l],tranvec[l]));
    }


    std::vector<Eigen::MatrixXd> mat1=doublesfmat(basm11,ystr11,
        basm02,ystr2,transtrmatm1,tranvec);
    for (int l = 0; l < sizedft; ++l)
    {
        Eigen::MatrixXd tempre=doublesfmatchange(basm11,bas1,
            basm02,bas2,
            changmat11,changmat2,
            mat1[l],tranvec[l]);
        matj1.push_back(tempre);
    }
    for (int tnum=0;tnum<sizedft;++tnum)
    {
        int t=tranvec[tnum];
        if (t!=0)
        {
            for (int i=0;i<bas1.size();++i)
            {
                int j1=bas1[i].sj.back();
                for (int j=0;j<bas2.size();++j)
                {
                    int j2=bas2[j].sj.back();
                    int num=(j1+j2+t)/2;
                    if ((num & 1) != 0)
                    {
                        matj0[tnum](i,j)=matj1[tnum](i,j);
                    }
                    matj0[tnum](i,j)=matj0[tnum](i,j);
                }
            }
        }

        matj0sch.push_back(schammit1*matj0[tnum]* schammit2.transpose());
        std::cout<<"tnum"<<" "<<tnum<<std::endl;
        std::cout<<matj0[tnum]<<std::endl;
    }
    return matj0sch;
}


std::vector<Eigen::MatrixXd> getdoublematresult(
    std::map<int, std::vector<CoupledBasis>> couplebasis1,
    std::map<int,std::vector<CoupledBasis>> couplebasis2,
    std::map<int, std::vector<std::vector<double>>> eigenre1,
    std::map<int, std::vector<std::vector<double>>> eigenre2,
    std::vector<int> jvecp1,
    std::vector<int> jvecp2,
    std::vector<int> jvecn1,
    std::vector<int> jvecn2,
    std::vector<int> tvec,
    std::vector<Eigen::MatrixXd> transmat,
    int norp,
    std::vector<basis> pbas1,
    std::vector<basism>pbasm01,
    std::vector<basism>pbasm11,
    std::vector<std::vector<Eigen::MatrixXd>> pystr1,
    std::vector<std::vector<Eigen::MatrixXd>> pystr11,
    Eigen::MatrixXd pchangmat1,
    Eigen::MatrixXd pchangmat11,
    Eigen::MatrixXd pschammit1,
    std::vector<basis> pbas2,
    std::vector<basism>pbasm02,
    std::vector<basism>pbasm12,
    std::vector<std::vector<Eigen::MatrixXd>> pystr2,
    std::vector<std::vector<Eigen::MatrixXd>> pystr21,
    Eigen::MatrixXd pchangmat2,
    Eigen::MatrixXd pchangmat21,
    Eigen::MatrixXd pschammit2
    )
{
    Eigen::MatrixXd overlapp=caloverlapmmatdif(pbasm01,pbasm02,pystr1,pystr2);
    Eigen::MatrixXd overlappch=overlapdifchange(pbasm01,pbas1,pbasm02,pbas2,
        pchangmat1,pchangmat2,overlapp,0);
    overlappch=pschammit1*overlappch* pschammit2.transpose();
    int sizedft=transmat.size();
    int totalStates1 = 0;
    for (auto &[key, vec] : eigenre1) {
        totalStates1 += vec.size();  // 每个key的本征态数
    }
    int totalStates2 = 0;
    for (auto &[key, vec] : eigenre2) {
        totalStates2 += vec.size();
    }
    std::vector<Eigen::MatrixXd> resultmat;
    resultmat.resize(sizedft);
    for (int tmatnum = 0; tmatnum < sizedft; ++tmatnum) {
        resultmat[tmatnum] = Eigen::MatrixXd::Zero(totalStates1, totalStates2);
    }
    int offset1 = 0;
    if (norp==2)
    {

        for (const auto &[key, value]: eigenre1)
        {
            int jisum = std::abs(key);
            int siz1 = couplebasis1[key].size();
            int vitnum = 0;
            int sizeeigen1 = eigenre1[key].size();
            for (int u = 0; u < sizeeigen1; u++)
            {
                int offset2 = 0;
                int keynum2 = 0;
                for (const auto &[key2, value2]: eigenre2)
                {
                    int jfsum = std::abs(key2);
                    int siz2 = couplebasis2[key2].size();
                    std::vector<std::vector<std::vector<double> > > transmatch
                    (sizedft, std::vector<std::vector<double> >
                     (siz1, std::vector<double>(siz2, 0.0)));
                    for (int i = 0; i < siz1; i++)
                    {
                        for (int j = 0; j < siz2; j++)
                        {
                            int n1 = couplebasis1[key][i].ni;
                            int n2 = couplebasis2[key2][j].ni;
                            int p1 = couplebasis1[key][i].pi;
                            int p2 = couplebasis2[key2][j].pi;
                            int jn = jvecn1[n1];
                            int jnp = jvecn2[n2];
                            int jp = jvecp1[p1];
                            int jpp = jvecp2[p2];


                            int deltann = deltatwo(jn, jnp);
                            int deltapp = deltatwo(jp, jpp);
                            for (int tmatnum = 0; tmatnum < sizedft; ++tmatnum)
                            {
                                int t = tvec[tmatnum];
                                Eigen::MatrixXd transmat1 = transmat[tmatnum];
                                double calnu = deltapp * overlappch(p1, p2) * std::pow(-1, (jn - jnp + jfsum - jisum) / 2)
                                               * ufrom6_j(jpp, jn, jfsum, t, jisum, jnp) * transmat1(n1, n2);
                                transmatch[tmatnum][i][j] = calnu;
                            }
                        }
                    }
                    if (key==0)
                    {
                        if (key2==4)
                        {
                            std::cout<<"_______________________"<<std::endl;
                            printMatrix(transmatch[1]);
                        }
                    }


                    int sizeeigen2 = eigenre2[key2].size();

                    for (int v = 0; v < sizeeigen2; v++)
                    {
                        std::vector<double> eigencal1 = eigenre1[key][u];
                        std::vector<double> eigencal2 = eigenre2[key2][v];

                        for (int tmatnum = 0; tmatnum < sizedft; ++tmatnum)
                        {
                            int t = tvec[tmatnum] / 2;
                            int fact = 1;
                            double sum = 0.0;
                            Eigen::MatrixXd tranmatcal(eigencal1.size(), eigencal2.size());

                            for (int i = 0; i < eigencal1.size(); ++i)
                            {
                                for (int j = 0; j < eigencal2.size(); ++j)
                                {
                                    double cal1 = eigencal1[i];
                                    double cal2 = eigencal2[j];
                                    sum += cal2 * transmatch[tmatnum][i][j] * cal1;
                                }
                            }
                            // sum = std::pow(sum, 2.0) * (jfsum + 1.0) / (jisum + 1.0);
                            resultmat[tmatnum](offset1 + u, offset2 + v) = sum;
                        }
                    }
                    offset2 += sizeeigen2;
                }
            }
            offset1 += sizeeigen1;
        }
    }
    return resultmat;
}

void writeResultMatricesToFile(
    const std::string &filename,
    const std::vector<Eigen::MatrixXd> &resultmat,
    const std::map<int, std::vector<std::vector<double>>> &eigenre1,
    const std::map<int, std::vector<std::vector<double>>> &eigenre2,
    const std::vector<int> &tvec,
    const std::vector<std::vector<std::vector<double> > > transtr
) {
    std::ofstream fout(filename);
    if (!fout.is_open()) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        return;
    }

    fout << std::fixed << std::setprecision(6); // 控制小数位数

    // 建立 row index 对应的 (key,u)
    std::vector<std::pair<int,int>> rowStates; // (key,u)
    for (const auto &[key, vec] : eigenre1) {
        for (int u = 0; u < vec.size(); ++u) {
            rowStates.emplace_back(key, u);
        }
    }

    // 建立 col index 对应的 (key2,v)
    std::vector<std::pair<int,int>> colStates; // (key2,v)
    for (const auto &[key2, vec] : eigenre2) {
        for (int v = 0; v < vec.size(); ++v) {
            colStates.emplace_back(key2, v);
        }
    }

    // 遍历每一个 tmatnum
    for (int tmatnum = 0; tmatnum < resultmat.size(); ++tmatnum) {
        fout << "==== Transition matrix for t = " << tvec[tmatnum] << " ====" << std::endl;
        fout << "Transition strengths (transtr) for t = " << tvec[tmatnum] << ":" << std::endl;

        for (int i = 0; i < transtr[tmatnum].size(); ++i) {
            for (int j = 0; j < transtr[tmatnum][i].size(); ++j) {
                fout << transtr[tmatnum][i][j] << " ";
            }
            fout << std::endl;
        }

        fout << std::endl;


        const Eigen::MatrixXd &mat = resultmat[tmatnum];

        for (int i = 0; i < mat.rows(); ++i) {
            auto [key1, u] = rowStates[i];
            for (int j = 0; j < mat.cols(); ++j) {
                auto [key2, v] = colStates[j];
                double val = mat(i, j);
                if (std::abs(val) > 1e-12) { // 筛掉太小的数
                    fout << "(key1=" << key1 << ", state=" << u
                         << ") -> (key2=" << key2 << ", state=" << v
                         << ") : " << val << std::endl;
                }
            }
        }
        fout << std::endl;
    }

    fout.close();
    std::cout << "Results written to " << filename << std::endl;
}




