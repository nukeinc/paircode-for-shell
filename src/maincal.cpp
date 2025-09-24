//
// Created by wang- on 25-7-29.
//


#define EIGEN_NO_DEBUG
#define EIGEN_UNROLLING_LIMIT 0
// #define EIGEN_USE_MKL_ALL
// #define EIGEN_USE_BLAS
// #define EIGEN_USE_LAPACKE
#include "maincal.h"
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
#include <sstream>
#include <string>
#include <map>
#include <utility>
#include "getbasis.h"



std::vector<std::map<int, Matrix4D>> buildVValue(const std::vector<DataRow>& data, const int num=1) {
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


std::vector<std::map<int, Matrix4D>> processFile(const std::string& filename) {
    std::vector<DataRow> V1val1;

    // 打开文件
    std::ifstream input_file(filename);
    if (!input_file.is_open()) {
        std::cerr << "无法打开文件！" << std::endl;
        return {}; // 返回一个空的 vector
    }

    std::string line;
    // 按行读取文件
    while (std::getline(input_file, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        DataRow row;

        // 解析每行数据
        if (!(iss >> row.t >> row.a >> row.b >> row.c >> row.d >> row.J >> row.value)) {
            std::cerr << "解析行出错: " << line << std::endl;
            continue;
        }

        // 根据需要进行值的处理
        DataRow newrow;
        newrow.t = row.t - 1;
        newrow.a = row.a - 1;
        newrow.b = row.b - 1;
        newrow.c = row.c - 1;
        newrow.d = row.d - 1;
        newrow.J = row.J;
        newrow.value = row.value*4;

        // 将处理后的 DataRow 添加到 V1val1 中
        V1val1.push_back(newrow);
    }

    // 使用 buildVValue 函数构建 V_value
    std::vector<std::map<int, Matrix4D>> V_value;
    V_value = buildVValue(V1val1); // 假设 buildVValue 已定义

    return V_value;  // 返回处理后的数据
}


std::vector<std::map<int, Matrix4D>> buildVValuepn(const std::vector<DataRow>& data) {
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


std::vector<std::map<int, Matrix4D>> processFilepn(const std::string& filename) {
    std::vector<DataRow> V1val1;

    // 打开文件
    std::ifstream input_file(filename);
    if (!input_file.is_open()) {
        std::cerr << "无法打开文件！" << std::endl;
        return {}; // 返回一个空的 vector
    }

    std::string line;
    // 按行读取文件
    while (std::getline(input_file, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        DataRow row;

        // 解析每行数据
        if (!(iss >> row.t >> row.a >> row.c >> row.b >> row.d >> row.J >> row.value)) {
            std::cerr << "解析行出错: " << line << std::endl;
            continue;
        }

        // 根据需要进行值的处理
        DataRow newrow;
        newrow.t = row.t - 1;
        newrow.a = row.a - 1;
        newrow.b = row.b - 1;
        newrow.c = row.c - 1;
        newrow.d = row.d - 1;
        newrow.J = row.J;
        newrow.value = row.value;

        // 将处理后的 DataRow 添加到 V1val1 中
        V1val1.push_back(newrow);
    }

    // 使用 buildVValue 函数构建 V_value
    std::vector<std::map<int, Matrix4D>> V_value;
    V_value = buildVValuepn(V1val1); // 假设 buildVValue 已定义

    return V_value;  // 返回处理后的数据
}

std::map<int,Matrix4D>getvpnval(const std::vector<std::map<int, Matrix4D>>& V_value,std::vector<double> strengthvec)
{
    int nucnu=nucleus.size();
    int nucnu2=nucleus2.size();
    int max_t = 0;
    int max_jnu1=0;
    int max_jnu2=0;
    std::map<int,Matrix4D> result;

    for (size_t i=0;i<nucleus.size();++i)
    {
        if (max_jnu1<nucleus[i].j)
        {
            max_jnu1=nucleus[i].j;
        }

    }
    for (size_t i=0;i<nucleus2.size();++i)
    {
        if (max_jnu1<nucleus2[i].j)
        {
            max_jnu2=nucleus2[i].j;
        }


    }
    max_t=max_jnu1+max_jnu2;
    for (int tt=0;tt<V_value.size();++tt)
    {
        std::map<int, Matrix4D>vvmatrix=V_value[tt];
        for (int t=0;t<=max_t;t=t+2)
        {
            Matrix4D vpnmatrix(nucnu, std::vector<Eigen::MatrixXd>(nucnu, Eigen::MatrixXd::Zero(nucnu2, nucnu2)));
            for (int a=0;a<nucnu;++a)
            {
                int ja=nucleus[a].j;
                for (int b=0;b<nucnu;++b)
                {
                    int jb=nucleus[b].j;
                    for (int c=0;c<nucnu2;++c)
                    {
                        int jc=nucleus2[c].j;
                        for (int d=0;d<nucnu2;++d)
                        {
                            int jd=nucleus2[d].j;
                            for (const auto& v_pair : vvmatrix)
                            {
                                int Jsum = v_pair.first;
                                const Matrix4D& v_matrix = v_pair.second;

                                if (v_matrix[a][b](c,d)!=0)
                                {
                                    double c1=std::pow(-1,(jc+jb+Jsum)/2);
                                    double c2=(Jsum+1);
                                    double c3=util::wigner_6j(ja,jc,Jsum,jd,jb,t);
                                    vpnmatrix[a][b](c,d)+=v_matrix[a][b](c,d)*c1*c2*c3*1;
                                }
                            }
                        }
                    }
                }
            }
            result[t]=vpnmatrix;

        }
    }
    return result;
}

// using namespace Eigen;


// 假设四维张量的维度为 (n, n, n, n)
void decompose_right_tensor(const Eigen::Tensor<double, 4>& right_tensor,
                          int rank,
                          std::vector<Eigen::MatrixXd>& q_pi,
                          Eigen::MatrixXd& V_it,
                          std::vector<Eigen::MatrixXd>& q_nu) {
    // 获取实际维度
    const int d0 = right_tensor.dimension(0); // j_alpha
    const int d1 = right_tensor.dimension(1); // j_beta
    const int d2 = right_tensor.dimension(2); // j_gamma
    const int d3 = right_tensor.dimension(3); // j_delta

    // 1. 展平张量为矩阵 (d0*d1) x (d2*d3)
    Eigen::MatrixXd mat(d0 * d1, d2 * d3);
    for (int ia = 0; ia < d0; ++ia) {
        for (int ib = 0; ib < d1; ++ib) {
            for (int ic = 0; ic < d2; ++ic) {
                for (int id = 0; id < d3; ++id) {
                    const int row = ia * d1 + ib;
                    const int col = ic * d3 + id;
                    mat(row, col) = right_tensor(ia, ib, ic, id);
                }
            }
        }
    }

    // 2. SVD分解
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // 3. 确定有效秩
    const int max_rank = std::min(mat.rows(), mat.cols());
    if (rank > max_rank) {
        std::cerr << "Warning: Reducing rank from " << rank
                  << " to max possible " << max_rank << std::endl;
        rank = max_rank;
    }

    // 4. 提取成分
    Eigen::MatrixXd U = svd.matrixU().leftCols(rank);
    Eigen::MatrixXd V = svd.matrixV().leftCols(rank);
    Eigen::VectorXd S = svd.singularValues().head(rank);


    // 5. 重组因子矩阵
    q_pi.resize(rank);
    q_nu.resize(rank);
    for (int i = 0; i < rank; ++i) {
        // 左因子重组为 d0 x d1
        q_pi[i] = Eigen::Map<Eigen::MatrixXd>(U.col(i).data(), d0, d1);
        // q_pi[i]= q_pi[i].transpose();
        q_pi[i].transposeInPlace();

        // 右因子重组为 d2 x d3
        q_nu[i] = Eigen::Map<Eigen::MatrixXd>(V.col(i).data(), d2, d3);
        q_nu[i].transposeInPlace();
    }

    // 6. 奇异值矩阵
    V_it = S.asDiagonal();
}
void getsvdresult(const std::map<int,Matrix4D>& vpnval,
                          std::vector<std::vector<Eigen::MatrixXd>>& q_pi,  // 输出: q_pi[i](j_alpha, j_beta)
                          std::vector<Eigen::MatrixXd>& V_it,               // 输出: V_it(i, t)
                          std::vector<std::vector<Eigen::MatrixXd>>& q_nu,std::vector<double>& strengthvec)
{
    for (const auto& vpn : vpnval) {
        int t = vpn.first;
        const Matrix4D& matrices = vpn.second;

        // 1. 将 Matrix4D 转换为四维张量
        Eigen::Tensor<double, 4> tensor = matrix4d_to_tensor(matrices);

        // 2. 使用 SVD 分解
        std::vector<Eigen::MatrixXd> q_pi_temp;
        Eigen::MatrixXd V_it_temp;
        std::vector<Eigen::MatrixXd> q_nu_temp;

        decompose_right_tensor(tensor, 100, q_pi_temp, V_it_temp, q_nu_temp); // 假设秩为10

        // 3. 存储结果
        q_pi.push_back(q_pi_temp);
        V_it.push_back(V_it_temp*strengthvec[0]);
        q_nu.push_back(q_nu_temp);
    }
}

void printDiagonals(const std::vector<Eigen::MatrixXd>& V_it) {
    for (size_t i = 0; i < V_it.size(); ++i) {
        // 提取当前矩阵的对角线元素（自动处理非方阵，取 min(rows, cols)）
        Eigen::VectorXd diag = V_it[i].diagonal();

        std::cout << "Matrix " << i << " diagonal: [";
        for (int j = 0; j < diag.size(); ++j) {
            std::cout << diag[j];
            if (j < diag.size() - 1) std::cout << ", ";
        }
        std::cout << "]" << std::endl;
    }
}



std::vector<MatrixIndex> findNonZeroDiagonal(const std::vector<Eigen::MatrixXd>& V_it) {
    std::vector<MatrixIndex> nonZeroIndices;

    // 遍历 V_it 中的每个矩阵
    for (int matrixIndex = 0; matrixIndex < V_it.size(); ++matrixIndex) {
        const Eigen::MatrixXd& matrix = V_it[matrixIndex];

        // 确保矩阵是方阵
        if (matrix.rows() != matrix.cols()) {
            std::cerr << "Matrix " << matrixIndex << " is not square!" << std::endl;
            continue;
        }

        // 遍历对角线元素
        for (int i = 0; i < matrix.rows(); ++i) {
            if (std::abs(matrix(i, i)) > 1e-7)
            { // 检查对角线元素是否非零
                // 如果对角线元素非零，记录它的矩阵索引和行索引
                nonZeroIndices.push_back(std::make_pair(matrixIndex, i));
            }
        }
    }

    return nonZeroIndices;
}


Eigen::MatrixXd transqjtoqm(Eigen::MatrixXd qj,int t, int mu)
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
            if (ma-mb==mu)
            {
                double cj=util::CG(ja, jb, t, ma, -mb, mu);
                qm(a,b)=cj*qj(numa,numb)*std::pow(-1,(jb+mb)/2.0);
            }
        }
    }
    return qm;
}



// std::map<int,std::map<int,Eigen::MatrixXd>> qmatjcal(const std::vector<std::vector<Eigen::MatrixXd>>& q_pi,Eigen::MatrixXd transformation_matrix,
//     std::vector<Eigen::MatrixXd>V_it,std::vector<std::vector<Eigen::MatrixXd>>allfab,std::vector<basis>bas,std::vector<basism>basm)
// {
//     std::map<int,std::map<int,Eigen::MatrixXd>> qmatall;
//     int nucnum=nucleusm.size();
//     int nucnu=nucleus.size();
//     std::vector<MatrixIndex> nonZeroIndices = findNonZeroDiagonal(V_it);
//     Eigen::MatrixXd qmatj=Eigen::MatrixXd::Zero(nucnum,nucnum);
//     int sizefab=allfab.size();
//     Eigen::MatrixXd resultqsig=Eigen::MatrixXd::Zero(nucnu,nucnu);
//     for (const auto& index : nonZeroIndices)
//     {
//         int matrixIndex = index.first;
//         int rowIndex = index.second;
//
//         int t= matrixIndex * 2; // 假设每个矩阵对应的 t 值是索引的两倍
//         std::vector<Eigen::MatrixXd> qt=q_pi[matrixIndex];
//         Eigen::MatrixXd qti=qt[rowIndex];
//         int mu=0;
//
//         Eigen::MatrixXd qmitmu=transqjtoqm(qti,t,mu);
//         Eigen::MatrixXd resultmat=Eigen::MatrixXd::Zero(sizefab,sizefab);
//         for (int l=0;l<sizefab;++l)
//         {
//             for (int m=0;m<sizefab;++m)
//             {
//                 Eigen::MatrixXd fab=allfab[l][m];
//                 Eigen::MatrixXd resultm = fab.cwiseProduct(qmitmu);
//                 double result = resultm.sum();
//                 resultmat(l,m)=result;
//             }
//         }
//         resultmat=hamchange(basm,bas,transformation_matrix, t, mu, resultmat);
//         qmatall[matrixIndex][rowIndex]=resultmat;
//
//     }
//
//     return qmatall;
// }

// 辅助函数：计算两个矩阵的逐元素乘积的和
double computeElementwiseProductSum(const Eigen::MatrixXd& mat1,
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

// 辅助函数：构建结果矩阵
// Eigen::MatrixXd buildResultMatrix(const Eigen::MatrixXd& qmitmu,
//                                 const std::vector<basism>& basm,
//                                 const std::vector<std::vector<Eigen::MatrixXd>>&ystrm) {
//     const int size = basm.size();
//     Eigen::MatrixXd result = Eigen::MatrixXd::Zero(size, size);
//
//
//     for (int l = 0; l < size; ++l) {
//         basism basl=basm[l];
//         std::vector<Eigen::MatrixXd> ystrml=ystrm[l];
//         for (int m = 0; m < size; ++m) {
//             basism basmm=basm[m];
//             std::vector<Eigen::MatrixXd> ystrmm=ystrm[m];
//             // std::cout << "_____________l = " << l<<","<<"m"<<m << std::endl;
//             Eigen::MatrixXd fab=calf(basl,basmm,ystrml,ystrmm);
//             result(l, m) = computeElementwiseProductSum(fab, qmitmu);
//         }
//     }
//     return result;
// }

Eigen::MatrixXd buildResultMatrix(const Eigen::MatrixXd& qmitmu,
                                 const std::vector<basism>& basm,
                                 const std::vector<std::vector<Eigen::MatrixXd>>& ystrm) {
    const int size = basm.size();
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(size, size);

    #pragma omp parallel for  // 并行化 l 循环
    for (int l = 0; l < size; ++l) {
        const basism& basl = basm[l];  // 用引用避免拷贝
        const auto& ystrml = ystrm[l]; // 用 const 引用
        for (int m = 0; m < size; ++m) {
            const basism& basmm = basm[m];
            const auto& ystrmm = ystrm[m];
            Eigen::MatrixXd fab = calf(basl, basmm, ystrml, ystrmm);
            result(l, m) = computeElementwiseProductSum(fab, qmitmu);
        }
    }
    return result;
}


Eigen::MatrixXd buildResultMatrix_1(const Eigen::MatrixXd& qmitmu,
                                 const std::vector<basism>& basm,
                                 const std::vector<std::vector<Eigen::MatrixXd>>& ystrm,
                                 const std::vector<basism>& basm_1,
                                 const std::vector<std::vector<Eigen::MatrixXd>>& ystrm_1) {
    const int size = basm.size();
    const int size1 = basm_1.size();
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(size1, size);

    #pragma omp parallel for  // 并行化 l 循环
    for (int l = 0; l < size1; ++l) {
        const basism& basl = basm_1[l];  // 用引用避免拷贝
        const auto& ystrml = ystrm_1[l]; // 用 const 引用
        for (int m = 0; m < size; ++m) {
            const basism& basmm = basm[m];
            const auto& ystrmm = ystrm[m];
            Eigen::MatrixXd fab = calf(basl, basmm, ystrml, ystrmm);
            result(l, m) = computeElementwiseProductSum(fab, qmitmu);
        }
    }
    return result;
}



std::map<int, std::map<int, Eigen::MatrixXd>> qmatjcal(
    const std::vector<std::vector<Eigen::MatrixXd>>& q_pi,
    const Eigen::MatrixXd& transformation_matrix,
    const std::vector<Eigen::MatrixXd>& V_it,
    const std::vector<basis>& bas,
    const std::vector<basism>& basm,
    const std::vector<std::vector<Eigen::MatrixXd>>&ystrm)
{
    std::map<int, std::map<int, Eigen::MatrixXd>> qmatall;
    const auto nonZeroIndices = findNonZeroDiagonal(V_it);

    for (const auto& [matrixIndex, rowIndex] : nonZeroIndices) {
        const int t = matrixIndex * 2;  // 假设t值是索引的两倍
        const int mu = 0;

        // 获取当前矩阵
        const Eigen::MatrixXd& qti = q_pi[matrixIndex][rowIndex];

        // 计算转换后的矩阵
        Eigen::MatrixXd qmitmu = transqjtoqm(qti, t, mu);
        std::cout  << " qmitmu"  << std::endl;
        std::cout << qmitmu << std::endl;

        // 构建结果矩阵
        Eigen::MatrixXd resultmat = buildResultMatrix(qmitmu,basm, ystrm);
        std::cout<<"resultmat = "<<"/n"<<resultmat<<std::endl;

        // 应用变换
        resultmat = hamchange(basm, bas, transformation_matrix, t, mu, resultmat);

        // 存储结果
        qmatall[matrixIndex][rowIndex] = resultmat;
    }

    return qmatall;
}


std::map<int, std::map<int, Eigen::MatrixXd>> qmatjcal_1(
    const std::vector<std::vector<Eigen::MatrixXd>>& q_pi,
    const Eigen::MatrixXd& transformation_matrix,
    const Eigen::MatrixXd& transformation_matrix_1,
    const std::vector<Eigen::MatrixXd>& V_it,
    const std::vector<basis>& bas,
    const std::vector<basism>& basm,
    const std::vector<basism>& basm_1,
    const std::vector<std::vector<Eigen::MatrixXd>>&ystrm,
    const std::vector<std::vector<Eigen::MatrixXd>>&ystrm_1)
{
    std::map<int, std::map<int, Eigen::MatrixXd>> qmatall;
    const auto nonZeroIndices = findNonZeroDiagonal(V_it);

    for (const auto& [matrixIndex, rowIndex] : nonZeroIndices) {
        const int t = matrixIndex * 2;  // 假设t值是索引的两倍
        const int mu = -2;

        // 获取当前矩阵
        const Eigen::MatrixXd& qti = q_pi[matrixIndex][rowIndex];

        // 计算转换后的矩阵
        Eigen::MatrixXd qmitmu = transqjtoqm(qti, t, mu);
        std::cout  << " qmitmu"  << std::endl;
        std::cout << qmitmu << std::endl;

        // 构建结果矩阵
        Eigen::MatrixXd resultmat = buildResultMatrix_1(qmitmu,basm, ystrm,basm_1,ystrm_1);
        std::cout<<"resultmat_1 = "<<"/n"<<resultmat<<std::endl;

        // 应用变换
        resultmat = hamchange_1(basm, bas, transformation_matrix, t, mu, resultmat,basm_1, transformation_matrix_1);

        // 存储结果
        qmatall[matrixIndex][rowIndex] = resultmat;
    }

    return qmatall;
}

Eigen::MatrixXd bemecal(
    const Eigen::MatrixXd& qti,
    const Eigen::MatrixXd& transformation_matrix,
    const std::vector<basis>& bas,
    const std::vector<basism>& basm,
    const std::vector<std::vector<Eigen::MatrixXd>>&ystrm,
    const int& t)
{

    const int mu = 0;
    // 计算转换后的矩阵
    Eigen::MatrixXd qmitmu = transqjtoqm(qti, t, mu);

    // 构建结果矩阵
    Eigen::MatrixXd resultmat = buildResultMatrix(qmitmu,basm, ystrm);
    std::cout<<"resultmatBE2 = "<<"/n"<<resultmat<<std::endl;

    // 应用变换
    resultmat = hamchange(basm, bas, transformation_matrix, t, mu, resultmat);
    return resultmat;
}


Eigen::MatrixXd bemecal_1(
    const Eigen::MatrixXd& qti,
    const Eigen::MatrixXd& transformation_matrix,
    const Eigen::MatrixXd& transformation_matrix_1,
    const std::vector<basis>& bas,
    const std::vector<basism>& basm,
    const std::vector<basism>& basm_1,
    const std::vector<std::vector<Eigen::MatrixXd>>&ystrm,
    const std::vector<std::vector<Eigen::MatrixXd>>&ystrm_1,
    const int& t)
{

    const int mu = -2;

    // 计算转换后的矩阵
    Eigen::MatrixXd qmitmu = transqjtoqm(qti, t, mu);

    // 构建结果矩阵
    Eigen::MatrixXd resultmat = buildResultMatrix_1(qmitmu,basm, ystrm,basm_1,ystrm_1);
    std::cout<<"resultmat_1 = "<<"/n"<<resultmat<<std::endl;

    // 应用变换
    resultmat = hamchange_1(basm, bas, transformation_matrix, t, mu, resultmat,basm_1, transformation_matrix_1);
    return resultmat;
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
            std::cout<< "delete basis no."<< a<<std::endl;
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

std::vector<Eigen::MatrixXd> qmatcalall(
    const std::vector<Eigen::MatrixXd>& qti,
    const Eigen::MatrixXd& transformation_matrix,
    const std::vector<basis>& bas,
    const std::vector<basism>& basm,
    const std::vector<std::vector<Eigen::MatrixXd>>&ystrm,
    const std::vector<int>& t)
{
    const int mu = 0;

    const int size = basm.size();
    std::vector<Eigen::MatrixXd> qmitmuvec={};
    std::vector<Eigen::MatrixXd> resultmatvec={};
    for (int a = 0; a < qti.size(); ++a)
    {
        Eigen::MatrixXd qmitmu = transqjtoqm(qti[a], t[a], mu);
        qmitmuvec.push_back(qmitmu);
        Eigen::MatrixXd resultmat = Eigen::MatrixXd::Zero(size, size);
        resultmatvec.push_back(resultmat);
    }
    #pragma omp parallel for  // 并行化 l 循环
    for (int l = 0; l < size; ++l) {
        const basism& basl = basm[l];  // 用引用避免拷贝
        const auto& ystrml = ystrm[l]; // 用 const 引用
        for (int m = 0; m < size; ++m) {
            const basism& basmm = basm[m];
            const auto& ystrmm = ystrm[m];
            Eigen::MatrixXd fab = calf(basl, basmm, ystrml, ystrmm);
            for (int n = 0; n < qmitmuvec.size(); ++n)
            {
                Eigen::MatrixXd qtimu=qmitmuvec[n];
                resultmatvec[n](l,m) = computeElementwiseProductSum(fab, qtimu);
            }
        }
    }
    for (int n = 0; n < resultmatvec.size(); ++n)
    {
        resultmatvec[n] =hamchange(basm, bas, transformation_matrix, t[n], mu, resultmatvec[n]);
    }
    return resultmatvec;
}


std::vector<Eigen::MatrixXd> qmatcalall_1(
    const std::vector<Eigen::MatrixXd>& qti,
    const Eigen::MatrixXd& transformation_matrix,
    const Eigen::MatrixXd& transformation_matrix_1,
    const std::vector<basis>& bas,
    const std::vector<basism>& basm,
    const std::vector<basism>& basm_1,
    const std::vector<std::vector<Eigen::MatrixXd>>&ystrm,
    const std::vector<std::vector<Eigen::MatrixXd>>&ystrm_1,
    const std::vector<int>& t)
{
    const int mu = -2;

    const int size = basm.size();
    const int size1=basm_1.size();
    std::vector<Eigen::MatrixXd> qmitmuvec={};
    std::vector<Eigen::MatrixXd> resultmatvec={};
    for (int a = 0; a < qti.size(); ++a)
    {
        Eigen::MatrixXd qmitmu = transqjtoqm(qti[a], t[a], mu);
        qmitmuvec.push_back(qmitmu);
        Eigen::MatrixXd resultmat = Eigen::MatrixXd::Zero(size1, size);
        resultmatvec.push_back(resultmat);
    }
#pragma omp parallel for  // 并行化 l 循环
    for (int l = 0; l < size1; ++l) {
        const basism& basl = basm_1[l];  // 用引用避免拷贝
        const auto& ystrml = ystrm_1[l]; // 用 const 引用
        for (int m = 0; m < size; ++m) {
            const basism& basmm = basm[m];
            const auto& ystrmm = ystrm[m];
            Eigen::MatrixXd fab = calf(basl, basmm, ystrml, ystrmm);
            for (int n = 0; n < qmitmuvec.size(); ++n)
            {
                if (t[n]==0)
                {
                    continue;
                }
                Eigen::MatrixXd qtimu=qmitmuvec[n];
                resultmatvec[n](l,m) = computeElementwiseProductSum(fab, qtimu);
            }
        }
    }
    for (int n = 0; n < resultmatvec.size(); ++n)
    {
        resultmatvec[n] =hamchange_1(basm, bas, transformation_matrix, t[n], mu, resultmatvec[n],basm_1,transformation_matrix_1);
    }
    return resultmatvec;
}


void lanczoscalham(
    const std::map<int, std::vector<CoupledBasis>>& coupleBasesall2,
    const std::map<int,Eigen::MatrixXd>& hamvec,
    const Eigen::MatrixXd& schmitmat1,
    const Eigen::MatrixXd& schmitmat2,
    std::map<int, std::vector<std::vector<double>>>& eigenre,
    std::map<int, std::vector<double>>& eigenvalue,
    std::map<int, std::vector<std::vector<double>>>& eigenreuncop)
{
    using lambda_lanczos::LambdaLanczos;
    using std::setprecision;
    for (const auto& [key, value] : hamvec)
    {
        size_t n = value.rows();
        Eigen::MatrixXd matrix=value;
        auto mv_mul=[&](const std::vector<double>& in,std::vector<double>&out)
        {
            auto eigen_in=Eigen::Map<const Eigen::VectorXd>(&in[0],in.size());
            auto eigen_out=Eigen::Map<Eigen::VectorXd>(&out[0],out.size());
            eigen_out=matrix*eigen_in;
        };
        LambdaLanczos<double> engine(mv_mul, n, false, 2);
        std::vector<double>eigenvalues1;
        std::vector<std::vector<double>>eigenvectors1;
        engine.run(eigenvalues1, eigenvectors1);
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
        eigenre[key]=eigenvectors1;
        eigenvalue[key]=eigenvalues1;
    }
    // for (const auto& [key,value]:eigenre)
    // {
    //     std::vector<CoupledBasis>coupledbasisvec= coupleBasesall2.at(key);
    //     int size = coupledbasisvec.size();
    //
    //         int sizeval=value.size();
    //         for (int k = 0; k < sizeval; ++k)
    //         {
    //             std::vector<double> eigenvecsig=value[k];
    //             int sizesch1=schmitmat1.rows();
    //             int sizesch2=schmitmat2.rows();
    //             std::vector<double>eigvecchange(sizesch1,0.0);
    //             for (int l = 0; l < sizesch1; ++l)
    //             {
    //                 std::vector<double>eigvecchangesum(sizesch1,0.0);
    //                 for (int m = 0; m < eigenvecsig.size(); ++m)
    //                 {
    //                     int pnu=coupledbasisvec[m].pi;
    //                     int val=eigenvecsig[m];
    //                     sum+=val*eigenvecsig[m];
    //                 }
    //             }
    //         }
    //
    // }
}



void calham(const int& nuculnum,
    const int& m1cal,
    const double& alpha1,
    const std::vector<double>& gvec,
    const std::vector<int>&qorderp1,
    const std::vector<double>& energyp,
    const std::vector<double>& strength,
    const std::vector<std::vector<double> >& ystrgetp,
    const std::vector<std::vector<int>> &rorderp,
    const std::vector<Eigen::MatrixXd>& V_it,
    const std::vector<std::map<int, Matrix4D>>& buildVValue2,
    const std::vector<std::vector<Eigen::MatrixXd>>& q_pi,
    std::vector<basism>& mtest,
    std::vector<basism>& mtest_1,
    Eigen::MatrixXd& schmitmat1,
    std::vector<basis>& allbasisp,
    std::vector<std::vector<double>>& hampcc,
    std::vector<Eigen::MatrixXd>& bemematcal,
    std::vector<int>& singleindex,
    std::vector<int>& tvec,
    std::map<int,std::map<int,Eigen::MatrixXd>>& qmatpicha,
    std::vector<int>& jvecp,
    std::vector<int>& parityvecp
    )
{
    int size1=nucleus.size();
    //构造正交基矢
    std::vector<std::vector<int>>rorderrp=rorderp;
    rorderrp.pop_back();
    rvecall=rorderp[0];
    rparityvec=rorderp.back();
    std::vector<std::pair<std::vector<int>, std::vector<int>>> catchpair1;
    catchpair1=generateValidPairs(nuculnum,rorderrp,51);
    allbasisp=calculateBasisWithRange(nuculnum, catchpair1);
    printBasisVectorOneLine(allbasisp);
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
    mtest=calculateBasisform(nuculnum,rorderrp);
    mtest_1=calculateBasisform_1(nuculnum,rorderrp);
    std::cout<<"mtestsize: "<<mtest.size()<<std::endl;
    std::cout<<"mtest1: "<<mtest_1.size()<<std::endl;
    Eigen::MatrixXd m1;
    Eigen::MatrixXd m1_1;
    m1=ComputeTransformationMatrix( mtest,allbasisp);
    m1_1=ComputeTransformationMatrix( mtest_1,allbasisp);
    std::vector<std::vector<Eigen::MatrixXd>>ystrm1= computeystrm( mtest, ystrgetp);
    removeZeroRows(m1,mtest,ystrm1);
    m1=ComputeTransformationMatrixeven( mtest,allbasisp);
    Eigen::MatrixXd overlapm=caloverlapmmat(mtest,mtest,ystrm1,ystrm1);
    std::cout << overlapm << std::endl;
    std::vector<std::vector<double> > overlapmch=overlapchange(overlapm,m1);
    printMatrix(overlapmch);
    updateMatrices(allbasisp, overlapmch,ystrallp);
    std::vector<std::vector<double>>schmitmatp;
    std::vector<std::pair<int, int>> blocksnum=findblocks(allbasisp);
    schmitmatp=gramSchmidtInOverlap(overlapmch,allbasisp,ystrallp,blocksnum);
    m1=ComputeTransformationMatrix( mtest,allbasisp);
    m1_1=ComputeTransformationMatrix( mtest_1,allbasisp);
    std::vector<std::vector<double> > overlapmch1=overlapchange(overlapm,m1);
    schmitmat1=convertToEigenMatrix(schmitmatp);
    std::cout<<" schmitmat1 = "<<"/n"<< schmitmat1<<std::endl;
    removeZeroRows(m1,mtest,ystrm1);
    std::vector<std::vector<Eigen::MatrixXd>>ystrm1_1= computeystrm( mtest_1, ystrgetp);
    removeZeroRows(m1_1,mtest_1,ystrm1_1);
    m1=ComputeTransformationMatrixeven( mtest,allbasisp);
    m1_1=ComputeTransformationMatrixeven( mtest_1,allbasisp);
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile<<"m1" << std::endl;
    outfile<<m1<<std::endl;
    outfile<<"m1_1" << std::endl;
    outfile<<m1_1<<std::endl;
    //正交基矢构造完成
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile << "allbasism01" << std::endl;
    writeBasismvecToFile(mtest);
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile << "allbasism11" << std::endl;
    writeBasismvecToFile(mtest_1);


    //计算单体矩阵元

    auto nonZeropos= findNonZeroDiagonal(V_it);
    std::vector<Eigen::MatrixXd> allqti;
    for (const auto& [matrixIndex, rowIndex] : nonZeropos) {
        const int t = matrixIndex * 2;  // 假设t值是索引的两倍
        Eigen::MatrixXd qti = q_pi[matrixIndex][rowIndex];
        allqti.push_back(qti);
        tvec.push_back(t);
        singleindex.push_back(0);
    }

    std::vector<std::vector<double> > qstrgetp1;
    qstrgetp1=qstrall(qorderp1,alpha1);
    int sizeq1=qorderp1.size();
    std::vector<std::vector<std::vector<double>>> qstrallp1(sizeq1, std::vector<std::vector<double>>(size1,
    std::vector<double>(size1, 0)));

    getystrbasis(qstrallp1, qstrgetp1, {0,1});

    for (int tnum=0;tnum<qstrallp1.size();++tnum)
    {
        Eigen::MatrixXd qti=convertToEigenMatrix(qstrallp1[tnum]);
        allqti.push_back(qti);
        tvec.push_back(qorderp1[tnum]);
        singleindex.push_back(1);

    }

    if (m1cal==1)
    {
        double gl=gvec[0];
        double gs=gvec[1];
        std::vector<std::vector<double> >m1qget=m1qstrall(gl,gs);
        std::vector<std::vector<std::vector<double>>> m1strallp1(1, std::vector<std::vector<double>>(size1,
        std::vector<double>(size1, 0)));
        getystrbasis(m1strallp1, m1qget, {0});
        Eigen::MatrixXd m1strallp1mat=convertToEigenMatrix(m1strallp1[0]);
        allqti.push_back(m1strallp1mat);
        tvec.push_back(2);
        singleindex.push_back(2);
    }
    std::vector<Eigen::MatrixXd>qmmat=qmatcalall(allqti,m1,allbasisp,mtest,
    ystrm1,tvec);
    std::vector<Eigen::MatrixXd>qmmat_1=qmatcalall_1(allqti,m1,m1_1,allbasisp,
        mtest,mtest_1,ystrm1,ystrm1_1,tvec);
    std::vector<Eigen::MatrixXd>qmat={};
    for (int tnum=0;tnum<qmmat.size();++tnum)
    {
        int t=tvec[tnum];
        if (t!=0)
        {
            for (int i=0;i<allbasisp.size();++i)
            {
                int j1=allbasisp[i].sj.back();
                for (int j=0;j<allbasisp.size();++j)
                {
                    int j2=allbasisp[j].sj.back();
                    int num=(j1+j2+t)/2;
                    if ((num & 1) != 0)
                    {
                        qmmat[tnum](i,j)=qmmat_1[tnum](i,j);
                    }
                }
            }
        }

        qmat.push_back(schmitmat1*qmmat[tnum]* schmitmat1.transpose());
    }
    int qmatpicha_size = 0;
    for (const auto& [matrixIndex, rowIndex] : nonZeropos) {
        qmatpicha[matrixIndex][rowIndex]= qmat[qmatpicha_size];
        ++qmatpicha_size;
    }
    for (int tnum=qmatpicha_size;tnum<qmat.size();++tnum)
    {
        bemematcal.push_back(qmat[tnum]);
    }
    Eigen::MatrixXd ham1cal=ham1( mtest,ystrm1,buildVValue2,strength);
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile<<"ham1cal"<<std::endl;
    outfile<<ham1cal<<std::endl;
    Eigen::MatrixXd changeham1=hamchange(mtest,allbasisp,m1, 0, 0, ham1cal);
    Eigen::MatrixXd hamsige1= hamsige(mtest,ystrm1,energyp);
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile<<"hamsige1"<<std::endl;
    outfile<<hamsige1<<std::endl;
    Eigen::MatrixXd hamchangesig1=hamchange(mtest,allbasisp,m1, 0, 0, hamsige1);
    Eigen::MatrixXd hamp=changeham1+ hamchangesig1;
    Eigen::MatrixXd hampchange=schmitmat1*hamp* schmitmat1.transpose();
    hampcc=eigenToNestedVector(hampchange);
    for (int i=0;i<allbasisp.size();++i)
    {
        jvecp.push_back(allbasisp[i].sj.back());
        parityvecp.push_back(allbasisp[i].parity);
    }
}

void calcouple(
    const std::vector<basis>& allbasisp1,
    const std::vector<basis>& allbasisp2,
    const std::vector<std::vector<double>>& hampcc1,
    const std::vector<std::vector<double>>& hampcc2,
    const std::vector<int>& jvecp1,
    const std::vector<int>& jvecp2,
    const std::vector<int>& parityvecp1,
    const std::vector<int>& parityvecp2,
    const std::vector<Eigen::MatrixXd>& bemematcal1,
    const std::vector<Eigen::MatrixXd>& bemematcal2,
    const std::map<int,std::map<int,Eigen::MatrixXd>>& qmatpicha1,
    const std::map<int,std::map<int,Eigen::MatrixXd>>& qmatpicha2,
    const std::vector<double>& eg1,
    const std::vector<double>& eg2,
    const std::vector<Eigen::MatrixXd>& V_it,
    const std::vector<int>& bej1,
    std::map<int, std::vector<CoupledBasis>>& coupleBasesall2,
    std::map<int, std::vector<std::vector<double>>>& eigenre,
    std::map<int,std::vector<double>>& eigenvalue,
    std::vector<Eigen::MatrixXd>& bematre
    )
{
    std::vector<int> bej=bej1;
    bej.push_back(2);
    int sizeq2=bej.size();
    coupleBasesall2.clear();
    coupleBasesall2=coupleBases2(jvecp2,jvecp1,parityvecp1,parityvecp2);

    #pragma omp parallel for
    for (int idx = 0; idx < (int)coupleBasesall2.size(); ++idx) {
        std::map<int, std::vector<CoupledBasis>>::iterator it;
        it = std::next(coupleBasesall2.begin(), idx); // 普通 iterator
        int key = it->first;
        auto& value = it->second;
        int siz1=coupleBasesall2[key].size();
        std::vector<std::vector<double>> ham(siz1, std::vector<double>(siz1, 0));
        std::vector<std::vector<double>> ham2(siz1, std::vector<double>(siz1, 0));
        std::vector<std::vector<double>> ham3(siz1, std::vector<double>(siz1, 0));
        std::vector<std::vector<double>> bemcal1(siz1, std::vector<double>(siz1, 0));
        std::vector<std::vector<double>> bemcal2(siz1, std::vector<double>(siz1, 0));
        for (int i=0;i< siz1;i++)
        {
            for (int j=0;j< siz1;j++)
            {
                int n1=coupleBasesall2[key][i].ni;
                int n2=coupleBasesall2[key][j].ni;
                int p1=coupleBasesall2[key][i].pi;
                int p2=coupleBasesall2[key][j].pi;
                int jn=jvecp2[n1];
                int jnp=jvecp2[n2];
                int jp=jvecp1[p1];
                int jpp=jvecp1[p2];
                if (n1==n2)
                {
                    ham[i][j]= ham[i][j]+hampcc1[p1][p2];
                }
                if (p1==p2)
                {
                    ham[i][j]=ham[i][j]+hampcc2[n1][n2];
                }
                int deltann=deltatwo(jn,jnp);
                int deltapp=deltatwo(jp,jpp);
                ham[i][j]=ham[i][j]*deltann*deltapp;
                for (auto qmatnuchat:qmatpicha1)
                {
                    int tnum=qmatnuchat.first;
                    for (auto qmatcalsec:qmatnuchat.second)
                    {
                        int rownum=qmatcalsec.first;
                        double Vit=V_it[tnum](rownum,rownum);
                        int tt=2*tnum;
                        std::vector<int> sigl={jpp,jn,std::abs(key),tt};
                        int sig=sign_func(sigl);
                        double c2=util::wigner_6j(jnp,jpp,std::abs(key),jp,jn,tt)*mysqrt(jnp+1)*mysqrt(jpp+1)*Vit;
                        Eigen::MatrixXd matsliqnch= qmatpicha2.at(tnum).at(rownum);
                        Eigen::MatrixXd matsliqpch= qmatpicha1.at(tnum).at(rownum);
                        double c3=matsliqnch(n1,n2)*matsliqpch(p1,p2);
                        ham2[i][j]+=c2*c3*sig;
                    }
                }

                ham3[i][j]=ham[i][j]+ham2[i][j];
            }
        }
        Eigen::MatrixXd hammat=convertToMatrixFast(ham3);
        int nmat=ham3.size();
//        auto mv_mul1 = [&](const std::vector<double>& in, std::vector<double>& out) {
//            for(int i = 0; i < nmat; ++i) {
//                for(int j = 0; j < nmat; ++j) {
//                    out[i] +=ham3[i][j]*in[j];
//                }
//            }
//        };
        auto mv_mul1 = [&](const vector<double>& in, vector<double>& out) {
            auto eigen_in = Eigen::Map<const Eigen::VectorXd>(&in[0], in.size());
            auto eigen_out = Eigen::Map<Eigen::VectorXd>(&out[0], out.size());

            eigen_out = hammat * eigen_in; // Easy version
            // eigen_out.noalias() += matrix * eigen_in; // Efficient version
        };
        LambdaLanczos<double> engine1(mv_mul1, nmat, false, 2); // true means to calculate the largest eigenvalue.
        std::vector<double> eigenvalues1;
        std::vector<std::vector<double>> eigenvectors1;
        engine1.run(eigenvalues1, eigenvectors1);
        #pragma omp critical
        {
            eigenre[key] = eigenvectors1;
            eigenvalue[key] = eigenvalues1;
        }

        //computeEigenvaluesAndEigenvectors(matchange1);
    }
    // 后处理输出
    outfile.open("basis_output.txt", std::ios::app);
    for (const auto& [key, values] : eigenvalue) {
        outfile << "J= " << key << std::endl;
        outfile << "Eigenvalues: " << std::endl;
        for (size_t i = 0; i < values.size(); ++i) {
            outfile << "Eigenvalue " << i + 1 << ": " << values[i] << std::endl;
        }

        const auto& vectors = eigenre[key];
        outfile << "Eigenvectors: " << std::endl;
        for (size_t i = 0; i < vectors.size(); ++i) {
            outfile << "Eigenvector " << i + 1 << ": [ ";
            for (size_t j = 0; j < vectors[i].size(); ++j) {
                outfile << vectors[i][j] << " ";
            }
            outfile << "]" << std::endl;
        }
    }
    outfile.close();


    //calculation transition
    int sizcpbsize=coupleBasesall2.size();
    int sizkey=eigenre.size();
    std::vector<Eigen::MatrixXd> bematrices(sizeq2, Eigen::MatrixXd::Zero(sizkey, sizkey));
    double ep=eg1[0];
    double en=eg2[0];
    int keynum=0;
    for (const auto& [key, value] : eigenre)
    {
        int jisum=std::abs(key);
        int siz1=coupleBasesall2[key].size();
        int vitnum=0;

        int keynum2=0;
        for (const auto& [key2, value2] : eigenre)
        {
            int jfsum=std::abs(key2);
            int siz2=coupleBasesall2[key2].size();
            std::vector<std::vector<std::vector<double>>> bmat1(
                sizeq2, std::vector<std::vector<double>>(siz1, std::vector<double>(siz2, 0.0)));
            std::vector<std::vector<std::vector<double>>> bmat2(
                sizeq2, std::vector<std::vector<double>>(siz1, std::vector<double>(siz2, 0.0)));
            for (int i=0;i< siz1;i++)
            {
                for (int j=0;j< siz2;j++)
                {
                    int n1=coupleBasesall2[key][i].ni;
                    int n2=coupleBasesall2[key2][j].ni;
                    int p1=coupleBasesall2[key][i].pi;
                    int p2=coupleBasesall2[key2][j].pi;
                    int jn=jvecp2[n1];
                    int jnp=jvecp2[n2];
                    int jp=jvecp1[p1];
                    int jpp=jvecp1[p2];


                    int deltann=deltatwo(jn,jnp);
                    int deltapp=deltatwo(jp,jpp);
                    for (int qmatnum=0; qmatnum<sizeq2;++qmatnum)
                    {
                        int t=bej[qmatnum];
                        Eigen::MatrixXd bmatcal1=bemematcal1[qmatnum];
                        Eigen::MatrixXd bmatcal2=bemematcal2[qmatnum];
                        double calpi=deltann*deltatwo(n1,n2)
                        *ufrom6_j(jnp,jp,jfsum,t,jisum,jpp)*bmatcal1(p1,p2);
                        double calnu=deltapp*deltatwo(p1,p2)*std::pow(-1,(jn-jnp+jfsum-jisum)/2)
                        *ufrom6_j(jpp,jn,jfsum,t,jisum,jnp)*bmatcal2(n1,n2);
                        bmat1[qmatnum][i][j]=calpi;
                        bmat2[qmatnum][i][j]=calnu;
                    }
                }
            }
            std::vector<double> eigencal1=eigenre[key][0];
            std::vector<double> eigencal2=eigenre[key2][0];

            for(int qmatnum=0; qmatnum<sizeq2;++qmatnum)
            {
                double epp=eg1[0];
                double enn=eg2[0];
                int t=bej[qmatnum]/2;
                int fact=1;
                if ((t%2)==1 )
                {
                    fact=-1;
                }
                if (qmatnum==sizeq2-1 )
                {
                    epp=1;
                    enn=1;
                }
                double sum=0.0;
                Eigen::MatrixXd matcal111(eigencal1.size(),eigencal2.size());
                Eigen::MatrixXd matcal222(eigencal1.size(),eigencal2.size());

                for (int i=0;i<eigencal1.size();++i)
                {
                    for (int j=0;j<eigencal2.size();++j)
                    {
                        double cal1=eigencal1[i];
                        double cal2=eigencal2[j];
                        sum+=cal1*bmat1[qmatnum][i][j]*cal2*epp;
                        sum+=cal2*bmat2[qmatnum][i][j]*cal1*enn*fact;


                    }
                }
                sum=std::pow(sum,2.0)*(jfsum+1.0)/(jisum+1.0);
                bematrices[qmatnum](keynum,keynum2)=sum;
            }
            keynum2++;
        }
        keynum++;
    }
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile << "bematrices" << std::endl;
    outfile<< bematrices[0] << std::endl;
    outfile << "bematrices" << std::endl;
    outfile<< bematrices[1] << std::endl;
    outfile << "bematrices" << std::endl;
    outfile<< bematrices[2] << std::endl;
    bematre=bematrices;
}




