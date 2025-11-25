//
// Created by wang- on 25-8-6.
//
#define EIGEN_NO_DEBUG
#define EIGEN_UNROLLING_LIMIT 0
// #define EIGEN_USE_MKL_ALL
// #define EIGEN_USE_BLAS
// #define EIGEN_USE_LAPACKE
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
struct DataRow {
    int t;       // 第1列
    int a, b, c, d; // 第2-5列
    int J;       // 第6列
    double value; // 第7列
};

using std::setprecision;
using lambda_lanczos::LambdaLanczos;





std::vector<std::map<int, Matrix4D>> buildVValue(const std::vector<DataRow>& data) {
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
        int j1=nucleus[row.a].j;
        int j2=nucleus[row.b].j;
        int j3=nucleus[row.c].j;
        int j4=nucleus[row.d].j;
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

        // 右因子重组为 d2 x d3
        q_nu[i] = Eigen::Map<Eigen::MatrixXd>(V.col(i).data(), d2, d3);
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

using MatrixIndex = std::pair<int, int>;

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
Eigen::MatrixXd buildResultMatrix(const Eigen::MatrixXd& qmitmu,
                                const std::vector<std::vector<Eigen::MatrixXd>>& allfab) {
    const int size = allfab.size();
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(size, size);


    for (int l = 0; l < size; ++l) {
        for (int m = 0; m < size; ++m) {
            // std::cout << "_____________l = " << l<<","<<"m"<<m << std::endl;
            result(l, m) = computeElementwiseProductSum(allfab[l][m], qmitmu);
        }
    }
    return result;
}

std::map<int, std::map<int, Eigen::MatrixXd>> qmatjcal(
    const std::vector<std::vector<Eigen::MatrixXd>>& q_pi,
    const Eigen::MatrixXd& transformation_matrix,
    const std::vector<Eigen::MatrixXd>& V_it,
    const std::vector<std::vector<Eigen::MatrixXd>>& allfab,
    const std::vector<basis>& bas,
    const std::vector<basism>& basm)
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

        // 构建结果矩阵
        Eigen::MatrixXd resultmat = buildResultMatrix(qmitmu, allfab);

        // 应用变换
        resultmat = hamchange(basm, bas, transformation_matrix, t, mu, resultmat);

        // 存储结果
        qmatall[matrixIndex][rowIndex] = resultmat;
    }

    return qmatall;
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

bool isSymmetricMatrix(const std::vector<std::vector<double>>& matrix) {
    // 检查是否为空矩阵
    if (matrix.empty() || matrix[0].empty()) {
        return false;
    }

    size_t n = matrix.size();

    // 检查是否为方阵
    for (const auto& row : matrix) {
        if (row.size() != n) {
            return false;
        }
    }

    // 检查对称性：matrix[i][j] == matrix[j][i]
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < i; ++j) {  // 只需检查下三角
            if (abs(matrix[i][j] - matrix[j][i])>1e-13) {
                return false;
            }
        }
    }

    return true;
}

bool isSymmetricWithTolerance(const Eigen::MatrixXd& matrix,
                             double absTol = 1e-13,
                             double relTol = 1e-13) {
    if (matrix.rows() != matrix.cols()) return false;

    return (matrix - matrix.transpose()).cwiseAbs().maxCoeff() <=
           std::max(absTol, relTol * matrix.cwiseAbs().maxCoeff());
}

bool isIdentityEigen(const Eigen::MatrixXd& matrix, double tol = 1e-10) {
    return matrix.isApprox(Eigen::MatrixXd::Identity(matrix.rows(), matrix.cols()), tol);
}

int main()
{


    nucleus={
        {0, 4, 3},
        {0, 4, 5},
        {2, 0, 1}
    };
    // nucleus={
    //     {0, 4, 3},
    //     // {0, 4, 5},
    //     // {2, 0, 1}
    // };
    int size1=nucleus.size();
    allnucleus={
        {0, 4, 3},
        {0, 4, 5},
        {2, 0, 1}
    };
    nucleus2=
    {
        {0, 4, 3},
        {0, 4, 5},
        {2, 0, 1}
    };
    int size2=nucleus.size();
    allnucleus2={
        {0, 4, 3},
        {0, 4, 5},
        {2, 0, 1}
    };
    // nucleus={
    //     {0, 4, 5},
    //     {0, 0, 1}
    // };
    // int size1=nucleus.size();
    // allnucleus={
    //     {0, 4, 5},
    //     {0, 0, 1}
    // };
    // nucleus={
    //     {0, 4, 5}
    // };
    // int size1=nucleus.size();
    // allnucleus={
    //     {0, 4, 5}
    // };
    std::vector<std::vector<double> > ystrgetp;
    // ystrgetp = {
    //     {4, 4, 0, -0.05916894},
    //     {0, 0, 0,  0.53584976},
    //     {1, 1, 0,  0.83771066},
    //     {2, 2, 0,  0.06838028},
    //     {3, 3, 0,  0.05412055},
    //     {4, 4, 1, -0.06584450},
    //     {0, 0, 1, -0.78415728},
    //     {0, 1, 1,  0.16363533},
    //     {0, 2, 1, -0.13233403},
    //     {1, 0, 1, -0.16363533},
    //     {1, 1, 1, -0.49331749},
    //     {1, 2, 1, -0.06878742},
    //     {1, 3, 1, -0.12271369},
    //     {2, 0, 1, -0.13233403},
    //     {2, 1, 1,  0.06878742},
    //     {2, 2, 1, -0.05967271},
    //     {2, 3, 1,  0.05329413},
    //     {3, 1, 1, -0.12271369},
    //     {3, 2, 1, -0.05329413}
    // };
    // ystrgetp = {
    //     {0, 0, 0,  1},
    //     {0, 0, 1, 1},
    // };
    ystrgetp = readstrfile("../strfilesd.txt");
    normalizeBySecondColumn(ystrgetp);
    nucleus=nucleus2;
    allnucleus=allnucleus2;
    size1=nucleus.size();


    std::vector<int> rorderp={0,    4 ,   2 ,   4  ,  6   , 8 ,   2  ,  4   , 0  ,  4    ,8   , 4  ,  6   , 0};
    // std::vector<int> rorderp={0,    4,   8 };
    std::vector<std::vector<int>>rorderrp={rorderp,{0,99,0,1,2,3,4,5,6,7,8,9,10,11,12,13}};
    // std::vector<std::vector<int>>rorderrp={rorderp,{0,99,0,1,2}};
    rvecall={0,    4 ,   2 ,   4  ,  6   , 8 ,   2  ,  4   , 0  ,  4    ,8   , 4  ,  6   , 0};
    rparityvec={1  ,  1   , 1   , 1  ,  1   , 1  ,  1   , 1    ,1   , 1   , 1 ,   1  ,  1  ,  1};
    // rvecall={0,    4,  8 };
    // rparityvec={1  ,  1,   1   };
    std::vector<std::pair<std::vector<int>, std::vector<int>>> catchpair1;
    catchpair1=generateValidPairs(4,rorderrp,51);
    std::vector<basis> allbasisp;
    allbasisp=calculateBasisWithRange(4, catchpair1);


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
    std::vector<basism> mtest;
    mtest=calculateBasisform(4,rorderrp);
    Eigen::MatrixXd m1;
    Eigen::MatrixXd mmm=ComputeTransformationMatrixeven({mtest[0]},allbasisp);

    m1=ComputeTransformationMatrixeven( mtest,allbasisp);


    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile<< "m1  " <<":\n"<< m1 << std::endl;
    std::vector<std::vector<Eigen::MatrixXd>>ystrm1= computeystrm( mtest, ystrgetp);
    std::cout<<"begin to get overlapm"<<std::endl;
    Eigen::MatrixXd overlapm=caloverlapmmat(mtest,mtest,ystrm1,ystrm1);
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }

    outfile<<"overlapm" << std::endl;
    outfile << overlapm << std::endl;
    std::cout << "对称矩阵测试: "
              << (isSymmetricWithTolerance(overlapm) ? "是" : "否") << "\n";
    std::cout << "按回车键继续...";
    std::cin.get();  // 等待用户输入回车

    std::vector<std::vector<double> > overlapmch=overlapchange(overlapm,m1);

    // printMatrix(overlapmch);
    updateMatrices(allbasisp, overlapmch,ystrallp);
    std::cout << "对称矩阵测试: "
              << (isSymmetricMatrix(overlapmch) ? "是" : "否") << "\n";
    std::cout << "按回车键继续...";
    std::cin.get();  // 等待用户输入回车
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile <<"overlapmch" << std::endl;
    writeMatrixToFile(overlapmch,outfile);
    std::vector<std::vector<double>>schmitmatp;
    std::vector<std::pair<int, int>> blocksnum=findblocks(allbasisp);


    schmitmatp=gramSchmidtInOverlap(overlapmch,allbasisp,ystrallp,blocksnum);
    std::cout << "对称矩阵测试:overlapch "
                  << (isSymmetricMatrix(overlapmch) ? "是" : "否") << "\n";
    std::cout << "按回车键继续...";
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile << "overlapmch" << std::endl;
    writeMatrixToFile(overlapmch,outfile);
    Eigen::MatrixXd schmitmat1=convertToEigenMatrix(schmitmatp);
    Eigen::MatrixXd overlap1=convertToEigenMatrix(overlapmch);
    Eigen::MatrixXd overlap1cha=schmitmat1 * overlap1 * schmitmat1.transpose();
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile << "overlap1cha" << std::endl;
    outfile << overlap1cha << std::endl;
    outfile << "schmitmat1" << std::endl;
    outfile << schmitmat1 << std::endl;

    std::cout << "单位矩阵测试: "
              << (isIdentityEigen(overlap1cha) ? "是" : "否") << "\n";
    std::cout << "按回车键继续...";
    std::cin.get();  // 等待用户输入回车
    m1=ComputeTransformationMatrixeven( mtest,allbasisp);

    outfile<<"overlapm" << std::endl;
    outfile << overlapm << std::endl;


    std::vector<std::vector<double> > overlapmch1=overlapchange(overlapm,m1);
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }

    outfile<<"overlapmch1" << std::endl;
    std::cout << "对称矩阵测试: "
          << (isSymmetricMatrix(overlapmch1) ? "是" : "否") << "\n";
    std::cout << "按回车键继续...";
    std::cin.get();  // 等待用户输入回车
    outfile<<"overlapmch1" << std::endl;
    writeMatrixToFile(overlapmch1,outfile);
    printBasisVectorOneLine(allbasisp);

    removeZeroRows(m1,mtest,ystrm1);
    m1=ComputeTransformationMatrixeven( mtest,allbasisp);
    std::vector<basism> mtest_1=calculateBasisform_1(4,rorderrp);
    Eigen::MatrixXd m1_1=ComputeTransformationMatrixeven( mtest_1,allbasisp);
    std::vector<std::vector<Eigen::MatrixXd>> ystrm1_1= computeystrm( mtest_1, ystrgetp);
    removeZeroRows(m1_1,mtest_1,ystrm1_1);
    m1_1=ComputeTransformationMatrixeven( mtest_1,allbasisp);


// std::cout<<ystrm1[0][1]<<std::endl;

    std::vector<std::vector<double> > V1val;
    std::vector<DataRow> V1val1;
    std::ifstream input_file("../efc1usda.txt"); // 替换为你的文件名

    if (!input_file.is_open()) {
        std::cerr << "无法打开文件！" << std::endl;
        return 1;
    }

    std::string line;
    while (std::getline(input_file, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        DataRow row;

        // 按新列名读取数据
        if (!(iss >> row.t >> row.a >> row.b >> row.c >> row.d >> row.J >> row.value)) {
            std::cerr << "解析行出错: " << line << std::endl;
            continue;
        }
        DataRow newrow;
        newrow.t = row.t-1;
        newrow.a = row.a-1;
        newrow.b = row.b-1;
        newrow.c = row.c-1;
        newrow.d = row.d-1;
        newrow.J = row.J;
        newrow.value = row.value;
        V1val1.push_back(newrow);
    }
    std::vector<std::map<int, Matrix4D>> V_value;
    std::vector<std::map<int, Matrix4D>> buildVValue1= buildVValue(V1val1);

    std::vector<std::map<int, Matrix4D>>buildVValue2= processFile("../efc1usda.txt");
    if (buildVValue1==buildVValue2)
    {
        std::cout<<"buildVValue1 and buildVValue2 are equal"<<std::endl;
    }
    else
    {
        std::cout<<"buildVValue1 and buildVValue2 are not equal"<<std::endl;
    }

    input_file.close();



    // std::ifstream input_file1("D:/paircalpro/mschemecode/efct3.txt"); // 替换为你的文件名
    //
    //
    // if (!input_file1.is_open()) {
    //     std::cerr << "无法打开文件！" << std::endl;
    //     return 1;
    // }


    // std::vector<std::map<int, Matrix4D>>buildVValue3= processFilepn("D:/paircalpro/mschemecode/efct3.txt");
    // std::map<int,Matrix4D>vpnmat= getvpnval(buildVValue3,{-0.06});
    // std::vector<std::vector<Eigen::MatrixXd>> q_pi;
    // std::vector<Eigen::MatrixXd> V_it;
    // std::vector<std::vector<Eigen::MatrixXd>> q_nu;
    // std::vector<double> strengthvec={-0.06};
    // getsvdresult( vpnmat,q_pi,  // 输出: q_pi[i](j_alpha, j_beta)
    //     V_it,               // 输出: V_it(i, t)
    //     q_nu,strengthvec);
    // printDiagonals(V_it);
    //
    // input_file1.close();
    // Eigen::MatrixXd term02=calf(mtest[1],mtest[2],ystrm1[1],ystrm1[2]);
    // Eigen::MatrixXd term20=calf( mtest[2],mtest[1],ystrm1[2],ystrm1[1]);





    // Matrix4D qcaltest= calg( mtest[0],mtest[0],ystrm1[0],ystrm1[0]);
    // Matrix4D qcaltest1= calg( mtest[2],mtest[0],ystrm1[2],ystrm1[0]);
    // printAsymmetricElements(qcaltest);
    // printNonZeroElements(qcaltest);
    // std::cout<<qcaltest[20][20]<< std::endl;
    // std::vector<std::vector<Matrix4D>> galltest= calgall( mtest, ystrm1);
    // checkInequality(galltest);

    std::vector<double>strength={0.9548078108};
    Eigen::MatrixXd ham1cal=ham1( mtest,ystrm1,buildVValue2,strength);
    std::cout << "ham1cal" << std::endl;
    std::cout << ham1cal << std::endl;
    Eigen::MatrixXd changeham1=hamchange(mtest,allbasisp,m1, 0,0 , ham1cal);
    std::cout << "hamchangecal" << std::endl;
    std::cout << changeham1 << std::endl;
    // std::vector<std::vector<Eigen::MatrixXd>> calfall1=calfall(mtest,ystrm1);
    // checkInequalityt(calfall1) ;
    // // std::vector<double> energyp={0.0 ,0.962  ,2.440,2.990, 2.763};
    // // std::vector<double> energyp={-0.09824028 ,0.84320220,2.25299294  ,2.7671830,2.68070966};
    std::vector<double> energyp={1.9798 ,    -3.9436 ,     -3.0612   };
    Eigen::MatrixXd hamsige1= hamsige(mtest,ystrm1,energyp);
    std::cout << "hamsige" << std::endl;
    std::cout << hamsige1 << std::endl;
    Eigen::MatrixXd hamchangesig1=hamchange(mtest,allbasisp,m1, 0, 0, hamsige1);
    std::cout << "hamsigechange" << std::endl;
    std::cout << hamchangesig1 << std::endl;
    Eigen::MatrixXd hamp=changeham1+ hamchangesig1;
    std::cout << "hamp" << std::endl;
    std::cout << hamp << std::endl;
    Eigen::MatrixXd hampchange=schmitmat1*hamp* schmitmat1.transpose();
    std::vector<std::vector<double>>hampcc=eigenToNestedVector(hampchange);
    std::vector<int> jvecp;
    std::vector<int> blockform1;
    std::vector<int> jjv;
    for (int i=0;i<allbasisp.size();++i)
    {
        int ii;
        if (i==0)
        {
            ii=0;
        }
        jvecp.push_back(allbasisp[i].sj.back()*allbasisp[i].parity);
        if (i>0)
        {
            if (jvecp[i]!=jvecp[i-1]||allbasisp[i].parity!=allbasisp[i-1].parity)
            {
                blockform1.push_back(i-ii);
                ii=i;
                jjv.push_back(allbasisp[i-1].sj.back()*allbasisp[i-1].parity);
            }
        }
        if (i==allbasisp.size()-1)
        {
            blockform1.push_back(i+1-ii);
            jjv.push_back(allbasisp[i].sj.back());

        }

    }
    std::vector<std::vector<std::vector<double>> > blockh;
    blockh=splitBlockDiagonalMatrix(hampcc,blockform1);
    std::map<int, std::vector<double>> eigenre1;
    std::map<std::pair<int, int>, std::vector<double>> engre;
    for (size_t num=0;num<blockh.size();++num)
    {
        int key=jjv[num];
        std::vector<std::vector<double> > ham3=blockh[num];
        int nmat=ham3.size();
        auto mv_mul1 = [&](const std::vector<double>& in, std::vector<double>& out) {
            for(int i = 0; i < nmat; ++i) {
                for(int j = 0; j < nmat; ++j) {
                    out[i] +=ham3[i][j]*in[j];
                }
            }
        };

        LambdaLanczos<double> engine1(mv_mul1, nmat, false, 2); // true means to calculate the largest eigenvalue.
        std::vector<double> eigenvalues1;
        std::vector<std::vector<double>> eigenvectors1;
        engine1.run(eigenvalues1, eigenvectors1);
        eigenre1[key]=eigenvectors1[0];
        for (size_t lever=0;lever<eigenvectors1.size();++lever)
        {
            engre[{key,lever}]=eigenvectors1[lever];
        }
        std::cout << "J= " <<key<< std::endl;
        std::cout << "Eigenvalues: " << std::endl;
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
    }
    std::map<std::pair<int, int>, int> myMap1;
    for (size_t num=0;num<blockh.size();++num)
    {
        int mapn;
        for(int i=0;i<blockh[num].size();++i)
        {
            if (num==0 && i==0)
            {
                mapn=0;
            }
            myMap1[{jjv[num],i}]=mapn;
            mapn=mapn+1;
        }
    }







    //2hezi
    std::vector<std::pair<std::vector<int>, std::vector<int>>> catchpair2;
    catchpair2=generateValidPairs(6,rorderrp,51);
    std::vector<basis> allbasisp2;
    allbasisp2=calculateBasisWithRange(6, catchpair2);


    sizep=allbasisp2.size();
    sizepr=allbasisp2[0].r.size();
    std::vector<std::vector<std::vector<std::vector<double>>>> ystrallp2(
        sizep,  // 第一维：大小为 sizep
        std::vector<std::vector<std::vector<double>>>(
            sizepr,  // 第二维：大小为 sizepr
            std::vector<std::vector<double>>(
                size1,  // 第三维：大小为 size1
                std::vector<double>(size1, 0.0)  // 第四维：大小为 size1，元素初始化为 0.0
            )
        )
    );
    calculateystr(ystrallp2,ystrgetp,allbasisp2);
    std::vector<basism> mtest2;
    mtest2=calculateBasisform(6,rorderrp);
    Eigen::MatrixXd m2;
    std::vector<std::vector<Eigen::MatrixXd>>ystrm2= computeystrm( mtest2, ystrgetp);
    m2=ComputeTransformationMatrixeven( mtest2,allbasisp2);

    Eigen::MatrixXd overlapm2=caloverlapmmat(mtest2,mtest2,ystrm2,ystrm2);
    std::vector<std::vector<double> > overlapmch2=overlapchange(overlapm2,m2);

    updateMatrices(allbasisp2, overlapmch2,ystrallp2);
    std::vector<std::pair<int, int>> blocksnum2=findblocks(allbasisp2);
    std::vector<std::vector<double>>schmitmatp2;
    schmitmatp2=gramSchmidtInOverlap(overlapmch2,allbasisp2,ystrallp2,blocksnum2);
    Eigen::MatrixXd schmitmat2=convertToEigenMatrix(schmitmatp2);
    Eigen::MatrixXd overlap2=convertToEigenMatrix(overlapmch2);
    Eigen::MatrixXd overlap2cha=schmitmat2 * overlap2 * schmitmat2.transpose();
    std::cout << overlap2cha << std::endl;
    m2=ComputeTransformationMatrixeven( mtest2,allbasisp2);
    removeZeroRows(m2,mtest2,ystrm2);
    m2=ComputeTransformationMatrixeven( mtest2,allbasisp2);

    // double overtext=caloverlapm(mtest2[0],mtest2[0],ystrm2[0],ystrm2[0]);
    // int bas2size=allbasisp2.size();
    // Eigen::MatrixXd overtextmat=Eigen::MatrixXd::Zero(bas2size,bas2size);
    // for (int i = 0; i < bas2size; ++i)
    // {
    //     for (int j = 0; j < bas2size; ++j)
    //     {
    //         overtextmat(i,j)=caloverlapm(mtest2[i],mtest2[j],ystrm2[i],ystrm2[j]);
    //     }
    // }
    // std::cout << overtextmat << std::endl;
    // double sftext=calsf(mtest[21],mtest2[0],ystrm1[21],ystrm2[0],4);
    std::vector<int>rget={4};
    std::vector<Eigen::MatrixXd>  qstrm0={};
    std::vector<Eigen::MatrixXd> qstrm1={};
    std::vector<int> rvec={};
    std::vector<std::vector<std::vector<double>>> transtrvec={};
    outfile << "allbasisp2" << std::endl;
    printBasisVectorOneLine(allbasisp2);
    transtrnocouple(rget, qstrm0, qstrm1, transtrvec,rvec);
    std::vector<Eigen::MatrixXd> doublemat=
                caldoublemat(
                allbasisp, mtest, mtest_1,
                ystrm1, ystrm1_1, m1, m1_1,
                schmitmat1,
                allbasisp2, mtest2, mtest2,
                ystrm2, ystrm2, m2, m2,
                schmitmat2,
                rvec, transtrvec
            );







    for (int i=0;i<rvec.size();++i)
    {
        Eigen::MatrixXd transtrt=convertToEigenMatrix(transtrvec[i]);
        if (!outfile.is_open())
        {
            // 如果文件未打开，则尝试打开文件
            outfile.open("basis_output.txt",std::ios::app);
        }

        outfile<< "dsfmatchsch for nucleus " << i << ":\n";
        outfile<<"sfmat"<<"\n"
        <<doublemat[i] << "\n"
        <<"rvec"<<"\n"
        <<rvec[i] << "\n"
        <<"transtrvec"<<"\n"
        << transtrt << "\n"
        << "------------------------" << std::endl;
    }



    Eigen::MatrixXd ham2cal=ham1( mtest2,ystrm2,buildVValue2,strength);
    std::cout << "ham2cal" << std::endl;
    std::cout << ham2cal << std::endl;
    Eigen::MatrixXd changeham2=hamchange(mtest2,allbasisp2,m2, 0, 0, ham2cal);
    std::cout << "hamchangecal2" << std::endl;
    std::cout << changeham2 << std::endl;
    Eigen::MatrixXd hamsige2= hamsige(mtest2,ystrm2,energyp);
    std::cout << "hamsige2" << std::endl;
    std::cout << hamsige2 << std::endl;
    Eigen::MatrixXd hamchangesig2=hamchange(mtest2,allbasisp2,m2, 0, 0, hamsige2);
    std::cout << "hamsigechange2" << std::endl;
    std::cout << hamchangesig2 << std::endl;
    Eigen::MatrixXd hamp2=changeham2+ hamchangesig2;
    std::cout << "hamp2" << std::endl;
    std::cout << hamp2 << std::endl;
    Eigen::MatrixXd hampchange2=schmitmat2*hamp2* schmitmat2.transpose();
    std::vector<std::vector<double>>hampcc2=eigenToNestedVector(hampchange2);




    std::vector<int> jvecp1;
    std::vector<int> blockform11;
    std::vector<int> jjv1;
    for (int i=0;i<allbasisp2.size();++i)
    {
        int ii;
        if (i==0)
        {
            ii=0;
        }
        jvecp1.push_back(allbasisp2[i].sj.back()*allbasisp2[i].parity);
        if (i>0)
        {
            if (jvecp1[i]!=jvecp1[i-1] || allbasisp2[i].parity!=allbasisp2[i-1].parity)
            {
                blockform11.push_back(i-ii);
                ii=i;
                jjv1.push_back(allbasisp2[i-1].sj.back()*allbasisp2[i].parity);
            }
        }
        if (i==allbasisp2.size()-1)
        {
            blockform11.push_back(i+1-ii);
            jjv1.push_back(allbasisp2[i].sj.back());

        }

    }
    std::vector<std::vector<std::vector<double>> > blockh1;
    blockh1=splitBlockDiagonalMatrix(hampcc2,blockform11);
    std::map<int, std::vector<double>> eigenre11;
    std::map<std::pair<int, int>, std::vector<double>> engre1;

    for (size_t num=0;num<blockh1.size();++num)
    {
        int key=jjv1[num];
        std::vector<std::vector<double> > ham3=blockh1[num];
        int nmat=ham3.size();
        auto mv_mul1 = [&](const std::vector<double>& in, std::vector<double>& out) {
            for(int i = 0; i < nmat; ++i) {
                for(int j = 0; j < nmat; ++j) {
                    out[i] +=ham3[i][j]*in[j];
                }
            }
        };

        LambdaLanczos<double> engine11(mv_mul1, nmat, false, 2); // true means to calculate the largest eigenvalue.
        std::vector<double> eigenvalues1;
        std::vector<std::vector<double>> eigenvectors1;
        engine11.run(eigenvalues1, eigenvectors1);
        eigenre11[key]=eigenvectors1[0];
        for (size_t lever=0;lever<eigenvectors1.size();++lever)
        {
            engre1[{key,lever}]=eigenvectors1[lever];
        }

        std::cout << "J= " <<key<< std::endl;
        std::cout << "Eigenvalues: " << std::endl;
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
    }

    std::map<std::pair<int, int>, int> myMap2;
    for (size_t num=0;num<blockh1.size();++num)
    {
        int mapn;
        for(int i=0;i<blockh1[num].size();++i)
        {
            if (num==0 && i==0)
            {
                mapn=0;
            }
            myMap2[{jjv1[num],i}]=mapn;
            mapn=mapn+1;
        }
    }
    double ree=0;

    for (int y=0;y<doublemat.size();++y)
    {
        std::vector<std::vector<double> > matrixa2ch=eigenToNestedVector(doublemat[y]);
        for (auto& entry : engre) {
            const std::pair<int, int>& key = entry.first;
            std::vector<double>& values = entry.second;
            for (auto& entry1:engre1)
            {
                const std::pair<int, int>& key1 = entry1.first;
                std::vector<double>& values1 = entry1.second;
                double result=0;
                for (int i=0;i<values.size();++i)
                {
                    int mnum=myMap1[{key.first,i}];
                    for (int j=0;j<values1.size();++j)
                    {
                        int mnum1=myMap2[{key1.first,j}];
                        result += values[i]*values1[j]*matrixa2ch[mnum][mnum1];
                        if ( isnan(result))
                        {
                            std::cout<< matrixa2ch[mnum][mnum1]<<" "<<values[i]<<" "<<values1[j]<<std::endl;
                        }
                    }
                }

                ree+=std::pow(result,2);
                if (!outfile.is_open())
                {
                    // 如果文件未打开，则尝试打开文件
                    outfile.open("basis_output.txt",std::ios::app);
                }
                if (outfile.is_open()) {
                    outfile <<"核子角动量"<<y
                            << "初态角动量 " << key.first
                            << " 初态位置 " << key.second
                            << " 末态角动量 " << key1.first
                            << " 末态位置 " << key1.second
                            << " 谱因子 S " << result
                            << std::endl;
                    outfile.close(); // 关闭文件
                } else {
                    std::cerr << "无法打开文件！" << std::endl;
                }




            }
            // 操作键（std::pair<int, int>）
        }
    }







    // std::map<int,std::map<int,Eigen::MatrixXd>>qmatpi= qmatjcal( q_pi, m1,V_it,calfall1,allbasisp,mtest);
    // std::cout<<"matre = "<<qmatpi[2][0]<<std::endl;







    return 0;



}