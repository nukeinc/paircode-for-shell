#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>

using namespace std;

// 矩阵与向量相乘（并行化）
vector<double> matVecMult(const vector<vector<double>>& A, const vector<double>& v) {
    size_t n = A.size();
    vector<double> result(n, 0.0);

    #pragma omp parallel for
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            result[i] += A[i][j] * v[j];
        }
    }
    return result;
}

// 计算向量的2范数
double norm(const vector<double>& v) {
    double sum = 0.0;
    for (double val : v) {
        sum += val * val;
    }
    return sqrt(sum);
}

// Lanczos 算法（并行化）
void lanczos(const vector<vector<double>>& A, size_t m) {
    size_t n = A.size();

    // 初始化向量
    vector<double> v_0(n, 0.0);
    v_0[0] = 1.0;  // 选择第一个分量为 1，其余为 0

    vector<vector<double>> V(m, vector<double>(n, 0.0));  // 存储正交化后的向量
    V[0] = v_0;

    vector<double> alpha(m, 0.0);  // 存储三对角矩阵 T 的对角元素
    vector<double> beta(m, 0.0);   // 存储三对角矩阵 T 的非对角元素

    for (size_t k = 1; k < m; ++k) {
        // 计算矩阵与向量的乘积
        vector<double> w = matVecMult(A, V[k-1]);

        // 计算 alpha_k
        alpha[k] = 0.0;
        for (size_t i = 0; i < n; ++i) {
            alpha[k] += V[k-1][i] * w[i];
        }

        // 更新 w
        for (size_t i = 0; i < n; ++i) {
            w[i] -= alpha[k] * V[k-1][i];
        }

        // 计算 beta_k
        beta[k] = norm(w);

        // 更新 V[k] 为正交化后的向量
        for (size_t i = 0; i < n; ++i) {
            V[k][i] = w[i] / beta[k];
        }

        // 正交化向量 V[k]
        #pragma omp parallel for
        for (size_t i = 0; i < k; ++i) {
            double dotProduct = 0.0;
            for (size_t j = 0; j < n; ++j) {
                dotProduct += V[k][j] * V[i][j];
            }
            for (size_t j = 0; j < n; ++j) {
                V[k][j] -= dotProduct * V[i][j];
            }
        }
    }

    // 输出三对角矩阵 T_k
    cout << "Tridiagonal matrix T_k:" << endl;
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < m; ++j) {
            cout << alpha[i] << " ";
        }
        cout << endl;
    }

    // 输出特征值
    cout << "\nEigenvalues of the tridiagonal matrix T_k:" << endl;
    for (size_t i = 0; i < m; ++i) {
        cout << "Eigenvalue " << i+1 << ": " << alpha[i] << endl;
    }
}