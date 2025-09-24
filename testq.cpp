#include <unsupported/Eigen/CXX11/Tensor>
#include <iostream>
#include <cassert>

int main() {
    // 定义两个矩阵 A (2×2) 和 B (2×2)
    Eigen::MatrixXd A(2, 2);
    Eigen::MatrixXd B(2, 2);

    A << 1, 2,
         3, 4;
    B << 5, 6,
         7, 8;

    // 方法1：广播法计算四维张量 T1
    Eigen::Tensor<double, 2> A_tensor = Eigen::TensorMap<Eigen::Tensor<double, 2>>(A.data(), A.rows(), A.cols());
    Eigen::Tensor<double, 2> B_tensor = Eigen::TensorMap<Eigen::Tensor<double, 2>>(B.data(), B.rows(), B.cols());

    Eigen::array<int, 4> reshape_A = {A.rows(), A.cols(), 1, 1};  // (a,b,1,1)
    Eigen::array<int, 4> reshape_B = {1, 1, B.rows(), B.cols()};  // (1,1,c,d)

    Eigen::Tensor<double, 4> A_reshaped = A_tensor.reshape(reshape_A);
    Eigen::Tensor<double, 4> B_reshaped = B_tensor.reshape(reshape_B);

    Eigen::Tensor<double, 4> T1 = A_reshaped * B_reshaped.broadcast(Eigen::array<int, 4>{1, 1, B.rows(), B.cols()});

    // 方法2：双循环法计算四维张量 T2
    Eigen::Tensor<double, 4> T2(2, 2, 2, 2);

    for (int i = 0; i < A.rows(); ++i) {
        for (int j = 0; j < A.cols(); ++j) {
            for (int k = 0; k < B.rows(); ++k) {
                for (int l = 0; l < B.cols(); ++l) {
                    T2(i, j, k, l) = A(i, j) * B(k, l);
                }
            }
        }
    }

    // 验证 T1 和 T2 是否完全一致
    bool is_equal = true;
    for (int i = 0; i < A.rows(); ++i) {
        for (int j = 0; j < A.cols(); ++j) {
            for (int k = 0; k < B.rows(); ++k) {
                for (int l = 0; l < B.cols(); ++l) {
                    if (T1(i, j, k, l) != T2(i, j, k, l)) {
                        is_equal = false;
                        break;
                    }
                }
                if (!is_equal) break;
            }
            if (!is_equal) break;
        }
        if (!is_equal) break;
    }

    if (is_equal) {
        std::cout << "✅ 方法1和方法2计算结果完全一致！" << std::endl;
    } else {
        std::cout << "❌ 方法1和方法2计算结果不一致！" << std::endl;
    }

    return 0;
}