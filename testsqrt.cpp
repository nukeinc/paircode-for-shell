#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <array>

constexpr int DIM1 = 5; // 第一维大小
constexpr int DIM2 = 5; // 第二维大小
constexpr int DIM3 = 5; // 第三维大小
constexpr int TOTAL_SIZE = DIM1 * DIM2 * DIM3; // 一维数组大小

// **1. 三维 `std::vector`**
std::vector<std::vector<std::vector<double>>> create3DVector() {
    return std::vector<std::vector<std::vector<double>>>(DIM1,
        std::vector<std::vector<double>>(DIM2, std::vector<double>(DIM3, 1.0)));
}

// **2. 一维 `std::vector`**
std::vector<double> create1DVector() {
    return std::vector<double>(TOTAL_SIZE, 1.0);
}

// **3. 一维 `std::array<double, TOTAL_SIZE>`**
std::array<double, TOTAL_SIZE> create1DArray() {
    std::array<double, TOTAL_SIZE> arr;
    arr.fill(1.0);
    return arr;
}

// **测试函数**
double test3DVector(const std::vector<std::vector<std::vector<double>>>& vec, const std::vector<int>& indices) {
    double sum = 0.0;
    for (const auto& idx : indices) {
        int i = (idx / (DIM2 * DIM3)) % DIM1;
        int j = (idx / DIM3) % DIM2;
        int k = idx % DIM3;
        sum += vec[i][j][k];
    }
    return sum;
}

double test1DVector(const std::vector<double>& vec, const std::vector<int>& indices) {
    double sum = 0.0;
    for (const auto& idx : indices) {
        sum += vec[idx];
    }
    return sum;
}

double test1DArray(const std::array<double, TOTAL_SIZE>& arr, const std::vector<int>& indices) {
    double sum = 0.0;
    for (const auto& idx : indices) {
        sum += arr[idx];
    }
    return sum;
}

int main() {
    // **初始化数据**
    auto ystrall3D = create3DVector();
    auto ystrall1DVec = create1DVector();
    auto ystrall1DArr = create1DArray();

    // **生成随机访问索引**
    std::vector<int> random_indices(1000000);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, TOTAL_SIZE - 1);

    for (auto& idx : random_indices) {
        idx = dist(gen);
    }

    // **测试 3D `std::vector`**
    auto start = std::chrono::high_resolution_clock::now();
    double sum3D = test3DVector(ystrall3D, random_indices);
    auto end = std::chrono::high_resolution_clock::now();
    double time3D = std::chrono::duration<double, std::milli>(end - start).count();

    // **测试 1D `std::vector`**
    start = std::chrono::high_resolution_clock::now();
    double sum1DVec = test1DVector(ystrall1DVec, random_indices);
    end = std::chrono::high_resolution_clock::now();
    double time1DVec = std::chrono::duration<double, std::milli>(end - start).count();

    // **测试 1D `std::array`**
    start = std::chrono::high_resolution_clock::now();
    double sum1DArr = test1DArray(ystrall1DArr, random_indices);
    end = std::chrono::high_resolution_clock::now();
    double time1DArr = std::chrono::duration<double, std::milli>(end - start).count();

    // **输出结果**
    std::cout << "读取 1,000,000 次数据的时间 (ms):\n";
    std::cout << "3D `std::vector`: " << time3D << " ms (sum=" << sum3D << ")\n";
    std::cout << "1D `std::vector`: " << time1DVec << " ms (sum=" << sum1DVec << ")\n";
    std::cout << "1D `std::array`: " << time1DArr << " ms (sum=" << sum1DArr << ")\n";

    return 0;
}
