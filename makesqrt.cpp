#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
    // 设置输出精度为小数点后16位
    std::cout << std::fixed << std::setprecision(16);

    std::cout << "const std::array<double, 18> sqrttab = {\n";
    for (int i = 0; i <= 17; ++i) {
        double sqrt_val = std::sqrt(i); // 计算平方根
        std::cout << "    " << sqrt_val;
        if (i < 17) {
            std::cout << ","; // 在最后一个值后不加逗号
        }
        std::cout << "\n";
    }
    std::cout << "};\n";

    return 0;
}
