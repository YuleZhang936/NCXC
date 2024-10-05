// Programmed by Xiaoyu Zhang at Peking University, Beijing, China 2024/10/05
// Ref: PHYSICAL REVIEW RESEARCH 5, 013036 (2023)
// This file tests the NCLibxc library.
// how to compile as an independent program:
/*
  module load libxc/5.2.3-icc17
  module load gcc/7.3.0-wzm
 g++ -std=c++17 -o my_program NCLibxc.cpp LebedevGrid.cpp interface_to_libxc.cpp NCLibxc.h test_NCLibxc.cpp -I/public1/soft/libxc/install/include -L/public1/soft/libxc/install/lib -lxc
*/

#include "NCLibxc.h"
#include <iostream>
#include <vector>
#include <array>
#include <complex>

int main() {
    try {
        // 示例输入数据
        std::vector<double> n = {1.0, 1.0, 1.0};
        std::vector<double> mx = {0.1, 0.1, 0.0};
        std::vector<double> my = {0.0, 0.1, 0.1414};
        std::vector<double> mz = {0.1, 0.0, 0.0};
        int xc_id = 1; // 例如，设置xc_id为1

        // 调用 lda_mc 函数
        auto [E_MC, V_MC] = NCLibxc::lda_mc(xc_id, n, mx, my, mz);

        // 输出结果
        std::cout << "Total E for each real-space grid point:" << std::endl;
        for (const auto &e : E_MC)
            std::cout << e << " ";
        std::cout << std::endl;

        std::cout << "Total V for each real-space grid point:" << std::endl;
        for (const auto &matrix : V_MC) {
            for (const auto &row : matrix) {
                for (const auto &elem : row) {
                    std::cout << elem << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

        auto [E_LC, V_LC] = NCLibxc::lda_lc(xc_id, n, mx, my, mz);

        std::cout << "Total E for each real-space grid point:" << std::endl;
        for (const auto &e : E_LC)
            std::cout << e << " ";
        std::cout << std::endl;

        std::cout << "Total V for each real-space grid point:" << std::endl;
        for (const auto &matrix : V_LC)
        {
            for (const auto &row : matrix)
            {
                for (const auto &elem : row)
                {
                    std::cout << elem << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
    }

    return 0;
}