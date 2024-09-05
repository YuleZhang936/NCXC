// Programmed by Xiaoyu Zhang at Peking University, Beijing, China 2024/08/26
// Ref: PHYSICAL REVIEW RESEARCH 5, 013036 (2023)
// This file implements multi-collinear approach for LDA,GGA,MGGA and locally collinear approach for LDA. Up to the first derivative.
// how to compile as an independent program:
/*
  module load libxc/5.2.3-icc17
  module load gcc/7.3.0-wzm
 g++ -std=c++17 -o my_program NCLibxc.cpp LebedevGrid.cpp interface_to_libxc.cpp -I/public1/soft/libxc/install/include -L/public1/soft/libxc/install/lib -lxc
*/
#include "interface_to_libxc.h"
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <complex>
#include <iomanip>
#include <xc.h>
#include <stdexcept>

std::vector<std::array<double, 4>> MakeAngularGrid(int grid_level);

///////////////////////////////////////////////////////////////////////////////////
//Matrix
// 定义 2x2 复数矩阵类型
using Matrix2x2 = std::array<std::array<std::complex<double>, 2>, 2>;

// 初始化 Pauli 矩阵
const Matrix2x2 sigma_x = {{{std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0)},
                            {std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0)}}};

const Matrix2x2 sigma_y = {{{std::complex<double>(0.0, 0.0), std::complex<double>(0.0, -1.0)},
                            {std::complex<double>(0.0, 1.0), std::complex<double>(0.0, 0.0)}}};

const Matrix2x2 sigma_z = {{{std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0)},
                            {std::complex<double>(0.0, 0.0), std::complex<double>(-1.0, 0.0)}}};

// 矩阵加法
Matrix2x2 add(const Matrix2x2 &a, const Matrix2x2 &b)
{
    Matrix2x2 result;
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            result[i][j] = a[i][j] + b[i][j];
        }
    }
    return result;
}

// 矩阵数乘
Matrix2x2 scalar_multiply(const std::complex<double> &scalar, const Matrix2x2 &matrix)
{
    Matrix2x2 result;
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            result[i][j] = scalar * matrix[i][j];
        }
    }
    return result;
}

// 定义 operator+ 来处理 Matrix2x2 类型的加法
Matrix2x2 operator+(const Matrix2x2& a, const Matrix2x2& b) {
    Matrix2x2 result;
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            result[i][j] = a[i][j] + b[i][j];
        }
    }
    return result;
}

// 根据直角坐标构建 Pauli 自旋矩阵
Matrix2x2 construct_pauli_matrix(double x, double y, double z)
{
    Matrix2x2 pauli_matrix = add(add(scalar_multiply(x, sigma_x), scalar_multiply(y, sigma_y)), scalar_multiply(z, sigma_z));
    return pauli_matrix;
}

// identity matrix
Matrix2x2 identity_matrix()
{
    Matrix2x2 result;
    result[0][0] = 1.0;
    result[0][1] = 0.0;
    result[1][0] = 0.0;
    result[1][1] = 1.0;
    return result;
}
///////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////

// post-processing of libxc. get partial derivatives from libxc and integrate them to get the 0th,1st,2nd derivatives of the functional
//Note that e is the exchange and correlation energy per electron per volume. You need to multiply by \rho before the integration.
/*
v1=0, v2=1, f1=0,0, f2=0,1, f3=1,1
*/
void postlibxc_lda(int xc_id, const std::vector<double>& rho_up, const std::vector<double>& rho_down, 
                   std::vector<double>& e, std::vector<double>& v1, std::vector<double>& v2, 
                   std::vector<double>& f1, std::vector<double>& f2, std::vector<double>& f3)
{
    LibxcInterface libxc(xc_id, true); // xc_id now passed from the caller

    std::vector<double> exc = libxc.lda_exc(rho_up, rho_down);
    std::vector<double> vrho_1(rho_up.size()), vrho_2(rho_down.size());
    libxc.lda_vxc(rho_up, rho_down, vrho_1, vrho_2);

    std::vector<double> v2rho2_1(rho_up.size()), v2rho2_2(rho_down.size()), v2rho2_3(rho_up.size());
    libxc.lda_fxc(rho_up, rho_down, v2rho2_1, v2rho2_2, v2rho2_3);

    
    e.resize(rho_up.size());
    v1.resize(rho_up.size());
    v2.resize(rho_down.size());
    f1.resize(rho_up.size());
    f2.resize(rho_down.size());
    f3.resize(rho_up.size());

    for (size_t i = 0; i < rho_up.size(); ++i) {
        e[i] = exc[i] ; 
        v1[i] = vrho_1[i] ; 
        v2[i] = vrho_2[i] ; 
        f1[i] = v2rho2_1[i] ; 
        f2[i] = v2rho2_2[i] ; 
        f3[i] = v2rho2_3[i] ; 
    }
}

///////////////////////////////////////////////////////////////////////////////////

// lda_mc函数，输入是xc_id, n, mx, my, mz，返回每个实空间格点的E和V
std::pair<std::vector<double>, std::vector<Matrix2x2>> lda_mc(int xc_id, const std::vector<double>& n, 
                                                              const std::vector<double>& mx, const std::vector<double>& my, const std::vector<double>& mz)
{
    std::vector<int> nums = {6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810}; // available grid levels

    int grid_level = nums[25];
    std::vector<std::array<double, 4>> grids = MakeAngularGrid(grid_level); // 假设 grid_level 为 14，这里可以调整

    size_t num_points = n.size();
    std::vector<double> E(num_points, 0.0);
    std::vector<Matrix2x2> V(num_points, {{{0.0, 0.0}, {0.0, 0.0}}});
    std::vector<double> m_omega(num_points, 0.0);

    for (const auto &coord : grids)
    {
        double x = coord[0];
        double y = coord[1];
        double z = coord[2];
        double w = coord[3];

        std::vector<double> rho_up(num_points, 0.0);
        std::vector<double> rho_down(num_points, 0.0);

        for (size_t i = 0; i < num_points; ++i)
        {
            m_omega[i] = mx[i] * x + my[i] * y + mz[i] * z;
            rho_up[i] = (n[i] + m_omega[i]) / 2.0;
            rho_down[i] = (n[i] - m_omega[i]) / 2.0;
        }

        std::vector<double> e(num_points, 0.0), v1(num_points, 0.0), v2(num_points, 0.0), f1(num_points, 0.0), f2(num_points, 0.0), f3(num_points, 0.0);
        postlibxc_lda(xc_id, rho_up, rho_down, e, v1, v2, f1, f2, f3);

        std::vector<double> Eeff(num_points, 0.0);
        std::vector<Matrix2x2> Veff(num_points, {{{0.0, 0.0}, {0.0, 0.0}}});

        for (size_t i = 0; i < num_points; ++i)
        {
            Eeff[i] = e[i] + 0.5 * (m_omega[i] / n[i]) * (v1[i] - v2[i]);

            Matrix2x2 pauli_matrix = construct_pauli_matrix(x, y, z);
            Matrix2x2 term1 = scalar_multiply(((v1[i] - v2[i]) + 0.25 * m_omega[i] * (f1[i] + f3[i] - 2 * f2[i])), pauli_matrix);
            Matrix2x2 term2 = scalar_multiply((0.5 * (v1[i] + v2[i]) + 0.25 * m_omega[i] * (f1[i] - f3[i])), identity_matrix());
            Veff[i] = add(term1, term2);

            // 将 Eeff 和 Veff 加到 E 和 V 上，乘权重
            E[i] += Eeff[i]*w;
            V[i] = add(V[i], scalar_multiply(w, Veff[i]));
        }
    }

    return {E, V};
}

// lda_lc函数，输入是xc_id, n, mx, my, mz，返回每个实空间格点的E和V
std::pair<std::vector<double>, std::vector<Matrix2x2>> lda_lc(int xc_id, const std::vector<double>& n, 
                                                              const std::vector<double>& mx, const std::vector<double>& my, const std::vector<double>& mz)
{
    size_t num_points = n.size();
    std::vector<double> E(num_points, 0.0);
    std::vector<Matrix2x2> V(num_points, {{{0.0, 0.0}, {0.0, 0.0}}});
    std::vector<double> m_mod(num_points, 0.0);

    // 计算m_mod，这是m的模长
    for (size_t i = 0; i < num_points; ++i)
    {
        m_mod[i] = std::sqrt(mx[i] * mx[i] + my[i] * my[i] + mz[i] * mz[i]);
    }

    // 计算rho_up和rho_down
    std::vector<double> rho_up(num_points, 0.0);
    std::vector<double> rho_down(num_points, 0.0);

    for (size_t i = 0; i < num_points; ++i)
    {
        rho_up[i] = (n[i] + m_mod[i]) / 2.0;
        rho_down[i] = (n[i] - m_mod[i]) / 2.0;
    }

    // 调用postlibxc_lda函数
    std::vector<double> e(num_points, 0.0), v1(num_points, 0.0), v2(num_points, 0.0), f1(num_points, 0.0), f2(num_points, 0.0), f3(num_points, 0.0);
    postlibxc_lda(xc_id, rho_up, rho_down, e, v1, v2, f1, f2, f3);

    // 计算E和V
    for (size_t i = 0; i < num_points; ++i)
    {
        E[i] = e[i];

        // 计算 x, y, z
        double x = mx[i] / m_mod[i];
        double y = my[i] / m_mod[i];
        double z = mz[i] / m_mod[i];

        // 构建 Pauli 矩阵
        Matrix2x2 pauli_matrix = construct_pauli_matrix(x, y, z);
        Matrix2x2 term1 = scalar_multiply(0.5 * (v1[i] - v2[i]), pauli_matrix);
        Matrix2x2 term2 = scalar_multiply(0.5 * (v1[i] + v2[i]), identity_matrix());
        V[i] = add(term1, term2);
    }

    return {E, V};
}

int main()
{
    try 
    {
        // 示例输入数据
        std::vector<double> n = {1.0, 1.0, 1.0};
        std::vector<double> mx = {0.1, 0.1, 0.0};
        std::vector<double> my = {0.0, 0.1, 0.1414};
        std::vector<double> mz = {0.1, 0.0, 0.0};
        int xc_id = 1; // 例如，设置xc_id为1

        auto [E_MC, V_MC] = lda_mc(xc_id, n, mx, my, mz);

        std::cout << "Total E for each real-space grid point:" << std::endl;
        for (const auto &e : E_MC)
            std::cout << e << " ";
        std::cout << std::endl;

        std::cout << "Total V for each real-space grid point:" << std::endl;
        for (const auto &matrix : V_MC)
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

        auto [E_LC, V_LC] = lda_lc(xc_id, n, mx, my, mz);

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
    } 
    catch (const std::exception& ex) 
    {
        std::cerr << "Error: " << ex.what() << std::endl;
    }

    return 0;
}

