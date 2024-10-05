#ifndef NCLIBXC_H
#define NCLIBXC_H

#include <vector>
#include <array>
#include <complex>
#include <utility>

// 定义 2x2 复数矩阵类型
using Matrix2x2 = std::array<std::array<std::complex<double>, 2>, 2>;

class NCLibxc {
public:
    // 初始化 Pauli 矩阵
    static const Matrix2x2 sigma_x;
    static const Matrix2x2 sigma_y;
    static const Matrix2x2 sigma_z;

    // 矩阵加法
    static Matrix2x2 add(const Matrix2x2 &a, const Matrix2x2 &b);

    // 矩阵数乘
    static Matrix2x2 scalar_multiply(const std::complex<double> &scalar, const Matrix2x2 &matrix);

    


    // 根据直角坐标构建 Pauli 自旋矩阵
    static Matrix2x2 construct_pauli_matrix(double x, double y, double z);

    // identity matrix
    static Matrix2x2 identity_matrix();

    // post-processing of libxc. get partial derivatives from libxc and integrate them to get the 0th,1st,2nd derivatives of the functional
    static void postlibxc_lda(int xc_id, const std::vector<double>& rho_up, const std::vector<double>& rho_down, 
                              std::vector<double>& e, std::vector<double>& v1, std::vector<double>& v2, 
                              std::vector<double>& f1, std::vector<double>& f2, std::vector<double>& f3);

    // lda_mc函数，输入是xc_id, n, mx, my, mz，返回每个实空间格点的E和V
    static std::pair<std::vector<double>, std::vector<Matrix2x2>> lda_mc(int xc_id, const std::vector<double>& n, 
                                                                  const std::vector<double>& mx, const std::vector<double>& my, const std::vector<double>& mz);

    // lda_lc函数，输入是xc_id, n, mx, my, mz，返回每个实空间格点的E和V
    static std::pair<std::vector<double>, std::vector<Matrix2x2>> lda_lc(int xc_id, const std::vector<double>& n, 
                                                                  const std::vector<double>& mx, const std::vector<double>& my, const std::vector<double>& mz);
};

// 声明 MakeAngularGrid 函数
std::vector<std::array<double, 4>> MakeAngularGrid(int grid_level);
// 定义 operator+ 来处理 Matrix2x2 类型的加法
Matrix2x2 operator+(const Matrix2x2& a, const Matrix2x2& b);

#endif // NCLIBXC_H