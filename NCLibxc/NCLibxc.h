#ifndef NCLIBXC_H
#define NCLIBXC_H

#include <vector>
#include <array>
#include <complex>
#include <utility>

namespace NCXC {

// define 2*2 complex matrix
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

    static void postlibxc_gga(int xc_id, 
                            const std::vector<double>& rho0, 
                            const std::vector<double>& rho1, 
                            const std::vector<double>& gradx_rho0, 
                            const std::vector<double>& grady_rho0, 
                            const std::vector<double>& gradz_rho0, 
                            const std::vector<double>& gradx_rho1, 
                            const std::vector<double>& grady_rho1, 
                            const std::vector<double>& gradz_rho1, 
                            const std::vector<double>& grad2xx_rho0, 
                            const std::vector<double>& grad2yy_rho0, 
                            const std::vector<double>& grad2zz_rho0, 
                            const std::vector<double>& grad2xy_rho0, 
                            const std::vector<double>& grad2yz_rho0, 
                            const std::vector<double>& grad2xz_rho0, 
                            const std::vector<double>& grad2xx_rho1, 
                            const std::vector<double>& grad2yy_rho1, 
                            const std::vector<double>& grad2zz_rho1, 
                            const std::vector<double>& grad2xy_rho1, 
                            const std::vector<double>& grad2yz_rho1, 
                            const std::vector<double>& grad2xz_rho1, 
                            const std::vector<double>& grad3xxx_rho0, 
                            const std::vector<double>& grad3xxy_rho0, 
                            const std::vector<double>& grad3xxz_rho0, 
                            const std::vector<double>& grad3xyy_rho0, 
                            const std::vector<double>& grad3xyz_rho0, 
                            const std::vector<double>& grad3xzz_rho0, 
                            const std::vector<double>& grad3yyy_rho0, 
                            const std::vector<double>& grad3yyz_rho0, 
                            const std::vector<double>& grad3yzz_rho0, 
                            const std::vector<double>& grad3zzz_rho0, 
                            const std::vector<double>& grad3xxx_rho1, 
                            const std::vector<double>& grad3xxy_rho1, 
                            const std::vector<double>& grad3xxz_rho1, 
                            const std::vector<double>& grad3xyy_rho1, 
                            const std::vector<double>& grad3xyz_rho1, 
                            const std::vector<double>& grad3xzz_rho1, 
                            const std::vector<double>& grad3yyy_rho1, 
                            const std::vector<double>& grad3yyz_rho1, 
                            const std::vector<double>& grad3yzz_rho1, 
                            const std::vector<double>& grad3zzz_rho1, 
                            std::vector<double>& e, 
                            std::vector<double>& v1, 
                            std::vector<double>& v2, 
                            std::vector<double>& f1, 
                            std::vector<double>& f2, 
                            std::vector<double>& f3);


    static void postlibxc_mgga(int xc_id, 
                                const std::vector<double>& rho0, 
                                const std::vector<double>& rho1, 
                                const std::vector<double>& gradx_rho0, 
                                const std::vector<double>& grady_rho0, 
                                const std::vector<double>& gradz_rho0, 
                                const std::vector<double>& gradx_rho1, 
                                const std::vector<double>& grady_rho1, 
                                const std::vector<double>& gradz_rho1, 
                                const std::vector<double>& grad2xx_rho0, 
                                const std::vector<double>& grad2yy_rho0, 
                                const std::vector<double>& grad2zz_rho0, 
                                const std::vector<double>& grad2xy_rho0, 
                                const std::vector<double>& grad2yz_rho0, 
                                const std::vector<double>& grad2xz_rho0, 
                                const std::vector<double>& grad2xx_rho1, 
                                const std::vector<double>& grad2yy_rho1, 
                                const std::vector<double>& grad2zz_rho1, 
                                const std::vector<double>& grad2xy_rho1, 
                                const std::vector<double>& grad2yz_rho1, 
                                const std::vector<double>& grad2xz_rho1, 
                                const std::vector<double>& grad3xxx_rho0, 
                                const std::vector<double>& grad3xxy_rho0, 
                                const std::vector<double>& grad3xxz_rho0, 
                                const std::vector<double>& grad3xyy_rho0, 
                                const std::vector<double>& grad3xyz_rho0, 
                                const std::vector<double>& grad3xzz_rho0, 
                                const std::vector<double>& grad3yyy_rho0, 
                                const std::vector<double>& grad3yyz_rho0, 
                                const std::vector<double>& grad3yzz_rho0, 
                                const std::vector<double>& grad3zzz_rho0, 
                                const std::vector<double>& grad3xxx_rho1, 
                                const std::vector<double>& grad3xxy_rho1, 
                                const std::vector<double>& grad3xxz_rho1, 
                                const std::vector<double>& grad3xyy_rho1, 
                                const std::vector<double>& grad3xyz_rho1, 
                                const std::vector<double>& grad3xzz_rho1, 
                                const std::vector<double>& grad3yyy_rho1, 
                                const std::vector<double>& grad3yyz_rho1, 
                                const std::vector<double>& grad3yzz_rho1, 
                                const std::vector<double>& grad3zzz_rho1, 
                                const std::vector<double>& tau0,
                                const std::vector<double>& tau1,
                                const std::vector<double>& gradx_tau0, 
                                const std::vector<double>& grady_tau0, 
                                const std::vector<double>& gradz_tau0, 
                                const std::vector<double>& gradx_tau1, 
                                const std::vector<double>& grady_tau1, 
                                const std::vector<double>& gradz_tau1,  
                                std::vector<double>& e, 
                                std::vector<double>& vn, 
                                std::vector<double>& vs,
                                std::vector<double>& vgn,
                                std::vector<double>& vgs 
                                );

    // lda_mc函数，输入是xc_id, n, mx, my, mz，返回每个实空间格点的E和V
    static std::pair<std::vector<double>, std::vector<Matrix2x2>> lda_mc(int xc_id, const std::vector<double>& n, 
                                                                  const std::vector<double>& mx, const std::vector<double>& my, const std::vector<double>& mz);

    // lda_lc函数，输入是xc_id, n, mx, my, mz，返回每个实空间格点的E和V
    static std::pair<std::vector<double>, std::vector<Matrix2x2>> lda_lc(int xc_id, const std::vector<double>& n, 
                                                                  const std::vector<double>& mx, const std::vector<double>& my, const std::vector<double>& mz);


    // gga_mc函数，输入是xc_id, n, mx, my, mz，grad,grad2，返回每个实空间格点的E和V
    static std::pair<std::vector<double>, std::vector<Matrix2x2>> gga_mc(int xc_id, const std::vector<double>& n, 
                                                              const std::vector<double>& mx, const std::vector<double>& my, const std::vector<double>& mz,const std::vector<double> gradx_n,const std::vector<double> grady_n,const std::vector<double> gradz_n, const std::vector<double>& gradx_mx, const std::vector<double>& grady_mx, const std::vector<double>& gradz_mx, const std::vector<double>& gradx_my, const std::vector<double>& grady_my, const std::vector<double>& gradz_my, const std::vector<double>& gradx_mz, const std::vector<double>& grady_mz, const std::vector<double>& gradz_mz, const std::vector<double> grad2xx_n, const std::vector<double> grad2yy_n, const std::vector<double> grad2zz_n, const std::vector<double> grad2xy_n, const std::vector<double> grad2yz_n, const std::vector<double> grad2xz_n, const std::vector<double> grad2xx_mx, const std::vector<double> grad2yy_mx, const std::vector<double> grad2zz_mx, const std::vector<double> grad2xy_mx, const std::vector<double> grad2yz_mx, const std::vector<double> grad2xz_mx, const std::vector<double> grad2xx_my, const std::vector<double> grad2yy_my, const std::vector<double> grad2zz_my, const std::vector<double> grad2xy_my, const std::vector<double> grad2yz_my, const std::vector<double> grad2xz_my, const std::vector<double> grad2xx_mz, const std::vector<double> grad2yy_mz, const std::vector<double> grad2zz_mz, const std::vector<double> grad2xy_mz, const std::vector<double> grad2yz_mz, const std::vector<double> grad2xz_mz, const std::vector<double> grad3xxx_n, const std::vector<double> grad3xxy_n, const std::vector<double> grad3xxz_n, const std::vector<double> grad3xyy_n, const std::vector<double> grad3xyz_n, const std::vector<double> grad3xzz_n, const std::vector<double> grad3yyy_n, const std::vector<double> grad3yyz_n, const std::vector<double> grad3yzz_n, const std::vector<double> grad3zzz_n, const std::vector<double> grad3xxx_mx, const std::vector<double> grad3xxy_mx, const std::vector<double> grad3xxz_mx, const std::vector<double> grad3xyy_mx, const std::vector<double> grad3xyz_mx, const std::vector<double> grad3xzz_mx, const std::vector<double> grad3yyy_mx, const std::vector<double> grad3yyz_mx, const std::vector<double> grad3yzz_mx, const std::vector<double> grad3zzz_mx, const std::vector<double> grad3xxx_my, const std::vector<double> grad3xxy_my, const std::vector<double> grad3xxz_my, const std::vector<double> grad3xyy_my, const std::vector<double> grad3xyz_my, const std::vector<double> grad3xzz_my, const std::vector<double> grad3yyy_my, const std::vector<double> grad3yyz_my, const std::vector<double> grad3yzz_my, const std::vector<double> grad3zzz_my, const std::vector<double> grad3xxx_mz, const std::vector<double> grad3xxy_mz, const std::vector<double> grad3xxz_mz, const std::vector<double> grad3xyy_mz, const std::vector<double> grad3xyz_mz, const std::vector<double> grad3xzz_mz, const std::vector<double> grad3yyy_mz, const std::vector<double> grad3yyz_mz, const std::vector<double> grad3yzz_mz, const std::vector<double> grad3zzz_mz);
// transform V and m to torque. torque is -m \times Bxc
static std::vector<double> gga_torque( 
    const std::vector<double>& mx, 
    const std::vector<double>& my, 
    const std::vector<double>& mz, 
    const std::vector<Matrix2x2>& V);

// transform V to Bxc
static std::vector<double> gga_Bxc(const std::vector<Matrix2x2> &V);


// mgga_mc，input: xc_id, n, mx, my, mz，grad rho, grad2rho, grad3rho, tau, ux, uy,uz,grad tau(u) return E and V for each real space grid point
    static std::tuple<std::vector<double>, std::vector<Matrix2x2>,std::vector<Matrix2x2>> mgga_mc(int xc_id, const std::vector<double>& n, 
    const std::vector<double>& mx, const std::vector<double>& my, const std::vector<double>& mz,const std::vector<double> gradx_n,const std::vector<double> grady_n,const std::vector<double> gradz_n, const std::vector<double>& gradx_mx, const std::vector<double>& grady_mx, const std::vector<double>& gradz_mx, const std::vector<double>& gradx_my, const std::vector<double>& grady_my, const std::vector<double>& gradz_my, const std::vector<double>& gradx_mz, const std::vector<double>& grady_mz, const std::vector<double>& gradz_mz, const std::vector<double> grad2xx_n, const std::vector<double> grad2yy_n, const std::vector<double> grad2zz_n, const std::vector<double> grad2xy_n, const std::vector<double> grad2yz_n, const std::vector<double> grad2xz_n, const std::vector<double> grad2xx_mx, const std::vector<double> grad2yy_mx, const std::vector<double> grad2zz_mx, const std::vector<double> grad2xy_mx, const std::vector<double> grad2yz_mx, const std::vector<double> grad2xz_mx, const std::vector<double> grad2xx_my, const std::vector<double> grad2yy_my, const std::vector<double> grad2zz_my, const std::vector<double> grad2xy_my, const std::vector<double> grad2yz_my, const std::vector<double> grad2xz_my, const std::vector<double> grad2xx_mz, const std::vector<double> grad2yy_mz, const std::vector<double> grad2zz_mz, const std::vector<double> grad2xy_mz, const std::vector<double> grad2yz_mz, const std::vector<double> grad2xz_mz, const std::vector<double> grad3xxx_n, const std::vector<double> grad3xxy_n, const std::vector<double> grad3xxz_n, const std::vector<double> grad3xyy_n, const std::vector<double> grad3xyz_n, const std::vector<double> grad3xzz_n, const std::vector<double> grad3yyy_n, const std::vector<double> grad3yyz_n, const std::vector<double> grad3yzz_n, const std::vector<double> grad3zzz_n, const std::vector<double> grad3xxx_mx, const std::vector<double> grad3xxy_mx, const std::vector<double> grad3xxz_mx, const std::vector<double> grad3xyy_mx, const std::vector<double> grad3xyz_mx, const std::vector<double> grad3xzz_mx, const std::vector<double> grad3yyy_mx, const std::vector<double> grad3yyz_mx, const std::vector<double> grad3yzz_mx, const std::vector<double> grad3zzz_mx, const std::vector<double> grad3xxx_my, const std::vector<double> grad3xxy_my, const std::vector<double> grad3xxz_my, const std::vector<double> grad3xyy_my, const std::vector<double> grad3xyz_my, const std::vector<double> grad3xzz_my, const std::vector<double> grad3yyy_my, const std::vector<double> grad3yyz_my, const std::vector<double> grad3yzz_my, const std::vector<double> grad3zzz_my, const std::vector<double> grad3xxx_mz, const std::vector<double> grad3xxy_mz, const std::vector<double> grad3xxz_mz, const std::vector<double> grad3xyy_mz, const std::vector<double> grad3xyz_mz, const std::vector<double> grad3xzz_mz, const std::vector<double> grad3yyy_mz, const std::vector<double> grad3yyz_mz, const std::vector<double> grad3yzz_mz, const std::vector<double> grad3zzz_mz, const std::vector<double>& tau, 
    const std::vector<double>& ux, const std::vector<double>& uy, const std::vector<double>& uz,
    const std::vector<double> gradx_tau, const std::vector<double> grady_tau, const std::vector<double> gradz_tau,
    const std::vector<double>& gradx_ux, const std::vector<double>& grady_ux, const std::vector<double>& gradz_ux,
    const std::vector<double>& gradx_uy, const std::vector<double>& grady_uy, const std::vector<double>& gradz_uy,
    const std::vector<double>& gradx_uz, const std::vector<double>& grady_uz, const std::vector<double>& gradz_uz);


    // print message and citation of the program
    static void print_NCLibxc();

    //print xc torque
    static void print_torque(const std::vector<double>& torque);

    // collinear test for gga
    static void gga_collinear_test();

    //consistency between lda_mc and lda_lc
    static void lda_mc_lc_test();

    // collinear test for lda_mc
    static void lda_mc_collinear_test();

    //xc local torque from gga
    static void gga_local_torque_test();

    // xc local torque and xc potential well-defined limit from gga
    static void gga_deri_limit();
};

// 声明 MakeAngularGrid 函数
std::vector<std::array<double, 4>> MakeAngularGrid(int grid_level);

}
#endif // NCLIBXC_H