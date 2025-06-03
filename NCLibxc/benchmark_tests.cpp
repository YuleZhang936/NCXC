#include "NCLibxc.h"
#include <iostream>
#include <iomanip>
#include <chrono>
namespace NCXC {
///////////////////////////////////////////////////////////////////////////////////
// collinear limit test for GGA
void NCLibxc::gga_collinear_test()
{
    try {
        NCLibxc::print_NCLibxc();
        
        // Sample 0-th order density
        std::vector<double> n = {1.0, 1.0, 2.0};
        std::vector<double> mx = {0.0, 0.0, 0.0};
        std::vector<double> my = {0.0, 0.0, 0.0};
        std::vector<double> mz = {0.1, -0.1, 0.1};

        // Sample gradients 
        std::vector<double> gradx_n = {0.1, 0.1, 0.1};
        std::vector<double> grady_n = {0.2, 0.2, 0.5};
        std::vector<double> gradz_n = {0.3, 0.3, 0.6};

        std::vector<double> gradx_mx = {0.0, 0.0, 0.0};
        std::vector<double> grady_mx = {0.0, 0.0, 0.0};
        std::vector<double> gradz_mx = {0.0, 0.0, 0.0};

        std::vector<double> gradx_my = {0.0, 0.0, 0.0};
        std::vector<double> grady_my = {0.0, 0.0, 0.0};
        std::vector<double> gradz_my = {0.0, 0.0, 0.0};

        std::vector<double> gradx_mz = {0.2, -0.2, 0.1};
        std::vector<double> grady_mz = {0.3, -0.3, 0.1};
        std::vector<double> gradz_mz = {0.1, -0.1, 0.1};

        // Second derivatives 
        std::vector<double> grad2xx_n = {0.5, 0.5, 0.0};
        std::vector<double> grad2yy_n = {0.2, 0.2, 0.0};
        std::vector<double> grad2zz_n = {0.3, 0.3, 0.0};
        std::vector<double> grad2xy_n = {0.7, 0.7, 0.0};
        std::vector<double> grad2yz_n = {0.8, 0.8, 0.0};
        std::vector<double> grad2xz_n = {0.9, 0.9, 0.0};
        std::vector<double> grad2xx_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad2yy_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad2zz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad2xy_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad2yz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad2xz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad2xx_my = {0.0, 0.0, 0.0};
        std::vector<double> grad2yy_my = {0.0, 0.0, 0.0};
        std::vector<double> grad2zz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad2xy_my = {0.0, 0.0, 0.0};
        std::vector<double> grad2yz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad2xz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad2xx_mz = {0.01, -0.01, 0.0};
        std::vector<double> grad2yy_mz = {0.02, -0.02, 0.0};
        std::vector<double> grad2zz_mz = {0.03, -0.03, 0.0};
        std::vector<double> grad2xy_mz = {0.06, -0.06, 0.0};
        std::vector<double> grad2yz_mz = {0.07, -0.07, 0.0};
        std::vector<double> grad2xz_mz = {0.08, -0.08, 0.0};

        // Third-order derivatives for n
        std::vector<double> grad3xxx_n = {0.1, 0.1, 0.2};
        std::vector<double> grad3xxy_n = {0.1, 0.1, 0.2};
        std::vector<double> grad3xxz_n = {0.1, 0.1, 0.2};
        std::vector<double> grad3xyy_n = {0.1, 0.1, 0.2};
        std::vector<double> grad3xyz_n = {0.1, 0.1, 0.2};
        std::vector<double> grad3xzz_n = {0.1, 0.1, 0.2};
        std::vector<double> grad3yyy_n = {0.1, 0.1, 0.2};
        std::vector<double> grad3yyz_n = {0.1, 0.1, 0.2};
        std::vector<double> grad3yzz_n = {0.1, 0.1, 0.2};
        std::vector<double> grad3zzz_n = {0.1, 0.1, 0.2};

        // Third-order derivatives for mx (all zeros)
        std::vector<double> grad3xxx_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad3xxy_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad3xxz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad3xyy_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad3xyz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad3xzz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad3yyy_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad3yyz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad3yzz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad3zzz_mx = {0.0, 0.0, 0.0};

        // Third-order derivatives for my (all zeros)
        std::vector<double> grad3xxx_my = {0.0, 0.0, 0.0};
        std::vector<double> grad3xxy_my = {0.0, 0.0, 0.0};
        std::vector<double> grad3xxz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad3xyy_my = {0.0, 0.0, 0.0};
        std::vector<double> grad3xyz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad3xzz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad3yyy_my = {0.0, 0.0, 0.0};
        std::vector<double> grad3yyz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad3yzz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad3zzz_my = {0.0, 0.0, 0.0};

        // Third-order derivatives for mz
        std::vector<double> grad3xxx_mz = {0.05, -0.05, 0.05};
        std::vector<double> grad3xxy_mz = {0.05, -0.05, 0.05};
        std::vector<double> grad3xxz_mz = {0.05, -0.05, 0.05};
        std::vector<double> grad3xyy_mz = {0.05, -0.05, 0.05};
        std::vector<double> grad3xyz_mz = {0.05, -0.05, 0.05};
        std::vector<double> grad3xzz_mz = {0.05, -0.05, 0.05};
        std::vector<double> grad3yyy_mz = {0.05, -0.05, 0.05};
        std::vector<double> grad3yyz_mz = {0.05, -0.05, 0.05};
        std::vector<double> grad3yzz_mz = {0.05, -0.05, 0.05};
        std::vector<double> grad3zzz_mz = {0.05, -0.05, 0.05};

        int xc_id = 106; // Example XC functional ID

        // Call gga_mc
        auto [E_GGA_MC, V_GGA_MC] = NCLibxc::gga_mc(
    xc_id, n, mx, my, mz,
    gradx_n, grady_n, gradz_n,
    gradx_mx, grady_mx, gradz_mx,
    gradx_my, grady_my, gradz_my,
    gradx_mz, grady_mz, gradz_mz,
    grad2xx_n, grad2yy_n, grad2zz_n, grad2xy_n, grad2yz_n, grad2xz_n,
    grad2xx_mx, grad2yy_mx, grad2zz_mx, grad2xy_mx, grad2yz_mx, grad2xz_mx,
    grad2xx_my, grad2yy_my, grad2zz_my, grad2xy_my, grad2yz_my, grad2xz_my,
    grad2xx_mz, grad2yy_mz, grad2zz_mz, grad2xy_mz, grad2yz_mz, grad2xz_mz,
    grad3xxx_n, grad3xxy_n, grad3xxz_n, grad3xyy_n, grad3xyz_n, grad3xzz_n,
    grad3yyy_n, grad3yyz_n, grad3yzz_n, grad3zzz_n,
    grad3xxx_mx, grad3xxy_mx, grad3xxz_mx, grad3xyy_mx, grad3xyz_mx, grad3xzz_mx,
    grad3yyy_mx, grad3yyz_mx, grad3yzz_mx, grad3zzz_mx,
    grad3xxx_my, grad3xxy_my, grad3xxz_my, grad3xyy_my, grad3xyz_my, grad3xzz_my,
    grad3yyy_my, grad3yyz_my, grad3yzz_my, grad3zzz_my,
    grad3xxx_mz, grad3xxy_mz, grad3xxz_mz, grad3xyy_mz, grad3xyz_mz, grad3xzz_mz,
    grad3yyy_mz, grad3yyz_mz, grad3yzz_mz, grad3zzz_mz
);

        // Print results from gga_mc
        std::cout << std::fixed << std::setprecision(8);
        std::cout << "Total E for each real-space grid point (GGA MC):" << std::endl;
        for (const auto &e : E_GGA_MC)
            std::cout << e << " ";
        std::cout << std::endl;
        std::cout << "Total V for each real-space grid point (GGA MC):" << std::endl;
        std::cout << "[ ";
        for (const auto &mat : V_GGA_MC) {
            std::cout << "[ ";
            std::cout << mat[0][0] << " " << mat[0][1] << " "
                      << mat[1][0] << " " << mat[1][1];
            std::cout << "] ";
        }
        std::cout << "]" << std::endl;
        // Prepare inputs for postlibxc_gga
        size_t size = n.size();
        std::vector<double> rho0(size), rho1(size);
        for (size_t i = 0; i < size; ++i) {
            rho0[i] = (n[i] + mz[i]) / 2.0;
            rho1[i] = (n[i] - mz[i]) / 2.0;
        }

        // Gradient calculations
        // First derivatives
        std::vector<double> gradx_rho0(size), grady_rho0(size), gradz_rho0(size);
        std::vector<double> gradx_rho1(size), grady_rho1(size), gradz_rho1(size);

        for (size_t i = 0; i < size; ++i) {
            gradx_rho0[i] = (gradx_n[i] + gradx_mz[i]) / 2.0;
            grady_rho0[i] = (grady_n[i] + grady_mz[i]) / 2.0;
            gradz_rho0[i] = (gradz_n[i] + gradz_mz[i]) / 2.0;

            gradx_rho1[i] = (gradx_n[i] - gradx_mz[i]) / 2.0;
            grady_rho1[i] = (grady_n[i] - grady_mz[i]) / 2.0;
            gradz_rho1[i] = (gradz_n[i] - gradz_mz[i]) / 2.0;
        }

        // Second derivatives
        std::vector<double> grad2xx_rho0(size), grad2yy_rho0(size), grad2zz_rho0(size);
        std::vector<double> grad2xy_rho0(size), grad2yz_rho0(size), grad2xz_rho0(size);
        std::vector<double> grad2xx_rho1(size), grad2yy_rho1(size), grad2zz_rho1(size);
        std::vector<double> grad2xy_rho1(size), grad2yz_rho1(size), grad2xz_rho1(size);

        for (size_t i = 0; i < size; ++i) {
            grad2xx_rho0[i] = (grad2xx_n[i] + grad2xx_mz[i]) / 2.0;
            grad2yy_rho0[i] = (grad2yy_n[i] + grad2yy_mz[i]) / 2.0;
            grad2zz_rho0[i] = (grad2zz_n[i] + grad2zz_mz[i]) / 2.0;
            grad2xy_rho0[i] = (grad2xy_n[i] + grad2xy_mz[i]) / 2.0;
            grad2yz_rho0[i] = (grad2yz_n[i] + grad2yz_mz[i]) / 2.0;
            grad2xz_rho0[i] = (grad2xz_n[i] + grad2xz_mz[i]) / 2.0;

            grad2xx_rho1[i] = (grad2xx_n[i] - grad2xx_mz[i]) / 2.0;
            grad2yy_rho1[i] = (grad2yy_n[i] - grad2yy_mz[i]) / 2.0;
            grad2zz_rho1[i] = (grad2zz_n[i] - grad2zz_mz[i]) / 2.0;
            grad2xy_rho1[i] = (grad2xy_n[i] - grad2xy_mz[i]) / 2.0;
            grad2yz_rho1[i] = (grad2yz_n[i] - grad2yz_mz[i]) / 2.0;
            grad2xz_rho1[i] = (grad2xz_n[i] - grad2xz_mz[i]) / 2.0;
        }

        // Third-order gradients
        std::vector<double> grad3xxx_rho0(size), grad3xxy_rho0(size), grad3xxz_rho0(size);
        std::vector<double> grad3xyy_rho0(size), grad3xyz_rho0(size), grad3xzz_rho0(size);
        std::vector<double> grad3yyy_rho0(size), grad3yyz_rho0(size), grad3yzz_rho0(size);
        std::vector<double> grad3zzz_rho0(size);

        std::vector<double> grad3xxx_rho1(size), grad3xxy_rho1(size), grad3xxz_rho1(size);
        std::vector<double> grad3xyy_rho1(size), grad3xyz_rho1(size), grad3xzz_rho1(size);
        std::vector<double> grad3yyy_rho1(size), grad3yyz_rho1(size), grad3yzz_rho1(size);
        std::vector<double> grad3zzz_rho1(size);

        for (size_t i = 0; i < size; ++i) {
            // Calculate third-order gradients for rho0
            grad3xxx_rho0[i] = (grad3xxx_n[i] + grad3xxx_mz[i]) / 2.0;
            grad3xxy_rho0[i] = (grad3xxy_n[i] + grad3xxy_mz[i]) / 2.0;
            grad3xxz_rho0[i] = (grad3xxz_n[i] + grad3xxz_mz[i]) / 2.0;
            grad3xyy_rho0[i] = (grad3xyy_n[i] + grad3xyy_mz[i]) / 2.0;
            grad3xyz_rho0[i] = (grad3xyz_n[i] + grad3xyz_mz[i]) / 2.0;
            grad3xzz_rho0[i] = (grad3xzz_n[i] + grad3xzz_mz[i]) / 2.0;
            grad3yyy_rho0[i] = (grad3yyy_n[i] + grad3yyy_mz[i]) / 2.0;
            grad3yyz_rho0[i] = (grad3yyz_n[i] + grad3yyz_mz[i]) / 2.0;
            grad3yzz_rho0[i] = (grad3yzz_n[i] + grad3yzz_mz[i]) / 2.0;
            grad3zzz_rho0[i] = (grad3zzz_n[i] + grad3zzz_mz[i]) / 2.0;

            // Calculate third-order gradients for rho1
            grad3xxx_rho1[i] = (grad3xxx_n[i] - grad3xxx_mz[i]) / 2.0;
            grad3xxy_rho1[i] = (grad3xxy_n[i] - grad3xxy_mz[i]) / 2.0;
            grad3xxz_rho1[i] = (grad3xxz_n[i] - grad3xxz_mz[i]) / 2.0;
            grad3xyy_rho1[i] = (grad3xyy_n[i] - grad3xyy_mz[i]) / 2.0;
            grad3xyz_rho1[i] = (grad3xyz_n[i] - grad3xyz_mz[i]) / 2.0;
            grad3xzz_rho1[i] = (grad3xzz_n[i] - grad3xzz_mz[i]) / 2.0;
            grad3yyy_rho1[i] = (grad3yyy_n[i] - grad3yyy_mz[i]) / 2.0;
            grad3yyz_rho1[i] = (grad3yyz_n[i] - grad3yyz_mz[i]) / 2.0;
            grad3yzz_rho1[i] = (grad3yzz_n[i] - grad3yzz_mz[i]) / 2.0;
            grad3zzz_rho1[i] = (grad3zzz_n[i] - grad3zzz_mz[i]) / 2.0;
        }

     

        // Call postlibxc_gga
        std::vector<double> e, v1, v2, f1, f2, f3;
        NCLibxc::postlibxc_gga(
    xc_id,
    rho0, rho1,
    gradx_rho0, grady_rho0, gradz_rho0,
    gradx_rho1, grady_rho1, gradz_rho1,
    grad2xx_rho0, grad2yy_rho0, grad2zz_rho0,
    grad2xy_rho0, grad2yz_rho0, grad2xz_rho0,
    grad2xx_rho1, grad2yy_rho1, grad2zz_rho1,
    grad2xy_rho1, grad2yz_rho1, grad2xz_rho1,
    grad3xxx_rho0, grad3xxy_rho0, grad3xxz_rho0,
    grad3xyy_rho0, grad3xyz_rho0, grad3xzz_rho0,
    grad3yyy_rho0, grad3yyz_rho0, grad3yzz_rho0,
    grad3zzz_rho0,
    grad3xxx_rho1, grad3xxy_rho1, grad3xxz_rho1,
    grad3xyy_rho1, grad3xyz_rho1, grad3xzz_rho1,
    grad3yyy_rho1, grad3yyz_rho1, grad3yzz_rho1,
    grad3zzz_rho1,
    e, v1, v2, f1, f2, f3
);

        // Print results from postlibxc_gga
       std::cout << "Total E for each real-space grid point (postlibxc_gga):" << std::endl;
        for (const auto &e : e)
            std::cout << e << " ";
        std::cout << std::endl;

        // Similarly, print v1, v2, f1, f2, f3 if needed
        std::cout << "Total V1 for each real-space grid point (postlibxc_gga):" << std::endl;
        for (const auto &v : v1)
            std::cout << v << " ";
        std::cout << std::endl;

        std::cout << "Total V2 for each real-space grid point (postlibxc_gga):" << std::endl;
        for (const auto &v : v2)
            std::cout << v << " ";
        std::cout << std::endl;

    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
    }
}


///////////////////////////////////////////////////////////////////////////////////
//consistency between lda_mc and lda_lc
void NCLibxc::lda_mc_lc_test()
{
try {
        NCLibxc::print_NCLibxc();
        // example input data
        std::vector<double> n = {1.0, 1.0, 1.0};
        std::vector<double> mx = {0.1, 0.1, 0.0};
        std::vector<double> my = {0.0, 0.1, 0.1414};
        std::vector<double> mz = {0.1, 0.0, 0.0};
        int xc_id = 1; // for example, set xc_id to 1(lda)

        auto [E_MC, V_MC] = NCLibxc::lda_mc(xc_id, n, mx, my, mz);

        std::cout << "Total E for each real-space grid point(LDA_MC):" << std::endl;
        for (const auto &e : E_MC)
            std::cout << e << " ";
        std::cout << std::endl;

        std::cout << "Total V for each real-space grid point(LDA_MC):" << std::endl;
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

        std::cout << "Total E for each real-space grid point(LDA_LC):" << std::endl;
        for (const auto &e : E_LC)
            std::cout << e << " ";
        std::cout << std::endl;

        std::cout << "Total V for each real-space grid point(LDA_LC):" << std::endl;
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
}

//consistency between lda_mc and postlibxc_lda(collinear)
void NCLibxc::lda_mc_collinear_test()
{
try {
        NCLibxc::print_NCLibxc();
        // example input data
        std::vector<double> n = {1.0, 1.0, 1.0};
        std::vector<double> mx = {0.0, 0.0, 0.0};
        std::vector<double> my = {0.0, 0.0, 0.0};
        std::vector<double> mz = {0.1, 0.2, 0.3};
        int xc_id = 1; // for example, set xc_id to 1(lda)

        auto [E_MC, V_MC] = NCLibxc::lda_mc(xc_id, n, mx, my, mz);
        
        std::cout << "Total E for each real-space grid point(LDA_MC):" << std::endl;
        for (const auto &e : E_MC)
            std::cout << e << " ";
        std::cout << std::endl;

        std::cout << "Total V for each real-space grid point(LDA_MC):" << std::endl;
        for (const auto &matrix : V_MC) {
            for (const auto &row : matrix) {
                for (const auto &elem : row) {
                    std::cout << elem << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

    
    std::vector<double> e;
    std::vector<double> v1, v2, f1, f2, f3;
    
    size_t np = n.size();

    std::vector<double> rho_up(np,0.0) ;     
    std::vector<double> rho_down(np,0.0);

    for(size_t i=0;i<np;i++){
        rho_up[i] = (n[i] + mz[i]) / 2.0;
        rho_down[i] = (n[i] - mz[i]) / 2.0;
    }   

    NCLibxc::postlibxc_lda(xc_id, rho_up, rho_down, e, v1, v2, f1, f2, f3);

    std::cout << "Results from postlibxc_lda:" << std::endl;

    std::cout << "e:" << std::endl;
    for (const auto &val : e)
        std::cout << val << " ";
    std::cout << std::endl;

    std::cout << "v1:" << std::endl;
    for (const auto &val : v1)
        std::cout << val << " ";
    std::cout << std::endl;

    std::cout << "v2:" << std::endl;
    for (const auto &val : v2)
        std::cout << val << " ";
    std::cout << std::endl;

        
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
    }
}

// xc local torque from gga, prove non-vanishing torque
void NCLibxc::gga_local_torque_test()
{
    
    std::vector<double> n  = {1.0, 1.0, 1.0};
    std::vector<double> mx = {0.1, 0.2, 0.3};
    std::vector<double> my = {0.2, 0.3, 0.4};
    std::vector<double> mz = {0.3, 0.4, 0.5};

    // 各阶梯度均采用3个元素的零初始数组（实际计算中请设置合适值）
    std::vector<double> gradx_n(3, 0.1), grady_n(3, 0.0), gradz_n(3, 0.0);
    std::vector<double> gradx_mx(3, 0.0), grady_mx(3, 0.3), gradz_mx(3, 0.0);
    std::vector<double> gradx_my(3, 0.1), grady_my(3, 0.0), gradz_my(3, 0.0);
    std::vector<double> gradx_mz(3, 0.0), grady_mz(3, 0.0), gradz_mz(3, 0.1);

    std::vector<double> grad2xx_n(3, 0.0), grad2yy_n(3, 0.0), grad2zz_n(3, 0.0);
    std::vector<double> grad2xy_n(3, 0.0), grad2yz_n(3, 0.0), grad2xz_n(3, 0.0);
    std::vector<double> grad2xx_mx(3, 0.0), grad2yy_mx(3, 0.0), grad2zz_mx(3, 0.0);
    std::vector<double> grad2xy_mx(3, 0.0), grad2yz_mx(3, 0.0), grad2xz_mx(3, 0.0);
    std::vector<double> grad2xx_my(3, 0.0), grad2yy_my(3, 0.0), grad2zz_my(3, 0.0);
    std::vector<double> grad2xy_my(3, 0.0), grad2yz_my(3, 0.0), grad2xz_my(3, 0.0);
    std::vector<double> grad2xx_mz(3, 0.0), grad2yy_mz(3, 0.0), grad2zz_mz(3, 0.0);
    std::vector<double> grad2xy_mz(3, 0.0), grad2yz_mz(3, 0.0), grad2xz_mz(3, 0.0);

    std::vector<double> grad3xxx_n(3, 0.0), grad3xxy_n(3, 0.0), grad3xxz_n(3, 0.0), 
                        grad3xyy_n(3, 0.0), grad3xyz_n(3, 0.0), grad3xzz_n(3, 0.0);
    std::vector<double> grad3yyy_n(3, 0.0), grad3yyz_n(3, 0.0), grad3yzz_n(3, 0.0), 
                        grad3zzz_n(3, 0.0);
    std::vector<double> grad3xxx_mx(3, 0.0), grad3xxy_mx(3, 0.0), grad3xxz_mx(3, 0.0),
                        grad3xyy_mx(3, 0.0), grad3xyz_mx(3, 0.0), grad3xzz_mx(3, 0.0);
    std::vector<double> grad3yyy_mx(3, 0.0), grad3yyz_mx(3, 0.0), grad3yzz_mx(3, 0.0),
                        grad3zzz_mx(3, 0.0);
    std::vector<double> grad3xxx_my(3, 0.0), grad3xxy_my(3, 0.0), grad3xxz_my(3, 0.0),
                        grad3xyy_my(3, 0.0), grad3xyz_my(3, 0.0), grad3xzz_my(3, 0.0);
    std::vector<double> grad3yyy_my(3, 0.0), grad3yyz_my(3, 0.0), grad3yzz_my(3, 0.0),
                        grad3zzz_my(3, 0.0);
    std::vector<double> grad3xxx_mz(3, 0.0), grad3xxy_mz(3, 0.0), grad3xxz_mz(3, 0.0),
                        grad3xyy_mz(3, 0.0), grad3xyz_mz(3, 0.0), grad3xzz_mz(3, 0.0);
    std::vector<double> grad3yyy_mz(3, 0.0), grad3yyz_mz(3, 0.0), grad3yzz_mz(3, 0.0),
                        grad3zzz_mz(3, 0.0);


    int xc_id = 106;  

    auto [E_MC, V_MC] = NCLibxc::gga_mc(
        xc_id, n, mx, my, mz,
        gradx_n, grady_n, gradz_n,
        gradx_mx, grady_mx, gradz_mx,
        gradx_my, grady_my, gradz_my,
        gradx_mz, grady_mz, gradz_mz,
        grad2xx_n, grad2yy_n, grad2zz_n, grad2xy_n, grad2yz_n, grad2xz_n,
        grad2xx_mx, grad2yy_mx, grad2zz_mx, grad2xy_mx, grad2yz_mx, grad2xz_mx,
        grad2xx_my, grad2yy_my, grad2zz_my, grad2xy_my, grad2yz_my, grad2xz_my,
        grad2xx_mz, grad2yy_mz, grad2zz_mz, grad2xy_mz, grad2yz_mz, grad2xz_mz,
        grad3xxx_n, grad3xxy_n, grad3xxz_n, grad3xyy_n, grad3xyz_n, grad3xzz_n,
        grad3yyy_n, grad3yyz_n, grad3yzz_n, grad3zzz_n,
        grad3xxx_mx, grad3xxy_mx, grad3xxz_mx, grad3xyy_mx, grad3xyz_mx, grad3xzz_mx,
        grad3yyy_mx, grad3yyz_mx, grad3yzz_mx, grad3zzz_mx,
        grad3xxx_my, grad3xxy_my, grad3xxz_my, grad3xyy_my, grad3xyz_my, grad3xzz_my,
        grad3yyy_my, grad3yyz_my, grad3yzz_my, grad3zzz_my,
        grad3xxx_mz, grad3xxy_mz, grad3xxz_mz, grad3xyy_mz, grad3xyz_mz, grad3xzz_mz,
        grad3yyy_mz, grad3yyz_mz, grad3yzz_mz, grad3zzz_mz
    );

    auto torque = NCLibxc::gga_torque(mx,my,mz,V_MC);

    NCLibxc::print_torque(torque);
}

// xc local torque and xc potential from gga, prove well-defined torque and potential
void NCLibxc::gga_deri_limit()
{
    // We now vary lambda from 1.0 down to 1e-20
    int xc_id = 106;
    std::complex<double> twoi(0.0, 2.0);
    std::complex<double> two(2.0, 0.0);
    for (int i = 0; i <= 20; ++i)
    {
        double lambda = std::pow(10.0, -static_cast<double>(i));

        // We only need one point
        std::vector<double> n (1, 2.0);
        std::vector<double> mx(1, lambda);
        std::vector<double> my(1, 0.1*lambda);
        std::vector<double> mz(1, 0.01 * lambda);

        // All gradients/derivatives set to 0.01
        std::vector<double> gradx_n(1, 0.01), grady_n(1, 0.01), gradz_n(1, 0.01);
        std::vector<double> gradx_mx(1, 0.01), grady_mx(1, 0.01), gradz_mx(1, 0.01);
        std::vector<double> gradx_my(1, 0.01), grady_my(1, 0.01), gradz_my(1, 0.01);
        std::vector<double> gradx_mz(1, 0.01), grady_mz(1, 0.01), gradz_mz(1, 0.01);

        std::vector<double> grad2xx_n(1, 0.01), grad2yy_n(1, 0.01), grad2zz_n(1, 0.01);
        std::vector<double> grad2xy_n(1, 0.01), grad2yz_n(1, 0.01), grad2xz_n(1, 0.01);
        std::vector<double> grad2xx_mx(1, 0.01), grad2yy_mx(1, 0.01), grad2zz_mx(1, 0.01);
        std::vector<double> grad2xy_mx(1, 0.01), grad2yz_mx(1, 0.01), grad2xz_mx(1, 0.01);
        std::vector<double> grad2xx_my(1, 0.01), grad2yy_my(1, 0.01), grad2zz_my(1, 0.01);
        std::vector<double> grad2xy_my(1, 0.01), grad2yz_my(1, 0.01), grad2xz_my(1, 0.01);
        std::vector<double> grad2xx_mz(1, 0.01), grad2yy_mz(1, 0.01), grad2zz_mz(1, 0.01);
        std::vector<double> grad2xy_mz(1, 0.01), grad2yz_mz(1, 0.01), grad2xz_mz(1, 0.01);

        std::vector<double> grad3xxx_n(1, 0.01), grad3xxy_n(1, 0.01), grad3xxz_n(1, 0.01);
        std::vector<double> grad3xyy_n(1, 0.01), grad3xyz_n(1, 0.01), grad3xzz_n(1, 0.01);
        std::vector<double> grad3yyy_n(1, 0.01), grad3yyz_n(1, 0.01), grad3yzz_n(1, 0.01);
        std::vector<double> grad3zzz_n(1, 0.01);

        std::vector<double> grad3xxx_mx(1, 0.01), grad3xxy_mx(1, 0.01), grad3xxz_mx(1, 0.01);
        std::vector<double> grad3xyy_mx(1, 0.01), grad3xyz_mx(1, 0.01), grad3xzz_mx(1, 0.01);
        std::vector<double> grad3yyy_mx(1, 0.01), grad3yyz_mx(1, 0.01), grad3yzz_mx(1, 0.01);
        std::vector<double> grad3zzz_mx(1, 0.01);

        std::vector<double> grad3xxx_my(1, 0.01), grad3xxy_my(1, 0.01), grad3xxz_my(1, 0.01);
        std::vector<double> grad3xyy_my(1, 0.01), grad3xyz_my(1, 0.01), grad3xzz_my(1, 0.01);
        std::vector<double> grad3yyy_my(1, 0.01), grad3yyz_my(1, 0.01), grad3yzz_my(1, 0.01);
        std::vector<double> grad3zzz_my(1, 0.01);

        std::vector<double> grad3xxx_mz(1, 0.01), grad3xxy_mz(1, 0.01), grad3xxz_mz(1, 0.01);
        std::vector<double> grad3xyy_mz(1, 0.01), grad3xyz_mz(1, 0.01), grad3xzz_mz(1, 0.01);
        std::vector<double> grad3yyy_mz(1, 0.01), grad3yyz_mz(1, 0.01), grad3yzz_mz(1, 0.01);
        std::vector<double> grad3zzz_mz(1, 0.01);


        // Compute E and V using gga_mc
        auto [energy, potential] = gga_mc(
            xc_id, n, mx, my, mz,
            gradx_n, grady_n, gradz_n,
            gradx_mx, grady_mx, gradz_mx,
            gradx_my, grady_my, gradz_my,
            gradx_mz, grady_mz, gradz_mz,
            grad2xx_n, grad2yy_n, grad2zz_n, grad2xy_n, grad2yz_n, grad2xz_n,
            grad2xx_mx, grad2yy_mx, grad2zz_mx, grad2xy_mx, grad2yz_mx, grad2xz_mx,
            grad2xx_my, grad2yy_my, grad2zz_my, grad2xy_my, grad2yz_my, grad2xz_my,
            grad2xx_mz, grad2yy_mz, grad2zz_mz, grad2xy_mz, grad2yz_mz, grad2xz_mz,
            grad3xxx_n, grad3xxy_n, grad3xxz_n, grad3xyy_n, grad3xyz_n, grad3xzz_n,
            grad3yyy_n, grad3yyz_n, grad3yzz_n, grad3zzz_n,
            grad3xxx_mx, grad3xxy_mx, grad3xxz_mx, grad3xyy_mx, grad3xyz_mx, grad3xzz_mx,
            grad3yyy_mx, grad3yyz_mx, grad3yzz_mx, grad3zzz_mx,
            grad3xxx_my, grad3xxy_my, grad3xxz_my, grad3xyy_my, grad3xyz_my, grad3xzz_my,
            grad3yyy_my, grad3yyz_my, grad3yzz_my, grad3zzz_my,
            grad3xxx_mz, grad3xxy_mz, grad3xxz_mz, grad3xyy_mz, grad3xyz_mz, grad3xzz_mz,
            grad3yyy_mz, grad3yyz_mz, grad3yzz_mz, grad3zzz_mz
        );

        auto torque = gga_torque(mx, my, mz, potential);

        // Print results
        std::cout << std::fixed << std::setprecision(30);
        std::cout << "For lambda=" << lambda << ":\n";

        // Torque
        std::cout << "Torque: ("
                  << torque[0] << ", " 
                  << torque[1] << ", "
                  << torque[2] << ")" << std::endl;

        // Energy
        std::cout << "E = " << n[0]*energy[0] << std::endl;
        

        // Potential V: 2x2 matrix
        const auto& mat = potential[0];
        double V0 = std::real((mat[0][0]+mat[1][1])/two);
        double V1 = std::real((mat[0][1]+mat[1][0])/two);
        double V2 = std::real((mat[1][0]-mat[0][1])/twoi); 
        double V3 = std::real((mat[0][0] - mat[1][1]) / two);
        std::cout << "DE/Dn =" << V0
                    << ", DE/Dmx = " << V1
                    << ", DE/Dmy = " << V2
                    << ", DE/Dmz = " << V3 << std::endl;
        

        std::cout << std::endl;
    }
}

// test for gga_mc with large scale data
void NCLibxc::gga_mc_large_scale_test()
{
    try
    {
        std::cout << "Running large-scale GGA MC test with 100,000 grid points..." << std::endl;
        
        // Set size to 100,000 elements
        const size_t size = 100000;
        
        // Sample 0-th order density (all n = 1.0, m = 0.1)
        std::vector<double> n(size, 1.0);
        std::vector<double> mx(size, 0.1);
        std::vector<double> my(size, 0.1);
        std::vector<double> mz(size, 0.1);

        // First derivatives (gradients) - all set to 0.1
        std::vector<double> gradx_n(size, 1);
        std::vector<double> grady_n(size, 1);
        std::vector<double> gradz_n(size, 1);

        std::vector<double> gradx_mx(size, 0.1);
        std::vector<double> grady_mx(size, 0.1);
        std::vector<double> gradz_mx(size, 0.1);

        std::vector<double> gradx_my(size, 0.1);
        std::vector<double> grady_my(size, 0.1);
        std::vector<double> gradz_my(size, 0.1);

        std::vector<double> gradx_mz(size, 0.1);
        std::vector<double> grady_mz(size, 0.1);
        std::vector<double> gradz_mz(size, 0.1);

        // Second derivatives - all set to 0.1
        std::vector<double> grad2xx_n(size, 1);
        std::vector<double> grad2yy_n(size, 1);
        std::vector<double> grad2zz_n(size, 1);
        std::vector<double> grad2xy_n(size, 1);
        std::vector<double> grad2yz_n(size, 1);
        std::vector<double> grad2xz_n(size, 1);
        
        std::vector<double> grad2xx_mx(size, 0.1);
        std::vector<double> grad2yy_mx(size, 0.1);
        std::vector<double> grad2zz_mx(size, 0.1);
        std::vector<double> grad2xy_mx(size, 0.1);
        std::vector<double> grad2yz_mx(size, 0.1);
        std::vector<double> grad2xz_mx(size, 0.1);
        
        std::vector<double> grad2xx_my(size, 0.1);
        std::vector<double> grad2yy_my(size, 0.1);
        std::vector<double> grad2zz_my(size, 0.1);
        std::vector<double> grad2xy_my(size, 0.1);
        std::vector<double> grad2yz_my(size, 0.1);
        std::vector<double> grad2xz_my(size, 0.1);
        
        std::vector<double> grad2xx_mz(size, 0.1);
        std::vector<double> grad2yy_mz(size, 0.1);
        std::vector<double> grad2zz_mz(size, 0.1);
        std::vector<double> grad2xy_mz(size, 0.1);
        std::vector<double> grad2yz_mz(size, 0.1);
        std::vector<double> grad2xz_mz(size, 0.1);

        // Third-order derivatives - all set to 0.1
        std::vector<double> grad3xxx_n(size, 1);
        std::vector<double> grad3xxy_n(size, 1);
        std::vector<double> grad3xxz_n(size, 1);
        std::vector<double> grad3xyy_n(size, 1);
        std::vector<double> grad3xyz_n(size, 1);
        std::vector<double> grad3xzz_n(size, 1);
        std::vector<double> grad3yyy_n(size, 1);
        std::vector<double> grad3yyz_n(size, 1);
        std::vector<double> grad3yzz_n(size, 1);
        std::vector<double> grad3zzz_n(size, 1);

        std::vector<double> grad3xxx_mx(size, 0.1);
        std::vector<double> grad3xxy_mx(size, 0.1);
        std::vector<double> grad3xxz_mx(size, 0.1);
        std::vector<double> grad3xyy_mx(size, 0.1);
        std::vector<double> grad3xyz_mx(size, 0.1);
        std::vector<double> grad3xzz_mx(size, 0.1);
        std::vector<double> grad3yyy_mx(size, 0.1);
        std::vector<double> grad3yyz_mx(size, 0.1);
        std::vector<double> grad3yzz_mx(size, 0.1);
        std::vector<double> grad3zzz_mx(size, 0.1);

        std::vector<double> grad3xxx_my(size, 0.1);
        std::vector<double> grad3xxy_my(size, 0.1);
        std::vector<double> grad3xxz_my(size, 0.1);
        std::vector<double> grad3xyy_my(size, 0.1);
        std::vector<double> grad3xyz_my(size, 0.1);
        std::vector<double> grad3xzz_my(size, 0.1);
        std::vector<double> grad3yyy_my(size, 0.1);
        std::vector<double> grad3yyz_my(size, 0.1);
        std::vector<double> grad3yzz_my(size, 0.1);
        std::vector<double> grad3zzz_my(size, 0.1);

        std::vector<double> grad3xxx_mz(size, 0.1);
        std::vector<double> grad3xxy_mz(size, 0.1);
        std::vector<double> grad3xxz_mz(size, 0.1);
        std::vector<double> grad3xyy_mz(size, 0.1);
        std::vector<double> grad3xyz_mz(size, 0.1);
        std::vector<double> grad3xzz_mz(size, 0.1);
        std::vector<double> grad3yyy_mz(size, 0.1);
        std::vector<double> grad3yyz_mz(size, 0.1);
        std::vector<double> grad3yzz_mz(size, 0.1);
        std::vector<double> grad3zzz_mz(size, 0.1);

        int xc_id = 106; // Example XC functional ID

        // Record start time
        auto start_time = std::chrono::high_resolution_clock::now();
        
        std::cout << "Starting gga_mc calculation with " << size << " points..." << std::endl;

        // Call gga_mc
        auto [E_GGA_MC, V_GGA_MC] = NCLibxc::gga_mc(
            xc_id, n, mx, my, mz,
            gradx_n, grady_n, gradz_n,
            gradx_mx, grady_mx, gradz_mx,
            gradx_my, grady_my, gradz_my,
            gradx_mz, grady_mz, gradz_mz,
            grad2xx_n, grad2yy_n, grad2zz_n, grad2xy_n, grad2yz_n, grad2xz_n,
            grad2xx_mx, grad2yy_mx, grad2zz_mx, grad2xy_mx, grad2yz_mx, grad2xz_mx,
            grad2xx_my, grad2yy_my, grad2zz_my, grad2xy_my, grad2yz_my, grad2xz_my,
            grad2xx_mz, grad2yy_mz, grad2zz_mz, grad2xy_mz, grad2yz_mz, grad2xz_mz,
            grad3xxx_n, grad3xxy_n, grad3xxz_n, grad3xyy_n, grad3xyz_n, grad3xzz_n,
            grad3yyy_n, grad3yyz_n, grad3yzz_n, grad3zzz_n,
            grad3xxx_mx, grad3xxy_mx, grad3xxz_mx, grad3xyy_mx, grad3xyz_mx, grad3xzz_mx,
            grad3yyy_mx, grad3yyz_mx, grad3yzz_mx, grad3zzz_mx,
            grad3xxx_my, grad3xxy_my, grad3xxz_my, grad3xyy_my, grad3xyz_my, grad3xzz_my,
            grad3yyy_my, grad3yyz_my, grad3yzz_my, grad3zzz_my,
            grad3xxx_mz, grad3xxy_mz, grad3xxz_mz, grad3xyy_mz, grad3xyz_mz, grad3xzz_mz,
            grad3yyy_mz, grad3yyz_mz, grad3yzz_mz, grad3zzz_mz
        );

        // Record end time and calculate duration
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        
        std::cout << "GGA MC calculation completed successfully!" << std::endl;
        std::cout << "Number of grid points processed: " << size << std::endl;
        std::cout << "Time taken: " << duration << " ms" << std::endl;
        std::cout << "Average processing time per point: " << static_cast<double>(duration) / size << " ms" << std::endl;

    } catch (const std::exception& ex) {
        std::cerr << "Error in gga_mc_large_scale_test: " << ex.what() << std::endl;
    }
}

void NCLibxc::gga_mc_ir_test()
{
    const int xc_id = 106;  // Functional ID
    const size_t size = 2;  // 两个样本点

    // 输入数据 - 一阶梯度
    std::vector<double> n        = {0.0206042, 0.020962};
    std::vector<double> mx       = {3.9212e-17, 3.93746e-17};
    std::vector<double> my       = {1.09293e-19, 1.08655e-19};
    std::vector<double> mz       = {0.0251652, 0.0262013};
    std::vector<double> gradx_n  = {0.0413746, 0.0416958};
    std::vector<double> grady_n  = {-8.30215e-18, 5.69725e-17};
    std::vector<double> gradz_n  = {0.00327891, 0.00164878};
    std::vector<double> gradx_mx = {8.20006e-17, 8.17718e-17};
    std::vector<double> grady_mx = {6.78554e-18, 6.14345e-18};
    std::vector<double> gradz_mx = {5.74577e-18, -3.4525e-18};
    std::vector<double> gradx_my = {3.39874e-19, 3.33983e-19};
    std::vector<double> grady_my = {6.58162e-23, 3.85928e-23};
    std::vector<double> gradz_my = {6.16462e-21, -1.41338e-20};
    std::vector<double> gradx_mz = {0.0652852, 0.0677117};
    std::vector<double> grady_mz = {2.77768e-17, -3.24048e-17};
    std::vector<double> gradz_mz = {0.00893192, 0.00511721};

    // 二阶梯度 (24个向量)
    std::vector<double> grad2xx_n   = {0.0814215, 0.0799918};
    std::vector<double> grad2yy_n   = {-0.0349828, -0.0352474};
    std::vector<double> grad2zz_n   = {-0.0128073, -0.0094527};
    std::vector<double> grad2xy_n   = {3.1792e-16, 4.54473e-16};
    std::vector<double> grad2yz_n   = {1.90782e-16, 3.12222e-16};
    std::vector<double> grad2xz_n   = {0.0037835, 0.00100217};
    
    std::vector<double> grad2xx_mx  = {1.57185e-16, 1.54371e-16};
    std::vector<double> grad2yy_mx  = {-6.91087e-17, -6.85987e-17};
    std::vector<double> grad2zz_mx  = {-6.35776e-17, -5.96898e-17};
    std::vector<double> grad2xy_mx  = {2.39442e-17, 2.1328e-17};
    std::vector<double> grad2yz_mx  = {-3.37759e-18, -5.18622e-18};
    std::vector<double> grad2xz_mx  = {1.16104e-17, -1.42126e-17};
    
    std::vector<double> grad2xx_my  = {9.45404e-19, 9.14004e-19};
    std::vector<double> grad2yy_my  = {-2.8729e-19, -2.82151e-19};
    std::vector<double> grad2zz_my  = {-1.5448e-19, -1.23271e-19};
    std::vector<double> grad2xy_my  = {3.38199e-22, 2.34028e-22};
    std::vector<double> grad2yz_my  = {-1.59456e-22, -2.05524e-22};
    std::vector<double> grad2xz_my  = {-6.93111e-21, -6.84909e-20};
    
    std::vector<double> grad2xx_mz  = {0.156385, 0.160828};
    std::vector<double> grad2yy_mz  = {-0.0551208, -0.0571468};
    std::vector<double> grad2zz_mz  = {-0.0263, -0.0254659};
    std::vector<double> grad2xy_mz  = {2.94412e-16, -1.72097e-16};
    std::vector<double> grad2yz_mz  = {-1.84602e-16, -4.82342e-16};
    std::vector<double> grad2xz_mz  = {0.0214891, 0.0117204};

    // 三阶梯度 (40个向量)
    std::vector<double> grad3xxx_n  = {0.137279, 0.12221};
    std::vector<double> grad3xxy_n  = {2.71305e-15, 7.13777e-16};
    std::vector<double> grad3xxz_n  = {-0.00824955, -0.00930039};
    std::vector<double> grad3xyy_n  = {-0.0984508, -0.0974493};
    std::vector<double> grad3xyz_n  = {5.15908e-16, 4.96079e-16};
    std::vector<double> grad3xzz_n  = {-0.0280179, -0.0105146};
    std::vector<double> grad3yyy_n  = {-3.35193e-15, -1.11822e-14};
    std::vector<double> grad3yyz_n  = {-0.00313887, -0.000807881};
    std::vector<double> grad3yzz_n  = {8.10856e-15, -5.91776e-15};
    std::vector<double> grad3zzz_n  = {0.0254741, 0.0185476};
    
    std::vector<double> grad3xxx_mx = {2.27888e-16, 2.04974e-16};
    std::vector<double> grad3xxy_mx = {7.4136e-17, 6.49588e-17};
    std::vector<double> grad3xxz_mx = {1.37304e-17, -4.99398e-17};
    std::vector<double> grad3xyy_mx = {-1.90707e-16, -1.86885e-16};
    std::vector<double> grad3xyz_mx = {-1.42255e-17, -2.06462e-17};
    std::vector<double> grad3xzz_mx = {-1.82635e-16, -1.62697e-16};
    std::vector<double> grad3yyy_mx = {-6.91616e-17, -6.00231e-17};
    std::vector<double> grad3yyz_mx = {-7.82111e-18, 1.42721e-17};
    std::vector<double> grad3yzz_mx = {-1.48981e-17, -9.45393e-18};
    std::vector<double> grad3zzz_mx = {-2.38755e-18, 5.21325e-17};
    
    std::vector<double> grad3xxx_my = {2.11451e-18, 1.92492e-18};
    std::vector<double> grad3xxy_my = {1.21908e-21, 8.72936e-22};
    std::vector<double> grad3xxz_my = {-1.33105e-19, -2.76356e-19};
    std::vector<double> grad3xyy_my = {-1.04135e-18, -1.01024e-18};
    std::vector<double> grad3xyz_my = {-6.06383e-22, -7.90474e-22};
    std::vector<double> grad3xzz_my = {-5.08053e-19, -3.32336e-19};
    std::vector<double> grad3yyy_my = {-1.02847e-21, -6.11762e-22};
    std::vector<double> grad3yyz_my = {7.10945e-21, 5.87502e-20};
    std::vector<double> grad3yzz_my = {-3.84482e-22, -2.30473e-22};
    std::vector<double> grad3zzz_my = {2.64306e-19, 1.64868e-19};
    
    std::vector<double> grad3xxx_mz = {0.30008, 0.29227};
    std::vector<double> grad3xxy_mz = {-3.91586e-15, -1.37685e-15};
    std::vector<double> grad3xxz_mz = {0.0399974, 0.0214668};
    std::vector<double> grad3xyy_mz = {-0.178534, -0.18398};
    std::vector<double> grad3xyz_mz = {-2.69306e-15, -1.82219e-15};
    std::vector<double> grad3xzz_mz = {-0.0731732, -0.059763};
    std::vector<double> grad3yyy_mz = {-1.21585e-15, 1.10801e-14};
    std::vector<double> grad3yyz_mz = {-0.0179568, -0.0097791};
    std::vector<double> grad3yzz_mz = {-3.89537e-15, 2.70304e-15};
    std::vector<double> grad3zzz_mz = {0.00700599, 0.00419113};

    // 调用 GGA_MC 泛函
    auto [E_ir, V_ir] = NCLibxc::gga_mc(
        xc_id, 
        n, mx, my, mz,
        gradx_n, grady_n, gradz_n,
        gradx_mx, grady_mx, gradz_mx,
        gradx_my, grady_my, gradz_my,
        gradx_mz, grady_mz, gradz_mz,
        // 二阶梯度
        grad2xx_n, grad2yy_n, grad2zz_n, grad2xy_n, grad2yz_n, grad2xz_n,
        grad2xx_mx, grad2yy_mx, grad2zz_mx, grad2xy_mx, grad2yz_mx, grad2xz_mx,
        grad2xx_my, grad2yy_my, grad2zz_my, grad2xy_my, grad2yz_my, grad2xz_my,
        grad2xx_mz, grad2yy_mz, grad2zz_mz, grad2xy_mz, grad2yz_mz, grad2xz_mz,
        // 三阶梯度
        grad3xxx_n, grad3xxy_n, grad3xxz_n, grad3xyy_n, grad3xyz_n, grad3xzz_n, 
        grad3yyy_n, grad3yyz_n, grad3yzz_n, grad3zzz_n,
        grad3xxx_mx, grad3xxy_mx, grad3xxz_mx, grad3xyy_mx, grad3xyz_mx, grad3xzz_mx,
        grad3yyy_mx, grad3yyz_mx, grad3yzz_mx, grad3zzz_mx,
        grad3xxx_my, grad3xxy_my, grad3xxz_my, grad3xyy_my, grad3xyz_my, grad3xzz_my,
        grad3yyy_my, grad3yyz_my, grad3yzz_my, grad3zzz_my,
        grad3xxx_mz, grad3xxy_mz, grad3xxz_mz, grad3xyy_mz, grad3xyz_mz, grad3xzz_mz,
        grad3yyy_mz, grad3yyz_mz, grad3yzz_mz, grad3zzz_mz
    );

    // 打印结果
    std::cout << "=== IR GGA_MC Test ===\n";
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "E: ";
    for (auto e : E_ir) {
        std::cout << e << " ";
    }
    std::cout << "\nV: [ ";
    for (auto &m : V_ir) {
        std::cout << "["
                  << m[0][0] << "," << m[0][1] << ";"
                  << m[1][0] << "," << m[1][1]
                  << "] ";
    }
    std::cout << "]\n";
}
}