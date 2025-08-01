// Programmed by Xiaoyu Zhang at Peking University, Beijing, China 2024/08/26
// Ref: PHYSICAL REVIEW RESEARCH 5, 013036 (2023)
// This file implements multi-collinear approach for LDA,GGA,MGGA and locally collinear approach for LDA. Up to the first derivative.
// how to compile as an independent program:
/*
  module load libxc/5.2.3-icc17
  module load gcc/7.3.0-wzm
 g++ -std=c++17 -o my_program NCLibxc.cpp LebedevGrid.cpp interface_to_libxc.cpp -I/public1/soft/libxc/install/include -L/public1/soft/libxc/install/lib -lxc
*/
#include "NCLibxc.h"
#include <vector>
#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>
#include "FibonacciGrid.h"
#include <tuple>

namespace NCXC{

std::vector<std::array<double, 4>> MakeAngularGrid(int grid_level);

// approximate a small collinear region
static bool is_collinear(double x, double y, double z)
{
    double r = std::sqrt(x*x + y*y + z*z);
    if (r < 1e-15) return true; //closed shell, consider it collinear
    return std::fabs(z / r) > 0.999; // mz dominates, consider it collinear
}


///////////////////////////////////////////////////////////////////////////////////

// input: xc_id, n, mx, my, mz，return E and V for each real space grid point
std::pair<std::vector<double>, std::vector<Matrix2x2>> NCLibxc::lda_mc(int xc_id, const std::vector<double>& n, 
                                                              const std::vector<double>& mx, const std::vector<double>& my, const std::vector<double>& mz)
{
    std::vector<int> nums = {6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810}; // available grid levels

    int grid_level = nums[20];// set number of grid points. more than 1202 is recommended. default is 1454.20
    std::vector<std::array<double, 4>> grids = MakeAngularGrid(grid_level); 

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
            if (n[i] == 0.0){
                continue;
            }
                
            Eeff[i] = e[i] + 0.5 * (m_omega[i] / n[i]) * (v1[i] - v2[i]);

            Matrix2x2 pauli_matrix = construct_pauli_matrix(x, y, z);
            Matrix2x2 term1 = scalar_multiply(((v1[i] - v2[i]) + 0.25 * m_omega[i] * (f1[i] + f3[i] - 2 * f2[i])), pauli_matrix);
            Matrix2x2 term2 = scalar_multiply((0.5 * (v1[i] + v2[i]) + 0.25 * m_omega[i] * (f1[i] - f3[i])), identity_matrix());
            Veff[i] = add(term1, term2);

            // integrate 
            E[i] += Eeff[i]*w;
            V[i] = add(V[i], scalar_multiply(w, Veff[i]));
        }
    }

    return {E, V};
}

// lda_lc, input: xc_id, n, mx, my, mz，return E and V for each real space grid point
std::pair<std::vector<double>, std::vector<Matrix2x2>> NCLibxc::lda_lc(int xc_id, const std::vector<double>& n, 
                                                              const std::vector<double>& mx, const std::vector<double>& my, const std::vector<double>& mz)
{
    size_t num_points = n.size();
    std::vector<double> E(num_points, 0.0);
    std::vector<Matrix2x2> V(num_points, {{{0.0, 0.0}, {0.0, 0.0}}});
    std::vector<double> m_mod(num_points, 0.0);

    // calculate the modulus of the magnetization |m|
    for (size_t i = 0; i < num_points; ++i)
    {
        m_mod[i] = std::sqrt(mx[i] * mx[i] + my[i] * my[i] + mz[i] * mz[i]);
    }

    std::vector<double> rho_up(num_points, 0.0);
    std::vector<double> rho_down(num_points, 0.0);

    for (size_t i = 0; i < num_points; ++i)
    {
        rho_up[i] = (n[i] + m_mod[i]) / 2.0;
        rho_down[i] = (n[i] - m_mod[i]) / 2.0;
    }


    std::vector<double> e(num_points, 0.0), v1(num_points, 0.0), v2(num_points, 0.0), f1(num_points, 0.0), f2(num_points, 0.0), f3(num_points, 0.0);
    postlibxc_lda(xc_id, rho_up, rho_down, e, v1, v2, f1, f2, f3);

    // calculate E and V
    for (size_t i = 0; i < num_points; ++i)
    {
        E[i] = e[i];

        
        double x = mx[i] / m_mod[i];
        double y = my[i] / m_mod[i];
        double z = mz[i] / m_mod[i];

        // construct the Pauli matrix
        Matrix2x2 pauli_matrix = construct_pauli_matrix(x, y, z);
        Matrix2x2 term1 = scalar_multiply(0.5 * (v1[i] - v2[i]), pauli_matrix);
        Matrix2x2 term2 = scalar_multiply(0.5 * (v1[i] + v2[i]), identity_matrix());
        V[i] = add(term1, term2);
    }

    return {E, V};
}

// gga_mc，input: xc_id, n, mx, my, mz，grad, grad2, grad3, return E and V for each real space grid point
std::pair<std::vector<double>, std::vector<Matrix2x2>> NCLibxc::gga_mc(int xc_id, const std::vector<double>& n, 
                                                              const std::vector<double>& mx, const std::vector<double>& my, const std::vector<double>& mz,const std::vector<double> gradx_n,const std::vector<double> grady_n,const std::vector<double> gradz_n, const std::vector<double>& gradx_mx, const std::vector<double>& grady_mx, const std::vector<double>& gradz_mx, const std::vector<double>& gradx_my, const std::vector<double>& grady_my, const std::vector<double>& gradz_my, const std::vector<double>& gradx_mz, const std::vector<double>& grady_mz, const std::vector<double>& gradz_mz, const std::vector<double> grad2xx_n, const std::vector<double> grad2yy_n, const std::vector<double> grad2zz_n, const std::vector<double> grad2xy_n, const std::vector<double> grad2yz_n, const std::vector<double> grad2xz_n, const std::vector<double> grad2xx_mx, const std::vector<double> grad2yy_mx, const std::vector<double> grad2zz_mx, const std::vector<double> grad2xy_mx, const std::vector<double> grad2yz_mx, const std::vector<double> grad2xz_mx, const std::vector<double> grad2xx_my, const std::vector<double> grad2yy_my, const std::vector<double> grad2zz_my, const std::vector<double> grad2xy_my, const std::vector<double> grad2yz_my, const std::vector<double> grad2xz_my, const std::vector<double> grad2xx_mz, const std::vector<double> grad2yy_mz, const std::vector<double> grad2zz_mz, const std::vector<double> grad2xy_mz, const std::vector<double> grad2yz_mz, const std::vector<double> grad2xz_mz, const std::vector<double> grad3xxx_n, const std::vector<double> grad3xxy_n, const std::vector<double> grad3xxz_n, const std::vector<double> grad3xyy_n, const std::vector<double> grad3xyz_n, const std::vector<double> grad3xzz_n, const std::vector<double> grad3yyy_n, const std::vector<double> grad3yyz_n, const std::vector<double> grad3yzz_n, const std::vector<double> grad3zzz_n, const std::vector<double> grad3xxx_mx, const std::vector<double> grad3xxy_mx, const std::vector<double> grad3xxz_mx, const std::vector<double> grad3xyy_mx, const std::vector<double> grad3xyz_mx, const std::vector<double> grad3xzz_mx, const std::vector<double> grad3yyy_mx, const std::vector<double> grad3yyz_mx, const std::vector<double> grad3yzz_mx, const std::vector<double> grad3zzz_mx, const std::vector<double> grad3xxx_my, const std::vector<double> grad3xxy_my, const std::vector<double> grad3xxz_my, const std::vector<double> grad3xyy_my, const std::vector<double> grad3xyz_my, const std::vector<double> grad3xzz_my, const std::vector<double> grad3yyy_my, const std::vector<double> grad3yyz_my, const std::vector<double> grad3yzz_my, const std::vector<double> grad3zzz_my, const std::vector<double> grad3xxx_mz, const std::vector<double> grad3xxy_mz, const std::vector<double> grad3xxz_mz, const std::vector<double> grad3xyy_mz, const std::vector<double> grad3xyz_mz, const std::vector<double> grad3xzz_mz, const std::vector<double> grad3yyy_mz, const std::vector<double> grad3yyz_mz, const std::vector<double> grad3yzz_mz, const std::vector<double> grad3zzz_mz)
{
    bool useFibonacciGrid = false; // Control flag for grid selection
    std::vector<std::array<double, 4>> grids;
    if(useFibonacciGrid){
        int samples = 100000; // Set the number of samples as needed
        std::vector<Point> fibonacciPoints = fibonacci_sphere(samples);

        // Convert Point structs to std::array<double, 4>
        for (const auto& p : fibonacciPoints)
        {
            grids.push_back({p.x, p.y, p.z, p.w});
        }
    }
    else{
        std::vector<int> nums = {6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810}; // available grid levels

        int grid_level = nums[20];// set number of grid points. more than 1202 is recommended. default is 1454(20-th)
        grids = MakeAngularGrid(grid_level);
    }
     
    size_t num_points = n.size();
    std::vector<double> E(num_points, 0.0);
    std::vector<Matrix2x2> Vpure(num_points, {{{0.0, 0.0}, {0.0, 0.0}}});
    std::vector<double> m_omega(num_points, 0.0);

    for (const auto &coord : grids)
    {
        double x = coord[0];
        double y = coord[1];
        double z = coord[2];
        double w = coord[3];

        std::vector<double> rho0(num_points, 0.0);
        std::vector<double> rho1(num_points, 0.0);
        std::vector<double> gradx_rho0(num_points, 0.0);
        std::vector<double> grady_rho0(num_points, 0.0);
        std::vector<double> gradz_rho0(num_points, 0.0);
        std::vector<double> gradx_rho1(num_points, 0.0);
        std::vector<double> grady_rho1(num_points, 0.0);
        std::vector<double> gradz_rho1(num_points, 0.0);
        std::vector<double> grad2xx_rho0(num_points, 0.0);
        std::vector<double> grad2yy_rho0(num_points, 0.0);
        std::vector<double> grad2zz_rho0(num_points, 0.0);
        std::vector<double> grad2xy_rho0(num_points, 0.0);
        std::vector<double> grad2yz_rho0(num_points, 0.0);
        std::vector<double> grad2xz_rho0(num_points, 0.0);
        std::vector<double> grad2xx_rho1(num_points, 0.0);
        std::vector<double> grad2yy_rho1(num_points, 0.0);
        std::vector<double> grad2zz_rho1(num_points, 0.0);
        std::vector<double> grad2xy_rho1(num_points, 0.0);
        std::vector<double> grad2yz_rho1(num_points, 0.0);
        std::vector<double> grad2xz_rho1(num_points, 0.0);

        // Initialize third-order gradients for rho0
        std::vector<double> grad3xxx_rho0(num_points, 0.0);
        std::vector<double> grad3xxy_rho0(num_points, 0.0);
        std::vector<double> grad3xxz_rho0(num_points, 0.0);
        std::vector<double> grad3xyy_rho0(num_points, 0.0);
        std::vector<double> grad3xyz_rho0(num_points, 0.0);
        std::vector<double> grad3xzz_rho0(num_points, 0.0);
        std::vector<double> grad3yyy_rho0(num_points, 0.0);
        std::vector<double> grad3yyz_rho0(num_points, 0.0);
        std::vector<double> grad3yzz_rho0(num_points, 0.0);
        std::vector<double> grad3zzz_rho0(num_points, 0.0);

        // Initialize third-order gradients for rho1
        std::vector<double> grad3xxx_rho1(num_points, 0.0);
        std::vector<double> grad3xxy_rho1(num_points, 0.0);
        std::vector<double> grad3xxz_rho1(num_points, 0.0);
        std::vector<double> grad3xyy_rho1(num_points, 0.0);
        std::vector<double> grad3xyz_rho1(num_points, 0.0);
        std::vector<double> grad3xzz_rho1(num_points, 0.0);
        std::vector<double> grad3yyy_rho1(num_points, 0.0);
        std::vector<double> grad3yyz_rho1(num_points, 0.0);
        std::vector<double> grad3yzz_rho1(num_points, 0.0);
        std::vector<double> grad3zzz_rho1(num_points, 0.0);

        // Initialize zero tau and tau gradients for GGA compatibility
        std::vector<double> tau0(num_points, 0.0);
        std::vector<double> tau1(num_points, 0.0);
        std::vector<double> gradx_tau0(num_points, 0.0);
        std::vector<double> grady_tau0(num_points, 0.0);
        std::vector<double> gradz_tau0(num_points, 0.0);
        std::vector<double> gradx_tau1(num_points, 0.0);
        std::vector<double> grady_tau1(num_points, 0.0);
        std::vector<double> gradz_tau1(num_points, 0.0);

        for (size_t i = 0; i < num_points; ++i)
        {
            m_omega[i] = mx[i] * x + my[i] * y + mz[i] * z;
            rho0[i] = (n[i] + m_omega[i]) / 2.0;
            rho1[i] = (n[i] - m_omega[i]) / 2.0;
            gradx_rho0[i] = (gradx_n[i] + gradx_mx[i] * x + gradx_my[i] * y + gradx_mz[i] * z) / 2.0;
            grady_rho0[i] = (grady_n[i] + grady_mx[i] * x + grady_my[i] * y + grady_mz[i] * z) / 2.0;
            gradz_rho0[i] = (gradz_n[i] + gradz_mx[i] * x + gradz_my[i] * y + gradz_mz[i] * z) / 2.0;
            gradx_rho1[i] = (gradx_n[i] - gradx_mx[i] * x - gradx_my[i] * y - gradx_mz[i] * z) / 2.0;
            grady_rho1[i] = (grady_n[i] - grady_mx[i] * x - grady_my[i] * y - grady_mz[i] * z) / 2.0;
            gradz_rho1[i] = (gradz_n[i] - gradz_mx[i] * x - gradz_my[i] * y - gradz_mz[i] * z) / 2.0;
            grad2xx_rho0[i] = (grad2xx_n[i] + grad2xx_mx[i] * x + grad2xx_my[i] * y + grad2xx_mz[i] * z) / 2.0;
            grad2yy_rho0[i] = (grad2yy_n[i] + grad2yy_mx[i] * x + grad2yy_my[i] * y + grad2yy_mz[i] * z) / 2.0;
            grad2zz_rho0[i] = (grad2zz_n[i] + grad2zz_mx[i] * x + grad2zz_my[i] * y + grad2zz_mz[i] * z) / 2.0;
            grad2xy_rho0[i] = (grad2xy_n[i] + grad2xy_mx[i] * x + grad2xy_my[i] * y + grad2xy_mz[i] * z) / 2.0;
            grad2yz_rho0[i] = (grad2yz_n[i] + grad2yz_mx[i] * x + grad2yz_my[i] * y + grad2yz_mz[i] * z) / 2.0;
            grad2xz_rho0[i] = (grad2xz_n[i] + grad2xz_mx[i] * x + grad2xz_my[i] * y + grad2xz_mz[i] * z) / 2.0;
            grad2xx_rho1[i] = (grad2xx_n[i] - grad2xx_mx[i] * x - grad2xx_my[i] * y - grad2xx_mz[i] * z) / 2.0;
            grad2yy_rho1[i] = (grad2yy_n[i] - grad2yy_mx[i] * x - grad2yy_my[i] * y - grad2yy_mz[i] * z) / 2.0;
            grad2zz_rho1[i] = (grad2zz_n[i] - grad2zz_mx[i] * x - grad2zz_my[i] * y - grad2zz_mz[i] * z) / 2.0;
            grad2xy_rho1[i] = (grad2xy_n[i] - grad2xy_mx[i] * x - grad2xy_my[i] * y - grad2xy_mz[i] * z) / 2.0;
            grad2yz_rho1[i] = (grad2yz_n[i] - grad2yz_mx[i] * x - grad2yz_my[i] * y - grad2yz_mz[i] * z) / 2.0;
            grad2xz_rho1[i] = (grad2xz_n[i] - grad2xz_mx[i] * x - grad2xz_my[i] * y - grad2xz_mz[i] * z) / 2.0;

            grad3xxx_rho0[i] = (grad3xxx_n[i] + grad3xxx_mx[i] * x + grad3xxx_my[i] * y + grad3xxx_mz[i] * z) / 2.0;
            grad3xxy_rho0[i] = (grad3xxy_n[i] + grad3xxy_mx[i] * x + grad3xxy_my[i] * y + grad3xxy_mz[i] * z) / 2.0;
            grad3xxz_rho0[i] = (grad3xxz_n[i] + grad3xxz_mx[i] * x + grad3xxz_my[i] * y + grad3xxz_mz[i] * z) / 2.0;
            grad3xyy_rho0[i] = (grad3xyy_n[i] + grad3xyy_mx[i] * x + grad3xyy_my[i] * y + grad3xyy_mz[i] * z) / 2.0;
            grad3xyz_rho0[i] = (grad3xyz_n[i] + grad3xyz_mx[i] * x + grad3xyz_my[i] * y + grad3xyz_mz[i] * z) / 2.0;
            grad3xzz_rho0[i] = (grad3xzz_n[i] + grad3xzz_mx[i] * x + grad3xzz_my[i] * y + grad3xzz_mz[i] * z) / 2.0;
            grad3yyy_rho0[i] = (grad3yyy_n[i] + grad3yyy_mx[i] * x + grad3yyy_my[i] * y + grad3yyy_mz[i] * z) / 2.0;
            grad3yyz_rho0[i] = (grad3yyz_n[i] + grad3yyz_mx[i] * x + grad3yyz_my[i] * y + grad3yyz_mz[i] * z) / 2.0;
            grad3yzz_rho0[i] = (grad3yzz_n[i] + grad3yzz_mx[i] * x + grad3yzz_my[i] * y + grad3yzz_mz[i] * z) / 2.0;
            grad3zzz_rho0[i] = (grad3zzz_n[i] + grad3zzz_mx[i] * x + grad3zzz_my[i] * y + grad3zzz_mz[i] * z) / 2.0;

            grad3xxx_rho1[i] = (grad3xxx_n[i] - grad3xxx_mx[i] * x - grad3xxx_my[i] * y - grad3xxx_mz[i] * z) / 2.0;
            grad3xxy_rho1[i] = (grad3xxy_n[i] - grad3xxy_mx[i] * x - grad3xxy_my[i] * y - grad3xxy_mz[i] * z) / 2.0;
            grad3xxz_rho1[i] = (grad3xxz_n[i] - grad3xxz_mx[i] * x - grad3xxz_my[i] * y - grad3xxz_mz[i] * z) / 2.0;
            grad3xyy_rho1[i] = (grad3xyy_n[i] - grad3xyy_mx[i] * x - grad3xyy_my[i] * y - grad3xyy_mz[i] * z) / 2.0;
            grad3xyz_rho1[i] = (grad3xyz_n[i] - grad3xyz_mx[i] * x - grad3xyz_my[i] * y - grad3xyz_mz[i] * z) / 2.0;
            grad3xzz_rho1[i] = (grad3xzz_n[i] - grad3xzz_mx[i] * x - grad3xzz_my[i] * y - grad3xzz_mz[i] * z) / 2.0;
            grad3yyy_rho1[i] = (grad3yyy_n[i] - grad3yyy_mx[i] * x - grad3yyy_my[i] * y - grad3yyy_mz[i] * z) / 2.0;
            grad3yyz_rho1[i] = (grad3yyz_n[i] - grad3yyz_mx[i] * x - grad3yyz_my[i] * y - grad3yyz_mz[i] * z) / 2.0;
            grad3yzz_rho1[i] = (grad3yzz_n[i] - grad3yzz_mx[i] * x - grad3yzz_my[i] * y - grad3yzz_mz[i] * z) / 2.0;
            grad3zzz_rho1[i] = (grad3zzz_n[i] - grad3zzz_mx[i] * x - grad3zzz_my[i] * y - grad3zzz_mz[i] * z) / 2.0;
        }

        std::vector<double> e(num_points, 0.0), vn(num_points, 0.0), vs(num_points, 0.0), vgn(num_points, 0.0), vgs(num_points, 0.0);

        NCLibxc::postlibxc_gga_eff(xc_id, 
            rho0, 
            rho1, 
            gradx_rho0, 
            grady_rho0, 
            gradz_rho0, 
            gradx_rho1, 
            grady_rho1, 
            gradz_rho1, 
            grad2xx_rho0, 
            grad2yy_rho0, 
            grad2zz_rho0, 
            grad2xy_rho0, 
            grad2yz_rho0, 
            grad2xz_rho0, 
            grad2xx_rho1, 
            grad2yy_rho1, 
            grad2zz_rho1, 
            grad2xy_rho1, 
            grad2yz_rho1, 
            grad2xz_rho1, 
            grad3xxx_rho0, 
            grad3xxy_rho0, 
            grad3xxz_rho0, 
            grad3xyy_rho0, 
            grad3xyz_rho0, 
            grad3xzz_rho0, 
            grad3yyy_rho0, 
            grad3yyz_rho0, 
            grad3yzz_rho0, 
            grad3zzz_rho0, 
            grad3xxx_rho1, 
            grad3xxy_rho1, 
            grad3xxz_rho1, 
            grad3xyy_rho1, 
            grad3xyz_rho1, 
            grad3xzz_rho1, 
            grad3yyy_rho1, 
            grad3yyz_rho1, 
            grad3yzz_rho1, 
            grad3zzz_rho1, 
            tau0,
            tau1,
            gradx_tau0, 
            grady_tau0, 
            gradz_tau0, 
            gradx_tau1, 
            grady_tau1, 
            gradz_tau1, 
            e, 
            vn, 
            vs,
            vgn,
            vgs);

        std::vector<double> Eeff(num_points, 0.0);
        std::vector<Matrix2x2> Vpure_eff(num_points, {{{0.0, 0.0}, {0.0, 0.0}}});

        for (size_t i = 0; i < num_points; ++i)
        {
            if (n[i] == 0.0) {
                continue;
            }
            
            Eeff[i] = e[i];

            Matrix2x2 pauli_matrix = construct_pauli_matrix(x, y, z);
            Matrix2x2 term1 = scalar_multiply(vs[i], pauli_matrix);
            Matrix2x2 term2 = scalar_multiply(vn[i], identity_matrix());
            Vpure_eff[i] = add(term1, term2);

            // integrate
            E[i] += Eeff[i]*w;
            Vpure[i] = add(Vpure[i], scalar_multiply(w, Vpure_eff[i]));
        }
    }

    return {E, Vpure};
}

// mgga_mc，input: xc_id, n, mx, my, mz，grad rho, grad2rho, grad3rho, tau, ux, uy,uz,grad tau(u), grad2 tau(u) return E, Vtau, Vpure for each real space grid point
std::tuple<std::vector<double>, std::vector<Matrix2x2>, std::vector<Matrix2x2>> NCLibxc::mgga_mc(int xc_id, const std::vector<double>& n, 
    const std::vector<double>& mx, const std::vector<double>& my, const std::vector<double>& mz,const std::vector<double> gradx_n,const std::vector<double> grady_n,const std::vector<double> gradz_n, const std::vector<double>& gradx_mx, const std::vector<double>& grady_mx, const std::vector<double>& gradz_mx, const std::vector<double>& gradx_my, const std::vector<double>& grady_my, const std::vector<double>& gradz_my, const std::vector<double>& gradx_mz, const std::vector<double>& grady_mz, const std::vector<double>& gradz_mz, const std::vector<double> grad2xx_n, const std::vector<double> grad2yy_n, const std::vector<double> grad2zz_n, const std::vector<double> grad2xy_n, const std::vector<double> grad2yz_n, const std::vector<double> grad2xz_n, const std::vector<double> grad2xx_mx, const std::vector<double> grad2yy_mx, const std::vector<double> grad2zz_mx, const std::vector<double> grad2xy_mx, const std::vector<double> grad2yz_mx, const std::vector<double> grad2xz_mx, const std::vector<double> grad2xx_my, const std::vector<double> grad2yy_my, const std::vector<double> grad2zz_my, const std::vector<double> grad2xy_my, const std::vector<double> grad2yz_my, const std::vector<double> grad2xz_my, const std::vector<double> grad2xx_mz, const std::vector<double> grad2yy_mz, const std::vector<double> grad2zz_mz, const std::vector<double> grad2xy_mz, const std::vector<double> grad2yz_mz, const std::vector<double> grad2xz_mz, const std::vector<double> grad3xxx_n, const std::vector<double> grad3xxy_n, const std::vector<double> grad3xxz_n, const std::vector<double> grad3xyy_n, const std::vector<double> grad3xyz_n, const std::vector<double> grad3xzz_n, const std::vector<double> grad3yyy_n, const std::vector<double> grad3yyz_n, const std::vector<double> grad3yzz_n, const std::vector<double> grad3zzz_n, const std::vector<double> grad3xxx_mx, const std::vector<double> grad3xxy_mx, const std::vector<double> grad3xxz_mx, const std::vector<double> grad3xyy_mx, const std::vector<double> grad3xyz_mx, const std::vector<double> grad3xzz_mx, const std::vector<double> grad3yyy_mx, const std::vector<double> grad3yyz_mx, const std::vector<double> grad3yzz_mx, const std::vector<double> grad3zzz_mx, const std::vector<double> grad3xxx_my, const std::vector<double> grad3xxy_my, const std::vector<double> grad3xxz_my, const std::vector<double> grad3xyy_my, const std::vector<double> grad3xyz_my, const std::vector<double> grad3xzz_my, const std::vector<double> grad3yyy_my, const std::vector<double> grad3yyz_my, const std::vector<double> grad3yzz_my, const std::vector<double> grad3zzz_my, const std::vector<double> grad3xxx_mz, const std::vector<double> grad3xxy_mz, const std::vector<double> grad3xxz_mz, const std::vector<double> grad3xyy_mz, const std::vector<double> grad3xyz_mz, const std::vector<double> grad3xzz_mz, const std::vector<double> grad3yyy_mz, const std::vector<double> grad3yyz_mz, const std::vector<double> grad3yzz_mz, const std::vector<double> grad3zzz_mz, const std::vector<double>& tau,const std::vector<double>& ux, const std::vector<double>& uy, const std::vector<double>& uz,
    const std::vector<double> gradx_tau, const std::vector<double> grady_tau, const std::vector<double> gradz_tau,
    const std::vector<double>& gradx_ux, const std::vector<double>& grady_ux, const std::vector<double>& gradz_ux,
    const std::vector<double>& gradx_uy, const std::vector<double>& grady_uy, const std::vector<double>& gradz_uy,
    const std::vector<double>& gradx_uz, const std::vector<double>& grady_uz, const std::vector<double>& gradz_uz)
{
    bool useFibonacciGrid = false; // Control flag for grid selection
    std::vector<std::array<double, 4>> grids;
    if(useFibonacciGrid){
        int samples = 100000; // Set the number of samples as needed
        std::vector<Point> fibonacciPoints = fibonacci_sphere(samples);

        // Convert Point structs to std::array<double, 4>
        for (const auto& p : fibonacciPoints)
        {
            grids.push_back({p.x, p.y, p.z, p.w});
        }
    }
    else{
        std::vector<int> nums = {6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810}; // available grid levels

        int grid_level = nums[20];// set number of grid points. more than 1202 is recommended. default is 1454(20-th)
        grids = MakeAngularGrid(grid_level);
    }
     
    size_t num_points = n.size();
    std::vector<double> E(num_points, 0.0);
    std::vector<Matrix2x2> Vpure(num_points, {{{0.0, 0.0}, {0.0, 0.0}}});
    std::vector<Matrix2x2> Vtau(num_points, {{{0.0, 0.0}, {0.0, 0.0}}});
    std::vector<double> m_omega(num_points, 0.0);
    std::vector<double> u_omega(num_points, 0.0);

    for (const auto &coord : grids)
    {
        double x = coord[0];
        double y = coord[1];
        double z = coord[2];
        double w = coord[3];

        std::vector<double> rho0(num_points, 0.0);
        std::vector<double> rho1(num_points, 0.0);
        std::vector<double> gradx_rho0(num_points, 0.0);
        std::vector<double> grady_rho0(num_points, 0.0);
        std::vector<double> gradz_rho0(num_points, 0.0);
        std::vector<double> gradx_rho1(num_points, 0.0);
        std::vector<double> grady_rho1(num_points, 0.0);
        std::vector<double> gradz_rho1(num_points, 0.0);
        std::vector<double> grad2xx_rho0(num_points, 0.0);
        std::vector<double> grad2yy_rho0(num_points, 0.0);
        std::vector<double> grad2zz_rho0(num_points, 0.0);
        std::vector<double> grad2xy_rho0(num_points, 0.0);
        std::vector<double> grad2yz_rho0(num_points, 0.0);
        std::vector<double> grad2xz_rho0(num_points, 0.0);
        std::vector<double> grad2xx_rho1(num_points, 0.0);
        std::vector<double> grad2yy_rho1(num_points, 0.0);
        std::vector<double> grad2zz_rho1(num_points, 0.0);
        std::vector<double> grad2xy_rho1(num_points, 0.0);
        std::vector<double> grad2yz_rho1(num_points, 0.0);
        std::vector<double> grad2xz_rho1(num_points, 0.0);

        // Initialize third-order gradients for rho0
        std::vector<double> grad3xxx_rho0(num_points, 0.0);
        std::vector<double> grad3xxy_rho0(num_points, 0.0);
        std::vector<double> grad3xxz_rho0(num_points, 0.0);
        std::vector<double> grad3xyy_rho0(num_points, 0.0);
        std::vector<double> grad3xyz_rho0(num_points, 0.0);
        std::vector<double> grad3xzz_rho0(num_points, 0.0);
        std::vector<double> grad3yyy_rho0(num_points, 0.0);
        std::vector<double> grad3yyz_rho0(num_points, 0.0);
        std::vector<double> grad3yzz_rho0(num_points, 0.0);
        std::vector<double> grad3zzz_rho0(num_points, 0.0);

        // Initialize third-order gradients for rho1
        std::vector<double> grad3xxx_rho1(num_points, 0.0);
        std::vector<double> grad3xxy_rho1(num_points, 0.0);
        std::vector<double> grad3xxz_rho1(num_points, 0.0);
        std::vector<double> grad3xyy_rho1(num_points, 0.0);
        std::vector<double> grad3xyz_rho1(num_points, 0.0);
        std::vector<double> grad3xzz_rho1(num_points, 0.0);
        std::vector<double> grad3yyy_rho1(num_points, 0.0);
        std::vector<double> grad3yyz_rho1(num_points, 0.0);
        std::vector<double> grad3yzz_rho1(num_points, 0.0);
        std::vector<double> grad3zzz_rho1(num_points, 0.0);

        std::vector<double> tau0(num_points, 0.0);
        std::vector<double> tau1(num_points, 0.0);
        std::vector<double> gradx_tau0(num_points, 0.0);
        std::vector<double> grady_tau0(num_points, 0.0);
        std::vector<double> gradz_tau0(num_points, 0.0);
        std::vector<double> gradx_tau1(num_points, 0.0);
        std::vector<double> grady_tau1(num_points, 0.0);
        std::vector<double> gradz_tau1(num_points, 0.0);


        for (size_t i = 0; i < num_points; ++i)
        {
            m_omega[i] = mx[i] * x + my[i] * y + mz[i] * z;
            u_omega[i] = ux[i] * x + uy[i] * y + uz[i] * z;
            rho0[i] = (n[i] + m_omega[i]) / 2.0;
            rho1[i] = (n[i] - m_omega[i]) / 2.0;
            gradx_rho0[i] = (gradx_n[i] + gradx_mx[i] * x + gradx_my[i] * y + gradx_mz[i] * z) / 2.0;
            grady_rho0[i] = (grady_n[i] + grady_mx[i] * x + grady_my[i] * y + grady_mz[i] * z) / 2.0;
            gradz_rho0[i] = (gradz_n[i] + gradz_mx[i] * x + gradz_my[i] * y + gradz_mz[i] * z) / 2.0;
            gradx_rho1[i] = (gradx_n[i] - gradx_mx[i] * x - gradx_my[i] * y - gradx_mz[i] * z) / 2.0;
            grady_rho1[i] = (grady_n[i] - grady_mx[i] * x - grady_my[i] * y - grady_mz[i] * z) / 2.0;
            gradz_rho1[i] = (gradz_n[i] - gradz_mx[i] * x - gradz_my[i] * y - gradz_mz[i] * z) / 2.0;
            grad2xx_rho0[i] = (grad2xx_n[i] + grad2xx_mx[i] * x + grad2xx_my[i] * y + grad2xx_mz[i] * z) / 2.0;
            grad2yy_rho0[i] = (grad2yy_n[i] + grad2yy_mx[i] * x + grad2yy_my[i] * y + grad2yy_mz[i] * z) / 2.0;
            grad2zz_rho0[i] = (grad2zz_n[i] + grad2zz_mx[i] * x + grad2zz_my[i] * y + grad2zz_mz[i] * z) / 2.0;
            grad2xy_rho0[i] = (grad2xy_n[i] + grad2xy_mx[i] * x + grad2xy_my[i] * y + grad2xy_mz[i] * z) / 2.0;
            grad2yz_rho0[i] = (grad2yz_n[i] + grad2yz_mx[i] * x + grad2yz_my[i] * y + grad2yz_mz[i] * z) / 2.0;
            grad2xz_rho0[i] = (grad2xz_n[i] + grad2xz_mx[i] * x + grad2xz_my[i] * y + grad2xz_mz[i] * z) / 2.0;
            grad2xx_rho1[i] = (grad2xx_n[i] - grad2xx_mx[i] * x - grad2xx_my[i] * y - grad2xx_mz[i] * z) / 2.0;
            grad2yy_rho1[i] = (grad2yy_n[i] - grad2yy_mx[i] * x - grad2yy_my[i] * y - grad2yy_mz[i] * z) / 2.0;
            grad2zz_rho1[i] = (grad2zz_n[i] - grad2zz_mx[i] * x - grad2zz_my[i] * y - grad2zz_mz[i] * z) / 2.0;
            grad2xy_rho1[i] = (grad2xy_n[i] - grad2xy_mx[i] * x - grad2xy_my[i] * y - grad2xy_mz[i] * z) / 2.0;
            grad2yz_rho1[i] = (grad2yz_n[i] - grad2yz_mx[i] * x - grad2yz_my[i] * y - grad2yz_mz[i] * z) / 2.0;
            grad2xz_rho1[i] = (grad2xz_n[i] - grad2xz_mx[i] * x - grad2xz_my[i] * y - grad2xz_mz[i] * z) / 2.0;

            grad3xxx_rho0[i] = (grad3xxx_n[i] + grad3xxx_mx[i] * x + grad3xxx_my[i] * y + grad3xxx_mz[i] * z) / 2.0;
            grad3xxy_rho0[i] = (grad3xxy_n[i] + grad3xxy_mx[i] * x + grad3xxy_my[i] * y + grad3xxy_mz[i] * z) / 2.0;
            grad3xxz_rho0[i] = (grad3xxz_n[i] + grad3xxz_mx[i] * x + grad3xxz_my[i] * y + grad3xxz_mz[i] * z) / 2.0;
            grad3xyy_rho0[i] = (grad3xyy_n[i] + grad3xyy_mx[i] * x + grad3xyy_my[i] * y + grad3xyy_mz[i] * z) / 2.0;
            grad3xyz_rho0[i] = (grad3xyz_n[i] + grad3xyz_mx[i] * x + grad3xyz_my[i] * y + grad3xyz_mz[i] * z) / 2.0;
            grad3xzz_rho0[i] = (grad3xzz_n[i] + grad3xzz_mx[i] * x + grad3xzz_my[i] * y + grad3xzz_mz[i] * z) / 2.0;
            grad3yyy_rho0[i] = (grad3yyy_n[i] + grad3yyy_mx[i] * x + grad3yyy_my[i] * y + grad3yyy_mz[i] * z) / 2.0;
            grad3yyz_rho0[i] = (grad3yyz_n[i] + grad3yyz_mx[i] * x + grad3yyz_my[i] * y + grad3yyz_mz[i] * z) / 2.0;
            grad3yzz_rho0[i] = (grad3yzz_n[i] + grad3yzz_mx[i] * x + grad3yzz_my[i] * y + grad3yzz_mz[i] * z) / 2.0;
            grad3zzz_rho0[i] = (grad3zzz_n[i] + grad3zzz_mx[i] * x + grad3zzz_my[i] * y + grad3zzz_mz[i] * z) / 2.0;

            grad3xxx_rho1[i] = (grad3xxx_n[i] - grad3xxx_mx[i] * x - grad3xxx_my[i] * y - grad3xxx_mz[i] * z) / 2.0;
            grad3xxy_rho1[i] = (grad3xxy_n[i] - grad3xxy_mx[i] * x - grad3xxy_my[i] * y - grad3xxy_mz[i] * z) / 2.0;
            grad3xxz_rho1[i] = (grad3xxz_n[i] - grad3xxz_mx[i] * x - grad3xxz_my[i] * y - grad3xxz_mz[i] * z) / 2.0;
            grad3xyy_rho1[i] = (grad3xyy_n[i] - grad3xyy_mx[i] * x - grad3xyy_my[i] * y - grad3xyy_mz[i] * z) / 2.0;
            grad3xyz_rho1[i] = (grad3xyz_n[i] - grad3xyz_mx[i] * x - grad3xyz_my[i] * y - grad3xyz_mz[i] * z) / 2.0;
            grad3xzz_rho1[i] = (grad3xzz_n[i] - grad3xzz_mx[i] * x - grad3xzz_my[i] * y - grad3xzz_mz[i] * z) / 2.0;
            grad3yyy_rho1[i] = (grad3yyy_n[i] - grad3yyy_mx[i] * x - grad3yyy_my[i] * y - grad3yyy_mz[i] * z) / 2.0;
            grad3yyz_rho1[i] = (grad3yyz_n[i] - grad3yyz_mx[i] * x - grad3yyz_my[i] * y - grad3yyz_mz[i] * z) / 2.0;
            grad3yzz_rho1[i] = (grad3yzz_n[i] - grad3yzz_mx[i] * x - grad3yzz_my[i] * y - grad3yzz_mz[i] * z) / 2.0;
            grad3zzz_rho1[i] = (grad3zzz_n[i] - grad3zzz_mx[i] * x - grad3zzz_my[i] * y - grad3zzz_mz[i] * z) / 2.0;

            tau0[i] = (tau[i] + u_omega[i]) / 2.0;
            tau1[i] = (tau[i] - u_omega[i]) / 2.0;
            gradx_tau0[i] = (gradx_tau[i] + gradx_ux[i] * x + gradx_uy[i] * y + gradx_uz[i] * z) / 2.0;
            grady_tau0[i] = (grady_tau[i] + grady_ux[i] * x + grady_uy[i] * y + grady_uz[i] * z) / 2.0;
            gradz_tau0[i] = (gradz_tau[i] + gradz_ux[i] * x + gradz_uy[i] * y + gradz_uz[i] * z) / 2.0;
            gradx_tau1[i] = (gradx_tau[i] - gradx_ux[i] * x - gradx_uy[i] * y - gradx_uz[i] * z) / 2.0;
            grady_tau1[i] = (grady_tau[i] - grady_ux[i] * x - grady_uy[i] * y - grady_uz[i] * z) / 2.0;
            gradz_tau1[i] = (gradz_tau[i] - gradz_ux[i] * x - gradz_uy[i] * y - gradz_uz[i] * z) / 2.0;
        }


        std::vector<double> e(num_points, 0.0), vn(num_points, 0.0), vgn(num_points, 0.0), vs(num_points, 0.0), vgs(num_points, 0.0);

        NCLibxc::postlibxc_mgga(xc_id, 
            rho0, 
            rho1, 
            gradx_rho0, 
            grady_rho0, 
            gradz_rho0, 
            gradx_rho1, 
            grady_rho1, 
            gradz_rho1, 
            grad2xx_rho0, 
            grad2yy_rho0, 
            grad2zz_rho0, 
            grad2xy_rho0, 
            grad2yz_rho0, 
            grad2xz_rho0, 
            grad2xx_rho1, 
            grad2yy_rho1, 
            grad2zz_rho1, 
            grad2xy_rho1, 
            grad2yz_rho1, 
            grad2xz_rho1, 
            grad3xxx_rho0, 
            grad3xxy_rho0, 
            grad3xxz_rho0, 
            grad3xyy_rho0, 
            grad3xyz_rho0, 
            grad3xzz_rho0, 
            grad3yyy_rho0, 
            grad3yyz_rho0, 
            grad3yzz_rho0, 
            grad3zzz_rho0, 
            grad3xxx_rho1, 
            grad3xxy_rho1, 
            grad3xxz_rho1, 
            grad3xyy_rho1, 
            grad3xyz_rho1, 
            grad3xzz_rho1, 
            grad3yyy_rho1, 
            grad3yyz_rho1, 
            grad3yzz_rho1, 
            grad3zzz_rho1, 
            tau0,
            tau1,
            gradx_tau0, 
            grady_tau0, 
            gradz_tau0, 
            gradx_tau1, 
            grady_tau1, 
            gradz_tau1, 
            e, 
            vn, 
            vs,
            vgn,
            vgs);

            std::vector<double> Eeff(num_points, 0.0);
            std::vector<Matrix2x2> Vpure_eff(num_points, {{{0.0, 0.0}, {0.0, 0.0}}});
            std::vector<Matrix2x2> Vtau_eff(num_points, {{{0.0, 0.0}, {0.0, 0.0}}});
    
            for (size_t i = 0; i < num_points; ++i)
            {
                if (n[i] == 0.0){
                    continue;
                }
                
                Eeff[i] = e[i];
    
                Matrix2x2 pauli_matrix = construct_pauli_matrix(x, y, z);
                Matrix2x2 term1 = scalar_multiply(vs[i], pauli_matrix);
                Matrix2x2 term2 = scalar_multiply(vn[i], identity_matrix());
                Matrix2x2 term3 = scalar_multiply(vgs[i], pauli_matrix);
                Matrix2x2 term4 = scalar_multiply(vgn[i], identity_matrix());
                Vpure_eff[i] = add(term1, term2);
                Vtau_eff[i] = add(term3, term4);
    
                // integrate
                E[i] += Eeff[i]*w;
                Vpure[i] = add(Vpure[i], scalar_multiply(w, Vpure_eff[i]));
                Vtau[i] = add(Vtau[i], scalar_multiply(w, Vtau_eff[i]));
            }

    }
    return std::make_tuple(E, Vpure, Vtau);
}


}// end of namespace