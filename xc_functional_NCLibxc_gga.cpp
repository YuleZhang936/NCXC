// Programmed by Xiaoyu Zhang at Peking University, Beijing, China 2024/10/14
// This file is for implementing multi-collinear appraoch for GGA functionals.
// NCLibxc package only completes multi-collienar approach for LDA functionals, because the GGA functionals need gradient that is entangled with ABACUS.
//Ref: PHYSICAL REVIEW RESEARCH 5, 013036 (2023)
#include "xc_functional.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "NCLibxc/NCLibxc.h"
#include <tuple>
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <complex>
#include <iomanip>
#include <xc.h>
#include <stdexcept>
#include "NCLibxc/interface_to_libxc.h"
#include <stdexcept>
#include <fstream>
#include <omp.h>

using namespace ModuleBase;

///////////////////////////////////////////////////////////////////////////////////

// post-processing of libxc. get partial derivatives from libxc and integrate them to get the 0th,1st,2nd derivatives of the functional
//Note that e is the exchange and correlation energy per electron per volume. You need to multiply by \rho before the integration.
/*
v1=0, v2=1, f1=0,0, f2=0,1, f3=1,1
*/
void XC_Functional::postlibxc_gga(int xc_id, const std::vector<double>& rho_up, const std::vector<double>& rho_down, 
                   std::vector<double>& e, std::vector<double>& v1, std::vector<double>& v2, 
                   std::vector<double>& f1, std::vector<double>& f2, std::vector<double>& f3,const Charge* const chr,const double tpiba)
{
    LibxcInterface libxc(xc_id, true); // xc_id now passed from the caller

    int nrxx = rho_up.size(); 
    std::vector<ModuleBase::Vector3<double>> grad_rhoup(nrxx);
    std::vector<ModuleBase::Vector3<double>> grad_rhodn(nrxx);

        

    std::vector<std::complex<double>> g_rhoup(chr->rhopw->npw); // temporary storage for a scalar in reciprocal space
    std::vector<std::complex<double>> g_rhodn(chr->rhopw->npw); // temporary storage for a scalar in reciprocal space
    

    chr->rhopw->real2recip(rho_up.data(), g_rhoup.data());
    XC_Functional::grad_rho(g_rhoup.data(), grad_rhoup.data(), chr->rhopw, tpiba);
    chr->rhopw->real2recip(rho_down.data(), g_rhodn.data());
    XC_Functional::grad_rho(g_rhodn.data(), grad_rhodn.data(), chr->rhopw, tpiba);

    std::vector<double> sigma_1(nrxx), sigma_2(nrxx), sigma_3(nrxx);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif 
    for (int ir = 0; ir < nrxx; ++ir)
    {
        sigma_1[ir] = grad_rhoup[ir] * grad_rhoup[ir];
        sigma_2[ir] = grad_rhoup[ir] * grad_rhodn[ir];
        sigma_3[ir] = grad_rhodn[ir] * grad_rhodn[ir];
    }
    std::vector<double> exc = libxc.gga_exc(rho_up, rho_down, sigma_1, sigma_2, sigma_3);

    std::vector<double> vrho_1(nrxx), vrho_2(nrxx);
    std::vector<double> vsigma_1(nrxx), vsigma_2(nrxx), vsigma_3(nrxx);
    libxc.gga_vxc(rho_up, rho_down, sigma_1, sigma_2, sigma_3,
                  vrho_1, vrho_2, vsigma_1, vsigma_2, vsigma_3);

    std::vector<double> v2rho2_1(nrxx), v2rho2_2(nrxx), v2rho2_3(nrxx);
    std::vector<double> v2rhosigma_1(nrxx), v2rhosigma_2(nrxx), v2rhosigma_3(nrxx);
    std::vector<double> v2rhosigma_4(nrxx), v2rhosigma_5(nrxx), v2rhosigma_6(nrxx);
    std::vector<double> v2sigma2_1(nrxx), v2sigma2_2(nrxx), v2sigma2_3(nrxx);
    std::vector<double> v2sigma2_4(nrxx), v2sigma2_5(nrxx), v2sigma2_6(nrxx);

    libxc.gga_fxc(rho_up, rho_down, sigma_1, sigma_2, sigma_3,
                  v2rho2_1, v2rho2_2, v2rho2_3,
                  v2rhosigma_1, v2rhosigma_2, v2rhosigma_3,
                  v2rhosigma_4, v2rhosigma_5, v2rhosigma_6,
                  v2sigma2_1, v2sigma2_2, v2sigma2_3,
                  v2sigma2_4, v2sigma2_5, v2sigma2_6);
/////////////////////////////////////////////////////////////////////////////////// For v
    std::vector<ModuleBase::Vector3<double>> hup(nrxx);// h is the additional part compared to lda
    std::vector<ModuleBase::Vector3<double>> hdn(nrxx);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif 
    for (int ir = 0; ir < nrxx; ++ir)
    {
        hup[ir] = 2.0*(grad_rhoup[ir]* vsigma_1[ir]) +grad_rhodn[ir]*vsigma_2[ir];
        hdn[ir] = 2.0*(grad_rhodn[ir]* vsigma_3[ir]) +grad_rhoup[ir]*vsigma_2[ir];
    }

    std::vector<double> div_hup(nrxx);// div_h is the divergence of h
    std::vector<double> div_hdn(nrxx);
    XC_Functional::grad_dot( hup.data(), div_hup.data(), chr->rhopw, tpiba);
    XC_Functional::grad_dot( hdn.data(), div_hdn.data(), chr->rhopw, tpiba);
/////////////////////////////////////////////////////////////////////////////////// For f
    std::vector<ModuleBase::Vector3<double>> h1(nrxx);
    std::vector<ModuleBase::Vector3<double>> h2(nrxx);
    std::vector<ModuleBase::Vector3<double>> h3(nrxx);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif 
    for (int ir = 0; ir < nrxx; ++ir)
    {
        h1[ir] = 2.0*(grad_rhoup[ir]* v2rhosigma_1[ir]) +grad_rhodn[ir]*v2rhosigma_2[ir];
        h2[ir] = 2.0*(grad_rhoup[ir]* v2rhosigma_4[ir]) +grad_rhodn[ir]*v2rhosigma_5[ir]+2.0*(grad_rhodn[ir]*v2rhosigma_3[ir]) +grad_rhoup[ir]*v2rhosigma_2[ir];
        h3[ir] = 2.0*(grad_rhodn[ir]* v2rhosigma_6[ir]) +grad_rhoup[ir]*v2rhosigma_5[ir];
    }

    std::vector<double> div_h1(nrxx);
    std::vector<double> div_h2(nrxx);
    std::vector<double> div_h3(nrxx);
    XC_Functional::grad_dot( h1.data(), div_h1.data(), chr->rhopw, tpiba);
    XC_Functional::grad_dot( h2.data(), div_h2.data(), chr->rhopw, tpiba);
    XC_Functional::grad_dot( h3.data(), div_h3.data(), chr->rhopw, tpiba);

    std::vector<std::vector<ModuleBase::Vector3<double>>> H1(3, std::vector<ModuleBase::Vector3<double>>(nrxx));
    std::vector<std::vector<ModuleBase::Vector3<double>>> H2(3, std::vector<ModuleBase::Vector3<double>>(nrxx));
    std::vector<std::vector<ModuleBase::Vector3<double>>> H3(3, std::vector<ModuleBase::Vector3<double>>(nrxx));
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif 
    for (int ir = 0; ir < nrxx; ++ir){
        H1[0][ir] = (2.0*grad_rhoup[ir][0]* v2sigma2_1[ir] +grad_rhodn[ir][0]*v2sigma2_2[ir])*(2.0*grad_rhoup[ir])+(2.0*grad_rhoup[ir][0]* v2sigma2_2[ir] +grad_rhodn[ir][0]*v2sigma2_4[ir])*grad_rhodn[ir];
        H1[1][ir] = (2.0*grad_rhoup[ir][1]* v2sigma2_1[ir] +grad_rhodn[ir][1]*v2sigma2_2[ir])*(2.0*grad_rhoup[ir])+(2.0*grad_rhoup[ir][1]* v2sigma2_2[ir] +grad_rhodn[ir][1]*v2sigma2_4[ir])*grad_rhodn[ir];
        H1[2][ir] = (2.0*grad_rhoup[ir][2]* v2sigma2_1[ir] +grad_rhodn[ir][2]*v2sigma2_2[ir])*(2.0*grad_rhoup[ir])+(2.0*grad_rhoup[ir][2]* v2sigma2_2[ir] +grad_rhodn[ir][2]*v2sigma2_4[ir])*grad_rhodn[ir];

        H2[0][ir] = (v2sigma2_3[ir]*2.0*grad_rhoup[ir][0]+v2sigma2_5[ir]*grad_rhodn[ir][0])*(2.0*grad_rhodn[ir])+(v2sigma2_2[ir]*2.0*grad_rhoup[ir][0]+v2sigma2_4[ir]*grad_rhodn[ir][0])*grad_rhoup[ir];
        H2[1][ir] = (v2sigma2_3[ir]*2.0*grad_rhoup[ir][1]+v2sigma2_5[ir]*grad_rhodn[ir][1])*(2.0*grad_rhodn[ir])+(v2sigma2_2[ir]*2.0*grad_rhoup[ir][1]+v2sigma2_4[ir]*grad_rhodn[ir][1])*grad_rhoup[ir];
        H2[2][ir] = (v2sigma2_3[ir]*2*grad_rhoup[ir][2]+v2sigma2_5[ir]*grad_rhodn[ir][2])*(2.0*grad_rhodn[ir])+(v2sigma2_2[ir]*2.0*grad_rhoup[ir][2]+v2sigma2_4[ir]*grad_rhodn[ir][2])*grad_rhoup[ir];

        H3[0][ir] = (2.0*grad_rhodn[ir][0]*v2sigma2_6[ir]+grad_rhoup[ir][0]*v2sigma2_5[ir])*2.0*grad_rhoup[ir]+(2.0*grad_rhodn[ir][0]*v2sigma2_5[ir]+grad_rhoup[ir][0]*v2sigma2_4[ir])*grad_rhodn[ir];
        H3[1][ir] = (2.0*grad_rhodn[ir][1]*v2sigma2_6[ir]+grad_rhoup[ir][1]*v2sigma2_5[ir])*2.0*grad_rhoup[ir]+(2.0*grad_rhodn[ir][1]*v2sigma2_5[ir]+grad_rhoup[ir][1]*v2sigma2_4[ir])*grad_rhodn[ir];
        H3[2][ir] = (2.0*grad_rhodn[ir][2]*v2sigma2_6[ir]+grad_rhoup[ir][2]*v2sigma2_5[ir])*2.0*grad_rhoup[ir]+(2.0*grad_rhodn[ir][2]*v2sigma2_5[ir]+grad_rhoup[ir][2]*v2sigma2_4[ir])*grad_rhodn[ir];
    }

    std::vector<std::vector<double>> div_H1(3, std::vector<double>(nrxx));
    std::vector<std::vector<double>> div_H2(3, std::vector<double>(nrxx));
    std::vector<std::vector<double>> div_H3(3, std::vector<double>(nrxx));

    for(int coords=0;coords<3;coords++){
        XC_Functional::grad_dot( H1[coords].data(), div_H1[coords].data(), chr->rhopw, tpiba);
        XC_Functional::grad_dot( H2[coords].data(), div_H2[coords].data(), chr->rhopw, tpiba);
        XC_Functional::grad_dot( H3[coords].data(), div_H3[coords].data(), chr->rhopw, tpiba);
    }

    std::vector<ModuleBase::Vector3<double>> vec_div_H1(nrxx);
    std::vector<ModuleBase::Vector3<double>> vec_div_H2(nrxx);
    std::vector<ModuleBase::Vector3<double>> vec_div_H3(nrxx);
    std::vector<double> div_div_H1(nrxx);
    std::vector<double> div_div_H2(nrxx);
    std::vector<double> div_div_H3(nrxx);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif 
    for(int ir=0; ir<nrxx;++ir){
        vec_div_H1[ir] = {div_H1[0][ir],div_H1[1][ir],div_H1[2][ir]};
        vec_div_H2[ir] = {div_H2[0][ir],div_H2[1][ir],div_H2[2][ir]};
        vec_div_H3[ir] = {div_H3[0][ir],div_H3[1][ir],div_H3[2][ir]};
    }

    
    XC_Functional::grad_dot( vec_div_H1.data(), div_div_H1.data(), chr->rhopw, tpiba);
    XC_Functional::grad_dot( vec_div_H2.data(), div_div_H2.data(), chr->rhopw, tpiba);
    XC_Functional::grad_dot( vec_div_H3.data(), div_div_H3.data(), chr->rhopw, tpiba);

    std::vector<std::complex<double>> g_vsigma_1(chr->rhopw->npw); // temporary storage for a vector in reciprocal space
    std::vector<std::complex<double>> g_vsigma_2(chr->rhopw->npw); // temporary storage for a vector in reciprocal space
    std::vector<std::complex<double>> g_vsigma_3(chr->rhopw->npw); // temporary storage for a vector in reciprocal space

    std::vector<ModuleBase::Vector3<double>> grad_vsigma_1(nrxx);
    std::vector<ModuleBase::Vector3<double>> grad_vsigma_2(nrxx);
    std::vector<ModuleBase::Vector3<double>> grad_vsigma_3(nrxx);

    chr->rhopw->real2recip(vsigma_1.data(), g_vsigma_1.data());
    XC_Functional::grad_rho(g_vsigma_1.data(), grad_vsigma_1.data(), chr->rhopw, tpiba);
    chr->rhopw->real2recip(vsigma_2.data(), g_vsigma_2.data());
    XC_Functional::grad_rho(g_vsigma_2.data(), grad_vsigma_2.data(), chr->rhopw, tpiba);
    chr->rhopw->real2recip(vsigma_3.data(), g_vsigma_3.data());
    XC_Functional::grad_rho(g_vsigma_3.data(), grad_vsigma_3.data(), chr->rhopw, tpiba);

    std::vector<double> div_grad_vsigma_1(nrxx);
    std::vector<double> div_grad_vsigma_2(nrxx);
    std::vector<double> div_grad_vsigma_3(nrxx);

    XC_Functional::grad_dot( grad_vsigma_1.data(), div_grad_vsigma_1.data(), chr->rhopw, tpiba);
    XC_Functional::grad_dot( grad_vsigma_2.data(), div_grad_vsigma_2.data(), chr->rhopw, tpiba);
    XC_Functional::grad_dot( grad_vsigma_3.data(), div_grad_vsigma_3.data(), chr->rhopw, tpiba);


    e.resize(rho_up.size());
    v1.resize(rho_up.size());
    v2.resize(rho_down.size());
    f1.resize(rho_up.size());
    f2.resize(rho_down.size());
    f3.resize(rho_up.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif 
    for (size_t i = 0; i < rho_up.size(); ++i) {
        e[i] = exc[i] ; 
        v1[i] = vrho_1[i]-div_hup[i] ; 
        v2[i] = vrho_2[i]-div_hdn[i] ; 
        f1[i] = v2rho2_1[i]-2.0*div_h1[i]+div_div_H1[i]-2.0*div_grad_vsigma_1[i] ; 
        f2[i] = v2rho2_2[i]-div_h2[i]+div_div_H2[i] - div_grad_vsigma_2[i] ; 
        f3[i] = v2rho2_3[i]-2.0*div_h3[i]+div_div_H3[i]-2.0*div_grad_vsigma_3[i] ; 
    }
}

///////////////////////////////////////////////////////////////////////////////////



std::pair<std::vector<double>, std::vector<Matrix2x2>> XC_Functional::gga_mc(int xc_id, const std::vector<double>& n, 
                                                              const std::vector<double>& mx, const std::vector<double>& my, const std::vector<double>& mz,const Charge* const chr,const double tpiba)
{
    std::vector<int> nums = {6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810}; // available grid levels

    int grid_level = nums[20]; // set number of grid points. more than 1202 is recommended. default is 1454（20-th）
    std::vector<std::array<double, 4>> grids = MakeAngularGrid(grid_level); 

    size_t num_points = n.size();
    std::vector<double> E(num_points, 0.0);
    std::vector<Matrix2x2> V(num_points, {{{0.0, 0.0}, {0.0, 0.0}}});
    

    for (size_t j = 0; j < grids.size(); ++j) // 遍历 grids 的索引
    {
        double x = grids[j][0];
        double y = grids[j][1];
        double z = grids[j][2];
        double w = grids[j][3];

        std::vector<double> rho_up(num_points, 0.0);
        std::vector<double> rho_down(num_points, 0.0);
        std::vector<double> m_omega(num_points, 0.0);

        for (size_t i = 0; i < num_points; ++i)
        {
            m_omega[i] = mx[i] * x + my[i] * y + mz[i] * z;
            rho_up[i] = (n[i] + m_omega[i]) / 2.0;
            rho_down[i] = (n[i] - m_omega[i]) / 2.0;
        }

        std::vector<double> e(num_points, 0.0), v1(num_points, 0.0), v2(num_points, 0.0), f1(num_points, 0.0), f2(num_points, 0.0), f3(num_points, 0.0);
        postlibxc_gga(xc_id, rho_up, rho_down, e, v1, v2, f1, f2, f3, chr, tpiba);

        std::vector<double> Eeff(num_points, 0.0);
        std::vector<Matrix2x2> Veff(num_points, {{{0.0, 0.0}, {0.0, 0.0}}});

        for (size_t i = 0; i < num_points; ++i)
        {
            Eeff[i] = e[i] + 0.5 * (m_omega[i] / n[i]) * (v1[i] - v2[i]);

            Matrix2x2 pauli_matrix = NCLibxc::construct_pauli_matrix(x, y, z);
            Matrix2x2 term1 = NCLibxc::scalar_multiply(((v1[i] - v2[i]) + 0.25 * m_omega[i] * (f1[i] + f3[i] - 2 * f2[i])), pauli_matrix);
            Matrix2x2 term2 = NCLibxc::scalar_multiply((0.5 * (v1[i] + v2[i]) + 0.25 * m_omega[i] * (f1[i] - f3[i])), NCLibxc::identity_matrix());
            Veff[i] = NCLibxc::add(term1, term2);
            E[i] += Eeff[i] * w;
            V[i] = NCLibxc::add(V[i], NCLibxc::scalar_multiply(w, Veff[i]));
            
        }
    }

    return {E, V};
}
