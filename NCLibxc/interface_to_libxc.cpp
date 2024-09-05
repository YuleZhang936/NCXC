
// This program is used to interface with libxc which is a library for exchange-correlation functionals in density functional theory.

// How to compile as an independent program:
//  module load libxc/5.2.3-icc17
//  g++ -std=c++11 -o my_program interface_to_libxc.cpp -I/public1/soft/libxc/install/include -L/public1/soft/libxc/install/lib -lxc
//  ./my_program

// func_id can be found in the libxc website https://libxc.gitlab.io/functionals/
#include "interface_to_libxc.h"
#include <iostream>

LibxcInterface::LibxcInterface(int xc_id, bool spin_polarized)
{
    int spin_option = spin_polarized ? XC_POLARIZED : XC_UNPOLARIZED;
    if (xc_func_init(&func, xc_id, spin_option) != 0)
    {
        std::cerr << "Failed to initialize functional." << std::endl;
        throw std::runtime_error("Libxc initialization error");
    }
}

LibxcInterface::~LibxcInterface()
{
    xc_func_end(&func);
}

    /////////////////////////////////////////////////////////////
    // LDA START, up to the fourth derivative
    // LDA Energy Density for spin-polarized systems
    std::vector<double> LibxcInterface::lda_exc(const std::vector<double> &rho_up, const std::vector<double> &rho_down)
    {
        int np = rho_up.size();
        std::vector<double> rho(2 * np);
        for (int i = 0; i < np; ++i)
        {
            rho[2 * i] = rho_up[i];
            rho[2 * i + 1] = rho_down[i];
        }
        std::vector<double> exc(np);
        xc_lda_exc(&func, np, rho.data(), exc.data());
        return exc;
    }

    // LDA Potential for spin-polarized systems
    void LibxcInterface::lda_vxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 std::vector<double> &vrho_1, std::vector<double> &vrho_2)
    {
        int np = rho_up.size();
        std::vector<double> rho(2 * np);
        std::vector<double> vrho(2 * np);
        for (int i = 0; i < np; ++i)
        {
            rho[2 * i] = rho_up[i];
            rho[2 * i + 1] = rho_down[i];
        }
        xc_lda_vxc(&func, np, rho.data(), vrho.data());
        for (int i = 0; i < np; ++i)
        {
            vrho_1[i] = vrho[2 * i];
            vrho_2[i] = vrho[2 * i + 1];
        }
    }

    // LDA the second derivative for spin-polarized systems
    void LibxcInterface::lda_fxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 std::vector<double> &v2rho2_1, std::vector<double> &v2rho2_2, std::vector<double> &v2rho2_3)
    {
        int np = rho_up.size();
        std::vector<double> rho(2 * np);
        std::vector<double> v2rho2(3 * np);
        for (int i = 0; i < np; ++i)
        {
            rho[2 * i] = rho_up[i];
            rho[2 * i + 1] = rho_down[i];
        }
        xc_lda_fxc(&func, np, rho.data(), v2rho2.data());
        for (int i = 0; i < np; ++i)
        {
            v2rho2_1[i] = v2rho2[3 * i];
            v2rho2_2[i] = v2rho2[3 * i + 1];
            v2rho2_3[i] = v2rho2[3 * i + 2];
        }
    }

    // LDA the third derivative for spin-polarized systems
    void LibxcInterface::lda_kxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 std::vector<double> &v3rho3_1, std::vector<double> &v3rho3_2, std::vector<double> &v3rho3_3, std::vector<double> &v3rho3_4)
    {
        int np = rho_up.size();
        std::vector<double> rho(2 * np);
        std::vector<double> v3rho3(4 * np);
        for (int i = 0; i < np; ++i)
        {
            rho[2 * i] = rho_up[i];
            rho[2 * i + 1] = rho_down[i];
        }
        xc_lda_kxc(&func, np, rho.data(), v3rho3.data());
        for (int i = 0; i < np; ++i)
        {
            v3rho3_1[i] = v3rho3[4 * i];
            v3rho3_2[i] = v3rho3[4 * i + 1];
            v3rho3_3[i] = v3rho3[4 * i + 2];
            v3rho3_4[i] = v3rho3[4 * i + 3];
        }
    }

    // LDA the fourth derivative for spin-polarized systems
    void LibxcInterface::lda_lxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 std::vector<double> &v4rho4_1, std::vector<double> &v4rho4_2, std::vector<double> &v4rho4_3, std::vector<double> &v4rho4_4, std::vector<double> &v4rho4_5)
    {
        int np = rho_up.size();
        std::vector<double> rho(2 * np);
        std::vector<double> v4rho4(5 * np);
        for (int i = 0; i < np; ++i)
        {
            rho[2 * i] = rho_up[i];
            rho[2 * i + 1] = rho_down[i];
        }
        xc_lda_lxc(&func, np, rho.data(), v4rho4.data());
        for (int i = 0; i < np; ++i)
        {
            v4rho4_1[i] = v4rho4[5 * i];
            v4rho4_2[i] = v4rho4[5 * i + 1];
            v4rho4_3[i] = v4rho4[5 * i + 2];
            v4rho4_4[i] = v4rho4[5 * i + 3];
            v4rho4_5[i] = v4rho4[5 * i + 4];
        }
    }
    //##################  LDA END ###############################################


    /////////////////////////////////////////////////////////////
    //GGA START, up to the second derivative
    // GGA Energy Density for spin-polarized systems
    std::vector<double> LibxcInterface::gga_exc(const std::vector<double> &rho_up, const std::vector<double> &rho_down, const std::vector<double> &sigma_1, const std::vector<double> &sigma_2, const std::vector<double> &sigma_3)
    {
        int np = rho_up.size();
        std::vector<double> rho(2 * np);
        std::vector<double> sigma(3 * np);
        for (int i = 0; i < np; ++i)
        {
            rho[2 * i] = rho_up[i];
            rho[2 * i + 1] = rho_down[i];
            sigma[3 * i] = sigma_1[i];
            sigma[3 * i + 1] = sigma_2[i];
            sigma[3 * i + 2] = sigma_3[i];
        }
        std::vector<double> exc(np);
        xc_gga_exc(&func, np, rho.data(), sigma.data(), exc.data());
        return exc;
    }

    // GGA Potential for spin-polarized systems
    void LibxcInterface::gga_vxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 const std::vector<double> &sigma_1, const std::vector<double> &sigma_2, const std::vector<double> &sigma_3, std::vector<double> &vrho_1, std::vector<double> &vrho_2,
                 std::vector<double> &vsigma_1, std::vector<double> &vsigma_2, std::vector<double> &vsigma_3)
    {
        int np = rho_up.size();
        std::vector<double> rho(2 * np);
        std::vector<double> sigma(3 * np);
        std::vector<double> vrho(2 * np);
        std::vector<double> vsigma(3 * np);
        for (int i = 0; i < np; ++i)
        {
            rho[2 * i] = rho_up[i];
            rho[2 * i + 1] = rho_down[i];
            sigma[3 * i] = sigma_1[i];
            sigma[3 * i + 1] = sigma_2[i];
            sigma[3 * i + 2] = sigma_3[i];
        }
        xc_gga_vxc(&func, np, rho.data(), sigma.data(), vrho.data(), vsigma.data());
        for (int i = 0; i < np; ++i)
        {
            vrho_1[i] = vrho[2 * i];
            vrho_2[i] = vrho[2 * i + 1];
            vsigma_1[i] = vsigma[3 * i];
            vsigma_2[i] = vsigma[3 * i + 1];
            vsigma_3[i] = vsigma[3 * i + 2];
        }
    }

    // GGA the second derivative for spin-polarized systems
    void LibxcInterface::gga_fxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 const std::vector<double> &sigma_1, const std::vector<double> &sigma_2, const std::vector<double> &sigma_3,
                 std::vector<double> &v2rho2_1, std::vector<double> &v2rho2_2, std::vector<double> &v2rho2_3,
                 std::vector<double> &v2rhosigma_1, std::vector<double> &v2rhosigma_2, std::vector<double> &v2rhosigma_3, std::vector<double> &v2rhosigma_4, std::vector<double> &v2rhosigma_5, std::vector<double> &v2rhosigma_6, std::vector<double> &v2sigma2_1, std::vector<double> &v2sigma2_2, std::vector<double> &v2sigma2_3, std::vector<double> &v2sigma2_4, std::vector<double> &v2sigma2_5, std::vector<double> &v2sigma2_6)
    {
        int np = rho_up.size();
        std::vector<double> rho(2 * np);
        std::vector<double> sigma(3 * np);
        std::vector<double> v2rho2(3 * np);
        std::vector<double> v2rhosigma(6 * np);
        std::vector<double> v2sigma2(6 * np);
        for (int i = 0; i < np; ++i)
        {
            rho[2 * i] = rho_up[i];
            rho[2 * i + 1] = rho_down[i];
            sigma[3 * i] = sigma_1[i];
            sigma[3 * i + 1] = sigma_2[i];
            sigma[3 * i + 2] = sigma_3[i];
        }
        xc_gga_fxc(&func, np, rho.data(), sigma.data(), v2rho2.data(), v2rhosigma.data(), v2sigma2.data());
        for (int i = 0; i < np; ++i)
        {
            v2rho2_1[i] = v2rho2[3 * i];
            v2rho2_2[i] = v2rho2[3 * i + 1];
            v2rho2_3[i] = v2rho2[3 * i + 2];
            v2rhosigma_1[i] = v2rhosigma[6 * i];
            v2rhosigma_2[i] = v2rhosigma[6 * i + 1];
            v2rhosigma_3[i] = v2rhosigma[6 * i + 2];
            v2rhosigma_4[i] = v2rhosigma[6 * i + 3];
            v2rhosigma_5[i] = v2rhosigma[6 * i + 4];
            v2rhosigma_6[i] = v2rhosigma[6 * i + 5];
            v2sigma2_1[i] = v2sigma2[6 * i];
            v2sigma2_2[i] = v2sigma2[6 * i + 1];
            v2sigma2_3[i] = v2sigma2[6 * i + 2];
            v2sigma2_4[i] = v2sigma2[6 * i + 3];
            v2sigma2_5[i] = v2sigma2[6 * i + 4];
            v2sigma2_6[i] = v2sigma2[6 * i + 5];
        }
    }
    // GGA END

    // Example usage for spin-polarized LDA
    // The input rho_up and rho_down are the densities for spin-up and spin-down electrons, respectively.
    // rhoup = (n+s)/2, rhodown = (n-s)/2
    //output:
    //exc[]: the energy density per unit particle
    //vrho[]: first derivative of the energy per unit volume
    // v2rho2[]: second derivative of the energy per unit volume
    //v3rho3[]: third derivative of the energy per unit volume
    // v4rho4[]: fourth derivative of the energy per unit volume
    // detailed expression can be found in https://libxc.gitlab.io/manual/libxc-5.1.x/
    //vrho [(0),(1)]
    //v2rho2 [(0,0),(0,1),(1,1)]
    //v3rho3 [(0,0,0),(0,0,1),(0,1,1),(1,1,1)]
    //v4rho4 [(0,0,0,0),(0,0,0,1),(0,0,1,1),(0,1,1,1),(1,1,1,1)]
    // 0:up 1:down (spin)

    void LibxcInterface::example_lda_spin()
    {
        std::vector<double> rho_up = {0.1, 0.2, 0.3};
        std::vector<double> rho_down = {0.2, 0.1, 0.3};

        auto exc = lda_exc(rho_up, rho_down);
        std::vector<double> vrho_1(rho_up.size()), vrho_2(rho_down.size());
        lda_vxc(rho_up, rho_down, vrho_1, vrho_2);

        std::vector<double> v2rho2_1(rho_up.size()), v2rho2_2(rho_down.size()), v2rho2_3(rho_up.size());
        lda_fxc(rho_up, rho_down, v2rho2_1, v2rho2_2, v2rho2_3);

        //        std::vector<double> v3rho3_1(rho_up.size()), v3rho3_2(rho_down.size()), v3rho3_3(rho_up.size()), v3rho3_4(rho_down.size());
        //        lda_kxc(rho_up, rho_down, v3rho3_1, v3rho3_2, v3rho3_3, v3rho3_4);

        //        std::vector<double> v4rho4_1(rho_up.size()), v4rho4_2(rho_down.size()), v4rho4_3(rho_up.size()), v4rho4_4(rho_down.size()), v4rho4_5(rho_up.size());
        //        lda_lxc(rho_up, rho_down, v4rho4_1, v4rho4_2, v4rho4_3, v4rho4_4, v4rho4_5);

        std::cout << "Exc: ";
        for (auto e : exc)
            std::cout << e << " ";
        std::cout << std::endl;

        std::cout << "VRho_1: ";
        for (auto v : vrho_1)
            std::cout << v << " ";
        std::cout << std::endl;

        std::cout << "VRho_2: ";
        for (auto v : vrho_2)
            std::cout << v << " ";
        std::cout << std::endl;

        std::cout << "V2Rho2_1: ";
        for (auto f : v2rho2_1)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "V2Rho2_2: ";
        for (auto f : v2rho2_2)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "V2Rho2_3: ";
        for (auto f : v2rho2_3)
            std::cout << f << " ";
        std::cout << std::endl;

        //        std::cout << "V3Rho3_1: ";
        //        for (auto k : v3rho3_1) std::cout << k << " ";
        //        std::cout << std::endl;

        //        std::cout << "V3Rho3_2: ";
        //        for (auto k : v3rho3_2) std::cout << k << " ";
        //        std::cout << std::endl;

        //        std::cout << "V3Rho3_3: ";
        //        for (auto k : v3rho3_3) std::cout << k << " ";
        //        std::cout << std::endl;

        //        std::cout << "V3Rho3_4: ";
        //        for (auto k : v3rho3_4) std::cout << k << " ";
        //        std::cout << std::endl;

        //        std::cout << "V4Rho4_1: ";
        //        for (auto l : v4rho4_1) std::cout << l << " ";
        //        std::cout << std::endl;

        //        std::cout << "V4Rho4_2: ";
        //        for (auto l : v4rho4_2) std::cout << l << " ";
        //        std::cout << std::endl;

        //        std::cout << "V4Rho4_3: ";
        //        for (auto l : v4rho4_3) std::cout << l << " ";
        //        std::cout << std::endl;

        //        std::cout << "V4Rho4_4: ";
        //        for (auto l : v4rho4_4) std::cout << l << " ";
        //        std::cout << std::endl;

        //        std::cout << "V4Rho4_5: ";
        //        for (auto l : v4rho4_5) std::cout << l << " ";
        //        std::cout << std::endl;
    }

    // // Example usage for spin-polarized GGA
    //    input:
    //      np: number of points
    //       rho[]: the density
    //      sigma[]: contracted gradients of the density
    //     output:
    //      exc[]: the energy per unit particle
    //      vrho[]: first partial derivative of the energy per unit volume in terms of the density
    //      vsigma[]: first partial derivative of the energy per unit volume in terms of sigma
    //      v2rho2[]: second partial derivative of the energy per unit volume in terms of the density
    //      v2rhosigma[]: second partial derivative of the energy per unit volume in terms of the density and sigma
    //      v2sigma2[]: second partial derivative of the energy per unit volume in terms of sigma
    //    // detailed expression can be found in https://libxc.gitlab.io/manual/libxc-5.1.x/

    void LibxcInterface::example_gga_spin()
    {
        // 定义自旋极化系统的示例密度和梯度数据
        std::vector<double> rho_up = {0.1, 0.2, 0.3};
        std::vector<double> rho_down = {0.1, 0.2, 0.3};
        std::vector<double> sigma_1 = {0.01, 0.02, 0.03};
        std::vector<double> sigma_2 = {0.02, 0.01, 0.03};
        std::vector<double> sigma_3 = {0.01, 0.02, 0.03};

        // 计算交换-相关能量密度
        auto exc = gga_exc(rho_up, rho_down, sigma_1, sigma_2, sigma_3);

        // 计算交换-相关势
        std::vector<double> vrho_1(rho_up.size()), vrho_2(rho_down.size());
        std::vector<double> vsigma_1(sigma_1.size()), vsigma_2(sigma_2.size()), vsigma_3(sigma_3.size());
        gga_vxc(rho_up, rho_down, sigma_1, sigma_2, sigma_3, vrho_1, vrho_2, vsigma_1, vsigma_2, vsigma_3);

        // 计算交换-相关势的二阶导数
        std::vector<double> v2rho2_1(rho_up.size()), v2rho2_2(rho_down.size()), v2rho2_3(rho_up.size());
        std::vector<double> v2rhosigma_1(rho_up.size()), v2rhosigma_2(rho_up.size()), v2rhosigma_3(rho_up.size()), v2rhosigma_4(rho_up.size()), v2rhosigma_5(rho_up.size()), v2rhosigma_6(rho_up.size());
        std::vector<double> v2sigma2_1(sigma_1.size()), v2sigma2_2(sigma_1.size()), v2sigma2_3(sigma_1.size()), v2sigma2_4(sigma_1.size()), v2sigma2_5(sigma_1.size()), v2sigma2_6(sigma_1.size());
        gga_fxc(rho_up, rho_down, sigma_1, sigma_2, sigma_3, v2rho2_1, v2rho2_2, v2rho2_3, v2rhosigma_1, v2rhosigma_2, v2rhosigma_3, v2rhosigma_4, v2rhosigma_5, v2rhosigma_6, v2sigma2_1, v2sigma2_2, v2sigma2_3, v2sigma2_4, v2sigma2_5, v2sigma2_6);

        // 输出结果
        std::cout << "GGA Exc: ";
        for (auto e : exc)
            std::cout << e << " ";
        std::cout << std::endl;

        std::cout << "GGA VRho_1: ";
        for (auto v : vrho_1)
            std::cout << v << " ";
        std::cout << std::endl;

        std::cout << "GGA VRho_2: ";
        for (auto v : vrho_2)
            std::cout << v << " ";
        std::cout << std::endl;

        std::cout << "GGA VSigma_1: ";
        for (auto v : vsigma_1)
            std::cout << v << " ";
        std::cout << std::endl;

        std::cout << "GGA VSigma_2: ";
        for (auto v : vsigma_2)
            std::cout << v << " ";
        std::cout << std::endl;

        std::cout << "GGA VSigma_3: ";
        for (auto v : vsigma_3)
            std::cout << v << " ";
        std::cout << std::endl;

        std::cout << "GGA V2Rho2_1: ";
        for (auto f : v2rho2_1)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2Rho2_2: ";
        for (auto f : v2rho2_2)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2Rho2_3: ";
        for (auto f : v2rho2_3)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2RhoSigma_1: ";
        for (auto f : v2rhosigma_1)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2RhoSigma_2: ";
        for (auto f : v2rhosigma_2)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2RhoSigma_3: ";
        for (auto f : v2rhosigma_3)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2RhoSigma_4: ";
        for (auto f : v2rhosigma_4)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2RhoSigma_5: ";
        for (auto f : v2rhosigma_5)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2RhoSigma_6: ";
        for (auto f : v2rhosigma_6)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2Sigma2_1: ";
        for (auto f : v2sigma2_1)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2Sigma2_2: ";
        for (auto f : v2sigma2_2)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2Sigma2_3: ";
        for (auto f : v2sigma2_3)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2Sigma2_4: ";
        for (auto f : v2sigma2_4)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2Sigma2_5: ";
        for (auto f : v2sigma2_5)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2Sigma2_6: ";
        for (auto f : v2sigma2_6)
            std::cout << f << " ";
        std::cout << std::endl;
    }




    /////////////////////////////////////////////////////////////
    //meta-GGA START, up to the second derivative
    // meta-GGA Energy Density for spin-polarized systems
    /*
    std::vector<double> gga_exc(const std::vector<double> &rho_up, const std::vector<double> &rho_down, const std::vector<double> &sigma_1, const std::vector<double> &sigma_2, const std::vector<double> &sigma_3)
    {
        int np = rho_up.size();
        std::vector<double> rho(2 * np);
        std::vector<double> sigma(3 * np);
        for (int i = 0; i < np; ++i)
        {
            rho[2 * i] = rho_up[i];
            rho[2 * i + 1] = rho_down[i];
            sigma[3 * i] = sigma_1[i];
            sigma[3 * i + 1] = sigma_2[i];
            sigma[3 * i + 2] = sigma_3[i];
        }
        std::vector<double> exc(np);
        xc_gga_exc(&func, np, rho.data(), sigma.data(), exc.data());
        return exc;
    }
*/










//////////////////////////////////////////////////////////////////////////////////////////////////

// MAIN FOR INDENPENDENT TEST
/*
int main() 
{
    // //LDA  test
   try 
   {
       LibxcInterface libxc(1, true); // Example with spin-polarized LDA (xc_id = 1)
       libxc.example_lda_spin();
   } 
   catch (const std::exception& ex) 
   {
       std::cerr << "Error: " << ex.what() << std::endl;
   }
    std::cout << "####################################################\n" << std::endl;

    // //GGA  test
    try
    {
        LibxcInterface libxc(101, true); // Example with spin-polarized GGA (xc_id = 101)
        libxc.example_gga_spin();
    }
    catch (const std::exception &ex)
    {
        std::cerr << "Error: " << ex.what() << std::endl;
    }
    std::cout << "####################################################\n" << std::endl;

//     // //meta-GGA test
//    try 
//    {
//        LibxcInterface libxc(202, true); // Example with spin-polarized meta-GGA (xc_id = 202)
//        libxc.example_mgga_spin();
//    } 
//    catch (const std::exception& ex) 
//    {
//        std::cerr << "Error: " << ex.what() << std::endl;
//    }

    return 0;
}

*/

//////////////////////////////////////////////////////////////////////////////////////////////////

