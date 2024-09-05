// Programmed by Xiaoyu Zhang at Peking University, Beijing, China 2024/08/26

#ifndef INTERFACE_TO_LIBXC_H
#define INTERFACE_TO_LIBXC_H

#include <xc.h>
#include <vector>
#include <stdexcept>

class LibxcInterface
{
private:
    xc_func_type func;

public:
    // Constructor with spin polarization
    LibxcInterface(int xc_id, bool spin_polarized = false);

    // Destructor
    ~LibxcInterface();

    std::vector<double> lda_exc(const std::vector<double> &rho_up, const std::vector<double> &rho_down);
    void lda_vxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 std::vector<double> &vrho_1, std::vector<double> &vrho_2);
    void lda_fxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 std::vector<double> &v2rho2_1, std::vector<double> &v2rho2_2, std::vector<double> &v2rho2_3);
    void lda_kxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 std::vector<double> &v3rho3_1, std::vector<double> &v3rho3_2, std::vector<double> &v3rho3_3, std::vector<double> &v3rho3_4);
    void lda_lxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 std::vector<double> &v4rho4_1, std::vector<double> &v4rho4_2, std::vector<double> &v4rho4_3, std::vector<double> &v4rho4_4, std::vector<double> &v4rho4_5);

    std::vector<double> gga_exc(const std::vector<double> &rho_up, const std::vector<double> &rho_down, const std::vector<double> &sigma_1, const std::vector<double> &sigma_2, const std::vector<double> &sigma_3);
    void gga_vxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 const std::vector<double> &sigma_1, const std::vector<double> &sigma_2, const std::vector<double> &sigma_3, std::vector<double> &vrho_1, std::vector<double> &vrho_2,
                 std::vector<double> &vsigma_1, std::vector<double> &vsigma_2, std::vector<double> &vsigma_3);
    void gga_fxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 const std::vector<double> &sigma_1, const std::vector<double> &sigma_2, const std::vector<double> &sigma_3,
                 std::vector<double> &v2rho2_1, std::vector<double> &v2rho2_2, std::vector<double> &v2rho2_3,
                 std::vector<double> &v2rhosigma_1, std::vector<double> &v2rhosigma_2, std::vector<double> &v2rhosigma_3, std::vector<double> &v2rhosigma_4, std::vector<double> &v2rhosigma_5, std::vector<double> &v2rhosigma_6, std::vector<double> &v2sigma2_1, std::vector<double> &v2sigma2_2, std::vector<double> &v2sigma2_3, std::vector<double> &v2sigma2_4, std::vector<double> &v2sigma2_5, std::vector<double> &v2sigma2_6);
    

    void example_lda_spin();
    void example_gga_spin();
};

#endif // INTERFACE_TO_LIBXC_H