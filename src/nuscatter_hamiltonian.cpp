#include <iostream>
#include <numeric>
#include <vector>
#include <cmath>
#include <complex>
#include "misc/matrixop.hpp"
#include "misc/vector.hpp"
#include "misc/crasher.hpp"
#include "nuscatter_hamiltonian.hpp"

namespace mqc {


    // --- ctor/dtor --- //


    Nuscatter_Hamiltonian::Nuscatter_Hamiltonian() 
    : Hamiltonian(2) { 
        m_doc = 
        "# Yanze's nuclear scattering Hamiltoinan \n"
        "# H00 = A * [1+tanh(eps*(m_theta-pi/2))] + B * [2 + tanh(alpha*(m_r-R)) - tanh(alpha*(m_r - r))] \n"
        "# H11 = A * [1-tanh(eps*(m_theta-pi/2))] + B * [2 + tanh(alpha*(m_r-R)) - tanh(alpha*(m_r - r))] \n"
        "# H01 = C * exp[-(eps*(m_theta-pi/2))^2 + i*W*m_r] * 1/2 * [tanh(alpha*(m_r-r)) - tanh(alpha*(m_r-R))] \n"
        "# m_r = sqrt((x+R/2)^2 + (y+R/2)^2) \n"
        "# m_theta = atan2(y+R/2, x+R/2) \n"
        "# parameters: "
        "# { A, C, R, W, alpha, eps, forbidden, r } \n"
        ;
        m_params["A"] = 0.02;
        m_params["C"] = 0.005;
        m_params["R"] = 8.0;
        m_params["W"] = 0.0;
        m_params["alpha"] = 40.0;
        m_params["eps"] = 5.0;
        m_params["forbidden"] = 10.0;
        m_params["r"] = 4.0;
    }


    // --- util --- //


    double Nuscatter_Hamiltonian::cal_m_r(const std::vector<double>& r) const {
        return sqrt( pow((r[0] + 0.5 * m_params.at("R")), 2) + pow((r[1] + 0.5 * m_params.at("R")), 2) );
    }

    std::vector<double> Nuscatter_Hamiltonian::cal_der_m_r(const std::vector<double>& r) const {
        return std::vector<double> { r[0] + 0.5 * m_params.at("R"), r[1] + 0.5 * m_params.at("R") } / cal_m_r(r);
    }

    double Nuscatter_Hamiltonian::cal_m_theta(const std::vector<double>& r) const {
        return atan2(r[1] + 0.5 * m_params.at("R"), r[0] + 0.5 * m_params.at("R"));
    }

    std::vector<double> Nuscatter_Hamiltonian::cal_der_m_theta(const std::vector<double>& r) const {
        const double x1 = r[0] + 0.5 * m_params.at("R");
        const double y1 = r[1] + 0.5 * m_params.at("R");
        return std::vector<double> { -y1, x1 } / (x1*x1 + y1*y1);
    }

    double Nuscatter_Hamiltonian::cal_phi(const std::vector<double>& r) const {
        return m_params.at("W") * cal_m_r(r);
    }

    std::vector<double> Nuscatter_Hamiltonian::cal_der_phi(const std::vector<double>& r) const {
        return m_params.at("W") * cal_der_m_r(r);
    }

    double Nuscatter_Hamiltonian::cal_param_coup() const {
        return m_params.at("eps") * m_params.at("eps");
    }

    // --- interfaces --- //

    std::vector<std::complex<double>> Nuscatter_Hamiltonian::cal_H(const std::vector<double>& r) const {
        // check 
        misc::confirm<misc::ValueError>(r.size() >= 2, "nablaH: r must have a size >= 2.");
        // setup
        const double m_r = cal_m_r(r);
        const double m_theta = cal_m_theta(r);
        const std::complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r));
        const double coup = cal_param_coup();
        std::vector<std::complex<double>> H(m_dim * m_dim);
        // H
        H.at(0+0*m_dim) = m_params.at("A") * (1.0 + tanh(m_params.at("eps") * (m_theta - M_PI/4))) 
                    + m_params.at("forbidden") * (2.0 + tanh(m_params.at("alpha") * (m_r - m_params.at("R"))) - tanh(m_params.at("alpha") * (m_r - m_params.at("r")))); 
        H.at(1+1*m_dim) = m_params.at("A") * (1.0 - tanh(m_params.at("eps") * (m_theta - M_PI/4))) 
                    + m_params.at("forbidden") * (2.0 + tanh(m_params.at("alpha") * (m_r - m_params.at("R"))) - tanh(m_params.at("alpha") * (m_r - m_params.at("r")))); 
        H.at(0+1*m_dim) = m_params.at("C") * exp(-coup* pow((m_theta - M_PI/4), 2)) * eip
                    * 0.5 * (tanh(m_params.at("alpha") * (m_r - m_params.at("r"))) - tanh(m_params.at("alpha") * (m_r - m_params.at("R"))));
        H.at(1+0*m_dim) = conj(H.at(0+1*m_dim));
        return H;
    }

    std::vector<std::vector<std::complex<double>>> Nuscatter_Hamiltonian::cal_nablaH(const std::vector<double>& r) const {
        /**
         * calculate nablaH for given r
         */
        // check 
        misc::confirm<misc::ValueError>(r.size() >= 2, "cal_nablaH: r must have a size >= 2.");
        // setup
        const double m_r = cal_m_r(r);
        const double m_theta = cal_m_theta(r);
        const std::vector<double> der_m_r = cal_der_m_r(r);
        const std::vector<double> der_m_theta = cal_der_m_theta(r);
        const std::complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r));
        const std::vector<double> der_phi = cal_der_phi(r);
        const double coup = cal_param_coup();
        std::vector<std::vector<std::complex<double>>> nablaH(r.size(), std::vector<std::complex<double>>(m_dim * m_dim, matrixop::ZEROZ));
        auto der_tanh = [](double x) { 
            return (1.0 - pow(tanh(x), 2)); 
        };
        // nablaHx
        std::vector<std::complex<double>>& nablaHx = nablaH.at(0);
        nablaHx.at(0+0*m_dim) = m_params.at("A") * der_tanh(m_params.at("eps") * (m_theta - M_PI/4)) * m_params.at("eps") * der_m_theta.at(0)
                            + m_params.at("forbidden") * ( 
                                    der_tanh(m_params.at("alpha") * (m_r - m_params.at("R"))) * m_params.at("alpha") * der_m_r.at(0) 
                                    - der_tanh(m_params.at("alpha") * (m_r - m_params.at("r"))) * m_params.at("alpha") * der_m_r.at(0) 
                                  );
        nablaHx.at(1+1*m_dim) = m_params.at("A") * -der_tanh(m_params.at("eps") * (m_theta - M_PI/4)) * m_params.at("eps") * der_m_theta.at(0)
                            + m_params.at("forbidden") * ( 
                                    der_tanh((m_r - m_params.at("R")) * m_params.at("alpha")) * m_params.at("alpha") * der_m_r.at(0) 
                                    - der_tanh((m_r - m_params.at("r")) * m_params.at("alpha")) * m_params.at("alpha") * der_m_r.at(0) 
                                  );
        // (e^{i*phi} * f)' = eip * (i*phi'*f + f')
        nablaHx.at(0+1*m_dim) = eip * (
                matrixop::IMAGIZ * der_phi.at(0) 
                    * m_params.at("C") * exp(-coup * pow((m_theta - M_PI/4), 2)) 
                    * 0.5 * (tanh(m_params.at("alpha") * (m_r - m_params.at("r"))) - tanh(m_params.at("alpha") * (m_r - m_params.at("R"))))
                + m_params.at("C") * exp(-coup * pow((m_theta - M_PI/4), 2)) * (-2.0) * coup * (m_theta - M_PI/4) * der_m_theta.at(0)
                    * 0.5 * (tanh(m_params.at("alpha") * (m_r - m_params.at("r"))) - tanh(m_params.at("alpha") * (m_r - m_params.at("R"))))
                + m_params.at("C") * exp(-coup * pow((m_theta - M_PI/4), 2))
                    * 0.5 * (der_tanh(m_params.at("alpha") * (m_r - m_params.at("r"))) * m_params.at("alpha") * der_m_r.at(0) - der_tanh(m_params.at("alpha") * (m_r - m_params.at("R"))) * m_params.at("alpha") * der_m_r.at(0))
                );
        nablaHx.at(1+0*2) = conj(nablaHx.at(0+1*2));
        // nablaHy
        std::vector<std::complex<double>>& nablaHy = nablaH.at(1);
        nablaHy.at(0+0*2) = m_params.at("A") * der_tanh(m_params.at("eps") * (m_theta - M_PI/4)) * m_params.at("eps") * der_m_theta.at(1)
                            + m_params.at("forbidden") * ( 
                                    der_tanh(m_params.at("alpha") * (m_r - m_params.at("R"))) * m_params.at("alpha") * der_m_r.at(1) 
                                    - der_tanh(m_params.at("alpha") * (m_r - m_params.at("r"))) * m_params.at("alpha") * der_m_r.at(1) 
                                  );
        nablaHy.at(1+1*2) = m_params.at("A") * -der_tanh(m_params.at("eps") * (m_theta - M_PI/4)) * m_params.at("eps") * der_m_theta.at(1)
                            + m_params.at("forbidden") * ( 
                                    der_tanh((m_r - m_params.at("R")) * m_params.at("alpha")) * m_params.at("alpha") * der_m_r.at(1) 
                                    - der_tanh((m_r - m_params.at("r")) * m_params.at("alpha")) * m_params.at("alpha") * der_m_r.at(1) 
                                  );
        // (e^{i*phi} * f)' = eip * (i*phi'*f + f')
        nablaHy.at(0+1*2) = eip * (
                matrixop::IMAGIZ * der_phi.at(1) 
                    * m_params.at("C") * exp(-coup * pow((m_theta - M_PI/4), 2)) 
                    * 0.5 * (tanh(m_params.at("alpha") * (m_r - m_params.at("r"))) - tanh(m_params.at("alpha") * (m_r - m_params.at("R"))))
                + m_params.at("C") * exp(-coup * pow((m_theta - M_PI/4), 2)) * (-2.0) * coup * (m_theta - M_PI/4) * der_m_theta.at(1)
                    * 0.5 * (tanh(m_params.at("alpha") * (m_r - m_params.at("r"))) - tanh(m_params.at("alpha") * (m_r - m_params.at("R"))))
                + m_params.at("C") * exp(-coup * pow((m_theta - M_PI/4), 2))
                    * 0.5 * (der_tanh(m_params.at("alpha") * (m_r - m_params.at("r"))) * m_params.at("alpha") * der_m_r.at(1) - der_tanh(m_params.at("alpha") * (m_r - m_params.at("R"))) * m_params.at("alpha") * der_m_r.at(1))
                );
        nablaHy.at(1+0*2) = conj(nablaHy.at(0+1*2));
        // return 
        return nablaH;
    }

} // namespace mqc