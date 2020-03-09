#include <iostream>
#include <vector>
#include <complex>
#include "misc/matrixop.hpp"
#include "misc/vector.hpp"
#include "tully1_hamiltonian.hpp"
#include "misc/crasher.hpp"

namespace mqc {


    // --- ctor/dtor/operator= --- //


    Tully1_Hamiltonian::Tully1_Hamiltonian() 
    : Hamiltonian(2) { 
        m_doc = 
        "# Classical 1D Tully1 Hamiltoinan \n"
        "# H_{00} = A * (1 - exp(-Bx)), for x >= 0 \n"
        "# H_{00} = -A * (1 - exp(Bx)), for x < 0 \n"
        "# H_{11} = -H_{00} \n"
        "# H_{01} = C * exp(-D * x^2) \n"
        "# H_{10} = conj(H_{01}) \n"
        "# paramters: "
        "# { A, B, C, D } \n"
        ;
        m_params["A"] = 0.01;
        m_params["B"] = 1.6;
        m_params["C"] = 0.005;
        m_params["D"] = 1.0;
    }


    // --- interfaces --- //


    std::vector<std::complex<double>> Tully1_Hamiltonian::cal_H(const std::vector<double>& r) const {
        /**
         * input : r
         * output : H(r)
         */
        // check
        misc::confirm(r.size() == 1, "Tully1_Hamiltonian: nuclear dim must be 1.");
        // params
        const double x = r.at(0);
        std::vector<std::complex<double>> H(m_dim * m_dim, matrixop::ZEROZ);
        // H
        if (x >= 0.0) {
            H.at(0+0*2) = m_params.at("A") * (1.0 - exp(-m_params.at("B") * x));
        }
        else {
            H.at(0+0*2) = -m_params.at("A") * (1.0 - exp(m_params.at("B") * x));
        }
        H.at(1+1*2) = -H.at(0+0*2);
        H.at(0+1*2) = m_params.at("C") * exp(-m_params.at("D") * x * x);
        H.at(1+0*2) = conj(H.at(0+1*2));
        return H;
    }

    std::vector<std::vector<std::complex<double>>> Tully1_Hamiltonian::cal_nablaH(const std::vector<double>& r) const {
        /**
         * input : r
         * output : nablaH(r)
         */
        // check
        misc::confirm(r.size() == 1, "Tully1_Hamiltonian: nuclear dim must be 1.");
        // params
        const double x = r.at(0);
        std::vector<std::vector<std::complex<double>>> nablaH(1);
        // dH/dx
        std::vector<std::complex<double>>& nablaHx = nablaH.at(0);
        nablaHx.resize(m_dim * m_dim);
        if (x >= 0.0 ) {
            nablaHx.at(0+0*2) = m_params.at("A") * m_params.at("B") * exp(-m_params.at("B") * x);
        }
        else {
            nablaHx.at(0+0*2) = m_params.at("A") * m_params.at("B") * exp(m_params.at("B") * x);
        }
        nablaHx.at(1+1*2) = -nablaHx.at(0+0*2);
        nablaHx.at(0+1*2) = -2.0 * m_params.at("C") * m_params.at("D") * x * exp(-m_params.at("D") * x * x);
        nablaHx.at(1+0*2) = conj(nablaHx.at(0+1*2));
        return nablaH;
    }

} // namespace mqc
