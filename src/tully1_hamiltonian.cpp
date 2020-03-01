#include <iostream>
#include <vector>
#include <complex>
#include "misc/matrixop.hpp"
#include "misc/vector.hpp"
#include "tully1_hamiltonian.hpp"
#include "misc/crasher.hpp"

namespace mqc {


    // --- ctor/dtor --- //


    Tully1_Hamiltonian::Tully1_Hamiltonian() 
    : Hamiltonian(2) { 
        m_doc = 
        "# 2D Tully1 Hamiltoinan \n"
        "# H_{00} = A * (1 - exp(-Bx)), for x >= 0 \n"
        "# H_{00} = -A * (1 - exp(Bx)), for x < 0 \n"
        "# H_{11} = -H_{00} \n"
        "# H_{01} = C * exp(-D * x^2) * exp(i * phi) \n"
        "# H_{10} = conj(H_{01}) \n"
        "# phi = Wx * x + Wy * y \n"
        "# paramters: "
        "# { A, B, C, D, Wx, Wy } \n"
        ;
        m_params.at("A") = 0.01;
        m_params.at("B") = 1.6;
        m_params.at("C") = 0.005;
        m_params.at("D") = 1.0;
        m_params.at("Wx") = 0.0;
        m_params.at("Wy") = 0.0;
    }


    // --- quantities --- //


    double Tully1_Hamiltonian::cal_phi(const std::vector<double>& r) const {
        return m_params.at("Wx") * r.at(0) + m_params.at("Wy") * r.at(1);
    }

    std::vector<double> Tully1_Hamiltonian::cal_nabla_phi(const std::vector<double>& r) const {
        return std::vector<double> { m_params.at("Wx"), m_params.at("Wy") };
    }

    std::vector<std::complex<double>> Tully1_Hamiltonian::cal_H(const std::vector<double>& r) const {
        /*
         * input : r
         * output : H(r)
         */
        // check
        misc::confirm(r.size() == 2, "Tully1_Hamiltonian: nuclear dim must be 2.");
        // params
        const double x = r.at(0);
        const double y = r.at(1);
        const std::complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r));
        std::vector<std::complex<double>> H(m_dim * m_dim, matrixop::ZEROZ);
        // H
        if (x >= 0.0) {
            H.at(0+0*2) = m_params.at("A") * (1.0 - exp(-m_params.at("B") * x));
        }
        else {
            H.at(0+0*2) = -m_params.at("A") * (1.0 - exp(m_params.at("B") * x));
        }
        H.at(1+1*2) = -H.at(0+0*2);
        H.at(0+1*2) = m_params.at("C") * exp(-m_params.at("D") * x * x) * eip;
        H.at(1+0*2) = conj(H.at(0+1*2));
        return H;
    }

    std::vector<std::vector<std::complex<double>>> Tully1_Hamiltonian::cal_nablaH(const std::vector<double>& r) const {
        /*
         * input : r
         * output : nablaH(r)
         */
        // check
        misc::confirm(r.size() == 2, "Tully1_Hamiltonian: nuclear dim must be 2.");
        // params
        const double x = r.at(0);
        const double y = r.at(1);
        const std::complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r));
        const std::vector<double> nabla_phi = cal_nabla_phi(r);
        std::vector<std::vector<std::complex<double>>> nablaH(2);
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
        nablaHx.at(0+1*2) = eip * (-2.0 * m_params.at("C") * m_params.at("D") * x * exp(-m_params.at("D") * x * x) 
                                + matrixop::IMAGIZ * m_params.at("C") * exp(-m_params.at("D") * x * x) * nabla_phi.at(0));
        nablaHx.at(1+0*2) = conj(nablaHx.at(0+1*2));
        // dH/dy
        std::vector<std::complex<double>>& nablaHy = nablaH.at(1);
        nablaHy.resize(m_dim * m_dim);
        nablaHy.at(0+0*2) = 0.0;
        nablaHy.at(1+1*2) = -nablaHy.at(0+0*2);
        nablaHy.at(0+1*2) = eip * (matrixop::IMAGIZ * m_params.at("C") * exp(-m_params.at("D") * x * x) * nabla_phi.at(1));
        nablaHy.at(1+0*2) = conj(nablaHy.at(0+1*2));

        return nablaH;
    }

};
