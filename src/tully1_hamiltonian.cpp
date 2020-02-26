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
    : Hamiltonian(2)
    { 
        m_params["A"] = 0.01;
        m_params["B"] = 1.6;
        m_params["C"] = 0.005;
        m_params["D"] = 1.0;
        m_params["Wx"] = 0.0;
        m_params["Wy"] = 0.0;
    }


    // --- quantities --- //


    double Tully1_Hamiltonian::cal_phi(const std::vector<double>& r) const {
        return m_params.at("Wx") * r[0] + m_params.at("Wy") * r[1];
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
        const double x = r[0];
        const double y = r[1];
        const std::complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r));
        std::vector<std::complex<double>> H(m_dim * m_dim, matrixop::ZEROZ);
        // H
        if (x >= 0.0) {
            H[0+0*2] = m_params.at("A") * (1.0 - exp(-m_params.at("B") * x));
        }
        else {
            H[0+0*2] = -m_params.at("A") * (1.0 - exp(m_params.at("B") * x));
        }
        H[1+1*2] = -H[0+0*2];
        H[0+1*2] = m_params.at("C") * exp(-m_params.at("D") * x * x) * eip;
        H[1+0*2] = conj(H[0+1*2]);
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
        std::vector<std::complex<double>>& nablaHx = nablaH[0];
        nablaHx.resize(m_dim * m_dim);
        if (x >= 0.0 ) {
            nablaHx[0+0*2] = m_params.at("A") * m_params.at("B") * exp(-m_params.at("B") * x);
        }
        else {
            nablaHx[0+0*2] = m_params.at("A") * m_params.at("B") * exp(m_params.at("B") * x);
        }
        nablaHx[1+1*2] = -nablaHx[0+0*2];
        nablaHx[0+1*2] = eip * (-2.0 * m_params.at("C") * m_params.at("D") * x * exp(-m_params.at("D") * x * x) 
                                + matrixop::IMAGIZ * m_params.at("C") * exp(-m_params.at("D") * x * x) * nabla_phi[0]);
        nablaHx[1+0*2] = conj(nablaHx[0+1*2]);
        // dH/dy
        std::vector<std::complex<double>>& nablaHy = nablaH[1];
        nablaHy.resize(m_dim * m_dim);
        nablaHy[0+0*2] = 0.0;
        nablaHy[1+1*2] = -nablaHy[0+0*2];
        nablaHy[0+1*2] = eip * (matrixop::IMAGIZ * m_params.at("C") * exp(-m_params.at("D") * x * x) * nabla_phi[1]);
        nablaHy[1+0*2] = conj(nablaHy[0+1*2]);

        return nablaH;
    }

};
