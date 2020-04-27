#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include "misc/matrixop.hpp"
#include "misc/vector.hpp"
#include "flat_2d_hamiltonian.hpp"
#include "misc/crasher.hpp"


namespace mqc {


    // --- ctor/dtor/operator= --- //


    Flat_2D_Hamiltonian::Flat_2D_Hamiltonian() 
    : Hamiltonian(2) { 
        m_doc = 
        "# 2D Flat Hamiltoinan \n"
        "# H_{00} = -A * cos(theta) \n"
        "# H_{11} = A * cos(theta) \n"
        "# H_{01} = A * sin(theta) * exp(i*phi)\n"
        "# H_{10} = conj(H_{01}) \n"
        "# theta = pi/2 * (erf(B*x) + 1) \n"
        "# phi = Wx * x + Wy * y \n"
        "# paramters: "
        "# { A, B, Wx, Wy } \n"
        ;
        m_params["A"] = 0.1;
        m_params["B"] = 3.0;
        m_params["Wx"] = 0.0;
        m_params["Wy"] = 0.0;
    }


    // --- utils --- //

    double Flat_2D_Hamiltonian::cal_theta(const std::vector<double>& r) const {
        return 0.5 * M_PI * (std::erf(m_params.at("B") * r.at(0)) + 1.0);
    }

    std::vector<double> Flat_2D_Hamiltonian::cal_nabla_theta(const std::vector<double>& r) const {
        std::vector<double> nabla_theta(r.size(), 0.0);
        nabla_theta.at(0) = std::sqrt(M_PI) * m_params.at("B") * std::exp(-std::pow(m_params.at("B") * r.at(0), 2));
        return nabla_theta;
    }

    double Flat_2D_Hamiltonian::cal_phi(const std::vector<double>& r) const {
        return m_params.at("Wx") * r.at(0) + m_params.at("Wy") * r.at(1);
    }

    std::vector<double> Flat_2D_Hamiltonian::cal_nabla_phi(const std::vector<double>& r) const {
        return std::vector<double> { m_params.at("Wx"), m_params.at("Wy") };
    }


    // --- interfaces --- //


    std::vector<std::complex<double>> Flat_2D_Hamiltonian::cal_H(const std::vector<double>& r) const {
        /**
         * input : r
         * output : H(r)
         */
        // check
        misc::confirm(r.size() == 2, "Flat_2D_Hamiltonian: nuclear dim must be 2.");
        // params
        const double theta = cal_theta(r);
        const std::complex<double> eip = std::exp(matrixop::IMAGIZ * cal_phi(r));
        std::vector<std::complex<double>> H(m_dim * m_dim, matrixop::ZEROZ);
        H.at(0+0*m_dim) = -m_params.at("A") * std::cos(theta);
        H.at(1+1*m_dim) = -H.at(0+0*m_dim);
        H.at(0+1*m_dim) = m_params.at("A") * std::sin(theta) * eip;
        H.at(1+0*m_dim) = std::conj(H.at(0+1*m_dim));
        return H;
    }

    std::vector<std::vector<std::complex<double>>> Flat_2D_Hamiltonian::cal_nablaH(const std::vector<double>& r) const {
        /**
         * input : r
         * output : nablaH(r)
         */
        // check
        misc::confirm(r.size() == 2, "Flat_2D_Hamiltonian: nuclear dim must be 2.");
        // params
        const double theta = cal_theta(r);
        const std::vector<double> nabla_theta = cal_nabla_theta(r);
        const std::vector<double> nabla_phi = cal_nabla_phi(r);
        const std::complex<double> eip = std::exp(matrixop::IMAGIZ * cal_phi(r));

        std::vector<std::vector<std::complex<double>>> nablaH(r.size());
        // dH/dx
        std::vector<std::complex<double>>& nablaHx = nablaH.at(0);
        nablaHx.resize(m_dim * m_dim);
        nablaHx.at(0+0*m_dim) = m_params.at("A") * std::sin(theta) * nabla_theta.at(0);
        nablaHx.at(1+1*m_dim) = -nablaHx.at(0+0*m_dim);
        nablaHx.at(0+1*m_dim) = m_params.at("A") * eip * (std::cos(theta) * nabla_theta.at(0) + matrixop::IMAGIZ * std::sin(theta) * nabla_phi.at(0));
        nablaHx.at(1+0*m_dim) = nablaHx.at(0+1*m_dim);
        // dH/dy
        std::vector<std::complex<double>>& nablaHy = nablaH.at(1);
        nablaHy.resize(m_dim * m_dim);
        nablaHy.at(0+0*m_dim) = m_params.at("A") * std::sin(theta) * nabla_theta.at(1);
        nablaHy.at(1+1*m_dim) = -nablaHy.at(0+0*m_dim);
        nablaHy.at(0+1*m_dim) = m_params.at("A") * eip * (std::cos(theta) * nabla_theta.at(1) + matrixop::IMAGIZ * std::sin(theta) * nabla_phi.at(1));
        nablaHy.at(1+0*m_dim) = std::conj(nablaHy.at(0+1*m_dim));

        return nablaH;
    }

};
