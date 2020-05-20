#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include "misc/matrixop.hpp"
#include "misc/vector.hpp"
#include "yanze_straight_hamiltonian.hpp"
#include "misc/crasher.hpp"


namespace mqc {


    // --- ctor/dtor/operator= --- //


    Yanze_Straight_Hamiltonian::Yanze_Straight_Hamiltonian() 
    : Hamiltonian(2) { 
        m_doc = 
        "# Yanze Straight Hamiltoinan \n"
        "# H_{00} = A * (1 + tanh(E*y))\n"
        "# H_{11} = A * (1 - tanh(E*y))\n"
        "# H_{01} = C * exp(-E^2 * y^2) * exp(i * W * phi)\n"
        "# H_{10} = conj(H_{01}) \n"
        "# phi = x"
        "# paramters: "
        "# { A, C, E, W } \n"
        ;
        m_params["A"] = 0.02;
        m_params["C"] = 0.005;
        m_params["E"] = 1.0;
        m_params["W"] = 0.2;
    }


    // --- utils --- //

    double Yanze_Straight_Hamiltonian::cal_phi(const std::vector<double>& r) const {
        return m_params.at("W") * r.at(0);
    }

    std::vector<double> Yanze_Straight_Hamiltonian::cal_nabla_phi(const std::vector<double>& r) const {
        return std::vector<double> { m_params.at("W"), 0.0 };
    }


    // --- interfaces --- //


    std::vector<std::complex<double>> Yanze_Straight_Hamiltonian::cal_H(const std::vector<double>& r) const {
        /**
         * input : r
         * output : H(r)
         */
        // check
        misc::confirm(r.size() == 2, "Yanze_Straight_Hamiltonian: nuclear dim must be 2.");
        // params
        const double A = m_params.at("A");
        const double C = m_params.at("C");
        const double E = m_params.at("E");
        const double x = r.at(0);
        const double y = r.at(1);
        const std::complex<double> eip = std::exp(matrixop::IMAGIZ * cal_phi(r));
        std::vector<std::complex<double>> H(m_dim * m_dim, matrixop::ZEROZ);
        H.at(0+0*m_dim) = A * (1.0 + std::tanh(E*y));
        H.at(1+1*m_dim) = A * (1.0 - std::tanh(E*y));
        H.at(0+1*m_dim) = C * std::exp(-std::pow(E*y, 2)) * eip;
        H.at(1+0*m_dim) = std::conj(H.at(0+1*m_dim));
        return H;
    }

    std::vector<std::vector<std::complex<double>>> Yanze_Straight_Hamiltonian::cal_nablaH(const std::vector<double>& r) const {
        /**
         * input : r
         * output : nablaH(r)
         */
        // check
        misc::confirm(r.size() == 2, "Yanze_Straight_Hamiltonian: nuclear dim must be 2.");
        // params
        const double A = m_params.at("A");
        const double C = m_params.at("C");
        const double E = m_params.at("E");
        const double x = r.at(0);
        const double y = r.at(1);
        const std::vector<double> nabla_phi = cal_nabla_phi(r);
        const std::complex<double> eip = std::exp(matrixop::IMAGIZ * cal_phi(r));

        std::vector<std::vector<std::complex<double>>> nablaH(r.size());
        // dH/dx
        std::vector<std::complex<double>>& nablaHx = nablaH.at(0);
        nablaHx.resize(m_dim * m_dim);
        nablaHx.at(0+0*m_dim) = 0.0;
        nablaHx.at(1+1*m_dim) = 0.0;
        nablaHx.at(0+1*m_dim) = C * std::exp(-std::pow(E*y, 2)) * eip * matrixop::IMAGIZ * nabla_phi.at(0);
        nablaHx.at(1+0*m_dim) = nablaHx.at(0+1*m_dim);
        // dH/dy
        std::vector<std::complex<double>>& nablaHy = nablaH.at(1);
        nablaHy.resize(m_dim * m_dim);
        nablaHy.at(0+0*m_dim) = A * (1.0 - std::pow(std::tanh(E*y), 2)) * E;
        nablaHy.at(1+1*m_dim) = -A * (1.0 - std::pow(std::tanh(E*y), 2)) * E;
        nablaHy.at(0+1*m_dim) = C * std::exp(-std::pow(E*y, 2)) * eip * (-2.0 * std::pow(E, 2) * y + matrixop::IMAGIZ * nabla_phi.at(1));
        nablaHy.at(1+0*m_dim) = std::conj(nablaHy.at(0+1*m_dim));

        return nablaH;
    }

};
