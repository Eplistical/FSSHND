#include <iostream>
#include <utility>
#include <numeric>
#include <vector>
#include <cmath>
#include <complex>
#include "misc/matrixop.hpp"
#include "misc/vector.hpp"
#include "three_state_1d_hamiltonian.hpp"
#include "misc/crasher.hpp"

namespace mqc {


    // --- ctor/dtor/operator= --- //


    Three_State_1d_Hamiltonian::Three_State_1d_Hamiltonian() 
    : Hamiltonian(3) { 
        m_doc = 
        "# Three_State_1d Hamiltoinan \n"
        "# H_00 = A * tanh(E * x) + A\n"
        "# H_11 = -A * tanh(E * x) + A\n"
        "# H_22 = -A * tanh(E * x) + B\n"

        "# H_01 = C * exp(-E^2 * x^2) * V01\n"
        "# H_02 = C * exp(-E^2 * x^2) * V02\n"
        "# H_12 = C * exp(-E^2 * x^2) * V12\n"

        "# H_10 = conj(H_01)\n"
        "# H_20 = conj(H_02)\n"
        "# H_21 = conj(H_12)\n"
        "# parameters: "
        "# { A, B, C, E} \n"
        ;
        m_params["A"] = 0.02;
        m_params["B"] = 0.03;
        m_params["C"] = 0.005;
        m_params["E"] = 1.0;

        // constants
        m_params["V01I"] = 0.1;
        m_params["V01R"] = 0.5;
        m_params["V02I"] = 0.5;
        m_params["V02R"] = 0.1;
        m_params["V12I"] = 0.1;
        m_params["V12R"] = 0.1;
    }


    // --- interfaces --- //


    std::vector<std::complex<double>> Three_State_1d_Hamiltonian::cal_H(const std::vector<double>& r) const {
        /**
         * input : r
         * output : H(r)
         */
        // check
        misc::confirm(r.size() == 1, "Three_State_1d_Hamiltonian: nuclear dim must be 1.");
        // tmps
        const double x = r.at(0);
        const double A = m_params.at("A");
        const double B = m_params.at("B");
        const double C = m_params.at("C");
        const double E = m_params.at("E");
        const std::complex<double> V01(m_params.at("V01R"), m_params.at("V01I"));
        const std::complex<double> V02(m_params.at("V02R"), m_params.at("V02I"));
        const std::complex<double> V12(m_params.at("V12R"), m_params.at("V12I"));
        // diagonal terms
        std::vector<std::complex<double>> H(m_dim * m_dim, matrixop::ZEROZ);
        H.at(0+0*m_dim) = A*std::tanh(E*x) + A;
        H.at(1+1*m_dim) = -A*std::tanh(E*x) + A;
        H.at(2+2*m_dim) = -A*std::tanh(E*x) + B;
        // off-diagonal terms
        H.at(0+1*m_dim) = C*exp(-std::pow(E*x, 2)) * V01;
        H.at(0+2*m_dim) = C*exp(-std::pow(E*x, 2)) * V02;
        H.at(1+2*m_dim) = C*exp(-std::pow(E*x, 2)) * V12;
        H.at(1+0*m_dim) = std::conj(H.at(0+1*m_dim));
        H.at(2+0*m_dim) = std::conj(H.at(0+2*m_dim));
        H.at(2+1*m_dim) = std::conj(H.at(1+2*m_dim));
        // return
        return H;
    }

    std::vector<std::vector<std::complex<double>>> Three_State_1d_Hamiltonian::cal_nablaH(const std::vector<double>& r) const {
        /**
         * input : r
         * output : nablaH(r)
         */
        // check
        misc::confirm(r.size() == 1, "Three_State_1d_Hamiltonian: nuclear dim must be 1.");
        // tmps
        const double x = r.at(0);
        const double A = m_params.at("A");
        const double B = m_params.at("B");
        const double C = m_params.at("C");
        const double E = m_params.at("E");
        const std::complex<double> V01(m_params.at("V01R"), m_params.at("V01I"));
        const std::complex<double> V02(m_params.at("V02R"), m_params.at("V02I"));
        const std::complex<double> V12(m_params.at("V12R"), m_params.at("V12I"));
        // nabla_H
        std::vector<std::vector<std::complex<double>>> nabla_H(m_dim, std::vector<std::complex<double>>(m_dim * m_dim, matrixop::ZEROZ));
        // dH/dx
        std::vector<std::complex<double>>& Hx = nabla_H.at(0);
        Hx.at(0+0*m_dim) = A * (1.0 - std::pow(tanh(E*x), 2)) * E;
        Hx.at(1+1*m_dim) = -A * (1.0 - std::pow(tanh(E*x), 2)) * E;
        Hx.at(2+2*m_dim) = -A * (1.0 - std::pow(tanh(E*x), 2)) * E;
        Hx.at(0+1*m_dim) = C * exp(-std::pow(E*x, 2)) * -2.0 * std::pow(E, 2) * x * V01;
        Hx.at(0+2*m_dim) = C * exp(-std::pow(E*x, 2)) * -2.0 * std::pow(E, 2) * x * V02;
        Hx.at(1+2*m_dim) = C * exp(-std::pow(E*x, 2)) * -2.0 * std::pow(E, 2) * x * V12;
        Hx.at(1+0*m_dim) = std::conj(Hx.at(0+1*m_dim));
        Hx.at(2+0*m_dim) = std::conj(Hx.at(0+2*m_dim));
        Hx.at(2+1*m_dim) = std::conj(Hx.at(1+2*m_dim));
        // return
        return nabla_H;
    }

};
