#include <iostream>
#include <utility>
#include <numeric>
#include <vector>
#include <cmath>
#include <complex>
#include "misc/matrixop.hpp"
#include "misc/vector.hpp"
#include "three_state_hamiltonian.hpp"
#include "misc/crasher.hpp"

namespace mqc {


    // --- ctor/dtor/operator= --- //


    Three_State_Hamiltonian::Three_State_Hamiltonian() 
    : Hamiltonian(3) { 
        m_doc = 
        "# Three_State Hamiltoinan \n"
        "# H_00 = 0.5 * MASS * OMEGA_A^2 * x^2 + A * tanh(E * y) + A\n"
        "# H_11 = 0.5 * MASS * OMEGA_A^2 * x^2 - A * tanh(E * y) + A\n"
        "# H_22 = 0.5 * MASS * OMEGA_B^2 * x^2 - A * tanh(E * y) + B\n"
        "# H_01 = C * exp(-E^2 * y^2) * V01\n"
        "# H_02 = C * exp(-E^2 * y^2) * V02\n"
        "# H_12 = C * exp(-E^2 * y^2) * V12\n"
        "# H_10 = conj(H_01)\n"
        "# H_20 = conj(H_02)\n"
        "# H_21 = conj(H_12)\n"
        "# parameters: "
        "# { A, B, C, E, MASS, OMEGA_A, OMEGA_B } \n"
        ;

        m_params["A"] = 0.02;
        m_params["B"] = 0.03;
        m_params["C"] = 0.005;
        m_params["E"] = 1.0;
        m_params["MASS"] = 1000.0;
        m_params["OMEGA_A"] = 0.02;
        m_params["OMEGA_B"] = 0.02;

        // constants
        m_params["V01I"] = 0.1;
        m_params["V01R"] = 0.5;
        m_params["V02I"] = 0.5;
        m_params["V02R"] = 0.1;
        m_params["V12I"] = 0.1;
        m_params["V12R"] = 0.1;
    }


    // --- interfaces --- //


    std::vector<std::complex<double>> Three_State_Hamiltonian::cal_H(const std::vector<double>& r) const {
        /**
         * input : r
         * output : H(r)
         */
        // check
        misc::confirm(r.size() == 2, "Three_State_Hamiltonian: nuclear dim must be 2.");
        // tmps
        const double x = r.at(0);
        const double y = r.at(1);
        const double A = m_params.at("A");
        const double B = m_params.at("B");
        const double C = m_params.at("C");
        const double E = m_params.at("E");
        const double M = m_params.at("MASS");
        const double Wa = m_params.at("OMEGA_A");
        const double Wb = m_params.at("OMEGA_B");
        const std::complex<double> V01(m_params.at("V01R"), m_params.at("V01I"));
        const std::complex<double> V02(m_params.at("V02R"), m_params.at("V02I"));
        const std::complex<double> V12(m_params.at("V12R"), m_params.at("V12I"));
        // diagonal terms
        std::vector<std::complex<double>> H(m_dim * m_dim, matrixop::ZEROZ);
        H.at(0+0*m_dim) = 0.5*M*std::pow(Wa*x, 2) + A*std::tanh(E*y) + A;
        H.at(1+1*m_dim) = 0.5*M*std::pow(Wa*x, 2) - A*std::tanh(E*y) + A;
        H.at(2+2*m_dim) = 0.5*M*std::pow(Wb*x, 2) - A*std::tanh(E*y) + B;
        // off-diagonal terms
        H.at(0+1*m_dim) = C*exp(-std::pow(E*y, 2)) * V01;
        H.at(0+2*m_dim) = C*exp(-std::pow(E*y, 2)) * V02;
        H.at(1+2*m_dim) = C*exp(-std::pow(E*y, 2)) * V12;
        H.at(1+0*m_dim) = std::conj(H.at(0+1*m_dim));
        H.at(2+0*m_dim) = std::conj(H.at(0+2*m_dim));
        H.at(2+1*m_dim) = std::conj(H.at(1+2*m_dim));
        // return
        return H;
    }

    std::vector<std::vector<std::complex<double>>> Three_State_Hamiltonian::cal_nablaH(const std::vector<double>& r) const {
        /**
         * input : r
         * output : nablaH(r)
         */
        // check
        misc::confirm(r.size() == 2, "Three_State_Hamiltonian: nuclear dim must be 2.");
        // tmps
        const double x = r.at(0);
        const double y = r.at(1);
        const double A = m_params.at("A");
        const double B = m_params.at("B");
        const double C = m_params.at("C");
        const double E = m_params.at("E");
        const double M = m_params.at("MASS");
        const double Wa = m_params.at("OMEGA_A");
        const double Wb = m_params.at("OMEGA_B");
        const std::complex<double> V01(m_params.at("V01R"), m_params.at("V01I"));
        const std::complex<double> V02(m_params.at("V02R"), m_params.at("V02I"));
        const std::complex<double> V12(m_params.at("V12R"), m_params.at("V12I"));
        // nabla_H
        std::vector<std::vector<std::complex<double>>> nabla_H(m_dim, std::vector<std::complex<double>>(m_dim * m_dim, matrixop::ZEROZ));
        // dH/dx
        std::vector<std::complex<double>>& Hx = nabla_H.at(0);
        Hx.at(0+0*m_dim) = M*Wa*Wa*x;
        Hx.at(1+1*m_dim) = M*Wa*Wa*x;
        Hx.at(2+2*m_dim) = M*Wb*Wb*x;
        // dH/dy
        std::vector<std::complex<double>>& Hy = nabla_H.at(1);
        Hy.at(0+0*m_dim) = A * (1.0 - std::pow(tanh(E*y), 2)) * E;
        Hy.at(1+1*m_dim) = -A * (1.0 - std::pow(tanh(E*y), 2)) * E;
        Hy.at(2+2*m_dim) = -A * (1.0 - std::pow(tanh(E*y), 2)) * E;
        Hy.at(0+1*m_dim) = C * exp(-std::pow(E*y, 2)) * -2.0 * std::pow(E, 2) * y * V01;
        Hy.at(0+2*m_dim) = C * exp(-std::pow(E*y, 2)) * -2.0 * std::pow(E, 2) * y * V02;
        Hy.at(1+2*m_dim) = C * exp(-std::pow(E*y, 2)) * -2.0 * std::pow(E, 2) * y * V12;
        Hy.at(1+0*m_dim) = std::conj(Hy.at(0+1*m_dim));
        Hy.at(2+0*m_dim) = std::conj(Hy.at(0+2*m_dim));
        Hy.at(2+1*m_dim) = std::conj(Hy.at(1+2*m_dim));
        // return
        return nabla_H;
    }

};
