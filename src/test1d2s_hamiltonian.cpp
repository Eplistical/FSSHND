#include <iostream>
#include <utility>
#include <numeric>
#include <vector>
#include <cmath>
#include <complex>
#include "misc/matrixop.hpp"
#include "misc/vector.hpp"
#include "test1d2s_hamiltonian.hpp"
#include "misc/crasher.hpp"

namespace mqc {


    // --- ctor/dtor/operator= --- //


    Test1d2s_Hamiltonian::Test1d2s_Hamiltonian() 
    : Hamiltonian(2) { 
        m_doc = 
        "# Test1d2s Hamiltoinan \n"
        "# 1-dim 2-state test Hamiltonian"
        ;
        /*
        m_params["A"] = 0.02;
        m_params["B"] = 0.3;
        m_params["C"] = 0.005;
        m_params["E"] = 1.0;

        // constants
        m_params["V01I"] = 0.1;
        m_params["V01R"] = 0.5;
        */

        m_params["A"] = 0.01;
        m_params["B"] = 1.6;
        m_params["C"] = 0.005;
        m_params["D"] = 1.0;
    }


    // --- interfaces --- //


    std::vector<std::complex<double>> Test1d2s_Hamiltonian::cal_H(const std::vector<double>& r) const {
        /**
         * input : r
         * output : H(r)
         */
        // check
        misc::confirm(r.size() == 1, "Test1d2s_Hamiltonian: nuclear dim must be 1.");
        /*
        // tmps
        const double x = r.at(0);
        const double A = m_params.at("A");
        const double B = m_params.at("B");
        const double C = m_params.at("C");
        const double E = m_params.at("E");
        const std::complex<double> V01(m_params.at("V01R"), m_params.at("V01I"));
        // diagonal terms
        std::vector<std::complex<double>> H(m_dim * m_dim, matrixop::ZEROZ);
        H.at(0+0*m_dim) = A*std::tanh(E*x) + A;
        H.at(1+1*m_dim) = -A*std::tanh(E*x) + A;
        // off-diagonal terms
        H.at(0+1*m_dim) = C*exp(-std::pow(E*x, 2)) * V01;
        H.at(1+0*m_dim) = std::conj(H.at(0+1*m_dim));
        // return
        return H;
        */

        // tmps
        const double x = r.at(0);
        const double A = m_params.at("A");
        const double B = m_params.at("B");
        const double C = m_params.at("C");
        const double D = m_params.at("D");
        std::vector<std::complex<double>> H(m_dim * m_dim, matrixop::ZEROZ);
        if (x > 0.0) {
            H[0+0*2] = A * (1.0 - exp(-B * x));
        }
        else {
            H[0+0*2] = -A * (1.0 - exp(B * x));
        }
        H[1+1*2] = -H[0+0*2];
        H[0+1*2] = C * exp(-D * x * x);
        H[1+0*2] = conj(H[0+1*2]);

        return H;
    }

    std::vector<std::vector<std::complex<double>>> Test1d2s_Hamiltonian::cal_nablaH(const std::vector<double>& r) const {
        /**
         * input : r
         * output : nablaH(r)
         */
        /*
        // check
        misc::confirm(r.size() == 1, "Test1d2s_Hamiltonian: nuclear dim must be 1.");
        // tmps
        const double x = r.at(0);
        const double A = m_params.at("A");
        const double B = m_params.at("B");
        const double C = m_params.at("C");
        const double E = m_params.at("E");
        const std::complex<double> V01(m_params.at("V01R"), m_params.at("V01I"));
        // nabla_H
        std::vector<std::vector<std::complex<double>>> nabla_H(m_dim, std::vector<std::complex<double>>(m_dim * m_dim, matrixop::ZEROZ));
        // dH/dx
        std::vector<std::complex<double>>& Hx = nabla_H.at(0);
        Hx.at(0+0*m_dim) = A * (1.0 - std::pow(tanh(E*x), 2)) * E;
        Hx.at(1+1*m_dim) = -A * (1.0 - std::pow(tanh(E*x), 2)) * E;
        Hx.at(0+1*m_dim) = C * exp(-std::pow(E*x, 2)) * -2.0 * std::pow(E, 2) * x * V01;
        Hx.at(1+0*m_dim) = std::conj(Hx.at(0+1*m_dim));
        // return
        return nabla_H;
        */

        const double x = r.at(0);
        const double A = m_params.at("A");
        const double B = m_params.at("B");
        const double C = m_params.at("C");
        const double D = m_params.at("D");
        std::vector<std::vector<std::complex<double>>> nabla_H(m_dim, std::vector<std::complex<double>>(m_dim * m_dim, matrixop::ZEROZ));
        std::vector<std::complex<double>>& Hx = nabla_H.at(0);
        Hx.resize(4);
        if (x > 0.0) {
            Hx[0+0*2] = A * B * exp(-B * x);
        }
        else {
            Hx[0+0*2] = A * B * exp(B * x);
        }
        Hx[1+1*2] = -Hx[0+0*2];
        Hx[0+1*2] = -2 * x * C * D * exp(-D * x * x);
        Hx[1+0*2] = conj(Hx[0+1*2]);
        return nabla_H;
    }
};
