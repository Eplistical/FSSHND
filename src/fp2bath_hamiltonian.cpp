#include <iostream>
#include <vector>
#include <complex>
#include "misc/matrixop.hpp"
#include "misc/vector.hpp"
#include "fp2bath_hamiltonian.hpp"
#include "misc/crasher.hpp"

namespace mqc {

    // --- ctor/dtor/operator= --- //


    FP2Bath_Hamiltonian::FP2Bath_Hamiltonian(double mu_l, double mu_r, double gamma_l, double gamma_r)
    : Hamiltonian(2), m_mu_l(mu_l), m_mu_r(mu_r), m_gamma_l(gamma_l), m_gamma_r(gamma_r) { 
        m_doc = 
        "# FP 2 Bath Hamiltoinan \n"
        "# H_{00} = 1/2 * MASS * OMEGA_X^2 * x^2 + 1/2 * MASS * OMEGA_Y^2 * y^2\n"
        "# H_{11} = 1/2 * MASS * OMEGA_X^2 * (x - X0)^2 + 1/2 * MASS * OMEGA_Y^2 * (y - Y0)^2 + E0\n"
        "# H_{01} = 0.0\n"
        "# H_{01} = 0.0\n"
        "# HERE H_{00} AND H_{11} ARE DIABATICAL SURFACES\n"
        "# paramters: "
        "# { E0, MASS, OMEGA_X, OMEGA_Y, X0, Y0} \n"
        ;
        m_params["E0"] = 0.0;
        m_params["MASS"] = 333.33;
        m_params["OMEGA_X"] = 0.003;
        m_params["OMEGA_Y"] = 0.0;
        m_params["X0"] = -9.428;
        m_params["Y0"] = 0.0;
    }


    // --- interfaces --- //


    std::vector<std::complex<double>> FP2Bath_Hamiltonian::cal_H(const std::vector<double>& r) const {
        /**
         * input : r
         * output : H(r)
         */
        // check
        misc::confirm(r.size() == 2, "FP2Bath_Hamiltonian: nuclear dim must be 2.");
        // params
        const double x = r.at(0);
        const double y = r.at(1);
        const double mass = m_params.at("MASS");
        const double wx = m_params.at("OMEGA_X");
        const double wy = m_params.at("OMEGA_Y");
        const double x0 = m_params.at("X0");
        const double y0 = m_params.at("Y0");
        const double E0 = m_params.at("E0");
        // H
        std::vector<std::complex<double>> H(m_dim * m_dim, matrixop::ZEROZ);
        H.at(0+0*m_dim) = 0.5 * mass * std::pow(wx*x, 2) + 0.5 * mass * std::pow(wy*y, 2);
        H.at(1+1*m_dim) = 0.5 * mass * std::pow(wx*(x-x0), 2) + 0.5 * mass * std::pow(wy*(y-y0), 2) + E0;
        H.at(0+1*m_dim) = 0.0;
        H.at(1+0*m_dim) = std::conj(H.at(1+0*m_dim));
        return H;
    }

    std::vector<std::vector<std::complex<double>>> FP2Bath_Hamiltonian::cal_nablaH(const std::vector<double>& r) const {
        /**
         * input : r
         * output : nablaH(r)
         */
        // check
        misc::confirm(r.size() == 2, "FP2Bath_Hamiltonian: nuclear dim must be 2.");
        // params
        const double x = r.at(0);
        const double y = r.at(1);
        const double mass = m_params.at("MASS");
        const double wx = m_params.at("OMEGA_X");
        const double wy = m_params.at("OMEGA_Y");
        const double x0 = m_params.at("X0");
        const double y0 = m_params.at("Y0");
        const double E0 = m_params.at("E0");
        // nablaH
        std::vector<std::vector<std::complex<double>>> nablaH(2);
        // dH/dx
        std::vector<std::complex<double>>& nablaHx = nablaH.at(0);
        nablaHx.resize(m_dim * m_dim);
        nablaHx.at(0+0*2) = mass * wx * wx * x;
        nablaHx.at(1+1*2) = mass * wx * wx * (x-x0);
        nablaHx.at(0+1*2) = 0.0;
        nablaHx.at(1+0*2) = conj(nablaHx.at(0+1*2));
        // dH/dy
        std::vector<std::complex<double>>& nablaHy = nablaH.at(1);
        nablaHy.resize(m_dim * m_dim);
        nablaHy.at(0+0*2) = mass * wy * wy * y;
        nablaHy.at(1+1*2) = mass * wy * wy * (y-y0);
        nablaHy.at(0+1*2) = 0.0;
        nablaHy.at(1+0*2) = conj(nablaHy.at(0+1*2));

        return nablaH;
    }

};
