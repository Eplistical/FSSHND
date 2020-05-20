#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include "misc/matrixop.hpp"
#include "misc/vector.hpp"
#include "ci_hamiltonian.hpp"
#include "misc/crasher.hpp"


namespace mqc {


    // --- ctor/dtor/operator= --- //


    CI_Hamiltonian::CI_Hamiltonian() 
    : Hamiltonian(2) { 
        m_doc = 
        "# CI Hamiltoinan \n"
        "# H_{00} = -x * f \n"
        "# H_{11} = -H_{00} \n"
        "# H_{01} = y * exp(i*phi) * f\n"
        "# H_{10} = conj(H_{01}) \n"
        "# phi = Wx * x + Wy * y \n"
        "# f = A*exp(-rnorm^2) + B/rnorm * (exp(C*rnorm)-1.0) / (exp(C*rnorm)+1.0) \n"
        "# rnorm = sqrt(x*x + y*y) \n"
        "# paramters: "
        "# { A, B, C, Wx, Wy } \n"
        ;
        m_params["A"] = 0.0025;
        m_params["B"] = 0.01;
        m_params["C"] = 0.8;
        m_params["Wx"] = 0.0;
        m_params["Wy"] = 0.0;
    }


    // --- utils --- //


    double CI_Hamiltonian::cal_phi(const std::vector<double>& r) const {
        return m_params.at("Wx") * r.at(0) + m_params.at("Wy") * r.at(1);
    }

    std::vector<double> CI_Hamiltonian::cal_nabla_phi(const std::vector<double>& r) const {
        return std::vector<double> { m_params.at("Wx"), m_params.at("Wy") };
    }


    // --- interfaces --- //


    std::vector<std::complex<double>> CI_Hamiltonian::cal_H(const std::vector<double>& r) const {
        /**
         * input : r
         * output : H(r)
         */
        // check
        misc::confirm(r.size() == 2, "CI_Hamiltonian: nuclear dim must be 2.");
        // params
        const double x = r.at(0);
        const double y = r.at(1);
        const double A = m_params.at("A");
        const double B = m_params.at("B");
        const double C = m_params.at("C");
        const double rnorm = sqrt(x*x + y*y);
        const std::complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r));
        const double P = A * exp(-rnorm*rnorm);
        const double Q = B / rnorm * (exp(C*rnorm) - 1.0) / (exp(C*rnorm) + 1.0);
        const double f = P + Q;
        std::vector< std::complex<double> > H(m_dim * m_dim, 0.0);
        if (rnorm > 1e-8) {
            H.at(0+0*m_dim) = -x * f;
            H.at(1+1*m_dim) = -H.at(0+0*m_dim);
            H.at(0+1*m_dim) = y * eip * f;
            H.at(1+0*m_dim) = std::conj(H.at(0+1*m_dim));
        }
        return H;
    }

    std::vector<std::vector<std::complex<double>>> CI_Hamiltonian::cal_nablaH(const std::vector<double>& r) const {
        /**
         * input : r
         * output : nablaH(r)
         */
        // check
        misc::confirm(r.size() == 2, "CI_Hamiltonian: nuclear dim must be 2.");
        // params
        const double x = r.at(0);
        const double y = r.at(1);
        const double A = m_params.at("A");
        const double B = m_params.at("B");
        const double C = m_params.at("C");
        const double rnorm = sqrt(x*x + y*y);
        const std::complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r));
        const double P = A * exp(-rnorm*rnorm);
        const double Q = B / rnorm * (exp(C*rnorm) - 1.0) / (exp(C*rnorm) + 1.0);
        const double f = P + Q;

        const double dr_dx = x / rnorm;
        const double dr_dy = y / rnorm;
        const double dP_dr = -2.0 * A * rnorm * exp(-rnorm*rnorm);
        const double dQ_dr = B / rnorm / rnorm * (2.0 * C * rnorm * exp(C*rnorm) - exp(2.0*C*rnorm) + 1.0) / pow(exp(C*rnorm) + 1.0, 2);
        const double df_dr = dP_dr + dQ_dr;
        const std::vector<double> nabla_phi = cal_nabla_phi(r);
        const double dphi_dx = nabla_phi.at(0);
        const double dphi_dy = nabla_phi.at(1);

        std::vector<std::vector<std::complex<double>>> nablaH(r.size());
        // dH/dx
        std::vector<std::complex<double>>& nablaHx = nablaH.at(0);
        nablaHx.assign(m_dim * m_dim, 0.0);
        if (rnorm > 1e-8) {
            nablaHx.at(0+0*m_dim) = -(f + x * df_dr * dr_dx);
            nablaHx.at(1+1*m_dim) = -nablaHx.at(0+0*m_dim);
            nablaHx.at(0+1*m_dim) = y * eip * matrixop::IMAGIZ * dphi_dx * f + y * eip * df_dr * dr_dx;
            nablaHx.at(1+0*m_dim) = conj(nablaHx.at(0+1*m_dim));
        }
        // dH/dy
        std::vector<std::complex<double>>& nablaHy = nablaH.at(1);
        nablaHy.assign(m_dim * m_dim, 0.0);
        if (rnorm > 1e-8) {
            nablaHy.at(0+0*m_dim) = -x * df_dr * dr_dy;
            nablaHy.at(1+1*m_dim) = -nablaHy.at(0+0*m_dim);
            nablaHy.at(0+1*m_dim) = (eip + y * eip * matrixop::IMAGIZ * dphi_dy) * f + y * eip * df_dr * dr_dy;
            nablaHy.at(1+0*m_dim) = conj(nablaHy.at(0+1*m_dim));
        }
        return nablaH;
    }

};
