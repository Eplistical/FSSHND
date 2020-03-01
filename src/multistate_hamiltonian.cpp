#include <iostream>
#include <numeric>
#include <vector>
#include <cmath>
#include <complex>
#include "misc/matrixop.hpp"
#include "misc/vector.hpp"
#include "multistate_hamiltonian.hpp"
#include "misc/crasher.hpp"

namespace mqc {


    // --- ctor/dtor --- //


    Multistate_Hamiltonian::Multistate_Hamiltonian(int dim) 
    : Hamiltonian(dim) { 
        m_doc = 
        "# Multistate Hamiltoinan \n"
        "# V_k = 1/2 * MASS * OMEGA_T^2 * (x_k - X0)^2 + \\sum_{j \\neq k} 1/2 * MASS * OMEGA^2 * x_j^2 \n"
        "# H = \\sum{k} p_k^2/2/MASS + \\sum_{k} V_k * |d_k><d_k| + \\sum_{k} (VC * exp(i * W * (x_k + x_{k+1})) |d_k><d_{k+1}| + c.c.) \n"
        "# parameters: "
        "# { MASS, OMEGA, OMEGA_T, VC, W, X0 } \n"
        ;
        m_params["MASS"] = 1.0;
        m_params["OMEGA"] = 1.0;
        m_params["OMEGA_T"] = 2.0;
        m_params["VC"] = 0.2;
        m_params["W"] = 0.0;
        m_params["X0"] = 2.0;
    }


    // --- quantities --- //

    double Multistate_Hamiltonian::cal_phi(const std::vector<double>& r, int k) const {
        return m_params.at("W") * (r[k] + r[k+1]);
    }

    std::vector<double> Multistate_Hamiltonian::cal_nabla_phi(const std::vector<double>& r, int k) const {
        std::vector<double> nabla_phi(m_dim, 0.0);
        nabla_phi[k] = m_params.at("W");
        nabla_phi[k+1] = m_params.at("W");
    }

    std::vector<std::complex<double>> Multistate_Hamiltonian::cal_H(const std::vector<double>& r) const {
        /*
         * input : r
         * output : H(r)
         */
        // check
        misc::confirm<misc::ValueError>(m_dim == r.size(), "cal_H: the size of r must equal to the dim of Hamiltonian.");
        // diagonal terms
        std::vector<std::complex<double>> H(m_dim * m_dim, matrixop::ZEROZ);
        const double r_norm = std::accumulate(r.begin(), r.end(), 0.0,  [](const double accu, const double x) {
            return accu + std::pow(x, 2);
        });
        for (int k(0); k < m_dim; ++k) {
            H[k+k*m_dim] = 0.5 * m_params.at("MASS") * (
                std::pow(m_params.at("OMEGA"), 2) * (r_norm - std::pow(r[k], 2)) 
                + std::pow(m_params.at("OMEGA_T"), 2) * pow(r[k] - m_params.at("X0"), 2)
                );
        }
        // off-diagonal terms
        for (int k(0); k < m_dim - 1; ++k) {
            const double phi = cal_phi(r, k);
            H[k+(k+1)*m_dim] = m_params.at("VC") * std::exp(matrixop::IMAGIZ * phi);
            H[(k+1)+k*m_dim] = std::conj(H[k+(k+1)*m_dim]);
        }
        return H;
    }

    std::vector<std::vector<std::complex<double>>> Multistate_Hamiltonian::cal_nablaH(const std::vector<double>& r) const {
        /*
         * input : r
         * output : nablaH(r)
         */
        // check
        misc::confirm<misc::ValueError>(m_dim == r.size(), "cal_H: the size of r must equal to the dim of Hamiltonian.");
        // diagonal terms
        std::vector<std::vector<std::complex<double>>> nabla_H(m_dim, std::vector<std::complex<double>>(m_dim * m_dim, matrixop::ZEROZ));
        for (int k(0); k < m_dim; ++k) {
            for (int j(0); j < m_dim; ++j) {
                if (j == k) {
                    nabla_H.at(j).at(k+k*m_dim) = m_params.at("MASS") * std::pow(m_params.at("OMEGA_T"), 2) * (r.at(j) - m_params.at("X0"));
                }
                else {
                    nabla_H.at(j).at(k+k*m_dim) = m_params.at("MASS") * std::pow(m_params.at("OMEGA"), 2) * r.at(j);
                }
            }
        }
        // off-diagonal terms
        std::vector<double> phi;
        std::vector<std::vector<double>> nabla_phi;
        for (int k(0); k < m_dim - 1; ++k) {
            phi.push_back(cal_phi(r, k));
            nabla_phi.push_back(cal_nabla_phi(r, k));
        }
        for (int k(0); k < m_dim - 1; ++k) {
            for (int j(0); j < m_dim; ++j) {
                nabla_H.at(j).at(k+(k+1)*m_dim) = m_params.at("VC") * exp(matrixop::IMAGIZ * phi.at(k)) * matrixop::IMAGIZ * nabla_phi.at(k).at(j);
                nabla_H.at(j).at((k+1)+k*m_dim) = std::conj(nabla_H.at(j).at(k+(k+1)*m_dim));
            }
        }
        return nabla_H;
    }

};
