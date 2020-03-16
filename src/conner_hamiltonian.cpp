#include <iostream>
#include <utility>
#include <numeric>
#include <vector>
#include <cmath>
#include <complex>
#include "misc/matrixop.hpp"
#include "misc/vector.hpp"
#include "conner_hamiltonian.hpp"
#include "misc/crasher.hpp"

namespace mqc {


    // --- ctor/dtor/operator= --- //


    Conner_Hamiltonian::Conner_Hamiltonian(int dim) 
    : Hamiltonian(dim) { 
        m_doc = 
        "# Conner Hamiltoinan \n"
        "# trap(x, lx, rx) = tanh(x-rx) - tanh(x-lx)\n"
        "# H_{kk} = D + \\sum_{j<k} trap(x_j, R-B, R+B) + \\sum_{j>k} trap(x_j, -R-B, -R+B)\n"
        "#              + trap(x_k, -R-A, R)   if k == 0\n"
        "#              + trap(x_k, -R, R+A)   if k == N-1\n"
        "#              + trap(x_k, -R, R)   else \n"
        "# H_{k,k+1} = C * exp(i * W * sqrt(x_k^2 + x_{k+1}^2))\n"
        "# parameters: "
        "# { A, B, C, D, W } \n"
        ;

        m_params["A"] = 10.0;
        m_params["B"] = 5.0;
        m_params["C"] = 0.15;
        m_params["D"] = 3.0;
        m_params["R"] = 10.0;
        m_params["W"] = 0.0;
    }


    // --- utils --- //


    double Conner_Hamiltonian::cal_phi(const std::vector<double>& r, int k) const {
        //return m_params.at("W") * (r.at(k) + r.at(k+1));
        const double R = std::sqrt(std::pow(r.at(k), 2) + std::pow(r.at(k+1), 2));
        return m_params.at("W") * R;
    }

    std::vector<double> Conner_Hamiltonian::cal_nabla_phi(const std::vector<double>& r, int k) const {
        /*
        std::vector<double> nabla_phi(m_dim, 0.0);
        nabla_phi.at(k) = m_params.at("W");
        nabla_phi.at(k+1) = m_params.at("W");
        */
        const double R = std::sqrt(std::pow(r.at(k), 2) + std::pow(r.at(k+1), 2));
        std::vector<double> nabla_phi(m_dim, 0.0);
        nabla_phi.at(k) = m_params.at("W") * r.at(k) / R;
        nabla_phi.at(k+1) = m_params.at("W") * r.at(k+1) / R;
        return nabla_phi;
    }

    double Conner_Hamiltonian::trap(double x, const std::pair<double, double>& lxrx) const {
        /**
         * trap function with tanh boundaries [lx, rx].
         */
        return std::tanh(x-lxrx.second) - std::tanh(x-lxrx.first);
    }

    double Conner_Hamiltonian::der_trap(double x, const std::pair<double, double>& lxrx) const {
        /**
         * dertivate of the trap function to x
         */
        return -std::pow(std::tanh(x-lxrx.second), 2) + std::pow(std::tanh(x-lxrx.first), 2);
    }
    
    std::pair<double, double> Conner_Hamiltonian::get_boundaries(int k, int j) const {
        /**
         * get boundaries (lx, rx) for the j^th term in H_{kk}
         */
        const double A = m_params.at("A");
        const double B = m_params.at("B");
        const double R = m_params.at("R");
        if (j == k) {
            if (k == 0) {
                return std::make_pair(-R - A, R);
            }
            else if (k == m_dim - 1) {
                return std::make_pair(-R, R + A);
            }
            else {
                return std::make_pair(-R, R);
            }
        }
        else if (j < k) {
            return std::make_pair(R - B, R + B);
        }
        else {
            return std::make_pair(-R - B, -R + B);
        }
    }


    // --- interfaces --- //


    std::vector<std::complex<double>> Conner_Hamiltonian::cal_H(const std::vector<double>& r) const {
        /**
         * input : r
         * output : H(r)
         */
        // check
        misc::confirm<misc::ValueError>(m_dim == r.size(), "cal_H: the size of r must equal to the dim of Hamiltonian.");
        // diagonal terms
        std::vector<double> aux_smallj(m_dim, 0.0);
        std::vector<double> aux_largej(m_dim, 0.0);
        std::vector<double> aux_equalj(m_dim, 0.0);
        for (int j(0); j < m_dim; ++j) {
            aux_smallj.at(j) = trap(r.at(j), get_boundaries(m_dim, j));
            aux_largej.at(j) = trap(r.at(j), get_boundaries(-1, j));
            aux_equalj.at(j) = trap(r.at(j), get_boundaries(j, j));
        }
        std::vector<std::complex<double>> H(m_dim * m_dim, matrixop::ZEROZ);
        for (int k(0); k < m_dim; ++k) {
            H.at(k+k*m_dim) = m_params.at("D") + aux_equalj.at(k);
            for (int j(0); j < m_dim; ++j) {
                if (j < k) {
                    H.at(k+k*m_dim) += aux_smallj.at(j);
                }
                else if (j > k) {
                    H.at(k+k*m_dim) += aux_largej.at(j);
                }
            }
        }
        // off-diagonal terms
        for (int k(0); k < m_dim - 1; ++k) {
            const double phi = cal_phi(r, k);
            H.at(k+(k+1)*m_dim) = m_params.at("C") * std::exp(matrixop::IMAGIZ * phi);
            H.at((k+1)+k*m_dim) = std::conj(H.at(k+(k+1)*m_dim));
        }
        return H;
    }

    std::vector<std::vector<std::complex<double>>> Conner_Hamiltonian::cal_nablaH(const std::vector<double>& r) const {
        /**
         * input : r
         * output : nablaH(r)
         */
        // check
        misc::confirm<misc::ValueError>(m_dim == r.size(), "cal_H: the size of r must equal to the dim of Hamiltonian.");
        // diagonal terms
        std::vector<double> aux_der_smallj(m_dim, 0.0);
        std::vector<double> aux_der_largej(m_dim, 0.0);
        std::vector<double> aux_der_equalj(m_dim, 0.0);
        for (int j(0); j < m_dim; ++j) {
            aux_der_smallj.at(j) = der_trap(r.at(j), get_boundaries(m_dim, j));
            aux_der_largej.at(j) = der_trap(r.at(j), get_boundaries(-1, j));
            aux_der_equalj.at(j) = der_trap(r.at(j), get_boundaries(j, j));
        }

        std::vector<std::vector<std::complex<double>>> nabla_H(m_dim, std::vector<std::complex<double>>(m_dim * m_dim, matrixop::ZEROZ));
        for (int k(0); k < m_dim; ++k) {
            for (int j(0); j < m_dim; ++j) {
                if (j < k) {
                    nabla_H.at(j).at(k+k*m_dim) = aux_der_smallj.at(j);
                }
                else if (j > k) {
                    nabla_H.at(j).at(k+k*m_dim) = aux_der_largej.at(j);
                }
                else {
                    nabla_H.at(j).at(k+k*m_dim) = aux_der_equalj.at(j);
                }
            }
        }
        // off-diagonal terms
        std::vector<double> phi;
        std::vector<std::vector<double>> nabla_phi;
        for (int k(0); k < m_dim - 1; ++k) {
            phi.push_back(cal_phi(r, k));
            nabla_phi.push_back(cal_nabla_phi(r, k)); // nabla_phi.at(k).at(j) = d phi_k / d x_j
        }
        for (int k(0); k < m_dim - 1; ++k) {
            for (int j(0); j < m_dim; ++j) {
                nabla_H.at(j).at(k+(k+1)*m_dim) = m_params.at("C") * exp(matrixop::IMAGIZ * phi.at(k)) * matrixop::IMAGIZ * nabla_phi.at(k).at(j);
                nabla_H.at(j).at((k+1)+k*m_dim) = std::conj(nabla_H.at(j).at(k+(k+1)*m_dim));
            }
        }
        return nabla_H;
    }

};
