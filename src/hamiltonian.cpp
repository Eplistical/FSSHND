#include <vector>
#include <iostream>
#include <complex>
#include <iterator>
#include <string>
#include "misc/vector.hpp"
#include "misc/matrixop.hpp"
#include "misc/crasher.hpp"
#include "adjust_evt_phase.hpp"
#include "hamiltonian.hpp"

namespace mqc {


    // --- ctor/dtor/operator= --- //


    Hamiltonian::Hamiltonian(int dim) 
    : m_dim(dim)
    { }


    // --- utils interfacs --- //


    void Hamiltonian::set_params(const std::vector<double>& params) {
        /**
         * set params with a vector
         *  params are ordered alphabetically
         */
        misc::confirm<misc::ValueError>(params.size() >= m_params.size(), "set_params: insufficient params size.");
        int i = 0;
        auto it = m_params.begin();
        while (it != m_params.end()) {
            it->second = params.at(i);
            std::advance(it, 1);
            i += 1;
        }
    }

    void Hamiltonian::output_params() const noexcept {
        /**
         * output parameters
         */
        std::cout << std::string(32, '#') << std::endl;
        std::cout << m_doc;
        std::cout << std::string(32, '#') << std::endl;
        std::cout << "# ";
        for (const auto& it : m_params) {
            std::cout << it.first << " = " << it.second << " ";
        }
        std::cout << std::endl;
        std::cout << std::string(32, '#') << std::endl;
    }


    // --- interfaces --- //


    void Hamiltonian::cal_info( const std::vector<double>& r,
        std::vector<std::vector<std::complex<double>>>& force,
        std::vector<std::vector<std::complex<double>>>& dc,
        std::vector<double>& eva,
        std::vector<std::complex<double>>& lastevt) const {
            /**
             * input r and lastevt (evt from the last step) 
             * output force, dc, eva and evt (repalces lastevt on exit)
             */

            // calculate eva & evt w/ correct phase
            std::vector<std::complex<double>> evt;
            matrixop::hdiag(cal_H(r), eva, evt);
            if (not lastevt.empty()) {
                /*
                zeyu::adjust_evt_phase(lastevt, evt, m_dim, 0.0);
                */
                auto tmp = matrixop::matCmat(lastevt, evt, m_dim);
                for (int j = 0; j < m_dim; ++j) {
                    std::complex<double> eip = tmp.at(j+j*m_dim) / abs(tmp.at(j+j*m_dim));
                    for (int k = 0; k < m_dim; ++k) {
                        evt.at(k+j*m_dim) /= eip;
                    }
                }
            }
            // calculate force & dc
            const int ndim = r.size();
            const std::vector<std::vector<std::complex<double>>> nablaH = cal_nablaH(r);
            dc.resize(ndim);
            force.resize(ndim);
            for (int i = 0; i < ndim; ++i) {
                dc.at(i) = matrixop::matCmatmat(evt, nablaH.at(i), evt, m_dim, m_dim);
                force.at(i).resize(m_dim * m_dim);
                for (int j = 0; j < m_dim; ++j) {
                    for (int k = 0; k < m_dim; ++k) {
                        force.at(i).at(j+k*m_dim) = -dc.at(i).at(j+k*m_dim);
                        if (j == k) {
                            dc.at(i).at(j+k*m_dim) = matrixop::ZEROZ;
                        }
                        else {
                            dc.at(i).at(j+k*m_dim) /= (eva.at(k) - eva.at(j));
                        }
                    }
                }
            }
            // store evt
            lastevt = std::move(evt);
    }


    void Hamiltonian::cal_info_Tmat( const std::vector<double>& r,
        std::vector<std::vector<std::complex<double>>>& force,
        std::vector<std::vector<std::complex<double>>>& dc,
        std::vector<double>& eva,
        std::vector<std::complex<double>>& lastevt,
        std::vector<std::complex<double>>& Tmat,
        double dt
        ) const {
            /**
             * input r and lastevt (evt from the last step) 
             * output force, dc, eva and evt (repalces lastevt on exit)
             */

            // calculate eva & evt w/ correct phase
            std::vector<std::complex<double>> evt;
            matrixop::hdiag(cal_H(r), eva, evt);
            if (not lastevt.empty()) {
                /*
                Tmat = zeyu::adjust_evt_phase(lastevt, evt, m_dim, dt);
                */
                auto tmp = matrixop::matCmat(lastevt, evt, m_dim);
                for (int j = 0; j < m_dim; ++j) {
                    std::complex<double> eip = tmp.at(j+j*m_dim) / abs(tmp.at(j+j*m_dim));
                    for (int k = 0; k < m_dim; ++k) {
                        evt.at(k+j*m_dim) /= eip;
                    }
                }
                auto U = matrixop::matCmat(lastevt, evt, m_dim);
                Tmat = matrixop::logmh(U) / dt;
            }
            // calculate force & dc
            const int ndim = r.size();
            const std::vector<std::vector<std::complex<double>>> nablaH = cal_nablaH(r);
            dc.resize(ndim);
            force.resize(ndim);
            for (int i = 0; i < ndim; ++i) {
                dc.at(i) = matrixop::matCmatmat(evt, nablaH.at(i), evt, m_dim, m_dim);
                force.at(i).resize(m_dim * m_dim);
                for (int j = 0; j < m_dim; ++j) {
                    for (int k = 0; k < m_dim; ++k) {
                        force.at(i).at(j+k*m_dim) = -dc.at(i).at(j+k*m_dim);
                        if (j == k) {
                            dc.at(i).at(j+k*m_dim) = matrixop::ZEROZ;
                        }
                        else {
                            dc.at(i).at(j+k*m_dim) /= (eva.at(k) - eva.at(j));
                        }
                    }
                }
            }
            // store evt
            lastevt = std::move(evt);
    }
}
