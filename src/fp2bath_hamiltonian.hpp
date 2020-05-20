#ifndef _FP2BATH_HAMILTONIAN_HPP
#define _FP2BATH_HAMILTONIAN_HPP

/* 
 * FP 2Bath Hamiltonian 
 */

#include <vector>
#include <complex>
#include "hamiltonian.hpp"

namespace mqc {

    class FP2Bath_Hamiltonian final : public Hamiltonian {
        public:
            // --- ctor/dtor/operator= --- //
            explicit FP2Bath_Hamiltonian(double /* mu_l */, double /* mu_r */, double /* gamma_l */, double /* gamma_r */);
            FP2Bath_Hamiltonian(const FP2Bath_Hamiltonian& /* other */) = default;
            FP2Bath_Hamiltonian& operator=(const FP2Bath_Hamiltonian& /* other */) = default;
            ~FP2Bath_Hamiltonian() = default;
        private:
            // --- utils --- //
        public:
            // --- interfaces --- //
            virtual std::vector<std::complex<double>> cal_H(const std::vector<double>& /* r */) const override;
            virtual std::vector<std::vector<std::complex<double>>> cal_nablaH(const std::vector<double>& /* r */) const override;
            double cal_gamma_l(const std::vector<double>& /* r */) const { return m_gamma_l; }
            double cal_gamma_r(const std::vector<double>& /* r */) const { return m_gamma_r; }
            double cal_mu_l() const { return m_mu_l; }
            double cal_mu_r() const { return m_mu_r; }

        private:
            double m_gamma_l, m_gamma_r;
            double m_mu_l, m_mu_r;

    };

} // namespace mqc

#endif // _FP2BATH_HAMILTONIAN_HPP
