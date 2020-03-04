#ifndef _NUSCATTER_HAMILTONIAN_HPP
#define _NUSCATTER_HAMILTONIAN_HPP

/* 
 * Yanze's nuclear scattering Hamiltonian
 */

#include <vector>
#include <complex>
#include "hamiltonian.hpp"

namespace mqc {

    class Nuscatter_Hamiltonian final : public Hamiltonian {
        public:
            // --- ctor/dtor --- //
            explicit Nuscatter_Hamiltonian();
            ~Nuscatter_Hamiltonian() = default;
        public:
            // --- quantities --- //
            double cal_m_r(const std::vector<double>& /* r */) const;
            std::vector<double> cal_der_m_r(const std::vector<double>& /* r */) const;
            double cal_m_theta(const std::vector<double>& /* r */) const;
            std::vector<double> cal_der_m_theta(const std::vector<double>& /* r */) const;
            double cal_phi(const std::vector<double>& /* r */) const;
            std::vector<double> cal_der_phi(const std::vector<double>& /* r */) const;
            double cal_param_coup() const;
        public:
            // --- interfaces --- //
            virtual std::vector<std::complex<double>> cal_H(const std::vector<double>& /* r */) const override;
            virtual std::vector<std::vector<std::complex<double>>> cal_nablaH(const std::vector<double>& /* r */) const override;
    };

} // namespace mqc

#endif // _NUSCATTER_HAMILTONIAN_HPP
