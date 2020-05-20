#ifndef _CI_HAMILTONIAN_HPP
#define _CI_HAMILTONIAN_HPP

/* 
 * 2D Tully1 Hamiltonian 
 */

#include <vector>
#include <complex>
#include "hamiltonian.hpp"

namespace mqc {

    class CI_Hamiltonian final : public Hamiltonian {
        public:
            // --- ctor/dtor/operator= --- //
            CI_Hamiltonian();
            CI_Hamiltonian(const CI_Hamiltonian& /* other */) = default;
            CI_Hamiltonian& operator=(const CI_Hamiltonian& /* other */) = default;
            ~CI_Hamiltonian() = default;
        private:
            // --- utils --- //
            double cal_phi(const std::vector<double>& /* r */) const;
            std::vector<double> cal_nabla_phi(const std::vector<double>& /* r */) const;
        public:
            // --- interfaces --- //
            virtual std::vector<std::complex<double>> cal_H(const std::vector<double>& /* r */) const override;
            virtual std::vector<std::vector<std::complex<double>>> cal_nablaH(const std::vector<double>& /* r */) const override;
    };

} // namespace mqc

#endif // _CI_HAMILTONIAN_HPP
