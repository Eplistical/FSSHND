#ifndef _TULLY1_HAMILTONIAN_HPP
#define _TULLY1_HAMILTONIAN_HPP

/* 
 * Classical 1D Tully1 Hamiltonian 
 */

#include <vector>
#include <complex>
#include "hamiltonian.hpp"

namespace mqc {

    class Tully1_Hamiltonian final : public Hamiltonian {
        public:
            // --- ctor/dtor/operator= --- //
            Tully1_Hamiltonian();
            Tully1_Hamiltonian(const Tully1_Hamiltonian& /* other */) = default;
            Tully1_Hamiltonian& operator=(const Tully1_Hamiltonian& /* other */) = default;
            ~Tully1_Hamiltonian() = default;
        public:
            // --- interfaces --- //
            virtual std::vector<std::complex<double>> cal_H(const std::vector<double>& /* r */) const override;
            virtual std::vector<std::vector<std::complex<double>>> cal_nablaH(const std::vector<double>& /* r */) const override;
    };

} // namespace mqc

#endif // _TULLY1_HAMILTONIAN_HPP
