#ifndef _TULLY1_2D_HAMILTONIAN_HPP
#define _TULLY1_2D_HAMILTONIAN_HPP

/* 
 * 2D Tully1 Hamiltonian 
 */

#include <vector>
#include <complex>
#include "hamiltonian.hpp"

namespace mqc {

    class Tully1_2D_Hamiltonian final : public Hamiltonian {
        public:
            // --- ctor/dtor/operator= --- //
            Tully1_2D_Hamiltonian();
            Tully1_2D_Hamiltonian(const Tully1_2D_Hamiltonian& /* other */) = default;
            Tully1_2D_Hamiltonian& operator=(const Tully1_2D_Hamiltonian& /* other */) = default;
            ~Tully1_2D_Hamiltonian() = default;
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

#endif // _TULLY1_2D_HAMILTONIAN_HPP
