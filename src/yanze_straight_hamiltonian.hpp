#ifndef _YANZE_STRAIGHT_HAMILTONIAN_HPP
#define _YANZE_STRAIGHT_HAMILTONIAN_HPP

/* 
 * Yanze straight Hamiltonian 
 */

#include <vector>
#include <complex>
#include "hamiltonian.hpp"

namespace mqc {

    class Yanze_Straight_Hamiltonian final : public Hamiltonian {
        public:
            // --- ctor/dtor/operator= --- //
            Yanze_Straight_Hamiltonian();
            Yanze_Straight_Hamiltonian(const Yanze_Straight_Hamiltonian& /* other */) = default;
            Yanze_Straight_Hamiltonian& operator=(const Yanze_Straight_Hamiltonian& /* other */) = default;
            ~Yanze_Straight_Hamiltonian() = default;
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

#endif // _YANZE_STRAIGHT_HAMILTONIAN_HPP
