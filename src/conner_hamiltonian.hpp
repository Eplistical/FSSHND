#ifndef _CONNER_HAMILTONIAN_HPP
#define _CONNER_HAMILTONIAN_HPP

/* 
 * Conner Hamiltonian
 */

#include <vector>
#include <utility>
#include <complex>
#include "hamiltonian.hpp"

namespace mqc {

    class Conner_Hamiltonian final : public Hamiltonian {
        public:
            // --- ctor/dtor/operator= --- //
            explicit Conner_Hamiltonian(int dim);
            Conner_Hamiltonian(const Conner_Hamiltonian& /* other */) = default;
            Conner_Hamiltonian& operator=(const Conner_Hamiltonian& /* other */) = default;
            ~Conner_Hamiltonian() = default;
        private:
            // --- utils --- //
            double cal_phi(const std::vector<double>& /* r */, int /* k */) const;
            std::vector<double> cal_nabla_phi(const std::vector<double>& /* r */, int /* k */) const;
            double trap(double /* x */, const std::pair<double, double>& /* lxrx */) const;
            double der_trap(double /* x */, const std::pair<double, double>& /* lxrx */) const;
            std::pair<double, double> get_boundaries(int /* k */, int /* j */) const;
        public:
            // --- interfaces --- //
            virtual std::vector<std::complex<double>> cal_H(const std::vector<double>& /* r */) const override;
            virtual std::vector<std::vector<std::complex<double>>> cal_nablaH(const std::vector<double>& /* r */) const override;
    };

} // namespace mqc

#endif // _CONNER_HAMILTONIAN_HPP
