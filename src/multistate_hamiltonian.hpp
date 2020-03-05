#ifndef _MULTISTATE_HAMILTONIAN_HPP
#define _MULTISTATE_HAMILTONIAN_HPP

/* 
 * Multistate Hamiltonian
 */

#include <vector>
#include <complex>
#include "hamiltonian.hpp"

namespace mqc {

    class Multistate_Hamiltonian final : public Hamiltonian {
        public:
            // --- ctor/dtor/operator= --- //
            explicit Multistate_Hamiltonian(int dim);
            Multistate_Hamiltonian(const Multistate_Hamiltonian& /* other */) = default;
            Multistate_Hamiltonian& operator=(const Multistate_Hamiltonian& /* other */) = default;
            ~Multistate_Hamiltonian() = default;
        private:
            // --- utils --- //
            double cal_phi(const std::vector<double>& /* r */, int /* k */) const;
            std::vector<double> cal_nabla_phi(const std::vector<double>& /* r */, int /* k */) const;
        public:
            // --- interfaces --- //
            virtual std::vector<std::complex<double>> cal_H(const std::vector<double>& /* r */) const override;
            virtual std::vector<std::vector<std::complex<double>>> cal_nablaH(const std::vector<double>& /* r */) const override;
    };

} // namespace mqc

#endif // _MULTISTATE_HAMILTONIAN_HPP
