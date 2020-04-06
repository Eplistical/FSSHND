#ifndef _THREE_STATE_HAMILTONIAN_HPP
#define _THREE_STATE_HAMILTONIAN_HPP

/* 
 * Three_State Hamiltonian
 */

#include <vector>
#include <utility>
#include <complex>
#include "hamiltonian.hpp"

namespace mqc {

    class Three_State_Hamiltonian final : public Hamiltonian {
        public:
            // --- ctor/dtor/operator= --- //
            explicit Three_State_Hamiltonian();
            Three_State_Hamiltonian(const Three_State_Hamiltonian& /* other */) = default;
            Three_State_Hamiltonian& operator=(const Three_State_Hamiltonian& /* other */) = default;
            ~Three_State_Hamiltonian() = default;
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

#endif // _THREE_STATE_HAMILTONIAN_HPP
