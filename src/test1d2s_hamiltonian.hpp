#ifndef _TEST1D2S_HAMILTONIAN_HPP
#define _TEST1D2S_HAMILTONIAN_HPP

/* 
 * Test1d2s Hamiltonian
 */

#include <vector>
#include <utility>
#include <complex>
#include "hamiltonian.hpp"

namespace mqc {

    class Test1d2s_Hamiltonian final : public Hamiltonian {
        public:
            // --- ctor/dtor/operator= --- //
            explicit Test1d2s_Hamiltonian();
            Test1d2s_Hamiltonian(const Test1d2s_Hamiltonian& /* other */) = default;
            Test1d2s_Hamiltonian& operator=(const Test1d2s_Hamiltonian& /* other */) = default;
            ~Test1d2s_Hamiltonian() = default;
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

#endif // _TEST1D2S_HAMILTONIAN_HPP
