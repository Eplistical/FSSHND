#ifndef _TULLY1_HAMILTONIAN_HPP
#define _TULLY1_HAMILTONIAN_HPP

/* 
 * Tully1 Hamiltonian 
 */

#include <vector>
#include <complex>
#include "hamiltonian.hpp"

namespace mqc {

    class Tully1_Hamiltonian final : public Hamiltonian {
        public:
            // --- ctor/dtor --- //
            Tully1_Hamiltonian();
            ~Tully1_Hamiltonian() = default;
        private:
            // --- quantities --- //
            double cal_phi(const std::vector<double>& /* r */) const;
            std::vector<double> cal_nabla_phi(const std::vector<double>& /* r */) const;
            virtual std::vector<std::complex<double>> cal_H(const std::vector<double>& /* r */) const override;
            virtual std::vector<std::vector<std::complex<double>>> cal_nablaH(const std::vector<double>& /* r */) const override;
    };

} // namespace mqc

#endif // _TULLY1_HAMILTONIAN_HPP