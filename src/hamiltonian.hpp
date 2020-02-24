#ifndef _HAMILTONIAN_HPP
#define _HAMILTONIAN_HPP

/* 
 * Hamiltonian base class
 */

#include <vector>
#include <complex>

namespace {

    class Hamiltonian {
        public:
            // --- ctor/dtor --- //
            Hamiltonian() = default;
            virtual ~Hamiltonian() = default;
        public:
            // --- parameters --- //
            virtual void set_params(const std::vector<double>& /* params */);
            virtual void output_params() const noexcept;
        public:
            // --- quantities --- //
            virtual std::vector<std::complex<double>> cal_H(const std::vector<double>& /* r */) const = 0;
            virtual std::vector<std::vector<std::complex<double>>> cal_nablaH(const std::vector<double>& /* r */) const = 0;
            virtual void cal_info( const std::vector<double>& /* r */,
                    std::vector<std::complex<double>>& /* force */,
                    std::vector<std::complex<double>>& /* dc */,
                    std::vector<double>& /* eva */,
                    std::vector<std::complex<double>>& /* lastevt */
                    ) const = 0;
    };

};

#endif
