#ifndef _HAMILTONIAN_HPP
#define _HAMILTONIAN_HPP

/* 
 * Hamiltonian base class
 */

#include <vector>
#include <complex>
#include <map>
#include <string>

namespace mqc {

    class Hamiltonian {
        public:
            // --- ctor/dtor --- //
            Hamiltonian(int dim);
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
                std::vector<std::vector<std::complex<double>>>& /* force */,
                std::vector<std::vector<std::complex<double>>>& /* dc */,
                std::vector<double>& /* eva */,
                std::vector<std::complex<double>>& /* lastevt */
                ) const;
        public:
            int get_dim() const noexcept { return m_dim; }

            std::string get_doc() const noexcept { return m_doc; }
            void set_doc(const std::string& doc) { m_doc = doc; }

            std::map<std::string, double> get_params() const noexcept { return m_params; }
            double get_param(const std::string& key) const { return m_params.at(key); }
            void set_params(const std::map<std::string, double>& params) { m_params = params; }

        protected:
            int m_dim;
            std::string m_doc;
            std::map<std::string, double> m_params;
    };

} // namespace mqc

#endif // _HAMILTONIAN_HPP