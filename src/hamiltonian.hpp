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
            // --- ctor/dtor/opetator= --- //
            explicit Hamiltonian(int /* dim */);
            Hamiltonian(const Hamiltonian& /* other */) = default;
            Hamiltonian& operator=(const Hamiltonian& /* other */) = default;
            virtual ~Hamiltonian() = default;
        public:
            // --- utils interfaces --- //
            virtual void set_params(const std::vector<double>& /* params */);
            virtual void output_params() const noexcept;
        public:
            // --- interfaces --- //
            virtual std::vector<std::complex<double>> cal_H(const std::vector<double>& /* r */) const = 0;
            virtual std::vector<std::vector<std::complex<double>>> cal_nablaH(const std::vector<double>& /* r */) const = 0;
            virtual void cal_info( const std::vector<double>& /* r */,
                std::vector<std::vector<std::complex<double>>>& /* force */,
                std::vector<std::vector<std::complex<double>>>& /* dc */,
                std::vector<double>& /* eva */,
                std::vector<std::complex<double>>& /* lastevt */
                ) const;
        public:
            inline int get_dim() const noexcept { return m_dim; }

            inline std::string get_doc() const noexcept { return m_doc; }
            inline void set_doc(const std::string& doc) { m_doc = doc; }

            inline std::map<std::string, double> get_params() const noexcept { return m_params; }
            inline void set_params(const std::map<std::string, double>& params) { m_params = params; }

            inline double get_param(const std::string& key) const { return m_params.at(key); }
            inline void set_param(const std::string& key, double val) { m_params[key] = val; }

        protected:
            int m_dim;
            std::string m_doc;
            std::map<std::string, double> m_params;
    };

} // namespace mqc

#endif // _HAMILTONIAN_HPP