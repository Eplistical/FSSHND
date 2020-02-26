#ifndef _FSSH_TRAJECTORY_HPP
#define _FSSH_TRAJECTORY_HPP

/* 
 * FSSH_Trajectory class
 */

#include <vector>
#include <complex>

namespace {

    template <typename HamiltonianType>
    class FSSH_Trajectory {
        public:
            // --- ctor/dtor --- //
            FSSH_Trajectory(const HamiltonianType& /* hamiltonian */);
            ~FSSH_Trajectory() = default;
        public:
            // --- simulation --- //
            void init_state( const std::vector<double>& /* r */, const std::vector<double>& /* v */,
                const std::vector<std::complex<double>>& /* c */, int /* s */);
            void integrator(double /* dt */);
            void hopper();
            bool check_end() const;
        public:
            // --- quantities --- //
            double cal_KE() const;
            double cal_PE() const;
        public:
            // --- getter/setter --- //
            double get_mass() const noexcept { return m_mass; }
            void set_mass(double mass) { m_mass = mass; }

            double get_kT() const noexcept { return m_kT; }
            void set_kT(double kT) { m_kT = kT; }

            std::vector<double> get_r() const noexcept { return m_r; }
            void set_r(const std::vector<double>& r) { m_r = r; }

            std::vector<double> get_v() const noexcept { return m_v; }
            void set_v(const std::vector<double>& v) { m_v = v; }

            std::vector<std::complex<double>> get_c() const noexcept { return m_c; }
            void set_c(const std::vector<std::complex<double>>& c) { m_c = c; }

            int get_s() const noexcept { return m_s; }
            void set_s(int s) { m_s = s; }

            double get_gamma() const noexcept { return m_gamma; }
            void set_gamma(double param_gamma) { m_gamma = param_gamma; }

            double get_enable_fric() const noexcept { return m_enable_fric; }
            void set_enable_fric(double enable_fric) { m_enable_fric = enable_fric; }
        private:
            // --- status --- //
            void update_status();
            std::vector<double> cal_berry_force();
        private:
            // --- basic --- //
            int m_ndim;
            int m_edim;
            double m_mass;
            double m_kT;
            std::vector<double> m_r;
            std::vector<double> m_v;
            std::vector<std::complex<double>> m_c;
            int m_s;
            // --- friction --- //
            double m_gamma;
            std::vector<double> m_randomforce;
            // --- Hamiltonian related --- //
            HamiltonianType m_hamiltonian;
            std::vector<std::vector<std::complex<double>>> m_force;
            std::vector<std::vector<std::complex<double>>> m_dc;
            std::vector<double> m_eva;
            std::vector<std::complex<double>> m_evt;
            // --- options/flags --- //
            bool m_enable_fric;
            bool m_initialized;
    };

};

#endif
