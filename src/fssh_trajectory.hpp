#ifndef _FSSH_TRAJECTORY_HPP
#define _FSSH_TRAJECTORY_HPP

/* 
 * FSSH_Trajectory class
 */

#include <vector>
#include <complex>
#include "misc/crasher.hpp"
#include "misc/matrixop.hpp"
#include "misc/randomer.hpp"
#include "misc/vector.hpp"


namespace mqc {


    // --- DECLARATION --- // 


    template <typename HamiltonianType>
    class FSSH_Trajectory {
        public:
            // --- ctor/dtor --- //
            FSSH_Trajectory(const HamiltonianType& /* hamiltonian */);
            ~FSSH_Trajectory() = default;
        private:
            // --- simualtion utils --- //
            double cal_hop_prob(int /* from */, int /* to */, double /* dt */) const;
        public:
            // --- simulation --- //
            void setup( double /* mass */, const std::vector<double>& /* r */, const std::vector<double>& /* v */,
                const std::vector<std::complex<double>>& /* c */, int /* s */);
            void die();
            void integrator(double /* dt */);
            void hopper(double /* dt */);
        public:
            // --- quantities --- //
            double cal_KE() const;
            double cal_PE() const;
        public:
            // --- getter/setter --- //
            int get_ndim() const noexcept { return m_ndim; }
            void set_ndim(int ndim) { m_ndim = ndim; }

            int get_edim() const noexcept { return m_ndim; }
            void set_edim(int ndim) { m_ndim = ndim; }

            double get_mass() const noexcept { return m_mass; }
            void set_mass(double mass) { m_mass = mass; }

            double get_kT() const noexcept { return m_kT; }
            void set_kT(double kT) { m_kT = kT; }

            double get_t() const noexcept { return m_t; }
            void set_t(double t) { m_t = t; }

            std::vector<double> get_r() const noexcept { return m_r; }
            void set_r(const std::vector<double>& r) { m_r = r; }

            std::vector<double> get_v() const noexcept { return m_v; }
            void set_v(const std::vector<double>& v) { m_v = v; }

            std::vector<std::complex<double>> get_c() const noexcept { return m_c; }
            void set_c(const std::vector<std::complex<double>>& c) { m_c = c; }

            int get_s() const noexcept { return m_s; }
            void set_s(int s) { m_s = s; }

            int get_Nhop_accepted() const noexcept { return m_Nhop_accepted; }
            void set_Nhop_accepted(int Nhop_accepted) { m_Nhop_accepted = Nhop_accepted; }

            int get_Nhop_frustrated() const noexcept { return m_Nhop_frustrated; }
            void set_Nhop_frustrated(int Nhop_frustrated) { m_Nhop_frustrated = Nhop_frustrated; }

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
            double m_t;
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
            // --- hop --- //
            int m_Nhop_accepted;
            int m_Nhop_frustrated;
            // --- options/flags --- //
            bool m_enable_fric;
            bool m_initialized;
    };


    // --- IMPLEMENTATION --- // 


    // --- ctor/dtor --- //


    template <typename HamiltonianType>
        FSSH_Trajectory<HamiltonianType>::FSSH_Trajectory(const HamiltonianType& hamiltonian) 
        : m_hamiltonian(hamiltonian), m_enable_fric(false), m_initialized(false)
        {  
            m_edim = hamiltonian.get_dim();
        }


    // --- simulation --- //


    template <typename HamiltonianType>
        void FSSH_Trajectory<HamiltonianType>::setup( double mass, const std::vector<double>& r, const std::vector<double>& v,
                const std::vector<std::complex<double>>& c, int s) { 
            /*
             * setup trajectory for simulation
             */
            // check
            misc::confirm<misc::ValueError>(r.size() == v.size(), "init_state: invalid r/v sizes.");
            misc::confirm<misc::ValueError>(c.size() == m_edim, "init_state: invalid c size.");
            misc::confirm<misc::ValueError>(s < m_edim, "init_state: invalid s number.");
            // setup
            m_ndim = r.size();
            set_mass(mass);
            set_t(0.0);
            set_r(r);
            set_v(v);
            set_c(c);
            set_s(s);
            set_Nhop_accepted(0);
            set_Nhop_frustrated(0);
            // update status
            update_status();
            m_initialized = true;
        }

    template <typename HamiltonianType>
        void FSSH_Trajectory<HamiltonianType>::die() {
            /*
             * kill the trajectory and save only simple configuration info
             */
            std::vector<double>().swap(m_randomforce);
            std::vector<std::complex<double>>().swap(m_evt);
            std::vector<std::vector<std::complex<double>>>().swap(m_force);
            std::vector<std::vector<std::complex<double>>>().swap(m_dc);
            m_initialized = false;
        }
            
    template <typename HamiltonianType>
        void FSSH_Trajectory<HamiltonianType>::integrator(double dt) {
            // propagate trajectory forward by dt
            misc::confirm<misc::ValueError>(m_initialized, "FSSH_Trajectory died / not initialzied.");

            // electronic part -- RK4
            std::vector<std::complex<double>> rk4_mat(m_edim * m_edim, 0.0);
            for (int j(0); j < m_edim; ++j) {
                for (int k(0); k < m_edim; ++k) {
                    for (int i(0); i < m_ndim; ++i) {
                        rk4_mat[j+k*m_edim] -= m_v[i] * m_dc[i][j+k*m_edim];
                    }
                }
            }
            for (int j(0); j < m_edim; ++j) {
                rk4_mat[j+j*m_edim] -= matrixop::IMAGIZ * m_eva[j];
            }
            std::vector<std::complex<double>> k1, k2, k3, k4;
            k1 = dt * matrixop::matmat(rk4_mat, m_c, m_edim);
            k2 = dt * matrixop::matmat(rk4_mat, m_c + 0.5 * k1, m_edim);
            k3 = dt * matrixop::matmat(rk4_mat, m_c + 0.5 * k2, m_edim);
            k4 = dt * matrixop::matmat(rk4_mat, m_c + k3, m_edim);
            m_c += 1.0 / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
            
            // TODO : ADD FRICTION & RANDOM FORCE
            // nuclear part -- velocity verlet
            std::vector<double> tot_force;
            // 1st force update
            tot_force = cal_berry_force(); // berry force
            for (int i(0); i < m_ndim; ++i) {
                tot_force[i] += m_force[i][m_s+m_s*m_edim].real(); // adiabatic force
            }
            // 1st conf update
            m_v += 0.5 * dt / m_mass * tot_force;
            m_r += dt * m_v;
            update_status();
            // 2nd force update 
            tot_force = cal_berry_force(); // berry force
            for (int i(0); i < m_ndim; ++i) {
                tot_force[i] += m_force[i][m_s+m_s*m_edim].real(); // adiabatic force
            }
            // 2nd conf update
            m_v += 0.5 * dt / m_mass * tot_force;
            m_t += dt;
            update_status();
        }

    template <typename HamiltonianType>
        double FSSH_Trajectory<HamiltonianType>::cal_hop_prob(int from, int to, double dt) const {
            /*
             * hopping probability for from->to
             */
            std::complex<double> v_dot_dc = matrixop::ZEROZ;
            for (int i(0); i < m_ndim; ++i) {
                v_dot_dc += m_v[i] * m_dc[i][to+from*m_edim];
            }
            return -2.0 * dt * (m_c[from] * conj(m_c[to]) * v_dot_dc).real() / (m_c[from] * conj(m_c[from])).real();
        }

    template <typename HamiltonianType>
        void FSSH_Trajectory<HamiltonianType>::hopper(double dt) {
            misc::confirm<misc::ValueError>(m_initialized, "FSSH_Trajectory died / not initialzied.");
            if (m_edim <= 1) {
                // one surface only, no hop
                return ;
            }
            // calculate hop probablities
            const int from = m_s;
            std::vector<double> prob(m_edim, 0.0);
            prob[from] = 1.0;
            for (int to(0); to < m_edim; ++to) {
                if (to != from) {
                    prob[to] = cal_hop_prob(from, to, dt);
                    prob[from] -= prob[to];
                }
            }
            misc::confirm<misc::ValueError>(prob[from] > 0.0, "hopper: tot prob > 1.0, dt too large.");
            // determine whether hop happens 
            const double randnum = randomer::rand(0.0, 1.0);
            double accu = 0.0;
            int to = 0;
            while (to < m_edim) {
                accu += prob[to];
                if (accu > randnum) {
                    break;
                }
                to += 1;
            }
            misc::confirm<misc::IndexError>(to < m_edim, "hopper: to >= edim.");
            // hop occurs
            if (to != from) {
                // momentum rescaling direction
                std::vector<double> n_rescale(m_ndim, 0.0);
                
                // x-direction rescaling
                n_rescale[0] = 1.0;

                /*
                // Re(dc_jk * (v \dot dc_kj)) rescaling
                std::complex<double> v_dot_dc = matrixop::ZEROZ;
                for (int i(0); i < m_ndim; ++i) {
                    v_dot_dc += m_v[i] * m_dc[i][to+from*m_edim];
                }
                for (int i(0); i < m_ndim; ++i) {
                    n_rescale[i] = (m_dc[i][from+to*m_edim] * v_dot_dc).real();
                }
                */

                // implement momemtum rescaling
                const double dE = m_eva[to] - m_eva[from];
                if (norm(n_rescale) > 1e-40) {
                    std::vector<double> vn = component(m_v, n_rescale);
                    double vn_norm = norm(vn);
                    double tmp = vn_norm * vn_norm - 2.0 * dE / m_mass; 
                    if (tmp > 0.0) {
                        // hop accepted
                        m_Nhop_accepted += 1;
                        double vn_norm_new = sqrt(tmp);
                        m_v += (vn_norm_new - vn_norm) / vn_norm * vn;
                        m_s = to;
                    }
                    else {
                        // hop frustrated
                        m_Nhop_frustrated += 1;
                    }
                }
            }
        }

    // --- quantities --- //


    template <typename HamiltonianType>
        double FSSH_Trajectory<HamiltonianType>::cal_KE() const {
            return 0.5 * m_mass * sum(m_v * m_v);
        }

    template <typename HamiltonianType>
        double FSSH_Trajectory<HamiltonianType>::cal_PE() const {
            misc::confirm<misc::ValueError>(not m_eva.empty(), "cal_PE: eigenvalues not yet calculated.");
            return m_eva.at(m_s);
        }


    // --- status --- //


    template <typename HamiltonianType>
        void FSSH_Trajectory<HamiltonianType>::update_status() {
            /*
             * update m_force, m_dc, m_eva, m_evt according to current m_r
             */
            m_hamiltonian.cal_info(m_r, m_force, m_dc, m_eva, m_evt);
        }
    
    template <typename HamiltonianType>
        std::vector<double> FSSH_Trajectory<HamiltonianType>::cal_berry_force() {
            /* 
             * calculate current berry force
             */
            std::vector<double> F_berry(m_ndim, 0.0);
            std::complex<double> v_dot_dc(0.0, 0.0);
            for (int k(0); k < m_edim; ++k) {
                if (k != m_s) {
                    v_dot_dc = 0.0;
                    for (int i(0); i < m_ndim; ++i) {
                        v_dot_dc += m_v[i] * m_dc[i][k+m_s*m_edim];
                    }
                    for (int i(0); i < m_ndim; ++i) {
                        F_berry[i] += (m_dc[i][m_s+k*m_edim] * v_dot_dc).imag();
                    }
                }
            }
            F_berry *= 2.0;
            return F_berry;
        }
} // namespace mqc

#endif // _FSSH_TRAJECTORY_HPP
