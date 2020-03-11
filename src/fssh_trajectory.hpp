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
#include "misc/fmtstring.hpp"


namespace mqc {


    // --- DECLARATION --- // 


    template <typename HamiltonianType>
    class FSSH_Trajectory {
        public:
            // --- ctor/dtor --- //
            FSSH_Trajectory(const HamiltonianType& /* hamiltonian */);
            FSSH_Trajectory(const FSSH_Trajectory<HamiltonianType>& /* other */) = default;
            FSSH_Trajectory<HamiltonianType>& operator=(const FSSH_Trajectory<HamiltonianType>& /* other */) = default;
            // --- simualtion utils --- //
            ~FSSH_Trajectory() = default;
        private:
            // --- simualtion utils --- //
            void cal_Tmat();
            void update_status(double /* dt */);
            std::vector<double> cal_berry_force();
            void nuclear_integrator(double /* dt */);
            void electronic_integrator(double /* dt */);
            double cal_hop_prob(int /* from */, int /* to */, double /* dt */) const;
            void hopper(double /* dt */);
        public:
            // --- simulation interfaces --- //
            void setup( double /* mass */, const std::vector<double>& /* r */, const std::vector<double>& /* v */,
                const std::vector<std::complex<double>>& /* c */, int /* s */);
            void die();
            void integrator(double /* dt */);
        public:
            // --- other interfaces --- //
            double cal_KE() const;
            double cal_PE() const;
            double cal_totE() const;
            std::vector<double> cal_diab_pop() const;
        public:
            // --- getter/setter --- //
            inline int get_ndim() const noexcept { return m_ndim; }
            inline void set_ndim(int ndim) { m_ndim = ndim; }

            inline int get_edim() const noexcept { return m_ndim; }
            inline void set_edim(int ndim) { m_ndim = ndim; }

            inline double get_mass() const noexcept { return m_mass; }
            inline void set_mass(double mass) { m_mass = mass; }

            inline double get_kT() const noexcept { return m_kT; }
            inline void set_kT(double kT) { m_kT = kT; }

            inline double get_t() const noexcept { return m_t; }
            inline void set_t(double t) { m_t = t; }

            inline std::vector<double> get_r() const noexcept { return m_r; }
            inline void set_r(const std::vector<double>& r) { m_r = r; }

            inline std::vector<double> get_v() const noexcept { return m_v; }
            inline void set_v(const std::vector<double>& v) { m_v = v; }

            inline std::vector<std::complex<double>> get_c() const noexcept { return m_c; }
            inline void set_c(const std::vector<std::complex<double>>& c) { m_c = c; }

            inline int get_s() const noexcept { return m_s; }
            inline void set_s(int s) { m_s = s; }

            inline int get_Nhop_accepted() const noexcept { return m_Nhop_accepted; }
            inline void set_Nhop_accepted(int Nhop_accepted) { m_Nhop_accepted = Nhop_accepted; }

            inline int get_Nhop_frustrated() const noexcept { return m_Nhop_frustrated; }
            inline void set_Nhop_frustrated(int Nhop_frustrated) { m_Nhop_frustrated = Nhop_frustrated; }

            inline double get_gamma() const noexcept { return m_gamma; }
            inline void set_gamma(double param_gamma) { m_gamma = param_gamma; }

            inline bool get_enable_hop() const noexcept { return m_enable_hop; }
            inline void set_enable_hop(bool enable_hop) { m_enable_hop = enable_hop; }

            inline bool get_enable_berry_force() const noexcept { return m_enable_berry_force; }
            inline void set_enable_berry_force(bool enable_berry_force) { m_enable_berry_force = enable_berry_force; }

            inline std::vector<double> get_eva() const noexcept { return m_eva; }
            inline std::vector<std::complex<double>> get_evt() const noexcept { return m_evt; }
            inline std::vector<std::vector<std::complex<double>>> get_force() const noexcept { return m_force; }
            inline std::vector<std::vector<std::complex<double>>> get_dc() const noexcept { return m_dc; }
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
            std::vector<std::complex<double>> m_Tmat;
            // --- hop --- //
            int m_Nhop_accepted;
            int m_Nhop_frustrated;
            // --- options/flags --- //
            bool m_initialized;
            bool m_enable_hop;
            bool m_enable_berry_force;
    };


    // --- IMPLEMENTATION --- // 


    // --- ctor/dtor --- //


    template <typename HamiltonianType>
        FSSH_Trajectory<HamiltonianType>::FSSH_Trajectory(const HamiltonianType& hamiltonian) 
        : m_kT(0.0), m_gamma(-1.0), m_hamiltonian(hamiltonian), m_initialized(false), m_enable_hop(true), m_enable_berry_force(true)
        {  
            m_edim = hamiltonian.get_dim();
        }


    // --- simulation --- //


    template <typename HamiltonianType>
        void FSSH_Trajectory<HamiltonianType>::cal_Tmat() {
            /**
             * calculate T_jk = v \dot dc_jk
             */
            m_Tmat.assign(m_edim * m_edim, matrixop::ZEROZ);
            for (int j(0); j < m_edim; ++j) {
                for (int k(0); k < j; ++k) {
                    for (int i(0); i < m_ndim; ++i) {
                        m_Tmat.at(j+k*m_edim) += m_v.at(i) * m_dc.at(i).at(j+k*m_edim);
                    }
                    m_Tmat.at(k+j*m_edim) = -std::conj(m_Tmat.at(j+k*m_edim));
                }
            }
        }

    template <typename HamiltonianType>
        void FSSH_Trajectory<HamiltonianType>::update_status(double dt) {
            /**
             * update status of the trajectory current configuration
             * if dt > 0, m_Tmat will be updated.
             */
            if (dt <= 0.0) {
                // dt <= 0.0 indicates this is the first step, calculate T = v \dot dc
                m_hamiltonian.cal_info(m_r, m_force, m_dc, m_eva, m_evt);
                cal_Tmat();
            }
            else {
                // dt > 0.0, calculate T = log[U(dt)] / dt, U_jk(dt) = <psi_j(t0) | psi_k(t0 + dt)>
                misc::confirm<misc::StateError> (not m_evt.empty(), "update_status: dt > 0.0 while m_evt is empty!");
                m_hamiltonian.cal_info_Tmat(m_r, m_force, m_dc, m_eva, m_evt, m_Tmat, dt);
            }
        }

    
    template <typename HamiltonianType>
        std::vector<double> FSSH_Trajectory<HamiltonianType>::cal_berry_force() {
            /** 
             * calculate current berry force
             */
            std::vector<double> F_berry(m_ndim, 0.0);
            if (m_enable_berry_force) {
                for (int i(0); i < m_ndim; ++i) {
                    for (int k(0); k < m_edim; ++k) {
                        if (k != m_s) {
                            F_berry.at(i) += (m_dc.at(i).at(m_s+k*m_edim) * m_Tmat.at(k+m_s*m_edim)).imag();
                        }
                    }
                }
                F_berry *= 2.0;
            }
            return F_berry;
        }

    template <typename HamiltonianType>
        void FSSH_Trajectory<HamiltonianType>::electronic_integrator(double dt) {
            /**
             * electronic integrator -- RK4
             */
            std::vector<std::complex<double>> rk4_mat = -m_Tmat;
            for (int j(0); j < m_edim; ++j) {
                rk4_mat.at(j+j*m_edim) -= matrixop::IMAGIZ * m_eva.at(j);
            }
            std::vector<std::complex<double>> k1, k2, k3, k4;
            k1 = dt * matrixop::matmat(rk4_mat, m_c, m_edim);
            k2 = dt * matrixop::matmat(rk4_mat, m_c + 0.5 * k1, m_edim);
            k3 = dt * matrixop::matmat(rk4_mat, m_c + 0.5 * k2, m_edim);
            k4 = dt * matrixop::matmat(rk4_mat, m_c + k3, m_edim);
            m_c += 1.0 / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        }

    template <typename HamiltonianType>
        void FSSH_Trajectory<HamiltonianType>::nuclear_integrator(double dt) {
            /**
             * nuclear integrator -- velocity verlet
             */
            std::vector<double> random_force;
            std::vector<double> tot_force;
            // 1st force update
            tot_force = cal_berry_force(); // berry force
            for (int i(0); i < m_ndim; ++i) {
                tot_force.at(i) += m_force.at(i).at(m_s+m_s*m_edim).real(); // adiabatic force
            }
            if (m_gamma > 0.0) {
                // friction & random force
                random_force = randomer::vnormal(m_ndim, 0.0, std::sqrt(2.0 * m_gamma * m_kT / dt));
                tot_force += -m_gamma * m_v + random_force;
            }
            // 1st conf update
            m_v += 0.5 * dt / m_mass * tot_force;
            m_r += dt * m_v;
            update_status(dt);
            // 2nd force update 
            tot_force = cal_berry_force(); // berry force
            for (int i(0); i < m_ndim; ++i) {
                tot_force.at(i) += m_force.at(i).at(m_s+m_s*m_edim).real(); // adiabatic force
            }
            if (m_gamma > 0.0) {
                // friction & random force
                tot_force += -m_gamma * m_v + random_force;
            }
            // 2nd conf update
            m_v += 0.5 * dt / m_mass * tot_force;
        }
            
    template <typename HamiltonianType>
        double FSSH_Trajectory<HamiltonianType>::cal_hop_prob(int j, int k, double dt) const {
            /**
             * calculate hopping probability for j->k
             */
            const std::complex<double> rho_jk = m_c.at(j) * conj(m_c.at(k));
            const std::complex<double> rho_jj = m_c.at(j) * conj(m_c.at(j));
            const double rst = -2.0 * dt * (rho_jk * m_Tmat.at(k+j*m_edim)).real() / rho_jj.real();
            return rst > 0.0 ? rst : 0.0;
        }

    template <typename HamiltonianType>
        void FSSH_Trajectory<HamiltonianType>::hopper(double dt) {
            /**
             * hopping function 
             */
            misc::confirm<misc::StateError>(m_initialized, "hopper: FSSH_Trajectory died / not initialzied.");
            if (m_edim <= 1) {
                // adiabatic dynamics, no hop
                return ;
            }
            // calculate hop probablities
            const int from = m_s;
            std::vector<double> prob(m_edim, 0.0);
            prob.at(from) = 1.0;
            for (int to(0); to < m_edim; ++to) {
                if (to != from) {
                    prob.at(to) = cal_hop_prob(from, to, dt);
                    prob.at(from) -= prob.at(to);
                }
            }
            misc::confirm<misc::ValueError>(prob.at(from) > 0.0, "hopper: hop prob too large, try a smaller dt.");
            // determine whether hop happens 
            const double randnum = randomer::rand(0.0, 1.0);
            double accu = 0.0;
            int to = 0;
            while (to < m_edim) {
                accu += prob.at(to);
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
                for (int i(0); i < m_ndim; ++i) {
                    n_rescale.at(i) = (m_dc.at(i).at(from+to*m_edim) * m_Tmat.at(to+from*m_edim)).real();
                }
                misc::confirm<misc::RuntimeError>(norm(n_rescale) > 1e-40, "hopper: n_rescale has zero norm.");
                // implement momemtum rescaling
                const double dE = m_eva.at(to) - m_eva.at(from);
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
                    // MOMENTUM REVERSAL ?
                }
            }
        }


    // --- simulation interfaces --- //


    template <typename HamiltonianType>
        void FSSH_Trajectory<HamiltonianType>::setup( double mass, const std::vector<double>& r, const std::vector<double>& v,
                const std::vector<std::complex<double>>& c, int s) { 
            /**
             * setup trajectory for simulation
             */
            // check
            misc::confirm<misc::ValueError>(r.size() == v.size(), "init_state: invalid r/v sizes.");
            misc::confirm<misc::ValueError>(c.size() == m_edim, "init_state: invalid c size.");
            misc::confirm<misc::ValueError>(s < m_edim, "init_state: invalid s number.");
            // setup
            m_ndim = r.size();
            set_mass(mass);
            set_r(r);
            set_v(v);
            set_c(c);
            set_s(s);
            set_t(0.0);
            set_Nhop_accepted(0);
            set_Nhop_frustrated(0);
            // update status
            update_status(0.0);
            m_initialized = true;
        }

    template <typename HamiltonianType>
        void FSSH_Trajectory<HamiltonianType>::die() {
            /**
             * kill the trajectory, cannot evolve any more
             */
            std::vector<double>().swap(m_randomforce);
            std::vector<std::vector<std::complex<double>>>().swap(m_force);
            std::vector<std::vector<std::complex<double>>>().swap(m_dc);
            m_initialized = false;
        }

    template <typename HamiltonianType>
        void FSSH_Trajectory<HamiltonianType>::integrator(double dt) {
            /**
             * INTERFACE METHOD
             * propagate trajectory forward by dt
             */
            misc::confirm<misc::StateError>(m_initialized, "integrator: FSSH_Trajectory died / not initialzied.");

            // --- nuclear part --- //

            auto lasteva = m_eva;
            nuclear_integrator(dt);

            // --- electronic part --- //

            // determine dtq and Ndq
            double max_abs_T = 0.0;
            for (auto& Tx : m_Tmat) {
                max_abs_T = std::max(std::abs(Tx), max_abs_T);
            } 
            double max_abs_dE = 0.0;
            auto diffE = m_eva - mean(m_eva);
            for (auto& diffEx : diffE) {
                max_abs_dE = std::max(std::abs(diffEx), max_abs_dE);
            }
            double dtq = std::min(dt, std::min(0.02 / max_abs_T, 0.02 / max_abs_dE));
            int Ndtq = std::round(dt / dtq);
            dtq = dt / Ndtq;

            // propagate: assuming const Tmat, and m_eva can be linearly interpolated
            auto deva = (m_eva - lasteva) / Ndtq;
            m_eva.swap(lasteva);
            for (int idtq(0); idtq < Ndtq; ++idtq) {
                electronic_integrator(dtq);
                misc::confirm<misc::ValueError>(norm(m_c) < 2.0, ("integrator: norm(c) diverges."));
                if (m_enable_hop) {
                    hopper(dtq);
                } 
                m_eva += deva;
            }
            m_eva.swap(lasteva);
            /*
            for (int idtq(0); idtq < Ndtq; ++idtq) {
                // propagate dtq forward
                double cur_dtq = dtq;
                double accu_dtq = 0.0;
                while (accu_dtq < dtq) {
                    if (dtq - accu_dtq < cur_dtq) {
                        cur_dtq = dtq - accu_dtq;
                    }
                    try {
                        hopper(cur_dtq);
                    } catch (const misc::ValueError& e) {
                        // catch exception: need to reduce cur_dtq
                        cur_dtq *= 0.5;
                        continue;
                    }
                    electronic_integrator(cur_dtq);
                    accu_dtq += cur_dtq;
                }
            */

            // --- time part --- //
            m_t += dt;
        }


    // --- other interfaces --- //


    template <typename HamiltonianType>
        double FSSH_Trajectory<HamiltonianType>::cal_KE() const {
            /**
             * calculate kinetic energy
             */
            return 0.5 * m_mass * sum(m_v * m_v);
        }

    template <typename HamiltonianType>
        double FSSH_Trajectory<HamiltonianType>::cal_PE() const {
            /**
             * calculate potential energy 
             */
            misc::confirm<misc::StateError>(not m_eva.empty(), "cal_PE: eigenvalues not yet calculated.");
            return m_eva.at(m_s);
        }

    template <typename HamiltonianType>
        double FSSH_Trajectory<HamiltonianType>::cal_totE() const {
            /**
             * calculate total energy 
             */
            return cal_KE() + cal_PE();
        }

    template <typename HamiltonianType>
        std::vector<double> FSSH_Trajectory<HamiltonianType>::cal_diab_pop() const {
            /**
             * calculate diab population using Landry Paper Method #3 (DOI: 10.1063/1.4837795)
             */
            std::vector<double> rst(m_edim);
            for (int a(0); a < m_edim; ++a) {
                rst.at(a) = std::pow(std::abs(m_evt.at(a+m_s*m_edim)), 2);
                for (int j(0); j < m_edim; ++j) {
                    for (int i(0); i < j; ++i) {
                        rst.at(a) += 2.0 * (m_evt.at(a+i*m_edim) * m_c.at(i) * std::conj(m_c.at(j)) * std::conj(m_evt.at(a+j*m_edim))).real();
                    }
                }
            }
            return rst;
        }


} // namespace mqc

#endif // _FSSH_TRAJECTORY_HPP
