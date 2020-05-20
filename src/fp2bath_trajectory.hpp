#ifndef _FP2BATH_TRAJECTORY_HPP
#define _FP2BATH_TRAJECTORY_HPP

/* 
 * FP2Bath_Trajectory class
 */

#include <vector>
#include <complex>
#include "misc/crasher.hpp"
#include "misc/matrixop.hpp"
#include "misc/randomer.hpp"
#include "misc/vector.hpp"
#include "misc/fmtstring.hpp"
#include "misc/fermi.hpp"


namespace mqc {


    // --- DECLARATION --- // 


    template <typename HamiltonianType>
        class FP2Bath_Trajectory {
            public:
                // --- ctor/dtor --- //
                FP2Bath_Trajectory(const HamiltonianType& /* hamiltonian */);
                FP2Bath_Trajectory(const FP2Bath_Trajectory<HamiltonianType>& /* other */) = default;
                FP2Bath_Trajectory<HamiltonianType>& operator=(const FP2Bath_Trajectory<HamiltonianType>& /* other */) = default;
                // --- simualtion utils --- //
                ~FP2Bath_Trajectory() = default;
            private:
                // --- simualtion utils --- //
                void update_status();
                void nuclear_integrator(double /* dt */);
                void cal_fric_noise_tensor(std::vector<double>& /* fric */, std::vector<double>& /* noise */) const;
            public:
                // --- simulation interfaces --- //
                void setup(double /* mass */, const std::vector<double>& /* r */, const std::vector<double>& /* v */, double /* kT */);
                void die();
                void integrator(double /* dt */);
            public:
                // --- other interfaces --- //
                double cal_KE() const { return 0.5 * m_mass * sum(m_v * m_v); }
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

                double get_nuclear_fric() const noexcept { return m_nuclear_fric; }
                void set_nuclear_fric(double nuclear_fric) { m_nuclear_fric = nuclear_fric; }

                std::vector<double> get_force() const noexcept { return m_force; }
            public:
                // --- NOT IMPLEMENTED THESE ARE DUMMY INTERFACES TO WORK WITH TRAJ_RECORDER --- //
                std::vector<std::complex<double>> get_c() const { 
                    misc::confirm(false, "get_c: NOT IMPLEMENTED!");
                    return std::vector<std::complex<double>>();
                }
                int get_s() const {
                    misc::confirm(false, "get_s: NOT IMPLEMENTED!");
                    return 0;
                }
                double cal_PE() const {
                    misc::confirm<>(false, "cal_PE: NOT IMPLEMENTED!");
                    return 0.0; 
                }
                double cal_totE() const {
                    misc::confirm<>(false, "cal_totE: NOT IMPLEMENTED!");
                    return 0.0; 
                }

            private:
                // --- basic --- //
                int m_ndim;
                int m_edim;
                double m_mass;
                double m_kT;
                double m_t;
                std::vector<double> m_r;
                std::vector<double> m_v;
                // --- friction --- //
                double m_nuclear_fric;
                // --- Hamiltonian related --- //
                HamiltonianType m_hamiltonian;
                std::vector<double> m_force;  // THIS IS THE DIABATIC FORCE!
                // --- options/flags --- //
                bool m_initialized;
        };


    // --- IMPLEMENTATION --- // 


    // --- ctor/dtor --- //


    template <typename HamiltonianType>
        FP2Bath_Trajectory<HamiltonianType>::FP2Bath_Trajectory(const HamiltonianType& hamiltonian) 
        :   m_ndim(0), m_edim(0), m_mass(0.0), m_kT(0.0), m_t(0.0),
        m_r(), m_v(), 
        m_nuclear_fric(-1.0), 
        m_hamiltonian(hamiltonian), 
        m_initialized(false)
    {  
        m_edim = hamiltonian.get_dim();
    }


    // --- simulation --- //


    template <typename HamiltonianType>
        void FP2Bath_Trajectory<HamiltonianType>::update_status() {
            /**
             * update status for the trajectory current configuration
             */
            // calculate useful quantities
            const std::vector<std::complex<double>> H = m_hamiltonian.cal_H(m_r);
            const std::vector<std::vector<std::complex<double>>> nablaH = m_hamiltonian.cal_nablaH(m_r);
            const double h = H.at(1+1*m_edim).real() - H.at(0+0*m_edim).real();
            const double gamma_l = m_hamiltonian.cal_gamma_l(m_r);
            const double gamma_r = m_hamiltonian.cal_gamma_r(m_r);
            const double gamma = gamma_l + gamma_r;
            const double mu_l = m_hamiltonian.cal_mu_l();
            const double mu_r = m_hamiltonian.cal_mu_r();
            const double f_l = misc::fermi((h - mu_l) / m_kT);
            const double f_r = misc::fermi((h - mu_r) / m_kT);
            const double f = (f_l * gamma_l + f_r * gamma_r) / gamma;
            // evaluate new force
            for (int i(0); i < m_ndim; ++i) {
                const double dhdri = nablaH.at(i).at(1+1*m_edim).real() - nablaH.at(i).at(0+0*m_edim).real();
                m_force.at(i) = -nablaH.at(i).at(0+0*m_edim).real() - dhdri * f;
            }
        }


    template <typename HamiltonianType>
        void FP2Bath_Trajectory<HamiltonianType>::nuclear_integrator(double dt) {
            /**
             * nuclear integrator -- velocity verlet
             */
            // evaluate friction & noise tensors
            std::vector<double> fric, noise;
            cal_fric_noise_tensor(fric, noise);
            // generate noise force vector
            std::vector<double> noise_force(m_ndim, 0.0);
            std::vector<double> eva, evt;
            matrixop::hdiag(noise, eva, evt);
            for (int i(0); i < m_ndim; ++i) {
                noise_force.at(i) = randomer::normal(0.0, sqrt(2.0 * eva.at(i) / dt));
            }
            noise_force = matrixop::matvec(evt, noise_force);
            // velocity verlet
            m_v += 0.5 * dt / m_mass * (m_force - matrixop::matvec(fric, m_v) + noise_force);
            m_r += m_v * dt;
            update_status();
            m_v += 0.5 * dt / m_mass * (m_force - matrixop::matvec(fric, m_v) + noise_force);
        }

    template <typename HamiltonianType>
        void FP2Bath_Trajectory<HamiltonianType>::cal_fric_noise_tensor(std::vector<double>& fric, std::vector<double>& noise) const {
            /**
             * calculate friction and noise tensors for the system
             *  ONLY WORK FOR TWO-STATE MODEL AND CONSTANT GAMMA
             */
            // check
            misc::confirm(m_edim == 2, "cal_fric_noise_tensor: only works for 2-state Hamiltonian!");
            // calculate useful quantities
            const std::vector<std::complex<double>> H = m_hamiltonian.cal_H(m_r);
            const std::vector<std::vector<std::complex<double>>> nablaH = m_hamiltonian.cal_nablaH(m_r);
            const double h = H.at(1+1*m_edim).real() - H.at(0+0*m_edim).real();
            std::vector<double> dhdr(m_ndim, 0.0);
            for (int i(0); i < m_ndim; ++i) {
                dhdr.at(i) = nablaH.at(i).at(1+1*m_edim).real() - nablaH.at(i).at(0+0*m_edim).real();
            }
            const double gamma_l = m_hamiltonian.cal_gamma_l(m_r);
            const double gamma_r = m_hamiltonian.cal_gamma_r(m_r);
            const double gamma = gamma_l + gamma_r;
            const double mu_l = m_hamiltonian.cal_mu_l();
            const double mu_r = m_hamiltonian.cal_mu_r();
            const double f_l = misc::fermi((h - mu_l) / m_kT);
            const double f_r = misc::fermi((h - mu_r) / m_kT);
            const double f = (f_l * gamma_l + f_r * gamma_r) / gamma;
            const std::vector<double> df_ldr = -1.0 / m_kT * f_l * (1.0 - f_l) * dhdr;
            const std::vector<double> df_rdr = -1.0 / m_kT * f_r * (1.0 - f_r) * dhdr;
            const std::vector<double> dfdr = (gamma_l * df_ldr + gamma_r * df_rdr) / gamma;
            // el part 
            fric.assign(m_ndim * m_ndim, 0.0);
            noise.assign(m_ndim * m_ndim, 0.0);
            for (int i(0); i < m_ndim; ++i) {
                for (int j(0); j < m_ndim; ++j) {
                    fric.at(i+j*m_ndim) = -1.0 / gamma * dfdr.at(j) * dhdr.at(i);
                    noise.at(i+j*m_ndim) = 1.0 / gamma * f * (1.0 - f) * dhdr.at(j) * dhdr.at(i);
                }
            }
            // nu part
            if (m_nuclear_fric > 0.0) {
                for (int i(0); i < m_ndim; ++i) {
                    fric.at(i+i*m_ndim) += m_nuclear_fric;
                    noise.at(i+i*m_ndim) += m_kT * m_nuclear_fric;
                }
            }
        }


    // --- simulation interfaces --- //


    template <typename HamiltonianType>
        void FP2Bath_Trajectory<HamiltonianType>::setup( double mass, const std::vector<double>& r, const std::vector<double>& v, double kT) { 
            /**
             * setup trajectory for simulation
             */
            // check
            misc::confirm<misc::ValueError>(r.size() == v.size(), "init_state: invalid r/v sizes.");
            // setup
            m_ndim = r.size();
            m_force.resize(m_ndim);
            set_mass(mass);
            set_r(r);
            set_v(v);
            set_kT(kT);
            set_t(0.0);
            // update status
            update_status();
            m_initialized = true;
        }

    template <typename HamiltonianType>
        void FP2Bath_Trajectory<HamiltonianType>::die() {
            /**
             * kill the trajectory, cannot evolve any more
             */
            std::vector<double>().swap(m_force);
            m_initialized = false;
        }

    template <typename HamiltonianType>
        void FP2Bath_Trajectory<HamiltonianType>::integrator(double dt) {
            /**
             * propagate trajectory forward by dt
             */
            misc::confirm<misc::StateError>(m_initialized, "integrator: FP2Bath_Trajectory died / not initialzied.");

            // --- nuclear part --- //
            nuclear_integrator(dt);

            // --- time part --- //
            m_t += dt;
        }

} // namespace mqc

#endif // _FP2BATH_TRAJECTORY_HPP
