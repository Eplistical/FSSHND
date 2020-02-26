#include <vector>
#include <complex>
#include <numeric>
#include <algorithm>
#include "misc/crasher.hpp"
#include "misc/matrixop.hpp"
#include "misc/vector.hpp"
#include "fssh_trajectory.hpp"

namespace {


    // --- ctor/dtor --- //


    template <typename HamiltonianType>
        FSSH_Trajectory<HamiltonianType>::FSSH_Trajectory(const HamiltonianType& hamiltonian) 
        : m_hamiltonian(hamiltonian), m_enable_fric(false), m_initialized(false)
        {  
            m_edim = hamiltonian.get_dim();
        }


    // --- simulation --- //


    template <typename HamiltonianType>
        void FSSH_Trajectory<HamiltonianType>::init_state( const std::vector<double>& r, const std::vector<double>& v,
                const std::vector<std::complex<double>>& c, int s) { 
            // setup trajectory for r, v, c, s
            // check
            misc::confirm<misc::ValueError>(r.size() == v.size(), "init_state: invalid r/v sizes.");
            misc::confirm<misc::ValueError>(c.size() == m_edim, "init_state: invalid c size.");
            misc::confirm<misc::ValueError>(s < m_edim, "init_state: invalid s number.");
            // setup
            m_ndim = m_r.size();
            set_r(r);
            set_v(v);
            set_c(c);
            set_s(s);
            // update status
            update_status();
        }
            
    template <typename HamiltonianType>
        void FSSH_Trajectory<HamiltonianType>::integrator(double dt) {
            // propagate trajectory forward by dt
            misc::confirm<misc::ValueError>(m_initialized, "FSSH_Trajectory not initialzied.");

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
            update_status();
        }

    template <typename HamiltonianType>
        void FSSH_Trajectory<HamiltonianType>::hopper() {
            misc::confirm<misc::ValueError>(m_initialized, "FSSH_Trajectory not initialzied.");
            // TODO : WRITE HOPPER
        }

    template <typename HamiltonianType>
        bool FSSH_Trajectory<HamiltonianType>::check_end() const {
            return false;
        }
    

    // --- quantities --- //


    template <typename HamiltonianType>
        double FSSH_Trajectory<HamiltonianType>::cal_KE() const {
            misc::confirm<misc::ValueError>(m_initialized, "FSSH_Trajectory not initialzied.");
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

};    
