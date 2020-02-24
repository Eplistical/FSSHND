#include <vector>
#include <complex>
#include <numeric>
#include <algorithm>
#include "exceptions.hpp"
#include "fssh_trajectory.hpp"

namespace {


    // --- ctor/dtor --- //


    template <typename HamiltonianType>
        FSSH_Trajectory<HamiltonianType>::FSSH_Trajectory(const HamiltonianType& hamiltonian) 
        : m_hamiltonian(hamiltonian), m_enable_fric(false)
        {  }


    // --- simulation --- //


    template <typename HamiltonianType>
        void FSSH_Trajectory<HamiltonianType>::init_state( const std::vector<double>& r, const std::vector<double>& v,
                const std::vector<std::complex<double>>& c, int s) { 
            set_r(r);
            set_v(v);
            set_c(c);
            set_s(s);
            m_hamiltonian.cal_info(m_r, m_force, m_dc, m_eva, m_evt);
        }

    template <typename HamiltonianType>
        void FSSH_Trajectory<HamiltonianType>::integrator(double dt) {
            // TODO
        }

    template <typename HamiltonianType>
        void FSSH_Trajectory<HamiltonianType>::hopper() {
            // TODO
        }

    template <typename HamiltonianType>
        bool FSSH_Trajectory<HamiltonianType>::check_end() const {
            return false;
        }
    

    // --- quantities --- //


    template <typename HamiltonianType>
        double FSSH_Trajectory<HamiltonianType>::cal_KE() const {
            return 0.5 * m_mass * std::accumulate(m_v.begin(), m_v.end(), 0.0, 
                    [](double accu, double vi) { return accu + vi*vi; });
        }

    template <typename HamiltonianType>
        double FSSH_Trajectory<HamiltonianType>::cal_PE() const {
            confirm(not m_eva.empty(), "cal_PE: eigenvalues not yet calculated.");
            return m_eva.at(m_s);
        }

};
