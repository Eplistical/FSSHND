#ifndef _ADJUST_EVT_PHASE_HPP
#define _ADJUST_EVT_PHASE_HPP

/* 
 * Function to calculate time derivative matrix
 *  Steal from Zeyu
 */

#include <vector>
#include <complex>

namespace mqc {
    namespace zeyu {

        std::vector<std::complex<double>> adjust_evt_phase(const std::vector<std::complex<double>>& /* curevt */, std::vector<std::complex<double>>& /* nextevt */, int /* NNN */, double /* timestep */);

    } // namespace zeyu
    
} // namespace mqc


#endif // _ADJUST_EVT_PHASE_HPP