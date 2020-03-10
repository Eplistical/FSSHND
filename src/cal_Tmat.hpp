#ifndef _CAL_TMAT_HPP
#define _CAL_TMAT_HPP

/* 
 * Function to calculate time derivative matrix
 *  Steal from Zeyu
 */

#include <vector>
#include <complex>

namespace mqc {
    namespace zeyu {

    std::vector<std::complex<double>> cal_Tmat(const std::vector<std::complex<double>>& /* curevt */, std::vector<std::complex<double>>& /* nextevt */, int /* NNN */, double /* timestep */);

    } // namespace zeyu
    
} // namespace mqc


#endif // _CAL_TMAT_HPP