#ifndef _MATRIXOP_LOGM_HPP
#define _MATRIXOP_LOGM_HPP

/*  
 *  matrixop::logmh(A) 
 *
 *  function for calculating matrix log(A)
 *
 *  param   A: Hermitian matrix A (N*N)
 *
 *  return logm(A)
 *
 */

#include <vector>
#include <cmath>
#include "matrixop_config.hpp"
#include "matmat.hpp"


namespace matrixop {
    using std::vector;

    // --- logmh --- //


    VOID _logmh(vector<ZTYPE>& A);


    // --- interfaces --- //


    template <typename T>
        VOID logmh_inplace(vector<T>& A) {
            _logmh(A);
        }

    template <typename T>
        vector<T> logmh(const vector<T>& A) {
            vector<T> rst(A);
            logmh_inplace(rst);
            return rst;
        }

} // namespace matrixop

#endif
