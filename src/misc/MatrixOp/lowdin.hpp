#ifndef _MATRIXOP_LOWDIN_HPP
#define _MATRIXOP_LOWDIN_HPP

/*  
 *  matrixop::lowdin(A, N) 
 *
 *  function for orthogonalizing matrix A (Lowdin algorithm)
 *  i.e. A => A * (A^H * A)^(-1/2)
 *
 *  param   A: matrix A (M*N)
 *  param   N: # columns in A
 *
 *  return orthogonalized matrix
 *  
 */

#include <vector>
#include <cmath>
#include "matrixop_config.hpp"
#include "matmat.hpp"
#include "mpower.hpp"


namespace matrixop {
    using std::vector;

    template <typename T>
        inline vector<T> lowdin(const vector<T>& A, CNST_ITYPE N) {
            return matmat(A, mpowerh(matCmat(A, A, N), -0.5), N);
        }

} // namespace matrixop

#endif
