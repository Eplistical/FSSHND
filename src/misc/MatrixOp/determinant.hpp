#ifndef _MATRIXOP_DETERMINANT_HPP
#define _MATRIXOP_DETERMINANT_HPP

/*  
 *  matrixop::determinant(A)  matrixop::det(A)
 *
 *  function for calculating determinant of a matrix A
 *
 *
 * 	param   A:  matrix A (N*N)
 *
 *  return  determinant of A
 *
 */


#include <vector>
#include <cmath>
#include "matrixop_config.hpp"

namespace matrixop {


    // --- helper --- //


    ITYPE _permutation_sign(CNST_ITYPE* permutation, CNST_ITYPE N);


    // --- determinant --- //
    

    STYPE _determinant(STYPE* A, CNST_ITYPE N);
    DTYPE _determinant(DTYPE* A, CNST_ITYPE N);
    CTYPE _determinant(CTYPE* A, CNST_ITYPE N);
    ZTYPE _determinant(ZTYPE* A, CNST_ITYPE N); 


    // --- interfaces --- //

    template <typename T>
        inline T determinant(const std::vector<T>& A) {
            MATRIXOP_STATIC std::vector<T> rst;
            rst.assign(A.begin(), A.end());
            return _determinant(&rst[0], static_cast<ITYPE>(std::sqrt(rst.size())));
        }

    template <typename T>
        inline T det(const std::vector<T>& A) {
            return determinant(A);
        }

} // namespace matrixop

#endif // _MATRIXOP_DETERMINANT_HPP
