#ifndef _MATRIXOP_MPOWER_HPP
#define _MATRIXOP_MPOWER_HPP

/*  
 *  matrixop::mpowerh(A, p) 
 *
 *  function for calculating matrix power A^p
 *
 *  param   A: Hermitian matrix A (N*N)
 *  param   p: power
 *
 *  return A^p
 *
 */

#include <vector>
#include <cmath>
#include "matrixop_config.hpp"
#include "matmat.hpp"
#include "eig.hpp"


namespace matrixop {
    using std::vector;

    template <typename T1, typename T2>
        vector<T1> mpowerh(const vector<T1>& A, const T2& p) {
            vector<DTYPE> eva;
            vector<T1> evamat;
            vector<T1> evt;
            eigh(A, eva, evt);
            CNST_ITYPE N = static_cast<ITYPE>(sqrt(A.size()));
            // construct evamat
            evamat.resize(N * N);
            for (ITYPE i(0); i < N; ++i) {
                for (ITYPE j(0); j < N; ++j) {
                    if (i == j) {
                        if (p < 1.0 and eva[i] < 0.0) {
                            throw std::runtime_error("matrixop::mpowerh : negative eigenvalues encountered while p < 1.0");
                        }
                        evamat[i+j*N] = std::pow(eva[i], p);
                    }
                    else {
                        evamat[i+j*N] = static_cast<T1>(0.0);
                    }
                }
            }
            return matmatmatC(evt, evamat, evt, N, N);
        }

} // namespace matrixop

#endif
