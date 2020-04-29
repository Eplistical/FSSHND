#ifndef _MATRIXOP_MATVEC_HPP
#define _MATRIXOP_MATVEC_HPP

/*  
 *  matrixop::matvec(A, V)
 *
 *  function for matrix-vector multiplication
 *
 * 	param   A:  matrix (M*N)
 * 	param   V:  vector (N)
 *
 * 	return  A * v
 *
 */


#include <vector>
#include "matrixop_config.hpp"

namespace matrixop {
	using std::vector;

    // --- matvec --- //

	VOID _matvec(CNST_STYPE* A, CNST_STYPE* V, CNST_ITYPE M, CNST_ITYPE N, STYPE* rst);
	VOID _matvec(CNST_DTYPE* A, CNST_DTYPE* V, CNST_ITYPE M, CNST_ITYPE N, DTYPE* rst);
	VOID _matvec(CNST_CTYPE* A, CNST_CTYPE* V, CNST_ITYPE M, CNST_ITYPE N, CTYPE* rst);
	VOID _matvec(CNST_ZTYPE* A, CNST_ZTYPE* V, CNST_ITYPE M, CNST_ITYPE N, ZTYPE* rst);

    // interfaces
    template<typename T>
        vector<T> matvec(const vector<T>& A, const vector<T>& V) {
            CNST_ITYPE N(static_cast<ITYPE>(V.size()));
            CNST_ITYPE M(static_cast<ITYPE>(A.size() / N));
            vector<T> rst(N);
            _matvec(&A[0], &V[0], M, N, &rst[0]);
            return rst;
        }
};

#endif // _MATRIXOP_MATVEC_HPP
