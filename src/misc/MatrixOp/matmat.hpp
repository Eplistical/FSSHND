#ifndef _MATRIXOP_MATMAT_HPP
#define _MATRIXOP_MATMAT_HPP

/*  
 *  matrixop::matmat(A, B, K)
 *
 *  function for matrix-matrix multiplication
 *
 * 	param   A:  matrix A (M*K)
 * 	param   B:  matrix B (K*N)
 * 	param   K:  matrix dimensions
 *
 * 	return A * B
 *
 */


#include <vector>
#include "matrixop_config.hpp"

namespace matrixop {
	using std::vector;

    // --- matmat --- //

	VOID _matmat(CNST_STYPE* A, CNST_STYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, STYPE* rst);
	VOID _matmat(CNST_DTYPE* A, CNST_DTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, DTYPE* rst);
	VOID _matmat(CNST_CTYPE* A, CNST_CTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, CTYPE* rst);
	VOID _matmat(CNST_ZTYPE* A, CNST_ZTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, ZTYPE* rst);

    // --- matCmat --- //

	VOID _matCmat(CNST_STYPE* A, CNST_STYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, STYPE* rst);
	VOID _matCmat(CNST_DTYPE* A, CNST_DTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, DTYPE* rst);
	VOID _matCmat(CNST_CTYPE* A, CNST_CTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, CTYPE* rst);
	VOID _matCmat(CNST_ZTYPE* A, CNST_ZTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, ZTYPE* rst);

    // --- matmatC --- //

	VOID _matmatC(CNST_STYPE* A, CNST_STYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, STYPE* rst);
	VOID _matmatC(CNST_DTYPE* A, CNST_DTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, DTYPE* rst);
	VOID _matmatC(CNST_CTYPE* A, CNST_CTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, CTYPE* rst);
	VOID _matmatC(CNST_ZTYPE* A, CNST_ZTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, ZTYPE* rst);

    // --- matCmatC --- //

	VOID _matCmatC(CNST_STYPE* A, CNST_STYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, STYPE* rst);
	VOID _matCmatC(CNST_DTYPE* A, CNST_DTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, DTYPE* rst);
	VOID _matCmatC(CNST_CTYPE* A, CNST_CTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, CTYPE* rst);
	VOID _matCmatC(CNST_ZTYPE* A, CNST_ZTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, ZTYPE* rst);


    // --- interfaces --- //


    template<typename T>
        vector<T> matmat(const vector<T>& A, const vector<T>& B, CNST_ITYPE K) {
            CNST_ITYPE M(static_cast<ITYPE>(A.size() / K));
            CNST_ITYPE N(static_cast<ITYPE>(B.size() / K));
            vector<T> rst(N * M);
            _matmat(&A[0], &B[0], M, K, N, &rst[0]);
            return rst;
        }

    template<typename T>
        vector<T> matCmat(const vector<T>& A, const vector<T>& B, CNST_ITYPE K) {
            CNST_ITYPE M(static_cast<ITYPE>(A.size() / K));
            CNST_ITYPE N(static_cast<ITYPE>(B.size() / K));
            vector<T> rst(N * M);
            _matCmat(&A[0], &B[0], M, K, N, &rst[0]);
            return rst;
        }

    template<typename T>
        vector<T> matmatC(const vector<T>& A, const vector<T>& B, CNST_ITYPE K) {
            CNST_ITYPE M(static_cast<ITYPE>(A.size() / K));
            CNST_ITYPE N(static_cast<ITYPE>(B.size() / K));
            vector<T> rst(N * M);
            _matmatC(&A[0], &B[0], M, K, N, &rst[0]);
            return rst;
        }

    template<typename T>
        vector<T> matCmatC(const vector<T>& A, const vector<T>& B, CNST_ITYPE K) {
            CNST_ITYPE M(static_cast<ITYPE>(A.size() / K));
            CNST_ITYPE N(static_cast<ITYPE>(B.size() / K));
            vector<T> rst(N * M);
            _matCmatC(&A[0], &B[0], M, K, N, &rst[0]);
            return rst;
        }

    template <typename T>
        vector<T> matCmatmat(const vector<T>& A, const vector<T>& B, const vector<T>& C, CNST_ITYPE L, CNST_ITYPE K) {
            return matCmat(A, matmat(B, C, K), L);
        }

    template <typename T>
        vector<T> matmatmatC(const vector<T>& A, const vector<T>& B, const vector<T>& C, CNST_ITYPE L, CNST_ITYPE K) {
            return matmat(A, matmatC(B, C, K), L);
        }

    template <typename T>
        vector<T> matmatCmat(const vector<T>& A, const vector<T>& B, const vector<T>& C, CNST_ITYPE L, CNST_ITYPE K) {
            return matmat(A, matCmat(B, C, K), L);
        }

    template <typename T>
        vector<T> matCmatmatC(const vector<T>& A, const vector<T>& B, const vector<T>& C, CNST_ITYPE L, CNST_ITYPE K) {
            return matCmat(A, matmatC(B, C, K), L);
        }

    template <typename T>
        vector<T> matCmatCmatC(const vector<T>& A, const vector<T>& B, const vector<T>& C, CNST_ITYPE L, CNST_ITYPE K) {
            return matCmat(A, matCmatC(B, C, K), L);
        }
};

#endif // _MATRIXOP_MATMAT_HPP
