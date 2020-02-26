#ifndef _MATRIXOP_EIG_HPP
#define _MATRIXOP_EIG_HPP

/*  
 *  matrixop::eigh(A, eva, [evt])  matrixop::hdiag(A, eva, [evt])
 *
 *  function for diaganalizing a Hermitian matrix A
 *
 *
 * 	param   A:  Hermitian matrix A (N*N)
 *
 * 	param eva:  result eigenvalues
 * 	param evt:  result eigenvectors
 *
 * the function detect the lower triangle of the matrix
 *
 ****************************************
 *
 *  matrixop::eigh_inplace(A, eva)  matrixop::hdiag_inplace(A, eva)
 *
 *  function for inplace diaganalizing a Hermitian matrix A
 *
 *
 * 	param   A:  Hermitian matrix A (N*N)
 * 	            on exit, A is result eigenvectors
 *
 * 	param eva:  result eigenvalues
 *
 * the function detect the lower triangle of the matrix
 *
 */


#include <vector>
#include <cmath>
#include <algorithm>
#include "matrixop_config.hpp"

namespace matrixop {
	using std::sqrt;
	using std::copy;
	using std::vector;


    // --- eigh --- //


	VOID _eigh(STYPE* A, STYPE* eva, CNST_ITYPE N, CNST_CHAR jobz);
	VOID _eigh(DTYPE* A, DTYPE* eva, CNST_ITYPE N, CNST_CHAR jobz);
	VOID _eigh(CTYPE* A, STYPE* eva, CNST_ITYPE N, CNST_CHAR jobz);
	VOID _eigh(ZTYPE* A, DTYPE* eva, CNST_ITYPE N, CNST_CHAR jobz);


    // --- eig --- //


	VOID _eig(STYPE* A, STYPE* eva, CNST_ITYPE N, CNST_CHAR jobz);
	VOID _eig(DTYPE* A, DTYPE* eva, CNST_ITYPE N, CNST_CHAR jobz);
	VOID _eig(CTYPE* A, CTYPE* eva, CNST_ITYPE N, CNST_CHAR jobz);
	VOID _eig(ZTYPE* A, ZTYPE* eva, CNST_ITYPE N, CNST_CHAR jobz);


    // --- interfaces --- //


    template <typename T1, typename T2>
        VOID eigh(const vector<T1>& A, vector<T2>& eva, vector<T1>& evt) {
            CNST_ITYPE N(static_cast<ITYPE>(sqrt(A.size())));
		    eva.resize(N);
            evt.resize(A.size());
            copy(A.begin(), A.end(), evt.begin());

            _eigh(&evt[0], &eva[0], N, CHARV);
        }

    template <typename T1, typename T2>
        VOID hdiag(const vector<T1>& A, vector<T2>& eva, vector<T1>& evt) {
            eigh(A, eva, evt);
        }

    template <typename T1, typename T2>
        VOID eigh(const vector<T1>& A, vector<T2>& eva) {
		    MATRIXOP_STATIC vector<T1> evt;

            CNST_ITYPE N(static_cast<ITYPE>(sqrt(A.size())));
		    eva.resize(N);
            evt.resize(A.size());
            copy(A.begin(), A.end(), evt.begin());

            _eigh(&evt[0], &eva[0], N, CHARN);
        }

    template <typename T1, typename T2>
        VOID hdiag(const vector<T1>& A, vector<T2>& eva) {
            eigh(A, eva);
        }

    template <typename T1, typename T2>
        VOID eigh_inplace(vector<T1>& A, vector<T2>& eva) {
            CNST_ITYPE N(static_cast<ITYPE>(sqrt(A.size())));
		    eva.resize(N);
            _eigh(&A[0], &eva[0], N, CHARV);
        }

    template <typename T1, typename T2>
        VOID hdiag_inplace(vector<T1>& A, vector<T2>& eva) {
            eigh_inplace(A, eva);
        }
};

#endif // _MATRIXOP_EIG_HPP
