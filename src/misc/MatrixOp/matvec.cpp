#include <vector>
#include "matrixop_config.hpp"
#include "matvec.hpp"

namespace matrixop {
    using std::vector;


    // --- matvec --- //
    

	VOID _matvec(CNST_STYPE* A, CNST_STYPE* V, CNST_ITYPE M, CNST_ITYPE N, STYPE* rst) {
        SGEMV(&CHARN, &M, &N, &ONES, A, &M, V, &ONEI, &ZEROS, rst, &ONEI);
	}

	VOID _matvec(CNST_DTYPE* A, CNST_DTYPE* V, CNST_ITYPE M, CNST_ITYPE N, DTYPE* rst) {
        DGEMV(&CHARN, &M, &N, &ONED, A, &M, V, &ONEI, &ZEROD, rst, &ONEI);
	}

	VOID _matvec(CNST_CTYPE* A, CNST_CTYPE* V, CNST_ITYPE M, CNST_ITYPE N, CTYPE* rst) {
        CGEMV(&CHARN, &M, &N, &ONEC, A, &M, V, &ONEI, &ZEROC, rst, &ONEI);
	}

	VOID _matvec(CNST_ZTYPE* A, CNST_ZTYPE* V, CNST_ITYPE M, CNST_ITYPE N, ZTYPE* rst) {
        ZGEMV(&CHARN, &M, &N, &ONEZ, A, &M, V, &ONEI, &ZEROZ, rst, &ONEI);
	}

};
