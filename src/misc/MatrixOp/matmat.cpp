#include <vector>
#include "matrixop_config.hpp"
#include "matmat.hpp"

namespace matrixop {
    using std::vector;


    // --- matmat --- //


	VOID _matmat(CNST_STYPE* A, CNST_STYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, STYPE* rst) {
        SGEMM(&CHARN, &CHARN, &M, &N, &K, &ONES, A, &M, B, &K, &ZEROS, rst, &M);
	}

	VOID _matmat(CNST_DTYPE* A, CNST_DTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, DTYPE* rst) {
        DGEMM(&CHARN, &CHARN, &M, &N, &K, &ONED, A, &M, B, &K, &ZEROD, rst, &M);
	}

	VOID _matmat(CNST_CTYPE* A, CNST_CTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, CTYPE* rst) {
        CGEMM(&CHARN, &CHARN, &M, &N, &K, &ONEC, A, &M, B, &K, &ZEROC, rst, &M);
	}

	VOID _matmat(CNST_ZTYPE* A, CNST_ZTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, ZTYPE* rst) {
        ZGEMM(&CHARN, &CHARN, &M, &N, &K, &ONEZ, A, &M, B, &K, &ZEROZ, rst, &M);
	}


    // --- matCmat --- //


	VOID _matCmat(CNST_STYPE* A, CNST_STYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, STYPE* rst) {
        SGEMM(&CHARC, &CHARN, &M, &N, &K, &ONES, A, &K, B, &K, &ZEROS, rst, &M);
	}

	VOID _matCmat(CNST_DTYPE* A, CNST_DTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, DTYPE* rst) {
        DGEMM(&CHARC, &CHARN, &M, &N, &K, &ONED, A, &K, B, &K, &ZEROD, rst, &M);
	}

	VOID _matCmat(CNST_CTYPE* A, CNST_CTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, CTYPE* rst) {
        CGEMM(&CHARC, &CHARN, &M, &N, &K, &ONEC, A, &K, B, &K, &ZEROC, rst, &M);
	}

	VOID _matCmat(CNST_ZTYPE* A, CNST_ZTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, ZTYPE* rst) {
        ZGEMM(&CHARC, &CHARN, &M, &N, &K, &ONEZ, A, &K, B, &K, &ZEROZ, rst, &M);
	}


    // --- matmatC --- //


	VOID _matmatC(CNST_STYPE* A, CNST_STYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, STYPE* rst) {
        SGEMM(&CHARN, &CHARC, &M, &N, &K, &ONES, A, &M, B, &N, &ZEROS, rst, &M);
	}

	VOID _matmatC(CNST_DTYPE* A, CNST_DTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, DTYPE* rst) {
        DGEMM(&CHARN, &CHARC, &M, &N, &K, &ONED, A, &M, B, &N, &ZEROD, rst, &M);
	}

	VOID _matmatC(CNST_CTYPE* A, CNST_CTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, CTYPE* rst) {
        CGEMM(&CHARN, &CHARC, &M, &N, &K, &ONEC, A, &M, B, &N, &ZEROC, rst, &M);
	}

	VOID _matmatC(CNST_ZTYPE* A, CNST_ZTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, ZTYPE* rst) {
        ZGEMM(&CHARN, &CHARC, &M, &N, &K, &ONEZ, A, &M, B, &N, &ZEROZ, rst, &M);
	}


    // --- matCmatC --- //
    

	VOID _matCmatC(CNST_STYPE* A, CNST_STYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, STYPE* rst) {
        SGEMM(&CHARC, &CHARC, &M, &N, &K, &ONES, A, &K, B, &N, &ZEROS, rst, &M);
	}

	VOID _matCmatC(CNST_DTYPE* A, CNST_DTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, DTYPE* rst) {
        DGEMM(&CHARC, &CHARC, &M, &N, &K, &ONED, A, &K, B, &N, &ZEROD, rst, &M);
	}

	VOID _matCmatC(CNST_CTYPE* A, CNST_CTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, CTYPE* rst) {
        CGEMM(&CHARC, &CHARC, &M, &N, &K, &ONEC, A, &K, B, &N, &ZEROC, rst, &M);
	}

	VOID _matCmatC(CNST_ZTYPE* A, CNST_ZTYPE* B, CNST_ITYPE M, CNST_ITYPE K, CNST_ITYPE N, ZTYPE* rst) {
        ZGEMM(&CHARC, &CHARC, &M, &N, &K, &ONEZ, A, &K, B, &N, &ZEROZ, rst, &M);
	}

};
