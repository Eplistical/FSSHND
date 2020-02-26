#include <vector>
#include <cmath>
#include <algorithm>
#include "matrixop_config.hpp"
#include "eig.hpp"

namespace matrixop {
	using std::vector;


    // --- eigh --- //


	VOID _eigh(STYPE* A, STYPE* eva, CNST_ITYPE N, CNST_CHAR jobz) {
		MATRIXOP_STATIC vector<STYPE> work;
		MATRIXOP_STATIC vector<ITYPE> iwork;

		ITYPE lwork(-1), liwork(-1), info(-1);
        work.resize(1);
        iwork.resize(1);

		SSYEVD(&jobz, &CHARL, &N, A, &N, eva, 
				&work[0], &lwork, &iwork[0], &liwork, &info);

		lwork = static_cast<ITYPE>(work[0]);
		liwork = iwork[0];
		work.resize(lwork);
		iwork.resize(liwork);

		SSYEVD(&jobz, &CHARL, &N, A, &N, eva, 
				&work[0], &lwork, &iwork[0], &liwork, &info);
	}

	VOID _eigh(DTYPE* A, DTYPE* eva, CNST_ITYPE N, CNST_CHAR jobz) {
		MATRIXOP_STATIC vector<DTYPE> work;
		MATRIXOP_STATIC vector<ITYPE> iwork;

		ITYPE lwork(-1), liwork(-1), info(-1);
        work.resize(1);
        iwork.resize(1);

		DSYEVD(&jobz, &CHARL, &N, A, &N, eva, 
				&work[0], &lwork, &iwork[0], &liwork, &info);

		lwork = static_cast<ITYPE>(work[0]);
		liwork = iwork[0];
		work.resize(lwork);
		iwork.resize(liwork);

		DSYEVD(&jobz, &CHARL, &N, A, &N, eva, 
				&work[0], &lwork, &iwork[0], &liwork, &info);
	}

	VOID _eigh(CTYPE* A, STYPE* eva, CNST_ITYPE N, CNST_CHAR jobz) {
		MATRIXOP_STATIC vector<CTYPE> work;
		MATRIXOP_STATIC vector<STYPE> rwork;
		MATRIXOP_STATIC vector<ITYPE> iwork;

		ITYPE lwork(-1), liwork(-1), lrwork(-1), info(-1);
        work.resize(1);
        rwork.resize(1);
        iwork.resize(1);

		CHEEVD(&jobz, &CHARL, &N, A, &N, eva, 
				&work[0], &lwork, &rwork[0], &lrwork, 
				&iwork[0], &liwork, &info);

		lwork = static_cast<ITYPE>(work[0].real());
		liwork = iwork[0];
		lrwork = rwork[0];
		work.resize(lwork);
		iwork.resize(liwork);
		rwork.resize(lrwork);

		CHEEVD(&jobz, &CHARL, &N, A, &N, eva, 
				&work[0], &lwork, &rwork[0], &lrwork, 
				&iwork[0], &liwork, &info);
	}

	VOID _eigh(ZTYPE* A, DTYPE* eva, CNST_ITYPE N, CNST_CHAR jobz) {
		MATRIXOP_STATIC vector<ZTYPE> work;
		MATRIXOP_STATIC vector<DTYPE> rwork;
		MATRIXOP_STATIC vector<ITYPE> iwork;

		ITYPE lwork(-1), liwork(-1), lrwork(-1), info(-1);
        work.resize(1);
        rwork.resize(1);
        iwork.resize(1);

		ZHEEVD(&jobz, &CHARL, &N, A, &N, eva, 
				&work[0], &lwork, &rwork[0], &lrwork, 
				&iwork[0], &liwork, &info);

		lwork = static_cast<ITYPE>(work[0].real());
		liwork = iwork[0];
		lrwork = rwork[0];
		work.resize(lwork);
		iwork.resize(liwork);
		rwork.resize(lrwork);

		ZHEEVD(&jobz, &CHARL, &N, A, &N, eva, 
				&work[0], &lwork, &rwork[0], &lrwork, 
				&iwork[0], &liwork, &info);
	}


    // --- eig --- //


	VOID _eig(STYPE* A, STYPE* eva, CNST_ITYPE N, CNST_CHAR jobz) {
        // TODO
    }

	VOID _eig(DTYPE* A, DTYPE* eva, CNST_ITYPE N, CNST_CHAR jobz) {
        // TODO
    }

	VOID _eig(CTYPE* A, CTYPE* eva, CNST_ITYPE N, CNST_CHAR jobz) {
        // TODO
    }

	VOID _eig(ZTYPE* A, ZTYPE* eva, CNST_ITYPE N, CNST_CHAR jobz) {
        // TODO
    }
};
