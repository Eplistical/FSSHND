#include <vector>
#include <cmath>
#include "matrixop_config.hpp"
#include "matmat.hpp"
#include "logm.hpp"


namespace matrixop {
    using std::vector;

    VOID _logmh(vector<ZTYPE>& A) {
        /**
         * calculate logm(A), where A must be an orthogonal complex-valued matrix.
         */
        CNST_ITYPE N = static_cast<ITYPE>(sqrt(A.size()));
        // check orthogonality
        vector<ZTYPE> AA(matCmat(A, A, N));
        for (int i(0); i < N; ++i) {
            for (int j(0); j < N; ++j) {
                if ((i == j and std::abs(AA[i+j*N] - 1.0) > 1e-10) or
                            (i != j and std::abs(AA[i+j*N]) > 1e-10))
                  throw std::runtime_error("matrixop::logmh : input matrix is not orthogonal!");
            }
        }
        // check positive elements
        for (int i(0); i < N; ++i) {
            if (std::abs(A[i+i*N].imag()) > 1e-10) {
                throw std::runtime_error("matrixop::logmh : complex-valued diagonal elements found!");
            }
            else if (A[i+i*N].real() <= 0.0) {
                throw std::runtime_error("matrixop::logmh : non-positive diagonal elements found!");
            }
        }
        // Schur decomposition A = Q * T * Q'
        vector<ZTYPE> w;
        vector<ZTYPE> work;
        vector<ZTYPE> vs;
        vector<DTYPE> rwork;
        ITYPE lwork, info, SDIM(0);
        w.resize(N);
        vs.resize(N * N);
        rwork.resize(N);

        lwork = -1;
        work.resize(1);
        ZGEES("V", "N", nullptr, &N, &A[0], &N, &SDIM, 
                    &w[0], &vs[0], &N, &work[0], &lwork, &rwork[0], nullptr, &info);
        lwork = static_cast<int>(work[0].real());
        work.resize(lwork);
        ZGEES("V", "N", nullptr, &N, &A[0], &N, &SDIM, 
                    &w[0], &vs[0], &N, &work[0], &lwork, &rwork[0], nullptr, &info);
        // construct diagonal unitary matrix logD (stored in A)
        for (ITYPE i(0); i < N; ++i) {
            for (ITYPE j(0); j < N; ++j) {
                if (i == j) 
                  A[i+j*N] = std::log(A[i+j*N] / std::abs(A[i+j*N]));
                else 
                  A[i+j*N] = 0.0;
            }
        }
        // get Q * log(D) * Q'
        A = matmatmatC(vs, A, vs, N, N);
    }


} // namespace matrixop
