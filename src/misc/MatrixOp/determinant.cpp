#include <vector>
#include <cmath>
#include "determinant.hpp"

namespace matrixop {


    // --- helper --- //


    ITYPE _permutation_sign(CNST_ITYPE* permutation, CNST_ITYPE N) {
        MATRIXOP_STATIC std::vector<BOOL> visited;
        // calculate sign of a permutation
        visited.assign(N, false);
        ITYPE rst(1), next(0), L(0);
        for (ITYPE k(1); k <= N; ++k) {
            if (!visited[k - 1]) {
                next = k;
                L = 0;
                while (!visited[next - 1]) {
                    ++L;
                    visited[next - 1] = true;
                    next = permutation[next - 1]; 
                }
                if (L % 2 == 0) {
                    rst *= -1; 
                }
            }
        }   
        return rst;
    }


    // --- determinant --- //
    

    // single
    STYPE _determinant(STYPE* A, CNST_ITYPE N) {
        MATRIXOP_STATIC std::vector<ITYPE> ipiv;
        MATRIXOP_STATIC std::vector<ITYPE> permutation;

        ITYPE info(0);

        ipiv.resize(N);
        SGETRF(&N, &N, A, &N, &ipiv[0], &info);

        permutation.resize(N);
        for (ITYPE i(0); i < N; ++i) {
            permutation[i] = i + 1;
        }
        for (ITYPE i(0); i < N; ++i) {
            std::swap(permutation[i], permutation[ipiv[i] - 1]);
        }

        STYPE rst(_permutation_sign(&permutation[0], N));

        CNST_ITYPE interval(N + 1);
        ITYPE idx(0);
        for (ITYPE i(0); i < N; ++i) {
            rst *= A[idx];
            idx += interval;
        }
        return rst;
    }

    // double
    DTYPE _determinant(DTYPE* A, CNST_ITYPE N) {
        MATRIXOP_STATIC std::vector<ITYPE> ipiv;
        MATRIXOP_STATIC std::vector<ITYPE> permutation;

        ITYPE info(0);

        ipiv.resize(N);
        DGETRF(&N, &N, A, &N, &ipiv[0], &info);

        permutation.resize(N);
        for (ITYPE i(0); i < N; ++i) {
            permutation[i] = i + 1;
        }
        for (ITYPE i(0); i < N; ++i) {
            std::swap(permutation[i], permutation[ipiv[i] - 1]);
        }

        DTYPE rst(_permutation_sign(&permutation[0], N));

        CNST_ITYPE interval(N + 1);
        ITYPE idx(0);
        for (ITYPE i(0); i < N; ++i) {
            rst *= A[idx];
            idx += interval;
        }
        return rst;
    }

    // complex
    CTYPE _determinant(CTYPE* A, CNST_ITYPE N) {
        MATRIXOP_STATIC std::vector<ITYPE> ipiv;
        MATRIXOP_STATIC std::vector<ITYPE> permutation;

        ITYPE info(0);

        ipiv.resize(N);
        CGETRF(&N, &N, A, &N, &ipiv[0], &info);

        permutation.resize(N);
        for (ITYPE i(0); i < N; ++i) {
            permutation[i] = i + 1;
        }
        for (ITYPE i(0); i < N; ++i) {
            std::swap(permutation[i], permutation[ipiv[i] - 1]);
        }

        CTYPE rst(_permutation_sign(&permutation[0], N));

        CNST_ITYPE interval(N + 1);
        ITYPE idx(0);
        for (ITYPE i(0); i < N; ++i) {
            rst *= A[idx];
            idx += interval;
        }
        return rst;
    }

    // complex*16
    ZTYPE _determinant(ZTYPE* A, CNST_ITYPE N) {
        MATRIXOP_STATIC std::vector<ITYPE> ipiv;
        MATRIXOP_STATIC std::vector<ITYPE> permutation;

        ITYPE info(0);

        ipiv.resize(N);
        ZGETRF(&N, &N, A, &N, &ipiv[0], &info);

        permutation.resize(N);
        for (ITYPE i(0); i < N; ++i) {
            permutation[i] = i + 1;
        }
        for (ITYPE i(0); i < N; ++i) {
            std::swap(permutation[i], permutation[ipiv[i] - 1]);
        }

        ZTYPE rst(_permutation_sign(&permutation[0], N));

        CNST_ITYPE interval(N + 1);
        ITYPE idx(0);
        for (ITYPE i(0); i < N; ++i) {
            rst *= A[idx];
            idx += interval;
        }
        return rst;
    }

} // namespace matrixop
