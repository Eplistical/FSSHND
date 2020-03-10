#include <vector>
#include <cmath>
#include <iostream>
#include <complex>
#include "cal_Tmat.hpp"
#include "misc/crasher.hpp"
#include "misc/vector.hpp"
#include "misc/matrixop.hpp"


namespace mqc {
    namespace zeyu {

        auto zgemm = matrixop::ZGEMM;
        auto dgemm = matrixop::DGEMM;
        auto zgees = matrixop::ZGEES;
        auto dgees = matrixop::DGEES;
        auto dgetrf = matrixop::DGETRF;
        auto zgetrf = matrixop::ZGETRF;
        auto dsyevd = matrixop::DSYEVD;
        auto zheevd = matrixop::ZHEEVD;

        void matprint(const std::vector<std::complex<double>>& A) {
            for (const auto& x : A) {
                std::cout << x << "\t";
            }
            std::cout << std::endl;
        }
        
        /******* Crucial function #1 *******/
        double rotation_angle(double AA, double BB, double CC, double DD) {
            /*This is to minimize a cos(2 \theta) + b sin(2 \theta) + c cos(\theta) + d sin(\theta) 
             *	    by calculating the roots of the first derivative function, i.e.
             *	    -2a sin(2 \theta) + 2b cos(2 \theta) - c sin(\theta) + d cos(\theta) = 0
             *	  This is equivalent to solving a quartic equation (4th degree algebraic equation)
             *	    x^4 + u x^3 + v x^2 + p x + q = 0
             *	  by diagonalizing its companion matrix and the eigenvalues (lapack routine "zgees") are the roots.
             *	    		|   0   0   0   -q   |
             *	    		|   		     |
             *	    		|   1   0   0   -p   |
             *	          C  =  |		     |
             *	    		|   0   1   0   -v   |
             *	    		|   		     |
             *	    		|   0   0   1   -u   |
             *	    Note that |x| = cos(\theta), and it must be real. So we
             *	    1. Calculate {u, v, p, q} by {a, b, c, d}
             *	    2. Diagonalize the companion matrix C and obtain 4 roots.
             *	    3. Drop the complex roots  and save the real ones.
             *	    4. Calculate arccos(roots) and -arccos(roots)
             *	    5. See if they are actually the root of the first derivative function.
             *	    6. See which of them has the lowest value.
             *
             * */
            double a = 2.0 * BB, b = - 2.0 * AA, c = DD, d = -CC;
            double x4 = 4.0 * (a * a + b * b);
            if(abs(x4) < 0.000000001) {
                double theta = atan(-c / d);
                if(c * sin(theta) - d * cos(theta) < 0.0) {
                    return theta;
                }
                else {
                    return theta + M_PI;
                }

            }
            double u, v, p, q;
            u = 4.0 * (a * c + b * d) / x4;
            v = (c * c - 4.0 * a * a + d * d - 4.0 * b * b) / x4;
            p = -(4.0 * b * d + 2.0 * a * c) / x4;
            q = (a * a - d * d) / x4;
            std::vector<double> companion_mat(16, 0.0);
            companion_mat[1] = 1.0;
            companion_mat[6] = 1.0;
            companion_mat[11] = 1.0;
            companion_mat[12] = -q;
            companion_mat[13] = -p;
            companion_mat[14] = -v;
            companion_mat[15] = -u;
            std::vector<std::complex<double> > A(16, 0.0);
            for (int i(0); i < A.size(); ++i) {
                A[i] = static_cast<std::complex<double> > (companion_mat[i]);
            }
            std::vector<std::complex<double> > w;
            std::vector<std::complex<double> > work;
            std::vector<std::complex<double> > vs;
            std::vector<double> rwork;
            int bwork[1], N = 4;
            int lwork, info, SDIM(0);
            w.resize(N);
            vs.resize(N * N);
            rwork.resize(N);
            lwork = -1;
            work.resize(1);
            zgees("V", "N", NULL, &N, &A[0], &N, &SDIM,
                        &w[0], &vs[0], &N, &work[0], &lwork, &rwork[0], &bwork[0], &info);

            lwork = static_cast<int>(work[0].real());
            work.resize(lwork);
            zgees("V", "N", NULL, &N, &A[0], &N, &SDIM,
                        &w[0], &vs[0], &N, &work[0], &lwork, &rwork[0], &bwork[0], &info);
            int numrealroot(0);
            //vector<double> realroots(0, 0.0);
            double attemptroot, func, mintheta, minval;
            bool flagfirstroot = true;
            int i;
            for(i=0; i < 4; i++) {
                if(abs(imag(A[i * 5])) > 0.00001) {
                    continue;
                }
                else {
                    if(real(A[i * 5]) > 1.0000010 || real(A[i * 5]) < -1.000000100) {
                        continue;
                    }
                    attemptroot = acos(std::max(std::min(real(A[i * 5]), 1.0), -1.0));
                    func = a * cos(2.0 * attemptroot) + b * sin(2.0 * attemptroot) + c * cos(attemptroot) + d * sin(attemptroot);
                    if(abs(func) < 0.001) {
                        if(flagfirstroot == true) {
                            mintheta = attemptroot;
                            minval = AA * cos(2.0 * attemptroot) + BB * sin(2.0 * attemptroot) + CC * cos(attemptroot) + DD * sin(attemptroot);
                            flagfirstroot = false;
                        }
                        else {
                            func = AA * cos(2.0 * attemptroot) + BB * sin(2.0 * attemptroot) + CC * cos(attemptroot) + DD * sin(attemptroot);
                            if(func < minval) {
                                mintheta = attemptroot;
                                minval = func;
                            }
                        }	
                    }
                    else {
                        attemptroot *= -1.0;
                        func = a * cos(2.0 * attemptroot) + b * sin(2.0 * attemptroot) + c * cos(attemptroot) + d * sin(attemptroot);
                        if(abs(func) < 0.000001) {
                            if(flagfirstroot == true) {
                                mintheta = attemptroot;
                                minval = AA * cos(2.0 * attemptroot) + BB * sin(2.0 * attemptroot) + CC * cos(attemptroot) + DD * sin(attemptroot);
                                flagfirstroot = false;
                            }
                            else {
                                func = AA * cos(2.0 * attemptroot) + BB * sin(2.0 * attemptroot) + CC * cos(attemptroot) + DD * sin(attemptroot);
                                if(func < minval) {
                                    mintheta = attemptroot;
                                    minval = func;
                                }
                            }
                        }
                        else {
                            std::cout << "There is a problem!" << std::endl;
                            std::cout << a << " " << b << " " << c << " " << d << std::endl;
                            matprint(A);
                            std::cout << "Ai" << real(A[0]) + 1.0 << "\nAttempt root" << attemptroot << "\nFunc:" << func<< std::endl;
                            break;
                        }
                    }
                }
            }
            if(flagfirstroot == false) {
                return mintheta;
            }
            else{
                std::cout << "There is a problem2!" << std::endl;
                flagfirstroot = true;
                std::cout << a << " " << b << " " << c << " " << d << std::endl;
                std::cout.precision(9);
                std::cout << "Ai" << real(A[0]) + 1.000000000 << A[i * 5]<< "\nAttempt root" << attemptroot << "\nFunc:" << func<< std::endl;
                matprint(A);
                abort();

                return 40.0;
            }
        }
        /**** Crucial function ****/
        double functiontominimize(const std::vector<std::complex<double>>& Mat) {
            /* This is the function to minimize */
            const int N = static_cast<int> (sqrt(Mat.size()));
            double functional = 0.0;
            for(int ii = 0; ii < N; ii++) {
                for(int jj = 0; jj < N; jj++) {
                    functional += real(Mat[ii + jj * N] * Mat[jj + ii * N]);
                }
                functional += - 16.0 / 3.0 * real(Mat[ii * (N + 1)]);
            }
            return functional;
        }
        /*****Calculate the determinant of a complex matrix *****/
        int per_sign(const std::vector<int>& nums)
        {
            /* calculate sign of a permutation*/
            const int N(nums.size());
            std::vector<bool> visited(N, false);
            int rst(1), next(0), L(0);
            for (int k(1); k <= N; ++k) {
                if (!visited[k - 1]) {
                    next = k;
                    L = 0;
                    while (!visited[next - 1]) {
                        ++L;
                        visited[next - 1] = true;
                        next = nums[next - 1];
                    }
                    if (L % 2 == 0) {
                        rst *= -1;
                    }
                }
            }
            return rst;
        }

        double det(const std::vector<double>& Mat) {
            /* calculate determinan of the matrix
             *  *          *
             *   *                   *  param Mat: N * N matrix
             *    *                            */
            const int N = static_cast<int>(sqrt(Mat.size()));
            std::vector<double> A(Mat);
            std::vector<int> ipiv(N);
            int info(0);

            dgetrf(&N, &N, &A[0], &N, &ipiv[0], &info);

            std::vector<int> per(N);
            for (int i(1); i <= N; ++i) {
                per[i - 1] = i;
            }

            for (int i(1); i <= N; ++i) {
                std::swap(per[i - 1], per[ipiv[i - 1] - 1]);
            }

            double rst( per_sign(per) );
            for (int i(0); i < N; ++i) {
                rst *= A[i + i * N];
            }
            return rst;
        }
        std::complex<double> cdet(const std::vector<std::complex<double> >& Mat) {
            /* calculate determinan of the matrix
             *  *  *          *
             *   *   *                   *  param Mat: N * N matrix
             *    *    *                            */
            const int N = static_cast<int>(sqrt(Mat.size()));
            std::vector<std::complex<double> > A(Mat);
            std::vector<int> ipiv(N);
            int info(0);

            zgetrf(&N, &N, &A[0], &N, &ipiv[0], &info);

            std::vector<int> per(N);
            for (int i(1); i <= N; ++i) {
                per[i - 1] = i;
            }

            for (int i(1); i <= N; ++i) {
                std::swap(per[i - 1], per[ipiv[i - 1] - 1]);
            }

            std::complex<double> rst( per_sign(per) );
            for (int i(0); i < N; ++i) {
                rst *= A[i + i * N];
            }
            return rst;
        }

        /***** multiply two complex matrices (std::vector<std::complex<double>>) *****/
        std::vector<std::complex<double> > zmatmat(
                    const std::vector<std::complex<double> >& Mat1,
                    const std::vector<std::complex<double> >& Mat2,
                    int K, const char* opa, const char* opb){
            /*  matrix-matrix multiplication for complex matrices,
             *
             *  param Mat1: M * K matrix
             *  param Mat2: K * N matrix
             *
             *  return a M * N matrix
             */
            std::complex<double> alpha(1.0);
            std::complex<double> beta(0.0);
            int M = static_cast<int>(Mat1.size() / K);
            int N = static_cast<int>(Mat2.size() / K);
            std::vector<std::complex<double> > rst(M * N);

            zgemm(opa, opb, &M, &N, &K, &alpha,
                        &Mat1[0], &M, &Mat2[0], &K,
                        &beta, &rst[0], &M);

            return rst;
        }
        /** Multiply two double matrices (std::vector<double>) **/
        std::vector<double> dmatmat(
                    const std::vector<double>& Mat1,
                    const std::vector<double>& Mat2,
                    int K, const char* opa, const char* opb){
            /*  matrix-matrix multiplication for double matrices,
             *
             *  param Mat1: M * K matrix
             *  param Mat2: K * N matrix
             *
             *  return a M * N matrix
             */
            double alpha(1.0);
            double beta(0.0);
            int M = static_cast<int>(Mat1.size() / K);
            int N = static_cast<int>(Mat2.size() / K);
            std::vector<double> rst(M * N);

            dgemm(opa, opb, &M, &N, &K, &alpha,
                        &Mat1[0], &M, &Mat2[0], &K,
                        &beta, &rst[0], &M);

            return rst;
        }

        void hdiag(const std::vector<double>& Mat,
                    std::vector<double>& eva, std::vector<double>& evt){
            /*  calc eva & evt for symmetric double Mat
             *
             *  param Mat:  N * N matrix
             *
             *  output eva, evt: eigenvalues & eigenvectors
             */

            int N = static_cast<int>(sqrt(Mat.size()));
            std::vector<double> evt_(Mat);

            int lwork(-1), liwork(-1), info(-1);
            std::vector<double> work(1);
            std::vector<int> iwork(1);
            eva.resize(N);

            dsyevd("V", "L", &N, &evt_[0], &N, &eva[0],
                        &work[0], &lwork, &iwork[0], &liwork, &info);

            lwork = (int)work[0];
            liwork = iwork[0];
            work.resize(lwork);
            iwork.resize(liwork);

            dsyevd("V", "L", &N, &evt_[0], &N, &eva[0],
                        &work[0], &lwork, &iwork[0], &liwork, &info);

            evt = evt_;
        }

        void zhediag(const std::vector<std::complex<double> >& Mat,
                    std::vector<double>& eva, std::vector<std::complex<double> >& evt){
            /*  calc eva & evt for Hermitian std::complex<double> Mat
             *
             *  param Mat:  N * N matrix
             *
             *  output eva, evt: eigenvalues & eigenvectors
             */

            int N = static_cast<int>(sqrt(Mat.size()));
            std::vector<std::complex<double> > evt_(Mat);

            int lwork(-1), lrwork(-1), liwork(-1), info(-1);
            std::vector<std::complex<double> > work(1);
            std::vector<double> rwork(1);
            std::vector<int> iwork(1);
            eva.resize(N);

            zheevd("V", "L", &N, &evt_[0], &N, &eva[0],
                        &work[0], &lwork, &rwork[0], &lrwork, &iwork[0], &liwork, &info);

            lwork = (int)real(work[0]);
            lrwork = (int)rwork[0];
            liwork = iwork[0];
            work.resize(lwork);
            rwork.resize(lrwork);
            iwork.resize(liwork);


            zheevd("V", "L", &N, &evt_[0], &N, &eva[0],
                        &work[0], &lwork, &rwork[0], &lrwork, &iwork[0], &liwork, &info);
            evt = evt_;
        }
        /***** Calculate the logarithm of a complex matrix *****/
        std::vector<std::complex<double> > mlog(const std::vector<std::complex<double> >& Mat, int N, bool& flagwrong){
            /*
             * calculate log(Mat)
             *
             * - Mat must be orthogonal
             * - All diagonal elements M[k,k] must be positive
             */

            // check orthogonality
            std::vector<std::complex<double>> MM(zmatmat(Mat, Mat, N, "C", "N"));
            for (int i(0); i < N; ++i) {
                for (int j(0); j < N; ++j) {
                    if ((i == j and std::abs(MM[i + j * N] - 1.0) > 0.00001) or
                                (i != j and std::abs(MM[i + j * N]) > 0.01)){
                        std::cout <<"matrixop::mlog : input matrix is not orthogonal!\n" << "i: " << i << "\nj: "<< j << "\nelement: " <<  MM[i + j * N] << std::endl;
                        //std::cout << "MM adjusted" << endl;
                        matprint(Mat);
                        flagwrong = true;
                        //abort();
                    }
                }
            }
            // convert from double to complex
            std::vector<std::complex<double> > A;
            A.resize(N * N);
            for (int i(0); i < A.size(); ++i) {
                A[i] = (std::complex<double>)Mat[i];
            }
            // Schur decomposition A = Q * T * Q'
            std::vector<std::complex<double> > w;
            std::vector<std::complex<double> > work;
            std::vector<std::complex<double> > vs;
            std::vector<double> rwork;
            int bwork[1];
            int lwork, info, SDIM(0);
            w.resize(N);
            vs.resize(N * N);
            rwork.resize(N);
            lwork = -1;
            work.resize(1);
            zgees("V", "N", NULL, &N, &A[0], &N, &SDIM,
                        &w[0], &vs[0], &N, &work[0], &lwork, &rwork[0], &bwork[0], &info);

            lwork = static_cast<int>(work[0].real());
            work.resize(lwork);
            zgees("V", "N", NULL, &N, &A[0], &N, &SDIM,
                        &w[0], &vs[0], &N, &work[0], &lwork, &rwork[0], &bwork[0], &info);
            // construct diagonal unitary matrix logD (stored in A)
            for (int i(0); i < N; ++i) {
                for (int j(0); j < N; ++j) {
                    if (i == j) A[i + j * N] = std::log(A[i + j * N] / std::abs(A[i + j * N]));
                    else A[i + j * N] = 0.0;
                }
            }
            // get Q * log(D) * Q'
            A = zmatmat(vs, zmatmat(A, vs, N, "N", "C"), N, "N", "N");
            /*//Check if A is pure real?
              bool flagimag = false;
              for (int i(0); i < N; ++i) {
              for (int j(0); j < N; ++j) {
              if (std::imag(A[i + j * N]) > 0.0001){
              flagimag = true;
              std::cout <<"matrixop::mlog: Dropping useful value in imaginary part!\n" << "i: " << i << "\nj: "<< j << "\nelement: " <<  A[i + j * N] << std::endl;
              }
              }
              }*/
            //Check if A is pure real?
            bool flagimag = false;
            for (int i(0); i < N; ++i) {
                if (abs(std::imag(A[i + i * N])) > 3.14){
                    flagimag = true;
                    //std::cout <<"matrixop::mlog: Dropping useful value in imaginary part!\n" << "i: " << i << "\nj: "<< j << "\nelement: " <<  A[i + j * N] << std::endl;
                }
            }
            if (flagimag == true) {
                std::cout << "Matrixop::mlog: I have big imaginary part along diagonal elements" << std::endl;
                matprint(Mat);
                flagwrong = true;
            }

            /*// convert from complex to double
              std::vector<double>rst(N * N);
              for (int i(0); i < rst.size(); ++i) {
              rst[i] = A[i].real();
              }*/
            return A;
        }


        std::vector<std::complex<double>> cal_Tmat(const std::vector<std::complex<double>>& curevt, std::vector<std::complex<double>>& nextevt, int NNN, double timestep) {
            /**
             * get the correct phases for next evt, and return correpsonding T matrix
             */
            const std::complex<double> complexnumberi(0.0, 1.0);
            std::vector<std::complex<double> > UUU, tempnextevt, logeigU;
            tempnextevt.resize(NNN * NNN, 0.0);
            std::vector<double> nexteva, absUUU(NNN * NNN, 0.0);
            int ii, jj, numtc, maximumvalue_rowindex;
            double maximumvalue_column = 0.0;
            bool flagtrivialcrossing = false;

            ///////Give me two eigenvectors here!
            //UUU = <\psi(t_0)|\psi(t_0+dt_c)>
            UUU = zmatmat(curevt, nextevt, NNN, "C", "N");
            numtc = 0;
            //absUUU = abs(UUU)
            for(ii = 0; ii < NNN * NNN; ii++) {
                absUUU[ii] = abs(UUU[ii]);
            }
            //Make sure Max(UUU(jj, iii)) element of column ii REAL and positive by changing the phase [exp(i\theta)] of eigenvectors at t_0 + dt_c
            for(ii = 0; ii < NNN; ii++) {
                //tempit = max_element(absUUU.begin() + ii * NNN, absUUU.begin() + ii * NNN + NNN - 1);
                maximumvalue_column = absUUU[ii * NNN];
                maximumvalue_rowindex = 0;
                for(jj = 0; jj < NNN; jj++) {
                    if(absUUU[ii * NNN + jj] > maximumvalue_column) {
                        maximumvalue_rowindex = jj;
                        maximumvalue_column = absUUU[ii * NNN + jj];
                    }
                }
                for(jj = 0; jj < NNN; jj++){
                    tempnextevt[ii * NNN + jj] = nextevt[ii * NNN + jj] * absUUU[maximumvalue_rowindex + ii * NNN] / UUU[maximumvalue_rowindex + ii * NNN];
                }
            }
            UUU = zmatmat(curevt, tempnextevt, NNN, "C", "N");
            std::vector<std::complex<double>> tempUUU(NNN * NNN, 0.0);
            std::complex<double> determinantofU = cdet(UUU); //////////I need a determinant function here
            if(abs(abs(determinantofU) - 1.0) > 0.0001){
                std::cout << "U is not unitary?!\ndet(U) = " << abs(determinantofU)<< std::endl;
                //matprint(UUU);
                abort();
            }
            //change first column to make det(UUU) = +1!//
            for(ii = 0; ii < NNN; ii++) {
                tempnextevt[ii] = tempnextevt[ii] / determinantofU;
                UUU[ii + 0 * NNN] = UUU[ii + 0 * NNN] / determinantofU;
            }
            /* Minimize f(theta) = a cos(theta) + b sin(theta), find the two coefficients {a, b} as denoted {coscoef, sincoef} */
            double coscoef, sincoef, theta, cos2coef, sin2coef;
            double minimumvalue = functiontominimize(UUU);
            double tempminvalue;
            tempUUU = UUU;
            while(1) {

                for(int inumtc = 0; inumtc < NNN; inumtc++) {
                    for(int jnumtc = inumtc + 1; jnumtc < NNN; jnumtc++) {
                        coscoef = 0.0;
                        sincoef = 0.0;
                        cos2coef = 0.0;
                        sin2coef = 0.0;
                        for(ii = 0; ii < NNN; ii++) {
                            coscoef += real(tempUUU[ii + inumtc * NNN] * tempUUU[ii * NNN + inumtc] + tempUUU[ii + jnumtc * NNN] * tempUUU[ii * NNN + jnumtc]);
                            sincoef += -imag(tempUUU[ii + inumtc * NNN] * tempUUU[ii * NNN + inumtc] - tempUUU[ii + jnumtc * NNN] * tempUUU[ii * NNN + jnumtc]);
                        }
                        coscoef -= real(tempUUU[jnumtc + inumtc * NNN] * tempUUU[inumtc + jnumtc * NNN]) * 2 + real(tempUUU[inumtc * (NNN + 1)] + tempUUU[jnumtc * (NNN + 1)]) * 8.0 / 3.0 + real(tempUUU[inumtc * (NNN + 1)] * tempUUU[inumtc * (NNN + 1)] + tempUUU[jnumtc * (NNN + 1)] * tempUUU[jnumtc * (NNN + 1)]);
                        sincoef -= 8.0 / 3.0 * imag(tempUUU[jnumtc * (NNN + 1)] - tempUUU[inumtc * (NNN + 1)]) + imag(tempUUU[jnumtc * (NNN + 1)] * tempUUU[jnumtc * (NNN + 1)] - tempUUU[inumtc * (NNN + 1)] * tempUUU[inumtc * (NNN + 1)]);
                        cos2coef = real(tempUUU[inumtc * (NNN + 1)] * tempUUU[inumtc * (NNN + 1)] + tempUUU[jnumtc * (NNN + 1)] * tempUUU[jnumtc * (NNN + 1)]) * 0.50;
                        sin2coef = imag(tempUUU[jnumtc * (NNN + 1)] * tempUUU[jnumtc * (NNN + 1)] - tempUUU[inumtc * (NNN + 1)] * tempUUU[inumtc * (NNN + 1)]) * 0.50;
                        theta = rotation_angle(cos2coef, sin2coef, coscoef, sincoef);
                        /*if(theta == 40.0) {
                          matprint(UUU);
                          matprint(tempUUU);

                          abort();
                          }*/

                        for(ii = 0; ii < NNN; ii++) {
                            tempUUU[ii + inumtc * NNN] = exp(complexnumberi * theta) * tempUUU[ii + inumtc * NNN];
                            tempUUU[ii + jnumtc * NNN] = exp(-complexnumberi * theta) * tempUUU[ii + jnumtc * NNN];
                            tempnextevt[ii + inumtc * NNN] = exp(complexnumberi * theta) * tempnextevt[ii + inumtc * NNN];
                            tempnextevt[ii + jnumtc * NNN] = exp(-complexnumberi * theta) * tempnextevt[ii + jnumtc * NNN];
                        }
                    }
                }
                tempminvalue = functiontominimize(tempUUU);
                //std::cout << "MIN: "<< minimumvalue << " -> " << tempminvalue << std::endl;
                if(abs(minimumvalue - tempminvalue) < 0.00001){
                    break;
                }
                else{
                    minimumvalue = tempminvalue;
                }
            }	
            UUU = zmatmat(curevt, tempnextevt, NNN, "C", "N");

            //Save eigenvectors at t_0+dt_c with correct signs
            nextevt = tempnextevt;
            //nextvals.usign = zmatmat(currentvals.usign, UUU, NNN, "N", "N");
            //T = logm(U)/dt_c
            std::fill(logeigU.begin(), logeigU.end(), 0.0);
            //logeigU = mlog(UUU, NNN, nextvals.flagcontinue);
            bool wrongflag = false;
            logeigU = mlog(UUU, NNN, wrongflag);
            //nextvals.TTT.resize(NNN * NNN, 0.0);
            std::vector<std::complex<double>> TTT(NNN * NNN, 0.0);
            for (ii = 0; ii < NNN * NNN; ii++) {
                TTT[ii] = logeigU[ii] / timestep;
            }
            //return nextvals;
            return TTT;
        }


    } // namespace zeyu
} // namespace mqc
