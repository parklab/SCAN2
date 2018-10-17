/*
 * Important note: using banding allows us to make an arbitrarily
 * good approximation to factorizing the full chromosome as the
 * number of bands approaches (n-1)/2.  HOWEVER, it is not true
 * in general that restricting a positive definite matrix to its
 * bands produces another positive definite matrix.  A non PSD
 * matrix cannot be a covariance matrix.  In many scenarios for
 * this application, the outer bands of the matrix are very close
 * to 0, so replacing them with 0 (banding) is less likely to
 * produce a non PSD matrix, but this is not guaranteed.  Further,
 * the outer bands for large b will NOT tend to 0, so in those
 * cases the probability of non-PSD is even more dangerous.
 * That said, in the cases where a banded approximation succeeds,
 * it can be compared against the chunk approximation as a sanity
 * check.
 */
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <stdlib.h>

#include "cpu.h"  /* for alloc */

#define USE_MKL 1
#undef SYMMETRIC

#ifdef USE_MKL
    #include <mkl.h>
    #include <mkl_lapacke.h>
#else
    #include <lapacke.h>
    #include <cblas.h>
#endif

#define IDX(_i, _j, _ld)   ((_i) + (_j)*(_ld))

/* define SYMMETRIC to use banding through Cholesky factorization */
/* the non-symmetric code was a trial of LU factorization to see if
 * the non-PSD results from Cholesky factorization could have been
 * due to numerical issues rather than real non-PSD.  But this did
 * not appear to be the case.
 */
#define SYMMETRIC


#ifdef SYMMETRIC
    #define IDXB(_i, _j, _k)  IDX((_k) + (_i) - (_j), (_j), (_k)+1)
#else
    /* special bit: because LU factorization needs more space, there
     * must  be _k blank rows on top of the stored matrix.  this is
     * why 2*(_k). */
    #define IDXB(_i, _j, _k)  IDX((_k) + (_i) - (_j), (_j), 2*(_k)+1)
    #define IDXLU(_i, _j, _k)  IDX(2*(_k) + (_i) - (_j), (_j), 3*(_k)+1)
#endif

lapack_int *piv;
double *U, *V, *K, *A, *B, *sqrtW;

/* CLOBBERS variables i, k */
#define PRINT_MATRIXB(_A, _n, _band) \
    for (k = 0; k < (_band) + 1; ++k) { \
        printf("[ "); \
        for (i = 0; i < _n; ++i) \
            printf("%0.5f ", _A[IDX(k, i, (_band)+1)]); \
        printf("]\n"); \
    } \

#define PRINT_MATRIXLU(_A, _n, _band) \
    for (k = 0; k < 3*(_band) + 1; ++k) { \
        printf("[ "); \
        for (i = 0; i < _n; ++i) \
            printf("%0.5f ", _A[IDX(k, i, 3*(_band)+1)]); \
        printf("]\n"); \
    } \



/* Allocating and deallocating these for each grid evaluation
 * can waste alot of time.  Instead, initialize all the necessary
 * memory just once.
 * */
void
init_cpu_band(int band, int len)
{
    int i, j;

    U = alloc_host_dvec(len, "init_memory_cpu_band: U");
    V = alloc_host_dvec(len, "init_memory_cpu_band: V");
    B = alloc_host_dvec(len, "init_memory_cpu_band: B");
    sqrtW = alloc_host_dvec(len, "init_memory_cpu_band: sqrtW");

    piv = malloc(len * sizeof(lapack_int));
    if (piv == NULL) {
        printf("init_cpu_band: couldn't allocate pivot array\n");
        exit(1);
    }

#ifdef SYMMETRIC
    K = alloc_host_dvec((band+1)*len, "init_memory_cpu_band: K");
    A = alloc_host_dvec((band+1)*len, "init_memory_cpu_band: A");
#else
    /* banded storage: for nonsymmetric banded matrices, can store
     * in a (ku + kl + 1) x n matrix, where ku is the number of
     * superdiagonals and kl the number of subdiagonals.  since
     * our banded matrices are symmetric, we need only (k + 1) x n.
     */
    /* see IDXB for why there's an additional k rows added on top,
     * hence 3*band */
    K = alloc_host_dvec((3*band+1)*len, "init_memory_cpu_band: K");
    A = alloc_host_dvec((3*band+1)*len, "init_memory_cpu_band: A");
#endif
}


void
finalize_cpu_band()
{
    free(U);
    free(V);
    free(B);
    free(K);
    free(A);
    free(sqrtW);
}



/*
 * Compute the Laplace approximation on the CPU
 * All memory has to be allocated by the caller.
 */
double
laplace_approx_band_cpu(int band, double a, double b, int *x, int *y, int *d, int len, int maxit, char verbose)
{
    int i, k;
    int iter = 0;
    double logp = -1.0, lastlogp = -INFINITY;
    double sqrt_eps = sqrt(2.2e-16);
    lapack_int info, m, n, lda, ldb, nrhs;

    /* Initialization */
    for (i = 0; i < len; ++i)
        B[i] = 0;

    /* K does not change during iteration */
    /* for banded matrix, makes more sense to iterate by band index */
    /* b=0: set the diagonal, b=1, set the first subdiagonal, ... */
    /* remember, working with lower triangle of sym banded mat */
    for (k = 0; k <= band; ++k) {
        for (i = k; i < len; ++i) {
            K[IDXB(i-k, i, band)] =
                exp(a - (double)abs(x[i] - x[i-k]) / (b*b));
#ifndef SYMMETRIC
            K[IDXB(i-k, i, band)] =
                exp(a - (double)abs(x[i] - x[i-k]) / (b*b));
#endif
        }
    }

if (verbose) {
    printf("matrix K initialized\n");
#ifdef SYMMETRIC
    PRINT_MATRIXB(K, len, band);
#endif
}

    do {
        iter += 1;

        /* O(n): sqrtW <- sqrt(D * exp(B) / (1 + exp(B))^2) */
        /* compute and upload */
        for (i = 0; i < len; ++i)
            sqrtW[i] = sqrt(d[i] * exp(B[i])) / (1.0 + exp(B[i]));

        /* O(n^2): I + outer(sqrtW, sqrtW) * K */
        for (k = 0; k <= band; ++k) {
            for (i = k; i < len; ++i) {
#ifdef SYMMETRIC
                A[IDXB(i-k, i, band)] =
                    sqrtW[i] * sqrtW[i-k] * K[IDXB(i-k, i, band)];
#else
                A[IDXLU(i, i-k, band)] =
                    sqrtW[i] * sqrtW[i-k] * K[IDXB(i, i-k, band)];
                A[IDXLU(i-k, i, band)] =
                    sqrtW[i] * sqrtW[i-k] * K[IDXB(i-k, i, band)];
#endif
            }
        }
        /* add I */
        for (i = 0; i < len; ++i)
#ifdef SYMMETRIC
            A[IDXB(i, i, band)] += 1.0;
#else
            A[IDXLU(i, i, band)] += 1.0;
#endif

if (verbose) {
    printf("matrix A initialized\n");
#ifdef SYMMETRIC
    PRINT_MATRIXB(A, len, band);
#endif
}
        
        /* O(n): W*B + (Y - D * exp(B) / (1 + exp(B))) */
        for (i = 0; i < len; ++i)
            V[i] = sqrtW[i]*sqrtW[i] * B[i] + y[i] - d[i] / (1.0 + exp(-B[i]));

        /* O(n^3): A <- t(chol(A)) */
#ifdef SYMMETRIC
        info = LAPACKE_dpbtrf(LAPACK_COL_MAJOR, 'U', len, band, A, band+1);
#else
        /* see IDXB for why 3*band+1 */
        info = LAPACKE_dgbtrf(LAPACK_COL_MAJOR, len, len, band, band, A, 3*band+1, piv);
#endif
        if (info != 0) {
#ifdef SYMMETRIC
            printf("laplace_approx(iter=%d): CPU Cholesky factorization failed: %d\n",
                iter, info);
#else
            printf("laplace_approx(iter=%d): CPU LU factorization failed: %d\n",
                iter, info);
#endif
            exit(1);
        }

        /* want to compute: A^{-1} sqrtW K V.
         * Set x = A^{-1} sqrtW K V, then it is clear that computing x
         * is equivalent to solving Ax = sqrtW K V, where A = LL^T was
         * the matrix just factorized.  Next code is to store sqrtW K V
         * into the scratch vector U.
         */
        /* U <- K V */
#ifdef SYMMETRIC
        cblas_dsbmv(CblasColMajor, CblasUpper, len, band, 1.0, K, band+1, V, 1, 0.0, U, 1);
#else
        cblas_dgbmv(CblasColMajor, CblasNoTrans,
            len, len, band, band, 1.0, K, 2*band+1, V, 1, 0.0, U, 1);
#endif

        /* U <- sqrtW U.  Note: sqrtW is a diagonal matrix */
        for (i = 0; i < len; ++i)
            U[i] = sqrtW[i] * U[i];

        /* Solve for x: Ax = U, gives x = LL^T sqrtW K V.  Store U <- x */
#ifdef SYMMETRIC
        info = LAPACKE_dpbtrs(LAPACK_COL_MAJOR, 'U', len, band, 1, A, band+1, U, len);
#else
        info = LAPACKE_dgbtrs(LAPACK_COL_MAJOR, 'N',
            len, band, band, 1, A, 3*band+1, piv, U, len);
#endif
        if (info != 0) {
#ifdef SYMMETRIC
            printf("laplace_approx: CPU Cholesky solve failed: %d\n", info);
#else
            printf("laplace_approx: CPU LU solve failed: %d\n", info);
#endif
            exit(1);
        }

        /* U <- sqrtW U and V <- V - U */
        for (i = 0; i < len; ++i) {
            U[i] = sqrtW[i] * U[i];
            V[i] = V[i] - U[i];
        }

        /* B <- K V */
#ifdef SYMMETRIC
        cblas_dsbmv(CblasColMajor, CblasUpper, len, band, 1.0, K, band+1, V, 1, 0.0, B, 1);
#else
        cblas_dgbmv(CblasColMajor, CblasNoTrans, len, len, band, band, 1.0, K, 2*band+1, V, 1, 0.0, B, 1);
#endif
        
if (verbose) {
printf("B: ");
for(i = 0; i < len; ++i) printf(" %0.4f", B[i]);
printf("\n");
printf("V: ");
for(i = 0; i < len; ++i) printf(" %0.4f", V[i]);
printf("\n");
}

#ifndef SYMMETRIC
        /* The diagonal stored in A is the diagonal of U.  The decomposition
         * is A = PLU, where P is a permutation matrix and L has unit diag.
         * det A = det P det L det U = det P det U = (-1)^k prod(diag(U)).
         * where k is the number of row exchanges.
         * To ensure the probability is valid, make sure the sign is +1.
         *   piv[i] = k means row i was exchanged with k, but note that
         *   k is 1-indexed while i is 0-indexed.
         */
        int sign = 1;
        double z;
        for (i = 0; i < len; ++i) {
            z = A[IDXLU(i,i,band)];
            sign = sign * (z > 0 ? 1 : (z < 0 ? -1 : 0));
            if (verbose) printf("sign=%d\n", sign);
            sign = sign * (i == piv[i] - 1 ? 1 : -1);
            if (verbose) printf("sign=%d\n", sign);
        }
        if (sign != 1) {
            printf("ERROR: non-positive determinant: %d\n", sign);
            if (verbose) {
                printf("diag(U):");
                for (k = 0; k < len; ++k) printf(" %0.4f", A[IDXLU(k,k,band)]);
                printf("\n");
                printf("piv:");
                for (k = 0; k < len; ++k) printf(" %d", piv[k]);
                printf("\n");
            }
            exit(1);
        }
#endif
        
        /* Compute the approximate log likelihood for this iteration */
        lastlogp = logp;
        logp = 0;
        for (i = 0; i < len; ++i) {
            logp += -V[i]*B[i]/2 + (B[i]*y[i] - 
                d[i]*log(1.0 + exp(B[i]))) -
#ifdef SYMMETRIC
                log(fabs(A[IDXB(i, i, band)]));
#else
                /* Using fabs() is justified as long as sign=1 above */
                log(fabs(A[IDXLU(i, i, band)]));
#endif
        }
        if (verbose) printf("    iter %2d: logp=%0.10f\n", iter, logp);
    } while (fabs((logp - lastlogp) / logp) > sqrt_eps && iter < maxit);


    if (verbose) printf("%d iterations\n", iter);
    if (iter > maxit) {
        printf("laplace_approx: reached maxit=%d iterations without converging\n",
               maxit);
        exit(1);
    }

    return logp;
}
