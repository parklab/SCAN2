#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <stdlib.h>

// Use Intel's MKL library for optimized BLAS and LAPACK
//#define USE_MKL 1
#undef USE_MKL
#define PACKED 1
//#undef PACKED

#ifdef USE_MKL
    #include <mkl.h>
    #include <mkl_lapacke.h>
#else
    #include <lapacke.h>
    #include <cblas.h>
#endif

#ifdef PACKED
    /* for packing symmetric matrices into lower triangular */
    //#define IDX(_i, _j, _ld) ((_i) - 1 + (2*(_ld) - j) * (j - 1)/2)
    // this is for one based indexing
    //#define IDX(_i, _j, _ld) ((_i) - 1 + (_j)*((_j) - 1)/2)
    /* Upper triangular, 0-based index */
    #define IDX(_i, _j, _ld) ((_i) + (_j)*((_j) + 1)/2)
#else
    /* For matrix indexing */
    #define IDX(_i, _j, _ld) ((_i) + (_j)*(_ld))
#endif


double *U, *V, *K, *A, *B, *sqrtW;

#define PVd(_v) \
printf("[ "); \
for(i = 0; i < 8; ++i) \
    printf("%d ", (_v)[i]); \
printf("]\n");

#define PV(_v) \
printf("[ "); \
for(i = 0; i < 8; ++i) \
    printf("%0.4f ", (_v)[i]); \
printf("]\n");

#define PM(A) \
for(i = 0; i < 6; ++i) { \
printf("[ "); \
for (j = i; j < 6; ++j) { \
    printf("%0.4f ", A[IDX(i,j,len)]); \
} \
printf("]\n"); \
} \


//PRINT_MATRIX(A)
//for(i = 0; i < len; ++i) {
//printf("[ ");
//for (j = i; j < len; ++j) {
    //printf("%0.4f ", A[IDX(i,j,-1)]);
//}
//printf("]\n");
//}


double *
alloc_host_dvec(int len, const char *message)
{
    double *x;

    if (posix_memalign((void **)&x, 32, len * sizeof(double)) != 0) {
        printf("alloc_host_dvec: failed: %s\n", message);
        exit(1);
    }

    return x;
}



/* Allocating and deallocating these for each grid evaluation
 * can waste alot of time.  Instead, initialize all the necessary
 * memory just once.
 * */
void
init_cpu(int n)
{
    U = alloc_host_dvec(n, "init_memory_cpu: U");
    V = alloc_host_dvec(n, "init_memory_cpu: V");
    B = alloc_host_dvec(n, "init_memory_cpu: B");
    sqrtW = alloc_host_dvec(n, "init_memory_cpu: sqrtW");
    K = alloc_host_dvec(n*n, "init_memory_cpu: K");
    A = alloc_host_dvec(n*n, "init_memory_cpu: A");
}


void
finalize_cpu()
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
/* Performance profile is roughly:
 *      10.1% - Computing K before the Newton-Raphson iterations
 *      89.8% - Newton-Raphson iterations
 *          80.0% - MKL-LAPACK routines for Cholesky factors and
 *                  matrix-vector products.
 *           8.8% - Computing A = I + sqrtW K sqrtW
 *           0.8% - Computing logp at the end of iteration
 *          <0.1% - Every other loop
 */
double
laplace_approx_cpu(double a, double b, double c, double paramd, int *x, int *y, int *d, int len, int maxit)
{
    int i, j;
    int iter = 0;
    double k;
    double logp = -1.0, lastlogp = -INFINITY;
    double sqrt_eps = sqrt(2.2e-16);
    lapack_int info, m, n, lda, ldb, nrhs;
    double one = 1.0, zero = 0.0;
    double alpha = 1.0, beta = 0.0;
    double r;


    /* Initialization */
    for (i = 0; i < len; ++i)
        B[i] = 0;

    /* K does not change during iteration */
    /* approx ~10% of program spent here */
    for (i = 0; i < len; ++i) {
        K[IDX(i, i, len)] = exp(a) + exp(c);
        for (j = i + 1; j < len; ++j) {
            r = x[i] - x[j];
            r *= r;
            k = exp(a - r / (b*b))
                + exp(c - r / (paramd*paramd));
            K[IDX(i, j, len)] = k;
#if !defined(PACKED)
            K[IDX(j, i, len)] = k;
#endif
        }
    }

    do {
        iter += 1;

        /* O(n): sqrtW <- sqrt(D * exp(B) / (1 + exp(B))^2) */
        /* compute and upload */
        for (i = 0; i < len; ++i)
            sqrtW[i] = sqrt(d[i] * exp(B[i])) / (1.0 + exp(B[i]));

        /* The CPU version of computing A, O(n^2) time */
        /* O(n^2): diag(n) + outer(sqrtW, sqrtW) * K */
        /* also approx. 10% of runtime */
        for (i = 0; i < len; ++i) {
            for (j = i; j < len; ++j) {
                k = sqrtW[i] * sqrtW[j] * K[IDX(i, j, len)];
                A[IDX(i, j, len)] = k;
#if !defined(PACKED)
                A[IDX(j, i, len)] = k;
#endif
            }
            A[IDX(i, i, len)] += 1.0;   /* Adding I */
        }
        
        /* O(n): W*B + (Y - D * exp(B) / (1 + exp(B))) */
        for (i = 0; i < len; ++i)
            V[i] = sqrtW[i]*sqrtW[i] * B[i] + y[i] - d[i] / (1.0 + exp(-B[i]));

        /* O(n^3): A <- t(chol(A)) */
        m = len;
        lda = m;
#ifdef PACKED
        info = LAPACKE_dpptrf(LAPACK_COL_MAJOR, 'U', m, A);
#else
        info = LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', m, A, lda);
#endif
        if (info != 0) {
            printf("laplace_approx: CPU Cholesky factorization failed: %d\n",
                info);
            exit(1);
        }

        /* want to compute: A^{-1} sqrtW K V.
         * Set x = A^{-1} sqrtW K V, then it is clear that computing x
         * is equivalent to solving Ax = sqrtW K V, where A = LL^T was
         * the matrix just factorized.  Next code is to store sqrtW K V
         * into the scratch vector U.
         */
        /* U <- K V */
#ifdef PACKED
        cblas_dspmv(CblasColMajor, CblasUpper, len, 1.0, K, V, 1, 0.0, U, 1);
#else
        cblas_dsymv(CblasColMajor, CblasLower, len, 1.0, K, len, V, 1, 0.0, U, 1);
#endif

        /* U <- sqrtW U.  Note: sqrtW is a diagonal matrix */
        for (i = 0; i < len; ++i)
            U[i] = sqrtW[i] * U[i];

        /* Solve for x: Ax = U, gives x = LL^T sqrtW K V.  Store U <- x */
        n = len;
        nrhs = 1;
        lda = m;
#ifdef PACKED
        info = LAPACKE_dpptrs(LAPACK_COL_MAJOR, 'U', n, nrhs, A, U, lda);
#else
        info = LAPACKE_dpotrs(LAPACK_COL_MAJOR, 'L', n, nrhs, A, lda, U, lda);
#endif
        if (info != 0) {
            printf("laplace_approx: CPU Cholesky solve failed\n");
            exit(1);
        }

        /* U <- sqrtW U and V <- V - U */
        for (i = 0; i < len; ++i) {
            U[i] = sqrtW[i] * U[i];
            V[i] = V[i] - U[i];
        }

        /* Compute the new B, storing in U for now: B <- K V */
#ifdef PACKED
        cblas_dspmv(CblasColMajor, CblasUpper, len, 1.0, K, V, 1, 0.0, B, 1);
#else
        cblas_dsymv(CblasColMajor, CblasLower, len, 1.0, K, len, V, 1, 0.0, B, 1);
#endif
        
        /* Compute the approximate log likelihood for this iteration */
        lastlogp = logp;
        logp = 0;
        for (i = 0; i < len; ++i) {
            logp += -V[i]*B[i]/2.0 + (B[i]*y[i] - 
                d[i]*log(1.0 + exp(B[i]))) - log(A[IDX(i, i, len)]);
        }
        //printf("    iter %2d: logp=%0.10f\n", iter, logp);
    } while (fabs((logp - lastlogp) / logp) > sqrt_eps && iter < maxit);


    //printf("%d iterations\n", iter);
    if (iter > maxit) {
        printf("laplace_approx: reached maxit=%d iterations without converging\n",
               maxit);
        exit(1);
    }

    return logp;
}



/*
 * Assume x, y and d have length len.  Instead of computing
 * the logp of the full (x,y,d), approximate by assuming each
 * chunk of size 'chunk' is independent.
 */
double
laplace_approx_chunk_cpu(int chunk,
    double a, double b, double c, double paramd, int *x, int *y, int *d, int len, int maxit, char verbose)
{
    int i;
    int *xchunk, *ychunk, *dchunk;
    int thislen;
    double logp = 0, thislogp;

    for (i = 0; i*chunk < len; ++i) {
        xchunk = x + i*chunk;
        ychunk = y + i*chunk;
        dchunk = d + i*chunk;
        thislen = (chunk < len - i*chunk) ? chunk : len - i*chunk;
        thislogp = laplace_approx_cpu(a, b, c, paramd, xchunk, ychunk, dchunk, thislen, maxit);
        logp += thislogp;
        if (verbose)
            printf("chunk %d, N=%d, thislogp=%0.5f, logp=%0.5f\n",
                i, thislen, thislogp, logp);
    }

    return logp;
}
