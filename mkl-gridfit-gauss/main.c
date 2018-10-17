#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <stdlib.h>

#include "cpu.h"
#include "cpu_band.h"

/* These relate to the parallelization grid */
/* Each job calculates SUBEDGE x SUBEDGE points, and jobs are
 * arranged into an interwoven SUPEREDGE x SUPEREDGE grid. */
#define SUB_EDGE     10
#define SUPER_N      1
#define SUPER_EDGE   SUB_EDGE * SUPER_N

double
runif(double min, double max)
{
    double u = drand48();
    return (u * (max-min) + min);
}

/*
 * Approximate logp for several (a,b) parameter sets.
 * Write the grid points and approx. logps in binary format
 * for consumption by R at a later time.
 */
void
laplace_approx_grid(FILE *outfile, int Nsamples, int chunk,
    int *x, int *y, int *d, int len,
    double amin, double amax, double bmin, double bmax,
    double cmin, double cmax, double dmin, double dmax, int strategy)
{
    int i;
    double a, b, c, paramd;
    double logp;

    printf("Computing Laplace approx grid, strategy=%s, %s=%d\n",
        strategy == 1 ? "CPU banded" : "CPU chunked",
        strategy == 1 ? "bands" : "chunk size", chunk);

    /* calculate and write the logp matrix in column major format */
    printf("Scheduled %d points for computation\n", Nsamples);
    fwrite(&Nsamples, sizeof(int), 1, outfile);
    for (i = 0; i < Nsamples; i += 1) {
        a = runif(amin, amax);
        b = pow(10.0, runif(bmin, bmax));
        c = runif(cmin, cmax);
        paramd = pow(10.0, runif(dmin, dmax));
        logp = laplace_approx_chunk_cpu(chunk, a, b, c, paramd, x, y, d, len, 50, 0);
        printf("a=%0.6f, b=%0.2f, c=%0.2f, d=%0.2f, logp=%0.3f\n",
            a, b, c, paramd, logp);
        fwrite(&a, sizeof(double), 1, outfile);
        fwrite(&b, sizeof(double), 1, outfile);
        fwrite(&c, sizeof(double), 1, outfile);
        fwrite(&paramd, sizeof(double), 1, outfile);
        fwrite(&logp, sizeof(double), 1, outfile);
    }
}


/*
 * input binary file format:
 * N <integer> number of rows
 * pos N*<integer> positions
 * hap1 N*<integer> reads supporting allele 1
 * dp N*<integer> total read depth at locus
 */
int
read_bin(const char *filename, int **x, int **y, int **d)
{
    int *xbuf;
    int *ybuf;
    int *dbuf;
    int N;
    FILE *f;

    printf("Opening binary input data file %s\n", filename);
    f = fopen(filename, "rb");
    if (f == NULL) {
        printf("fopen: %s", strerror(errno));
        exit(1);
    }

    fread((void *)&N, sizeof(int), 1, f);
    printf("   binary data contains %d rows\n", N);
    xbuf = (int *)malloc(N*sizeof(int));
    ybuf = (int *)malloc(N*sizeof(int));
    dbuf = (int *)malloc(N*sizeof(int));
    if (xbuf == NULL || ybuf == NULL || dbuf == NULL) {
        printf("read_bin: couldn't allocate memory to read input file\n");
        exit(1);
    }
    fread((void *)xbuf, sizeof(int), N, f);
    fread((void *)ybuf, sizeof(int), N, f);
    fread((void *)dbuf, sizeof(int), N, f);
    fclose(f);

    *x = xbuf;
    *y = ybuf;
    *d = dbuf;

    return N;
}



int
main(int argc, char **argv)
{
    /* some hard coded input data */
    //int x[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    //int y[] = { 4, 3, 0, 0, 4, 6, 12, 6, 13, 9 };
    //int d[] = { 46, 14, 7, 12, 23, 31, 37, 11, 21, 16 };
    int i;
    int *x, *y, *d, N;
    int chunksize = 250; //100;  // for strategy=1, number of subdiagonals
    int strategy = 0;

    if (argc != 13) {
        printf("usage: %s bindata output.bin amin amax bmin bmax cmin cmax dmin dmax Nsamples seed\n", argv[0]);
        printf("  output - write binary format results to this file.\n");
        printf("  seed   - random seed.\n");
        exit(1);
    }
    const char *filename = argv[1];
    const char *output = argv[2];
    double amin = strtod(argv[3], NULL);
    double amax = strtod(argv[4], NULL);
    double bmin = strtod(argv[5], NULL);
    double bmax = strtod(argv[6], NULL);
    double cmin = strtod(argv[7], NULL);
    double cmax = strtod(argv[8], NULL);
    double dmin = strtod(argv[9], NULL);
    double dmax = strtod(argv[10], NULL);
    int Nsamples = strtol(argv[11], NULL, 10);
    int seed = strtol(argv[12], NULL, 10);
    FILE *outfile;

    srand48(seed);

    printf("   bindata = %s\n", filename);
    printf("    output = %s\n", output);
    printf("      amin = %0.6f\n", amin);
    printf("      amax = %0.6f\n", amax);
    printf("      bmin = %0.6f\n", bmin);
    printf("      bmax = %0.6f\n", bmax);
    printf("      cmin = %0.6f\n", cmin);
    printf("      cmax = %0.6f\n", cmax);
    printf("      dmin = %0.6f\n", dmin);
    printf("      dmax = %0.6f\n", dmax);
    printf("  Nsamples = %d\n", Nsamples);
    printf("      seed = %d\n", seed);

    outfile = fopen(output, "wb");
    if (outfile == NULL) {
        printf("couldn't open output file %s\n", output);
        exit(1);
    }

    N = read_bin(filename, &x, &y, &d);
    printf("read x[%d] = [ %d, %d, %d, %d, %d, ... ]\n",
           N, x[0], x[1], x[2], x[3], x[4]);
    printf("read y[%d] = [ %d, %d, %d, %d, %d ... ]\n",
        N, y[0], y[1], y[2], y[3], y[4]);
    printf("read d[%d] = [ %d, %d, %d, %d, %d ... ]\n",
        N, d[0], d[1], d[2], d[3], d[4]);

    switch (strategy) {
    case 0:
        init_cpu(chunksize);
        laplace_approx_grid(outfile, Nsamples, chunksize,
            x, y, d, N, amin, amax, bmin, bmax, cmin, cmax, dmin, dmax, strategy);
        finalize_cpu();
        break;
    case 1:
        init_cpu_band(chunksize, N);
        printf("*** WARNING: approximating a full matrix by a banded matrix\n"
               "*** does not generally produce a positive semidefinite matrix,\n"
               "*** which can lead to a nonpositive determinant and therefore\n"
               "*** a NaN log probability.\n"
               "*** For large values of b, non-PSD banded matrices are frequent!\n"
               "*** USE AT YOUR OWN RISK.\n");
        //laplace_approx_grid(chunksize, x, y, d, 10, strategy);
        //laplace_approx_grid(outfile, gridoffset, chunksize,
            //x, y, d, N, amin, amax, bmin, bmax, c, paramd, strategy);
        finalize_cpu_band();
        break;
    default:
        printf("invalid strategy: %d\n", strategy);
        exit(1);
    }

    fclose(outfile);

    return 0;
}
