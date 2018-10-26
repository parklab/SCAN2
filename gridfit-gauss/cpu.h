void init_cpu(int n);
void finalize_cpu();

double *alloc_host_dvec(int len, const char *message);

double
laplace_approx_cpu(double a, double b, int *x, int *y, int *d, int intlen, int maxit);

double
laplace_approx_chunk_cpu(int chunk,
    double a, double b, double c, double paramd, int *x, int *y, int *d, int len, int maxit, char verbose);
