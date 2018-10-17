void init_cpu_band(int band, int len);
void finalize_cpu_band();

double
laplace_approx_band_cpu(int band,
    double a, double b, int *x, int *y, int *d, int len, int maxit, char verbose);
