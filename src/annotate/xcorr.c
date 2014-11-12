#include <annotate/xcorr.h>

/*
 * gnoise
 *
 * @see src/include/annotate/xcorr.h
 */
double gnoise(const double signal[], const double mean, const double variance)
{
  double x = 0;
  int i;

  for (i = 0; i < MAX_GNOISE_N; i++) {
    double u = (double) rand() / (double)RAND_MAX;
    x += u;
  }
	
  // for uniform randoms in [0,1], mu = 0.5 and var = 1/12
  x = x - MAX_GNOISE_N / 2;        // set mean to 0
  x = x * sqrt(12 / MAX_GNOISE_N); // adjust variance to 1

  // modify x in order to have a particular mean and variance
  return(mean + sqrt(variance) * x);
}

/*
 * nxcorr
 *
 * @see src/include/annotate/xcorr.h
 */
int nxcorr(double corr[], const double x[], size_t n_x, const double y[], size_t n_y)
{
  double *ext_x, *ext_y;
  int index;
  int lag;
  double length;
  double rxx = 0, ryy = 0, rnm;

  // Extend shorter signal with gaussian white noise
  if (n_x > n_y) {
    length = n_x;
    double mean = gsl_stats_mean(y, 1, n_y);
    double variance = gsl_stats_variance(y, 1, n_y);
    ext_x = (double*) malloc(n_x * sizeof(double));
    ext_y = (double*) malloc(n_x * sizeof(double));
    for (index = 0; index < n_x; index++) {
      ext_x[index] = x[index];
      if (index < n_y)
        ext_y[index] = y[index];
      else {
        ext_y[index] = gnoise(y, mean, variance) * (-1);
      }
    }
  }
  else if (n_x < n_y) {
    length = n_y;
    double mean = gsl_stats_mean(x, 1, n_x);
    double variance = gsl_stats_variance(x, 1, n_x);
    ext_x = (double*) malloc(n_y * sizeof(double));
    ext_y = (double*) malloc(n_y * sizeof(double));
    for (index = 0; index < n_y; index++) {
      ext_y[index] = y[index];
      if (index < n_x)
        ext_x[index] = x[index];
      else {
        ext_x[index] = gnoise(x, mean, variance) * (-1);
      }
    }
  }
  else {
    length = n_x;
    ext_x = (double*) malloc(n_x * sizeof(double));
    ext_y = (double*) malloc(n_y * sizeof(double));
    for (index = 0; index < n_x; index++) {
      ext_y[index] = y[index];
      ext_x[index] = x[index];
    }
  }

  // Calculate rxx and ryy at lag 0 for normalization
  for (index = 0; index < length; index++) {
    rxx += ext_x[index] * ext_x[index];
    ryy += ext_y[index] * ext_y[index];
  }
  rnm = sqrt(rxx * ryy);

  // For lag 0 : (max(n_x, n_y)-1) calculate x-corr
  index = 0;
  for (lag = 0; lag <= 0; lag++) {//(length - 1); lag++) {
    int i;
    double rxy = 0;
    for (i = 0; i < (length - lag); i++)
      rxy += ext_x[i+lag] * ext_y[i];
    corr[lag] = rxy / rnm;
    if (corr[lag] > corr[index])
      index = lag;
  }

  free(ext_x);
  free(ext_y);

  return index;
}
