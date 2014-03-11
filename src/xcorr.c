#include <xcorr.h>

int nxcorr(double corr[], const double x[], size_t n_x, const double y[], size_t n_y)
{
  double *ext_x, *ext_y;
  int index;
  int lag;
  double length;
  double rxx = 0, ryy = 0, rnm;

  // Extend shorter signal with zeroes
  if (n_x > n_y) {
    length = n_x;
    ext_x = (double*) malloc(n_x * sizeof(double));
    ext_y = (double*) malloc(n_x * sizeof(double));
    for (index = 0; index < n_x; index++) {
      ext_x[index] = x[index];
      if (index < n_y)
        ext_y[index] = y[index];
      else
        ext_y[index] = 0;
    }
  }
  else if (n_x < n_y) {
    length = n_y;
    ext_x = (double*) malloc(n_y * sizeof(double));
    ext_y = (double*) malloc(n_y * sizeof(double));
    for (index = 0; index < n_y; index++) {
      ext_y[index] = y[index];
      if (index < n_x)
        ext_x[index] = x[index];
      else
        ext_x[index] = 0;
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
  for (lag = 0; lag <= (length - 1); lag++) {
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
