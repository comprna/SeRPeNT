#include <annotate/dclust.h>

int cmpd(const void *x, const void *y)
{
  double xx = *(double*)x, yy = *(double*)y;
  if (xx < yy) return -1;
  if (xx > yy) return  1;
  return 0;
}


/*
 * dcoptimize
 *
 * @see include/annotate/dclust.h
 */
double dcoptimize(double** dist, int n, double* max)
{
  double dc, z, h, hmin, lower, upper, sigma;
  double* allds;
  int i, j, counter;
  double* potentials;

  // Initialize variables
  potentials = (double*) malloc(n * sizeof(double));
  allds = (double*) malloc(sizeof(double) * (n * (n - 1) / 2));
  hmin = DBL_MAX;

  // Convert dissimilarity matrix into a vector
  // Sort dissimilarity vector ascendently
  counter = 0;
  for (i = 0; i < n; i++) {
    for(j = i + 1; j < n; j++) {
      allds[counter] = dist[i][j];
      counter++;
    }
  }
  qsort(allds, n * (n - 1) / 2, sizeof(double), cmpd);
  *max = allds[counter - 1];

  // Calculate lower and upper boundaries for sigma (impact factor)
  //   lower : minimum distance > 0
  //   upper : distance at 10 percentile
  counter = 0;
  while(allds[counter] == 0) counter++;
  lower = allds[counter];
  upper = gsl_stats_quantile_from_sorted_data(allds, 1, (n * (n - 1) / 2), 0.10);

  // Calculate entropy for values of sigma between lower and upper in increments of 0.005
  for (sigma = lower; sigma <= upper; sigma += 0.005f) {

    // Calculate potential for each point
    for(i = 0; i < n; i++) {
      potentials[i] = 0.0f;
      for(j = 0; j < n; j++) {
        double sq = (dist[i][j] / sigma) * (dist[i][j] / sigma);
        double expsq = exp(-sq);
        if (j != i) potentials[i] += expsq;
      }
    }

    // Calculate Z (normalization factor)
    z = 0;
    for (i = 0; i < n; i++)
      z += potentials[i];

    // Calculate H (entropy)
    h = 0;
    for (i = 0; i < n; i++) {
      double a = potentials[i] / z;
      double b = log(a);
      h += a*b;
    }
    h = h * (-1);

    // Store dc and hmin
    if (h < hmin) {
      hmin = h;
      dc = sigma;
    }
  }

  free(potentials);
  free(allds);

  return((3/sqrt(2)) * dc);
}
