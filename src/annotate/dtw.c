#include <annotate/dtw.h>

/*
 * min
 *   Compute the minimum of three values
 *
 * @arg double x, double y, double z
 *   Three double values
 *
 * @return
 *   The minimum of x, y and z
 */
double min(double x, double y, double z)
{
  if( ( x <= y ) && ( x <= z ) ) return x;
  if( ( y <= x ) && ( y <= z ) ) return y;
  if( ( z <= x ) && ( z <= y ) ) return z;
  return 0;
}

/*
 * sqrtsqr
 *   Compute square root of square of a double
 *
 * @arg double x
 *
 * @return
 *   The square root of the square of x
 */
double sqrtsqr(double x)
{
  return (sqrt(x * x));
}


/*
 * ndtw
 *
 * @see include/annotate/ndtw.h
 */
double ndtw(double* s, int n, double max_s, double* q, int m, double max_q) {
  double** warping;
  double cf_s, cf_q, dtwscore;
  int i, j;

  // Calculate peak normalization correction factors
  if (max_s > max_q) {
    cf_s = 1;
    cf_q = ((double) max_s) / ((double) max_q);
  }
  else {
    cf_s = ((double) max_q) / ((double) max_s);
    cf_q = 1;
  }

  // Initialization
  warping = (double**) malloc(n * sizeof(double*));
  for (i = 0; i < n; i++) warping[i] = (double*) malloc(m * sizeof(double));

  // Dynamic algorithm
  warping[0][0] = sqrtsqr((s[0] * cf_s) - (q[0] * cf_q)) * 2;
  for (i = 1; i < n; i++) warping[i][0] = sqrtsqr((s[i] * cf_s) - (q[0] * cf_q)) + warping[i - 1][0];
  for (i = 1; i < m; i++) warping[0][i] = sqrtsqr((s[0] * cf_s) - (q[i] * cf_q)) + warping[0][i - 1];
  for(i = 1; i < n; i++) 
   for(j = 1; j < m; j++) 
     warping[i][j] = min (sqrtsqr((s[i] * cf_s) - (q[j] * cf_q)) + warping[i - 1][j],
                          sqrtsqr((s[i] * cf_s) - (q[j] * cf_q)) * 2 + warping[i - 1][j - 1],
                          sqrtsqr((s[i] * cf_s) - (q[j] * cf_q)) + warping[i][j - 1]);

  dtwscore = warping[n - 1][m - 1] / (m + n);

  for (i = 0; i < n; i++) free(warping[i]);
  free(warping);

  return (dtwscore);
}


/*
 * dtw
 *
 * @see include/annotate/dtw.h
 */
double dtw(double* s, int n, double* q, int m) {
  double** warping;
  double dtwscore;
  int i, j;

  // Initialization
  warping = (double**) malloc(n * sizeof(double*));
  for (i = 0; i < n; i++) warping[i] = (double*) malloc(m * sizeof(double));

  // Dynamic algorithm
  warping[0][0] = abs(s[0] - q[0]) * 2;
  for (i = 1; i < n; i++) warping[i][0] = sqrtsqr(s[i] - q[0]) + warping[i - 1][0];
  for (i = 1; i < m; i++) warping[0][i] = sqrtsqr(s[0] - q[i]) + warping[0][i - 1];
  for(i = 1; i < n; i++)
   for(j = 1; j < m; j++)
     warping[i][j] = min(sqrtsqr(s[i] - q[j]) + warping[i - 1][j],
                         sqrtsqr(s[i] - q[j]) * 2 + warping[i - 1][j - 1],
                         sqrtsqr(s[i] - q[j]) + warping[i][j - 1]);

  dtwscore = warping[n - 1][m - 1];

  for (i = 0; i < n; i++) free(warping[i]);
  free(warping);

  return dtwscore;
}
