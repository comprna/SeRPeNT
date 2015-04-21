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
double tmin(double x, double y, double z)
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
  warping[0][0] = sqrtsqr(s[0] - q[0]);
  for (i = 1; i < n; i++) warping[i][0] = sqrtsqr(s[i] - q[0]) + warping[i - 1][0];
  for (i = 1; i < m; i++) warping[0][i] = sqrtsqr(s[0] - q[i]) + warping[0][i - 1];
  for(i = 1; i < n; i++)
    for(j = 1; j < m; j++)
      warping[i][j] = tmin(sqrtsqr(s[i] - q[j]) + warping[i - 1][j],
                           sqrtsqr(s[i] - q[j]) * 2 + warping[i - 1][j - 1],
                           sqrtsqr(s[i] - q[j]) + warping[i][j - 1]);

  dtwscore = warping[n - 1][m - 1];

  for (i = 0; i < n; i++) free(warping[i]);
  free(warping);

  return dtwscore;
}


/*
 * xdtw
 *
 * @see include/annotate/dtw.h
 */
double xdtw(profile_struct_annotation* p1, profile_struct_annotation* p2) {
  double warping[MAX_PROFILE_LENGTH][MAX_PROFILE_LENGTH][3];
  int i, j, n, m, w;
  double score;

  // Initialization
  n = p1->length;
  m = p2->length;
  w = abs(n - m);
  double* s = p1->profile;
  double* q = p2->profile;
  srand(time(NULL));

  // Initial condition
  warping[0][0][0] = s[0] * q[0];
  warping[0][0][1] = s[0] * s[0];
  warping[0][0][2] = q[0] * q[0];

  // First row and column
  for (i = 1; i < n; i++) {
    double noise = p2->noise[rand() % MAX_PROFILE_LENGTH];
    warping[i][0][0] = s[i] * noise + warping[i - 1][0][0];
    warping[i][0][1] = s[i] * s[i] + warping[i - 1][0][1];
    warping[i][0][2] = noise * noise + warping[i - 1][0][2];
  }

  for (j = 1; j < m; j++) {
    double noise = p1->noise[rand() % MAX_PROFILE_LENGTH];
    warping[0][j][0] = noise * q[j] + warping[0][j - 1][0];
    warping[0][j][1] = noise * noise + warping[0][j - 1][1];
    warping[0][j][2] = q[j] * q[j] + warping[0][j - 1][2];
  }

  // Fill matrix
  for(i = 1; i < n; i++) {
    int start, stop, x, y;

    start= MAX(1, i - w);
    stop = MIN(i + w, m - 1);
    x = MIN(i + w, m - 1);
    y = MIN(i + w, n - 1);

    for(j = start; j <= stop; j++) {
      double c1, c2, c3;
      double noise1 = p1->noise[rand() % MAX_PROFILE_LENGTH];
      double noise2 = p2->noise[rand() % MAX_PROFILE_LENGTH];

      c2 = (s[i] * q[j] + warping[i - 1][j - 1][0]) / sqrt((s[i]*s[i] + warping[i - 1][j - 1][1]) * (q[j]*q[j] + warping[i - 1][j - 1][2]));
      
      if (j == start && j != 1)
        c3 = INFINITY * (-1);
      else
        c3 = (noise1 * q[j] + warping[i][j - 1][0]) / sqrt((noise1 * noise1 + warping[i][j - 1][1]) * (q[j] * q[j] + warping[i][j - 1][2]));

      if (j == stop && i != 1)
        c1 = INFINITY * (-1);
      else
        c1 = (s[i] * noise2 + warping[i - 1][j][0]) / sqrt((s[i] * s[i] + warping[i - 1][j][1]) * (noise2 * noise2 + warping[i - 1][j][2]));

      if (c2 >= c1 && c2 >= c3) {
        warping[i][j][0] = warping[i - 1][j - 1][0] + s[i] * q[j];
        warping[i][j][1] = warping[i - 1][j - 1][1] + s[i] * s[i];
        warping[i][j][2] = warping[i - 1][j - 1][2] + q[j] * q[j];
      }
      else if (c1 >= c2 && c1 >= c3) {
        warping[i][j][0] = warping[i - 1][j][0] + s[i] * noise2;
        warping[i][j][1] = warping[i - 1][j][1] + s[i] * s[i];
        warping[i][j][2] = warping[i - 1][j][2] + noise2 * noise2;
      }
      else {
        warping[i][j][0] = warping[i][j - 1][0] + noise1 * q[j];
        warping[i][j][1] = warping[i][j - 1][1] + noise1 * noise1;
        warping[i][j][2] = warping[i][j - 1][2] + q[j] * q[j];
      }
    }
  }

  score = warping[n - 1][m - 1][0] / sqrt(warping[n - 1][m - 1][1] * warping[n - 1][m - 1][2]);

  return score;
}


/*
 * adtw
 *
 * @see include/annotate/dtw.h
 */
double adtw(profile_struct_annotation* p1, profile_struct_annotation* p2) {
  double warping[MAX_PROFILE_LENGTH][MAX_PROFILE_LENGTH][3];
  double xcorr[MAX_PROFILE_LENGTH][MAX_PROFILE_LENGTH][3];
  int i, j, n, m, w;
  double score;

  // Initialization
  n = p1->length;
  m = p2->length;
  w = abs(n - m);
  double* s = (double*) malloc(n * sizeof(double));
  double* q = (double*) malloc(m * sizeof(double));
  srand(time(NULL));

  for (i = 0; i < n; i++) s[i] = p1->profile[i] + 1;
  for (j = 0; j < m; j++) q[j] = p2->profile[j] + 1;
  int max_s = p1->max_height + 1;
  int max_q = p2->max_height + 1;

  // Dynamic algorithm
  warping[0][0][0] = sqrtsqr(s[0] / max_s - q[0] / max_q);
  warping[0][0][1] = -1;
  warping[0][0][2] = -1;
  xcorr[0][0][0] = s[0] * q[0];
  xcorr[0][0][1] = s[0] * s[0];
  xcorr[0][0][2] = q[0] * q[0];

  for (i = 1; i < MIN(w + 1, n); i++) {
    double noise = p2->noise[rand() % MAX_PROFILE_LENGTH];
    warping[i][0][0] = sqrtsqr(s[i] / max_s - q[0] / max_q) + warping[i - 1][0][0];
    warping[i][0][1] = i - 1;
    warping[i][0][2] = 0;
    xcorr[i][0][0] = s[i] * noise + xcorr[i - 1][0][0];
    xcorr[i][0][1] = s[i] * s[i] + xcorr[i - 1][0][1];
    xcorr[i][0][2] = noise * noise + xcorr[i - 1][0][2];
  }
  if (i < n - 1) warping[w + 1][0][0] = INFINITY;

  for (j = 1; j < MIN(w + 1, m); j++) {
    double noise = p1->noise[rand() % MAX_PROFILE_LENGTH];
    warping[0][j][0] = sqrtsqr(s[0] / max_s - q[j] / max_q) + warping[0][j - 1][0];
    warping[0][j][1] = 0;
    warping[0][j][2] = j - 1;
    xcorr[0][j][0] = noise * q[j] + xcorr[0][j - 1][0];
    xcorr[0][j][1] = noise * noise + xcorr[0][j - 1][1];
    xcorr[0][j][2] = q[j] * q[j] + xcorr[0][j - 1][2];
  }
  if (j < m - 1) warping[0][w + 1][0] = INFINITY;

  for(i = 1; i < n; i++) {
    int start, stop, x, y;

    start= MAX(1, i - w);
    stop = MIN(i + w, m - 1);
    x = MIN(i + w, m - 1);
    y = MIN(i + w, n - 1);

    if (x < m - 1) warping[i][x + 1][0] = INFINITY;
    if (y < n - 1) warping[y + 1][i][0] = INFINITY;

    for(j = start; j <= stop; j++) {
      double c1 = sqrtsqr(s[i] / max_s - q[j] / max_q) + warping[i - 1][j][0];
      double c2 = sqrtsqr(s[i] / max_s - q[j] / max_q) * 0.25 + warping[i - 1][j - 1][0];
      double c3 = sqrtsqr(s[i] / max_s - q[j] / max_q) + warping[i][j - 1][0];

      if (c2 <= c1 && c2 <= c3) {
        warping[i][j][0] = c2;
        warping[i][j][1] = i - 1;
        warping[i][j][2] = j - 1;
        xcorr[i][j][0] = s[i] * q[j] + xcorr[i - 1][j - 1][0];
        xcorr[i][j][1] = s[i] * s[i] + xcorr[i - 1][j - 1][1];
        xcorr[i][j][2] = q[j] * q[j] + xcorr[i - 1][j - 1][2];
      }
      else if (c1 <= c2 && c1 <= c3) {
        double noise = p2->noise[rand() % MAX_PROFILE_LENGTH];
        warping[i][j][0] = c1;
        warping[i][j][1] = i - 1;
        warping[i][j][2] = j;
        xcorr[i][j][0] = s[i] * noise + xcorr[i - 1][j][0];
        xcorr[i][j][1] = s[i] * s[i] + xcorr[i - 1][j][1];
        xcorr[i][j][2] = noise * noise + xcorr[i - 1][j][2];
      }
      else {
        double noise = p1->noise[rand() % MAX_PROFILE_LENGTH];
        warping[i][j][0] = c3;
        warping[i][j][1] = i;
        warping[i][j][2] = j - 1;
        xcorr[i][j][0] = noise * q[j] + xcorr[i][j - 1][0];
        xcorr[i][j][1] = s[i] * s[i] + xcorr[i][j - 1][1];
        xcorr[i][j][1] = noise * noise + xcorr[i][j - 1][1];
      }
    }
  }

  score = xcorr[n - 1][m - 1][0] / sqrt(xcorr[n - 1][m - 1][1] * xcorr[n - 1][m - 1][2]);

  // Backtracking
  double* a1 = (double*) malloc(sizeof(double) * (n + m));
  double* b1 = (double*) malloc(sizeof(double) * (n + m));
  int index = 0;
  i = n - 1;
  j = m - 1;

  while (i >= 0 || j >= 0) {
    a1[index] = s[i];
    b1[index] = q[j];
    int i1 = warping[i][j][1];
    int j1 = warping[i][j][2];
    i = i1;
    j = j1;
    index++;
  }

  free(a1);
  free(b1);

  // Finalization
  free(s);
  free(q);
  return(score);
}
