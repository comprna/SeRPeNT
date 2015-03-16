#include <gsl/gsl_statistics_double.h>
#include <core/structs.h>

/*
 * dtw
 *   Normalized Standard Dynamic Time Warping algorithm.
 *   Peak normalization is performed prior to applying DTW algorithm
 *   DTW score is normalized respect to the sum of the lengths of the time series (m + n)
 *
 * @see Sakoe and Chiba. 1978.
 * @see Naymat, Chawla and Taheri. 2012.
 *
 * @arg double* s
 *   Time series S = s1, s2, ..., sn
 * @arg n
 *   Length of time series s
 * @arg max_s
 *   Height of the highest peak from time series S
 * @arg double* q
 *   Time series Q = q1, q2, ..., qm
 * @arg double m
 *   Length of time series Q
 * @arg max_q
 *   Height of the highest peak from time series Q
 *
 * @return
 *   Optimal warping distance between time series S and Q
 */
double ndtw(double* s, int n, double max_s, double* q, int m, double max_q);

/*
 * dtw
 *   Standard Dynamic Time Warping algorithm
 *
 * @see Sakoe and Chiba. 1978.
 * @see Naymat, Chawla and Taheri. 2012.
 *
 * @arg double* s
 *   Time series S = s1, s2, ..., sn
 * @arg n
 *   Length of time series s
 * @arg double* q
 *   Time series Q = q1, q2, ..., qm
 * @arg double m
 *   Length of time series Q
 *
 * @return
 *   Optimal warping distance between time series S and Q
 */
double dtw(double* s, int n, double* q, int m);
