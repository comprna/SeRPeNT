#include <core/structs.h>

/*
 * dtw
 *   Sakoe-Chiba Dynamic Time Warping algorithm
 *
 * @see Sakoe and Chiba. 1978.
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

/*
 * xdtw
 *   X-Correlation-based dynamic time warping with Sakoe Chiba band as global constraint.
 *
 *   Calculate the optimal deterministic normalized cross-correlation between two
 *   deterministic signals.
 *   
 *   Gaps are penalized by introducing gaussian white noise.
 *
 * @see Sakoe and Chiba. 1978.
 * @see Simulation of communication systems. 1992.
 * @see Orfanidis. 1988.
 *
 * @arg profile_struct_annotation* p1
 *   Profile handler struct containing the first time series
 * @arg profile_struct_annotation* p2
 *   Profile handler struct containing the second time series
 *
 * @return
 *   Optimal normalized X-Correlation between signals in p1 and p2
 */
double xdtw(profile_struct_annotation* p1, profile_struct_annotation* p2);

/*
 * adtw
 *   Normalized alignment for Standard Dynamic Time Warping algorithm.
 *   Auto Peak normalization is performed prior to applying DTW algorithm
 *   X-Correlation from the alignment is reported
 *
 * @see Sakoe and Chiba. 1978.
 *
 * @arg profile_struct_annotation* p1
 *   Profile handler struct containing the first time series
 * @arg profile_struct_annotation* p2
 *   Profile handler struct containing the second time series
 *
 * @return
 *   Nomrmalized x-correlation from the warping alignment
 */
double adtw(profile_struct_annotation* p1, profile_struct_annotation* p2);
