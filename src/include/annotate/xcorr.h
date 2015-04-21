#include <core/structs.h>

/*
 * Calculate the deterministic cross-correlation between two deterministic
 * signals x and y over the range of lags : 0 to max(n_x-1, n_y-1).
 *
 * The correlation is normalized so that the auto-correlations at 
 * 0 lag are identically 1.0
 *
 * @reference: Orfanidis, "Optimum Signal Processing. An Introduction"
 *             2nd Ed. Macmillan, 1988.
 *
 * @arg profile_struct p1 
 *   Profile struct containing the first deterministic signal
 * @arg profile_struct p2
 *   Profile struct containing the second deterministic signal
 *
 * @return The maximum x-corr score over the range of lags
 */
double nxcorr(profile_struct_annotation* p1, profile_struct_annotation* p2);
