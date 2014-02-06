# include <stdlib.h>
# include <math.h>
# include <stdio.h>
# include <debug.h>

/*
 * Calculate the deterministic correlation between two deterministic
 * signals x and y, of length n_x and n_y respectively, over the range
 * of lags : 0 to max(n_x-1, n_y-1).
 *
 * The correlation is normalized so that the auto-correlations at 
 * 0 lag are identically 1.0
 *
 * @reference: Orfanidis, "Optimum Signal Processing. An Introduction"
 *             2nd Ed. Macmillan, 1988.
 *
 * @arg double corr[] 
 *   Normalized correlation for the range of lags
 * @arg double x[]
 *   Deterministic signal x
 * @arg n_x
 *   Length of the signal x
 * @arg double y[]
 *   Deterministic signal y
 * @arg n_y
 *   Length of the signal y
 *
 * @return The index position on the corr vector with the maximum x-corr score
 */
int nxcorr(double corr[], const double x[], size_t n_x, const double y[], size_t n_y);
