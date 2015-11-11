#include <core/structs.h>

/*
 * Perform Mann-Whitney U test on two independent samples using direct method
 *
 * For each observation in one sample, count the number of times this first value is higher than any other observation in the other sample (win). Count 0.5 if observations are equal (tie).
 * The sum of wins and ties is U for the first sample. U for the other sample is the converse.
 *
 * Being n and m the size of the first and second samples respectively, and U1 and U2 the U statitics of the first and the second samples respectively:
 *   U1 = n * m - U2
 *   U2 = n * m - U1
 *   U = min(U1, U2)
 *
 * If U <= Ucritical then null hypothesis is rejected at a chosen alpha level.
 *
 * @reference https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test
 *
 * @arg double* p
 *   All the observations in the first sample
 * @arg int n
 *   Number of observations in the first sample
 * @arg double* q
 *   All the observations in the second sample
 * @arg int m
 *   Number of observations in the second sample
 * @arg double pval
 *   Desired p-value of the given test
 *
 * @return
 *   1 if the null hipothesis is rejected with a p-value <= pval. 0 otherwise.
 */
int mannwhitney_d (double* p, int n, double* q, int m, double pval);

/*
 * Perform Mann-Whitney U test on two independent samples using indirect method
 *
 * Being n and m the size of the first and second samples respectively:
 *   1. Assign numeric ranks to all the observations in both samples, beginning with 1 for the smallest value.
 *      Where there are groups of tied values, assign a rank equal to the midpoint of unadjusted rankings [e.g., the ranks of (3, 5, 5, 9) are (1, 2.5, 2.5, 4)].
 *   2. Add up the ranks for the observations which came from sample 1 (R1).
 *   3. R2 is now determinate, since the sum of all the ranks equals N(N + 1)/2 where N is the total number of observations.
 *   4. U1 and U2 (statistics of sample 1 and respectively) are given by:
 *      U1 = R1 - n * (n+1) / 2
 *      U2 = R2 - m * (m+1) / 2
 *
 * p-value is computed for U using the normal distribution approximation:
 *   Z = (U - mu) / sd
 * where mu and sd are the mean and the standard deviation of U. Are given by:
 *   mu = n*m / 2
 *   sd = sqrt(n*m*(n+m+1) / 12)
 *
 * Since U1+U2 = n*m, mu is the mean of the two values of U. Therefore, the absolute value of the z statistic calculated will be same whichever value of U is used.
 *
 * @reference https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test
 *
 * @arg double* p
 *   All the observations in the first sample
 * @arg int n
 *   Number of observations in the first sample
 * @arg double* q
 *   All the observations in the second sample
 * @arg int m
 *   Number of observations in the second sample
 * @arg double pval
 *   Desired p-value of the given test
 *
 * @return
 *   Mann-Whitney U test p-value
 */
int mannwhitney_i (double* p, int n, double* q, int m, double pval);

/*
 * Tests applicability of Mann-Whitney U Test based on the size of the samples and the alpha level (p-value)
 * 
 * @reference http://math.usask.ca/~laverty/S245/Tables/wmw.pdf
 *
 * @arg int n
 *   Number of observations in the first sample
 * @arg int m
 *   Number of observation in the second sample
 * @arg double pvalue
 *   alpha level
 *
 * @return
 *   1 if Mann-Whitney U test is applicable. 0 otherwise.
 */
int mannwhitney_a (int n, int m, double pval);
