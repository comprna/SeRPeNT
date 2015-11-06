#include <core/structs.h>
#include <float.h>

/*
 * Optimization of the distance cutoff (dc) for the dclust clustering algorithm.
 *
 * dc is calculated by finding the impact factor (sigma) of a gaussian density function
 * that minimizes the value of the entropy of the potentials (H).
 *
 * @reference: Rodriguez A. and Laio A. "Clustering by fast search and find of density peaks".
 *             Science 27 June 2014
 * @reference: Wang et al. "Data Field for Hierarchical Clustering".
 *             International Journal of Data Warehousing and Mining 7(4) October 2011
 * @reference: Wang et al. "Comment on <Clustering by fast search and find of density peaks>".
 *             aRxiv:1501.04267
 *
 * @arg double** dist
 *   Distance/dissimilarity matrix
 * @arg int n
 *   Number of elements in dist
 * @arg double* max
 *   Pointer to a double variable where the maximum distance in dist will be stored
 *
 * @return
 *   Optimal distance cutoff (dc) for the dclust clustering algorithm
 */
double dcoptimize(double** dist, int n, double* max);
