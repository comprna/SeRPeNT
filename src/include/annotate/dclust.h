#include <core/structs.h>

/*
 * Calculate a clustering by fast search and find of density peaks
 *
 * @reference: Rodriguez A. and Laio A. "Clustering by fast search and find of density peaks".
 *             Science 27 June 2014.
 *
 * @arg double** dist
 *   Distance matrix
 * @arg int n
 *   Number of elements in the distance matrix
 * @arg profile_struct_annotation* profiles
 *   An array of profiles
 * @arg double cutoff
 *   Distance cutoff (dc)
 * @arg int gaussian
 *   0 if no gaussian kernel for density calculation. 1 otherwise.
 *
 * @return
 *   Number of clusters
 */
int dclust(double** dist, int n, profile_struct_annotation* profiles, double cutoff, int gaussian);
