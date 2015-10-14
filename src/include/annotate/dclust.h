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
 * @arg int k
 *   Number of clusters. If k = 0 then number of clusters is computed.
 * @arg profile_struct_annotation* profiles
 *   An array of profiles
 *
 * @return
 *   Number of clusters
 */
int dclust(double** dist, int n, int k, profile_struct_annotation* profiles);
