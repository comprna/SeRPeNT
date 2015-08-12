#include <core/structs.h>
#include <annotate/cluster.h>

/*
 * hc_cluster
 *   Performs pairwise maximum- (or complete)- linkage clustering
 *
 * @reference
 *    De Hoon et al. "Open Source Clustering Software". Bioinformatics(2004).
 *
 * @arg double** correlation
 *   Correlation matrix containing distances between profiles
 * @arg int nprofiles
 *   Total number of profiles
 *
 * @return
 *   A pointer to a newly allocated array of hcnode_struct structures that describes
 *   the calculated hierarchical clustering solution. NULL if hc_cluster fails due to
 *   memory allocation error
 */
hcnode_struct* hc_cluster(double** correlation, int nprofiles);

/*
 * hc_print
 *   Prints a hierarchical clustering solution to a file in neWick format
 *
 * @arg FILE* fp
 *   The file descriptor where the hierarchical clustering solution will be printed
 * @arg hcnode_struct* hc
 *   A pointer to the hierarchical clustering solution
 * @arg int nprofiles
 *   Total number of profiles
 * @arg profile_struct_annotation* profiles
 *   An array of profiles
 */
void hc_print(FILE* fp, hcnode_struct* hc, int nprofiles, profile_struct_annotation* profiles);

/*
 * hc_branch
 *   Branches the tree given a cutoff value. All the leafs under a branch belong to the parent node
 *   that has the greatest distance that is lower or equal to the cutoff value
 *
 * @arg hcnode_struct* hc
 *   A pointer to the hierarchical clustering solution
 * @arg int nprofiles
 *   Total number of profiles
 * @arg profile_struct_annotation* profiles
 *   An array of profiles
 * @arg double cutoff
 *   The cutoff value to branch the tree
 *
 * @return
 *   The number of clusters
 */
int hc_branch(hcnode_struct* hc, int nprofiles, profile_struct_annotation* profiles, double cutoff);

/*
 * hc_eval
 *   Evaluates the branched hierarchical solution in terms of the external information-based measure V
 *
 * @reference
 *   Rosenberg and Hirchsberg.
 *   V-Measure : a conditional entropy-based external cluster evaluation measure
 *   EMNLP (2007)
 *
 * @arg int nprofiles
 *   Total number of profiles
 * @arg profile_struct_annotation* profiles
 *   An array of profiles
 * 
 * @return
 *   The V-measure score of the branched hierarchical solution
 */
double hc_eval(int nprofiles, profile_struct_annotation* profiles);
