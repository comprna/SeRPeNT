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
 *   memory allocation error.
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
 * @arg profile_struct* profiles
 *   An array of profiles for labeling the hierarchical clustering solution
 */
void hc_print(FILE* fp, hcnode_struct* hc, int nprofiles, profile_struct_annotation* profiles);

/*
 * hc_annotate
 *   Annotates the unannotated leafs in the hierarchical solution by using a leaf-weighted average-based algorithm
 *
 * @reference
 *   Malik and Kender. "Classification by Pattern-based hierarchical clustering".
 *
 * @arg hcnode_struct* hc
 *   A pointer to the hierarchical clustering solution
 * @arg int nprofiles
 *   Total number of profiles
 * @arg profile_struct* profiles
 *   An array of profiles for labeling the hierarchical clustering solution
 */
void hc_annotate(hcnode_struct* hc, int nprofiles, profile_struct_annotation* profiles, double** correlations, double cutoff);

void xcorr_annotate(annotation_struct** xcorr, int nprofiles, profile_struct_annotation* profiles);
