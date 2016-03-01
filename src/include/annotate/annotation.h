#include <core/structs.h>
#include <annotate/strmap.h>

/*
 * xcorr_annotate
 * 
 * Annotation of profiles based on top 3 nearest neighbours
 *
 * @arg annotation_struct** xcorr
 *   Cross-correlation matrix. Annotation structs are triplets <profile_index:profile_index:x-correlation>
 * @arg int nprofiles
 *   Number of profiles
 * @arg profile_struct_annotation* profiles
 *   Array of profile structs
 */
void xcorr_annotate(annotation_struct** xcorr, int nprofiles, profile_struct_annotation* profiles);

/*
 * cluster_annotate
 *
 * Annotation of profiles based on majority vote within clusters
 *
 * @arg int nclusters
 *   Total number of clusters
 * @arg int nprofiles
 *   Total number of profiles
 * @arg profile_struct_annotation* profiles
 *   Array of profile structs
 */
void cluster_annotate(int nclusters, int nprofiles, profile_struct_annotation* profiles);
