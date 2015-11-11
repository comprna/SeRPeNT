#include <core/structs.h>

/*
 * next_diffproc_feature
 *   Reads a line of the BED cluster file and stores the contents in a given pointer
 *
 * @arg
 *
 * @return
 *   1 if more lines available. 0 if no more lines. -1 if file is ill-formatted.
 */
int next_diffproc_feature(FILE* bedf, feature_struct_diffproc* feature);

/*
 * next_profile
 *   Reads a line of the profile file and stores the profile in a given pointer
 *
 * @arg 
 *
 * @return
 *   1 if more lines available. 0 if no more lines. -1 if file is ill-formatted.
 */
int next_diffproc_profile(FILE* fp, profile_struct_diffproc* profile);

/*
 * find_clusters
 *   Reads all the lines in a BED cluster file and returns the numbe of clusters
 *
 * @arg
 *
 * @return
 *   The number of clusters
 */
int find_clusters(FILE* bedf);

void allocate_clusters(FILE* bedf, int* profiles_per_cluster);
