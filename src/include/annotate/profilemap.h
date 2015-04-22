#include <annotate/itvltree.h>

/*
 * map_init
 *   Initializes the profile map
 *
 * @arg map_struct* map
 *   Pointer to a profile map
 */
void map_init(map_struct* map);

/*
 * map_sort
 *   Sorts the profile map
 *
 * @arg map_struct* map
 *   Pointer to a profile map
 */
void map_sort(map_struct* map);

/*
 * map_search
 *   Searches a profile within the map
 *
 * @arg map_struct* map
 *   Pointer to a profile map
 * @arg profile_struct* profile
 *   The profile to be found
 *
 * @return
 *   The index of the profile within the map.
 *   -1 if profile is not in the map.
 */
int map_search(map_struct* map, char* chrom, int strand);

/*
 * map_profile
 *   Adds a profile to the map
 *
 * @args map_struct* map
 *   Pointer to a profile map
 * @args profile_struct* profile
 *   Profile to be mapped
 */
void map_add_profile(map_struct* map, profile_struct_annotation* profile);

/*
 * map_annotate
 *   Annotates a profile in the map
 *
 * @args args_a_struct* arguments
 *   Pointer to a struct handling the command line parameters
 * @args map_struct* map
 *   Pointer to a profile map
 * @args char* chromosome
 *   Chromosome of the feature for annotation
 * @args int start
 *   Start genomic coordinate of the feature
 * @args int end
 *   End genomic coordinate of the feature
 * @args int strand
 *   Genomic strand of the feature
 */
void map_annotate(args_a_struct* arguments, map_struct* map, char* chrom, int start, int end, int strand, char* feature);

/*
 * map_destroy
 *   Destroys a map
 *
 * @arg map_struct* map
 *   Pointer to a profile map 
 */
void map_destroy(map_struct* map);
