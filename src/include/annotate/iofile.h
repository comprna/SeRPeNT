#include <core/structs.h>
     
/*
 * wcl
 *   Line count
 *
 * @arg *fp
 *   File descriptor
 *
 * @return
 *   The number of lines in a file
 */
int wcl (FILE *fp);

/*
 * next_profile
 *   Reads a line of the profile file and stores the profile in a given pointer
 *
 * @arg 
 *
 * @return
 *   1 if more lines available. 0 if no more lines. -1 if file is ill-formatted.
 */
int next_profile(FILE* fp, profile_struct* profile);

/*
 * next_profile
 *   Reads a line of the additional profile file and stores the profile in a given pointer
 *   ID of the profile is preceeded by "species_"
 *
 * @arg 
 *
 * @return
 *   1 if more lines available. 0 if no more lines. -1 if file is ill-formatted.
 */
int next_additional_profile(FILE* fp, profile_struct* profile, char* species);

/*
 * next_feature
 *   Reads a line of the BED annotation file and stores the contents in a given pointer
 *
 * @arg
 *
 * @return
 *   1 if more lines available. 0 if no more lines. -1 if file is ill-formatted.
 */
int next_feature(FILE* bedf, feature_struct* feature);
