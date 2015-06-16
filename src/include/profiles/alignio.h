#include <core/structs.h>
#include <samtools/sam.h>

/*
 * next_alignment
 *   Read and process line in a BAM file and store data in an alignment handler struct
 *
 * @arg samfile_t *bam_file
 *   Pointer to a BAM/SAM file descriptor
 * @arg alignment_struct alignment
 *   Pointer to an alignment handler struct
 * @arg int replica
 *   Replicate number
 * @arg args_p_struct *arguments
 *   Pointer to an arguments handler struct
 *
 * @return
 *   -2 if error ocurred. -1 if EOF. 0 othwerwise.
 */
int next_alignment(samfile_t *bam_file, alignment_struct *alignment, int replica, args_p_struct *arguments);

/*
 * parse_alignment
 *   Parse processed alignment in an alignment handler struct and store date in a contig handler struct
 *   Store processed contigs in a temporal profiles file
 *
 * @arg args_p_struct *arguments
 *   Pointer to an arguments handler struct
 * @arg alignment_struct *alignment
 *   Pointer to an alignment handler struct
 * @arg contig_struct *primary
 *   Pointer to a contig handler struct
 * @arg FILE* tmprofiles_file
 *   Pointer to a profile temporal file
 *
 * @return
 *   -1 if memory corruption. 0 otherwise.
 */
int parse_alignment(const args_p_struct *arguments, const alignment_struct *alignment, contig_struct *primary, FILE* tmprofiles_file);

/*
 * next_tmprofile
 *   Read and process stored alignment in the alignment temporal file and store data in a profile handler struct
 *
 * @arg FILE* fp
 *   Pointer to a profile temporal file
 * @arg profile_struct* profile
 *   Pointer to a profile handler struct
 *
 * @return
 *   -1 if file corrupted. 0 otherwise.
 */
int next_tmprofile(FILE* fp, profile_struct* profile);
