#include <core/structs.h>
#include <samtools/sam.h>

/*
 * next_alignment
 *   Read and process line in a BAM file
 *
 * @arg samfile_t *bam_file
 *   Pointer to a BAM/SAM file descriptor
 * @arg alignment_struct alignment
 *   Pointer to an alignment handler struct
 * @arg int replica
 *   Number of replicate
 * @arg args_p_struct *arguments
 *   Pointer to an arguments handler struct
 *
 * @return
 *   -2 if error ocurred. -1 if EOF. 0 othwerwise.
 */
int next_alignment(samfile_t *bam_file, alignment_struct *alignment, int replica, args_p_struct *arguments);

/*
 * parse_alignment
 *   Parse processed alignment from BAM file
 *
 * @arg args_p_struct *arguments
 *   Pointer to an arguments handler struct
 * @arg alignment_struct *alignment
 *   Pointer to an alignment handler struct
 * @arg profile_struct *current_profile
 *   Pointer to the current profile handler struct
 * @arg contig_struct *primary
 *   Pointer to a contig handler struct
 * @arg int *index
 *   Pointer to an integer containing the profile index
 *
 * @return
 *   -1 if memory corruption. 0 otherwise.
 */
int parse_alignment(const args_p_struct *arguments, const alignment_struct *alignment, profile_struct *current_profile, contig_struct *primary, int *index);
