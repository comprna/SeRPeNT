#include <core/structs.h>

/*
 * create_sere
 *   Create struct that handles data to calculate SERE scores
 *
 * @args int** reads_per_contig
 *   Array containing number of reads per each contig and replicate
 * @args int n_contigs
 *   Number of contigs
 * @args int n_replicates
 *   Number of replicates
 *
 * @return
 *   A SERE data handler struct. NULL if an error ocurred.
 */
sere_struct* create_sere(int** reads_per_contig, int n_contigs, int n_replicates);

/*
 * calculate_sere_scores
 *   Filter contigs according their Irreproducibility score and a cutoff value
 *   Irreproducibility scores are calculated by using the SERE method
 *
 * @references
 *   Schulze et al. "SERE : Single-parameter quality control and sample comparison for RNA-Seq". BMC Genomics (2012).
 *
 * @args profile_struct* profiles
 *   Array of contigs
 * @args int n_contigs
 *   Total number of contigs
 * @args int n_replicates
 *   Total number of replicates
 * @args double cutoff
 *   Threshold for irreproducibility acceptance
 */
void calculate_sere_score(profile_struct* profile, sere_struct* sere_s, int n_replicates, double cutoff);

/*
 * destroy_sere
 *   Destroys a struct that handles data to calculate SERE scores
 *
 * @args sere_struct* sere
 *   A SERE data handler struct
 */
void destroy_sere(sere_struct* sere);

/*
 * calculate_common_scores
 *   Filter contigs according by checking if they are common to all the replicates
 *
 * @args profile_struct* profiles
 *   Array of contigs
 * @args int n_contigs
 *   Total number of contigs
 * @args int n_replicates
 *   Total number of replicates
 */
void calculate_common_score(profile_struct* profile, int n_replicates);

/*
 * create_npidr
 *   Create struct that handles data to calculate npIDR scores
 *
 * @args int** reads_per_contig
 *   Array containing number of reads per each contig and replicate
 * @args int n_contigs
 *   Number of contigs
 * @args int n_replicates
 *   Number of replicates
 *
 * @return
 *   A npIDR data handler struct. NULL if an error ocurred.
 */
npidr_struct* create_npidr(int** reads_per_contig, int n_contigs, int n_replicates);

/*
 * calculate_npidr_scores
 *   Filter contigs according their Irreproducibility score and a cutoff value
 *   Irreproducibility scores are calculated by using the npIDR method
 *
 * @reference
 *   Dobin et al. "STAR : ultrafast universal RNA-Seq aligner". Bioinformatics (2013).
 *
 * @args profile_struct* profiles
 *   Array of contigs
 * @args int n_contigs
 *   Total number of contigs
 * @args int n_replicates
 *   Total number of replicates
 * @args double cutoff
 *   Threshold for irreproducibility acceptance
 */
void calculate_npidr_score(profile_struct* profiles, npidr_struct* npidr, int index, int n_contigs, int n_replicates, double cutoff);

/*
 * destroy_npidr
 *   Destroys a struct that handles data to calculate npIDR scores
 *
 * @args sere_npidr* sere
 *   A npIDR data handler struct
 */
void destroy_npidr(npidr_struct* npidr);
