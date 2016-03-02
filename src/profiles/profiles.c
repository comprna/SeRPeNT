#include <profiles/profiles.h>


/*
 * compare_double:
 *   Function for double comparison
 *
 * @arg const void *a
 *   Pointer to the first value to be compared
 * @arg const void *b
 *   Pointer to the second value to be compared
 *
 * @return -1 if a < b. 0 if a = b. 1 if a > b.
 */
int compare_double (const void *a, const void *b)
{
  double aa = *(double*) a;
  double bb = *(double*) b;

  if (aa < bb) return -1;
  if (aa > bb) return  1;
  return 0;
}


/*
 * deepcpy
 *  Function to copy values from one alignment handler struct to another
 *
 * @arg alignment_struct* destiny
 *   Pointer to the alignment handler struct that is updated
 *
 * @arg alignment_struct* source
 *   Pointer to the alignment handler struct that will be copied
 */
void deepcpy(alignment_struct* destiny, alignment_struct* source)
{
  strncpy(destiny->chromosome, source->chromosome, MAX_FEATURE);
  destiny->start = source->start;
  destiny->end = source->end;
  destiny->strand = source->strand;
  destiny->valid = source->valid;
  destiny->replicate = source->replicate;
  destiny->nreads = source->nreads;
}


/*
 * morcgez
 *   Multiple OR C Greater or Equal to Zero
 *
 * @arg int[] operators
 *   Each position in the operators array will be compared >= 0.
 * @args pos
 *   Number of positions to be compared
 *
 * @return
 *   1 if any position in the operators is >= 0. 0 otherwise.
 */
int morcgez(const int operators[], const int pos)
{
  int i, result = 0;

  if (pos <= 0)
    return 0;

  for (i = 0; i < pos; i++)
    result = result || (operators[i] >= 0);

  return result;
}


/*
 * nmstrcmp
 *   Negative Multiple STRCMP
 *
 * @args alignment_struct operators[]
 *   Each chromosome field in the array will be strcmp
 * @args pos
 *   Number of positions to be compared
 * @args str
 *   String to be compared
 *
 * @return
 *   1 if all the positions are different than str. 0 otherwise.
 */
int nmstrcmp(const alignment_struct operators[], const int pos, const char* str)
{
  int i, result = 1;

  if (pos <= 0)
    return 0;

  for (i = 0; i < pos; i++)
    result = result && (strcmp(operators[i].chromosome, str) != 0);

  return(result);
}


/*
 * Application entry point
 */
int profiles_sc(int argc, char **argv)
{
  // Define and declare variables
  args_p_struct arguments;                             // Struct for handling command line parameters
  FILE* tmprofiles_file;                               // Temporal output file descriptor
  FILE* profiles_file;                                 // Profiles output file desciptor
  FILE* contigs_file;                                  // Contigs output file descriptor
  char* tmprofiles_file_name;                          // Absolute path of the temporal output file
  char* profiles_file_name;                            // Absolute path of the profiles output file
  char* contigs_file_name;                             // Absolute path of the contigs output file
  samfile_t* replicate_file[MAX_REPLICATES];           // Array of BAM file descriptors
  alignment_struct current_alignments[MAX_REPLICATES]; // Array of alignment_struct struct
  int result;                                          // Result of any operation
  int results[MAX_REPLICATES];                         // Replicate-specific results
  char* error_message;                                 // Error message to display in case of abnormal termination
  contig_struct contig_fwd, contig_rev;                // Struct for building profiles
  int i, j, index, chridx, blkidx, maxstart;           // Multi-purpose indexes and checkpoint variables
  char curr_chrom[MAX_FEATURE];                        // Current chromosome
  int curr_len;                                        // Current length
  int ncontigs;                                        // Number of contigs
  int** reads_per_contig;                              // Data matrix containing number of reads per contig and replicate
  profile_struct profile;                              // Profile struct
  sere_struct* sere_s;                                 // SERE fast calculation struct pointer
  npidr_struct* npidr_s;                               // npIDR fast calculation struct pointer

  // Initialize options with default values
  arguments.number_replicates = MAX_REPLICATES;
  arguments.min_len = MIN_LEN;
  arguments.max_len = MAX_LEN;
  arguments.spacing = SPACING;
  arguments.min_reads = MIN_READS;
  arguments.replicate_treat = REPLICATE_POOL;
  arguments.replicate_number = REPLICATE_NUMBER;
  arguments.idr_method = IDR_COMMON;
  arguments.idr_cutoff = CUTOFF;
  arguments.read_length = MAX_READ_LENGTH;
  arguments.min_read_len = MIN_READ_LEN;
  arguments.trim_threshold = TRIM_THRESHOLD;
  arguments.trim_min = TRIM_MIN;
  arguments.trim_max = TRIM_MAX;

  // Parse command line
  // Exit if command is not well-formed
  if (parse_command_line_p(argc, argv, &error_message, &arguments) < 0) {
    fprintf(stderr, "%s\n", error_message);
    if ((strcmp(error_message, PROFILES_HELP_MSG) == 0) || (strcmp(error_message, VERSION_MSG) == 0))
      return(0);
    fprintf(stderr, "%s\n", ERR_PROFILES_HELP_MSG);
    return(1);
  }

  // Open replicate file for reading and read BAM header
  // Exit if replicate BAM files do not exist, are not readable or do not have header
  for (index = 0; index < arguments.number_replicates; index++) {
    if ((replicate_file[index] = samopen(arguments.replicate_f_path[index], "rb", 0)) == 0) {
      fprintf(stderr, "%s - %s\n", ERR_BAM_F_NOT_READABLE, arguments.replicate_f_path[index]);
      return(1);
    }
    if (replicate_file[index]->header == 0) {
      fprintf(stderr, "%s - %s\n", ERR_BAM_F_NOT_HEADER, arguments.replicate_f_path[index]);
      return(1);
    }
  }

  // Check that all the replicates have reads in the same chromosomes
  // TODO => replicates can have different number of chromosomes
  result = replicate_file[0]->header->n_targets;
  for (i = 1; i < arguments.number_replicates; i++) {
    if (replicate_file[0]->header->n_targets != result) {
      fprintf(stderr, "%s - %s\n", ERR_BAM_F_TRUNCATED, arguments.replicate_f_path[i]);
      return(1);
    }
  }
  for (i = 0; i < replicate_file[0]->header->n_targets; i++) {
    for (j = 1; j < arguments.number_replicates; j++) {
      if (strcmp(replicate_file[0]->header->target_name[i], replicate_file[j]->header->target_name[i]) != 0) {
        fprintf(stderr, "%s - %s\n", ERR_BAM_F_TRUNCATED, arguments.replicate_f_path[j]);
        return(1);
      }
    }
  }

  // Build absolute paths for temporal and output files
  tmprofiles_file_name = malloc((MAX_PATH + strlen(TMPROFILES_SUFFIX) + 2) * sizeof(char));
  strncpy(tmprofiles_file_name, arguments.output_f_path, MAX_PATH);
  strcat(tmprofiles_file_name, PATH_SEPARATOR);
  strcat(tmprofiles_file_name, TMPROFILES_SUFFIX);
  profiles_file_name = malloc((MAX_PATH + strlen(PROFILES_SUFFIX) + 2) * sizeof(char));
  strncpy(profiles_file_name, arguments.output_f_path, MAX_PATH);
  strcat(profiles_file_name, PATH_SEPARATOR);
  strcat(profiles_file_name, PROFILES_SUFFIX);
  contigs_file_name = malloc((MAX_PATH + strlen(CONTIGS_SUFFIX) + 2) * sizeof(char));
  strncpy(contigs_file_name, arguments.output_f_path, MAX_PATH);
  strcat(contigs_file_name, PATH_SEPARATOR);
  strcat(contigs_file_name, CONTIGS_SUFFIX);

  // Open temporal output file for writing temporal results
  tmprofiles_file = fopen(tmprofiles_file_name, "w");
  if (!tmprofiles_file) {
    fprintf(stderr, "%s\n", ERR_OUTPUT_F_NOT_WRITABLE);
    return (1);
  }

  // Initialize variables for reading BAM files
  // Exit if BAM files are truncated or ill-formed, or if not enough memory
  fprintf(stderr, "[LOG] GENERATING PROFILES\n");
  contig_fwd.start = 0; contig_fwd.end = 0; contig_rev.start = 0; contig_rev.end = 0;
  index = 0;
  chridx = 0;
  blkidx = 1;
  strncpy(curr_chrom, replicate_file[0]->header->target_name[chridx], MAX_FEATURE);
  curr_len = replicate_file[0]->header->target_len[chridx];
  fprintf(stderr, "[LOG]   Parsing chromosome %s\n", curr_chrom);

  if (curr_len > (MAX_BLOCK * blkidx))
    maxstart = MAX_BLOCK * blkidx;
  else
    maxstart = curr_len;

  for (i = 0; i < arguments.number_replicates; i++) {
    do {
      results[i] = next_alignment(replicate_file[i], &current_alignments[i], i, &arguments);
    } while(results[i] > -1 && !current_alignments[i].valid);
    if (results[i] < 0) {
      fprintf(stderr, "%s\n", ERR_BAM_F_TRUNCATED);
      return(1);
    }
  }

  // Read BAM file and build sRNA profiles
  // Exit if BAM files are truncated or ill-formed, or if not enough memory
  while (morcgez(results, arguments.number_replicates)) {
    alignment_struct* alignments;
    int alignment_counter = 0;
    heap_struct sorter; 

    alignments = (alignment_struct*) malloc (MAX_BLOCK * arguments.read_length * 2 * arguments.number_replicates * sizeof(alignment_struct));//TODO -> code read length in arguments

    // Add all the reads in replicates to the alignments vector.
    // Collapse all the identical reads into one single alignment.
    for (i = 0; i < arguments.number_replicates; i++) {
      while ((results[i] > -1) &&
             (strcmp(current_alignments[i].chromosome, curr_chrom) == 0) &&
             (current_alignments[i].start <= maxstart))
      {
        int add_heap = 1;

        if(alignment_counter != 0) {
          int pointer = alignment_counter - 1;

          while((add_heap == 1) && (pointer >= 0) && (current_alignments[i].start == alignments[pointer].start) &&
                (strcmp(current_alignments[i].chromosome, alignments[pointer].chromosome) == 0)) {
            if ((current_alignments[i].end == alignments[pointer].end) &&
                (current_alignments[i].strand == alignments[pointer].strand) &&
                (i == alignments[pointer].replicate)) {
              add_heap = 0;
              alignments[pointer].nreads++;
            }
            pointer--;
          }
        }

        if (add_heap) {
          deepcpy(&alignments[alignment_counter], &current_alignments[i]);
          alignment_counter++;
        }

        do {
          results[i] = next_alignment(replicate_file[i], &current_alignments[i], i, &arguments);
        } while(results[i] > -1 && !current_alignments[i].valid);
      }
    }

    // Insert all alignments into heap
    initbh(&sorter, alignment_counter);
    for (i = 0; i < alignment_counter; i++) insertbh(&sorter, &alignments[i]);

    // Process all elements in the heap
    for (i = 0; i < alignment_counter; i++) {
      alignment_struct* algn = deletebh(&sorter);

      if (algn->strand == FWD_STRAND) {
        if (parse_alignment(&arguments, algn, &contig_fwd, tmprofiles_file) < 0) {
          fprintf(stderr, "%s\n", ERR_REALLOC_FAILED);
          return(1);
        }
      }
      else {
        if (parse_alignment(&arguments, algn, &contig_rev, tmprofiles_file) < 0) {
          fprintf(stderr, "%s\n", ERR_REALLOC_FAILED);
          return(1);
        }
      }
    }

    // Update curr_len
    blkidx++;
    if (curr_len > (MAX_BLOCK * blkidx))
      maxstart = MAX_BLOCK * blkidx;
    else
      maxstart = curr_len;

    // Check chromosome change
    if (nmstrcmp(current_alignments, arguments.number_replicates, curr_chrom)) {
      chridx++;
      blkidx = 1;
      strncpy(curr_chrom, replicate_file[0]->header->target_name[chridx], MAX_FEATURE);
      curr_len = replicate_file[0]->header->target_len[chridx];

      if (curr_len > (MAX_BLOCK * blkidx))
        maxstart = MAX_BLOCK * blkidx;
      else
        maxstart = curr_len;

      fprintf(stderr, "[LOG]   Parsing chromosome %s\n", curr_chrom);
    }

    free(alignments);
    destroybh(&sorter);
  }

  // Flush the contig contents in the forward contig struct
  fprintf(tmprofiles_file, "%s", contig_fwd.chromosome);
  fprintf(tmprofiles_file, "\t%d", contig_fwd.start);
  fprintf(tmprofiles_file, "\t%d", contig_fwd.end);
  fprintf(tmprofiles_file, "\t%d", FWD_STRAND);
  fprintf(tmprofiles_file, "\t%d", arguments.number_replicates);
  for (i = 0; i < arguments.number_replicates; i++) fprintf(tmprofiles_file, "\t%d", contig_fwd.nreads[i]);
  free(contig_fwd.nreads);

  // Flush the profile contents in the forward contig struct if allowed by parameters
  if ((gsl_stats_max(contig_fwd.profile, 1, (contig_fwd.end - contig_fwd.start + 1)) >= arguments.min_reads) &&  // Contig has more than r reads
      ((contig_fwd.end - contig_fwd.start + 1) >= arguments.min_len))                                            // Contig has, at most, M nucleotides
  {
    fprintf(tmprofiles_file, "\t%d", contig_fwd.end - contig_fwd.start + 1);
    for (i = 0; i < (contig_fwd.end - contig_fwd.start + 1); i++) {
      if (arguments.replicate_treat == REPLICATE_MEAN)
        fprintf(tmprofiles_file, "\t%f", contig_fwd.profile[i] / ((double) arguments.number_replicates));
      else
        fprintf(tmprofiles_file, "\t%f", contig_fwd.profile[i]);
    }
    fprintf(tmprofiles_file, "\n");
  }
  else
    fprintf(tmprofiles_file, "\t%d\n", 0);
  free(contig_fwd.profile);

  // Flush the contents in the reverse contig struct
  fprintf(tmprofiles_file, "%s", contig_rev.chromosome);
  fprintf(tmprofiles_file, "\t%d", contig_rev.start);
  fprintf(tmprofiles_file, "\t%d", contig_rev.end);
  fprintf(tmprofiles_file, "\t%d", REV_STRAND);
  fprintf(tmprofiles_file, "\t%d", arguments.number_replicates);
  for (i = 0; i < arguments.number_replicates; i++) fprintf(tmprofiles_file, "\t%d", contig_rev.nreads[i]);
  free(contig_rev.nreads);

  // Flush the profile contents in the reverse contig struct if allowed by parameters
  if ((gsl_stats_max(contig_rev.profile, 1, (contig_rev.end - contig_rev.start + 1)) >= arguments.min_reads) &&  // Contig has more than r reads
      ((contig_rev.end - contig_rev.start + 1) >= arguments.min_len))                                            // Contig has, at most, M nucleotides
  {
    fprintf(tmprofiles_file, "\t%d", contig_rev.end - contig_rev.start + 1);
    for (i = 0; i < (contig_rev.end - contig_rev.start + 1); i++) {
      if (arguments.replicate_treat == REPLICATE_MEAN)
        fprintf(tmprofiles_file, "\t%f", contig_rev.profile[contig_rev.end - contig_rev.start - i] / ((double) arguments.number_replicates));
      else
        fprintf(tmprofiles_file, "\t%f", contig_rev.profile[contig_rev.end - contig_rev.start - i]);
    }
    fprintf(tmprofiles_file, "\n");
  }
  else
    fprintf(tmprofiles_file, "\t%d\n", 0);
  free(contig_rev.profile);
  fclose(tmprofiles_file);

  // Count lines in profiles file
  tmprofiles_file = fopen(tmprofiles_file_name, "r");
  if (!tmprofiles_file) {
    fprintf(stderr, "%s\n", ERR_INPUT_F_NOT_READABLE);
    return (1);
  }
  ncontigs = wcl(tmprofiles_file);
  fclose(tmprofiles_file);
  
  // Allocate memory for storing irreproducibility data
  reads_per_contig = (int**) malloc(ncontigs * sizeof(int*));
  for (i = 0; i < ncontigs; i++) reads_per_contig[i] = (int*) malloc(arguments.number_replicates * sizeof(int));
  
  // Read profiles and store irreproducibility and trimming data
  tmprofiles_file = fopen(tmprofiles_file_name, "r");
  if (!tmprofiles_file) {
    fprintf(stderr, "%s\n", ERR_INPUT_F_NOT_READABLE);
    return (1);
  }
  index = 0;
  while((result = next_tmprofile(tmprofiles_file, &profile) > 0)) {
    for (i = 0; i < arguments.number_replicates; i++)
      reads_per_contig[index][i] = profile.nreads[i];
    free(profile.nreads);
    if (profile.free) free(profile.profile);
    index++;
  }
  if (result < 0) {
    fprintf(stderr, "%s\n", ERR_INPUT_F_NOT_READABLE);
    return(1);
  }
  fclose(tmprofiles_file);

  // Open temporal file for reading profiles
  // Open profiles and contigs output files for writing results
  fprintf(stderr, "[LOG] CALCULATING IRREPRODUCIBILITY SCORES\n");
  tmprofiles_file = fopen(tmprofiles_file_name, "r");
  profiles_file = fopen(profiles_file_name, "w");
  contigs_file = fopen(contigs_file_name, "w");
  if (!tmprofiles_file) {
    fprintf(stderr, "%s\n", ERR_INPUT_F_NOT_READABLE);
    return (1);
  }
  if (!profiles_file || !contigs_file) {
    fprintf(stderr, "%s\n", ERR_OUTPUT_F_NOT_WRITABLE);
    return (1);
  }

  // Generate data structures for ID
  sere_s = create_sere(reads_per_contig, ncontigs, arguments.number_replicates);
  npidr_s = create_npidr(reads_per_contig, ncontigs, arguments.number_replicates);

  // Read profiles and print results
  index = 0;
  while((result = next_tmprofile(tmprofiles_file, &profile) > 0)) {

    // Calculate irreproducibility scores
    if (arguments.idr_method == IDR_SERE) 
      calculate_sere_score(&profile, sere_s, arguments.number_replicates, arguments.idr_cutoff);
    else if (arguments.idr_method == IDR_COMMON)
      calculate_common_score(&profile, arguments.number_replicates);
    else if (arguments.idr_method == IDR_IDR)
      calculate_npidr_score(&profile, npidr_s, index, ncontigs, arguments.number_replicates, arguments.idr_cutoff);
    else
      profile.idr_score = 0;

    // Print contig
    fprintf(contigs_file, "%s\t%d\t%d\t%s", profile.chromosome, profile.start, profile.end, STR(profile.strand));
    for (i = 0; i < arguments.number_replicates; i++)
      fprintf(contigs_file, "\t%d", profile.nreads[i]);
    fprintf(contigs_file, "\t%f\n", profile.idr_score);

    // Trim
    if (profile.valid)
      trim(&profile, &arguments);

    // Internal trim if profile is too long
    if ((profile.valid) && (profile.length > arguments.max_len)) {
      double max_height = gsl_stats_max(profile.profile, 1, profile.olength);
      double threshold = max_height * arguments.trim_threshold;
      int ix, jx;
      int pstart, pend;
      int sstart, send;

      // Initialize
      pstart = -1; pend = -1;
      sstart = -1; send = -1;

      // Iterate through profile
      for (ix = profile.tstart; ix <= profile.tend; ix++) {
 
        // Profile
        if ((profile.profile[ix] >= arguments.trim_max) || ((profile.profile[ix] > arguments.trim_min && profile.profile[ix] > threshold))) {
          if (pstart == -1) {
            pstart = ix;
          }
          else if (sstart != -1) {
            send = ix - 1;
            if ((send - sstart + 1) > arguments.spacing) { // print profile part if proceeds
              max_height = -1;
              for (jx = pstart; jx <= pend; jx++)
                if (max_height < profile.profile[jx]) max_height = profile.profile[jx];
              if (((pend - pstart + 1) >= arguments.min_len) && ((pend - pstart + 1) <= arguments.max_len) && (max_height >= arguments.min_reads)) {
                int prfps, prfpe;
                if (profile.strand == FWD_STRAND) {
                  prfps = profile.start + (pstart - profile.tstart);
                  prfpe = profile.start + (pend - profile.tstart);
                }
                else {
                  prfpe = profile.end - (pstart - profile.tstart);
                  prfps = profile.end - (pend - profile.tstart);
                }
                fprintf(profiles_file, "%s:%d-%d:%s", profile.chromosome, prfps, prfpe, STR(profile.strand));
                for (jx = pstart; jx <= pend; jx++)
                  fprintf(profiles_file, "\t%f", profile.profile[jx]);
                fprintf(profiles_file, "\n");
              }
              sstart = -1; send = -1;
              pstart = ix; pend = -1;
            }
            else {
              sstart = -1; send = -1;
              pend = -1;
            }
          }
        }

        // Empty space
        else {
          if (sstart == -1) {
            sstart = ix;
            pend = ix - 1;
          }
        }
      }

      max_height = -1;
      pend = ix - 1;
      for (jx = pstart; jx <= pend; jx++)
        if (max_height < profile.profile[jx]) max_height = profile.profile[jx];
      if (((pend - pstart + 1) >= arguments.min_len) && ((pend - pstart + 1) <= arguments.max_len) && (max_height >= arguments.min_reads)) {
        int prfps, prfpe;
        if (profile.strand == FWD_STRAND) {
          prfps = profile.start + (pstart - profile.tstart);
          prfpe = profile.start + (pend - profile.tstart);
        }
        else {
          prfpe = profile.end - (pstart - profile.tstart);
          prfps = profile.end - (pend - profile.tstart);
        }
        fprintf(profiles_file, "%s:%d-%d:%s", profile.chromosome, prfps, prfpe, STR(profile.strand));
        for (jx = pstart; jx <= pend; jx++)
          fprintf(profiles_file, "\t%f", profile.profile[jx]);
        fprintf(profiles_file, "\n");
      }
    }

    // Print profile
    else if ((profile.valid) && (profile.length <= arguments.max_len)) {
      fprintf(profiles_file, "%s:%d-%d:%s", profile.chromosome, profile.start, profile.end, STR(profile.strand));
      for (i = profile.tstart; i <= profile.tend; i++)
        fprintf(profiles_file, "\t%f", profile.profile[i]);
      fprintf(profiles_file, "\n");
    }

    // Free structures in profile
    free(profile.nreads);
    if (profile.free) free(profile.profile);

    index++;
  }

  // Destroy data structures for ID
  destroy_sere(sere_s);
  destroy_npidr(npidr_s);

  // Close file descriptors
  fclose(tmprofiles_file);
  fclose(profiles_file);
  fclose(contigs_file);
  for(i = 0; i < arguments.number_replicates; i++)
    samclose(replicate_file[i]);

  // Delete temporary files
  result = unlink(tmprofiles_file_name);

  // Free pointers
  free(tmprofiles_file_name);
  free(profiles_file_name);
  free(contigs_file_name);
  for (i = 0; i < index; i++)
    free(reads_per_contig[i]);
  free(reads_per_contig);

  return(0);
}
