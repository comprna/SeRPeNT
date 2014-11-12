#include <profiles/profiles.h>

/*
 * next_alignment
 *   Retrieves one alignment from a BAM file
 *
 * @arg const samfile_t *bam_file
 *   Pointer to a BAM file descriptor
 * @arg const bam1_t *bam_alignment
 *   Pointer to a BAM alignment information structure
 * @arg alignment *alignment
 *   Pointer to an alignment handler
 *
 * @return -2 if error ocurred. -1 if EOF. 0 othwerwise.
 */
int next_alignment(samfile_t *bam_file, alignment_struct *alignment, int replica)
{
  int r;
  bam1_t *bam_alignment = bam_init1();
  
  if ((r = samread(bam_file, bam_alignment)) >= 0) {
    int32_t pos = bam_alignment->core.pos + 1;
    int32_t end = bam_alignment->core.pos;
    int32_t flag = bam_alignment->core.flag;
    char *chr = bam_file->header->target_name[bam_alignment->core.tid];
    uint32_t *cigar = bam1_cigar(bam_alignment);
    int i;
    char spliced = 0;

    // Determine if alignment is spliced by checking cigar string
    for (i = 0; i < bam_alignment->core.n_cigar; ++i) {
      char operation = bam_cigar_opchr(cigar[i]);
      
      spliced = spliced || (operation == 'N');

      if ((operation == 'D') ||
          (operation == 'M') ||
          (operation == 'X') ||
          (operation == '='))
        end += bam_cigar_oplen(cigar[i]);
    }

    // Mark the alignment as valid if:
    // - it is not spliced, and
    // - it is not paired, and
    // - it is not unmapped, and
    // - it passed qc check by aligner, and
    // - it is not a PCR or optical duplicate
    if (!(spliced)            &&
        !(flag & BAM_FPAIRED) &&
        !(flag & BAM_FUNMAP)  &&
        !(flag & BAM_FQCFAIL) &&
        !(flag & BAM_FDUP))
    {
      alignment->valid = VALID_ALIGNMENT;
      alignment->replicate = replica;
      alignment->nreads = 1;
      strncpy(alignment->chromosome, chr, MAX_FEATURE);
      alignment->start = pos;
      alignment->end = end;
      if (flag & BAM_FREVERSE)
        alignment->strand = REV_STRAND;
      else 
        alignment->strand = FWD_STRAND;
    }
    else
      alignment->valid = INVALID_ALIGNMENT;
  }

  bam_destroy1(bam_alignment);
  return(r);
}


/*
 * parse_alignment:
 *   Parse alignment and build contig accordingly
 *
 * @arg 
 */
int parse_alignment(const args_p_struct *arguments, const alignment_struct *alignment, profile_struct *current_profile, contig_struct *primary, int *index)
{
  // FIRST alignment
  if (primary->start == 0) {
    primary->start = alignment->start;
    primary->end = alignment->end;
    strncpy(primary->chromosome, alignment->chromosome, MAX_FEATURE);
    primary->profile = (double*) malloc(sizeof(double) * (primary->end - primary->start + 1));
    primary->nreads = (int*) malloc(sizeof(int) * (arguments->number_replicates));
    int i;
    if ((arguments->replicate_treat != REPLICATE_REPLICATE) || ((arguments->replicate_treat == REPLICATE_REPLICATE) && ((arguments->replicate_number - 1) == alignment->replicate)))
      for (i = 0; i < (primary->end - primary->start + 1); i++) primary->profile[i] = alignment->nreads;
    for (i = 0; i < arguments->number_replicates; i++) primary->nreads[i] = 0;
    primary->nreads[alignment->replicate] += alignment->nreads;
  }
  // NOT FIRST alignment
  else {
    // Cases:
    //       current contig    chr A  |----------------|
    //
    //       alignments        chr A        |------|
    //                         chr A  |----------|
    //                         chr A      |------------|
    //                         chr A  |----------------|
    if (strcmp(alignment->chromosome, primary->chromosome) == 0 && alignment->end <= primary->end && alignment->end > primary->start) {
      int i;
      if ((arguments->replicate_treat != REPLICATE_REPLICATE) || ((arguments->replicate_treat == REPLICATE_REPLICATE) && ((arguments->replicate_number - 1) == alignment->replicate)))
        for (i = alignment->start - primary->start; i < (alignment->end - primary->start + 1); i++) primary->profile[i] += alignment->nreads;
      primary->nreads[alignment->replicate] += alignment->nreads;
    }

    // Cases:
    //       current contig    chr A  |----------------|
    //
    //       alignments        chr A               |--------|
    //                         chr A  |---------------------|
    //                         chr A                   |----|
    //                         chr A                    spacing |-------|
    else if (strcmp(alignment->chromosome, primary->chromosome) == 0 && alignment->start <= primary->end + arguments->spacing && alignment->end > primary->start) {
      double *profile_realloc = (double*) realloc(primary->profile, sizeof(double) * (alignment->end - primary->start + 1));
      if (profile_realloc == NULL)
        return(-1);
      primary->profile = profile_realloc;
      int i;
      for (i = (primary->end - primary->start + 1); i < (alignment->end - primary->start + 1); i++) primary->profile[i] = 0;
      if ((arguments->replicate_treat != REPLICATE_REPLICATE) || ((arguments->replicate_treat == REPLICATE_REPLICATE) && ((arguments->replicate_number - 1) == alignment->replicate)))
        for (i = (alignment->start - primary->start); i < (alignment->end - primary->start + 1); i++) primary->profile[i] += alignment->nreads;
      primary->end = alignment->end;
      primary->nreads[alignment->replicate] += alignment->nreads;
    }

    // Cases:
    //       current contig    chr A             |----------------|
    //
    //       alignments        chr A                                > spacing   |------|
    //                         chr B  |--------|
    else {
      int i;

      // Store contig data
      current_profile->nreads = (int*) malloc(sizeof(int) * arguments->number_replicates);
      for (i = 0; i < arguments->number_replicates; i++) current_profile->nreads[i] = primary->nreads[i];
      strncpy(current_profile->chromosome, primary->chromosome, MAX_FEATURE);
      current_profile->start = primary->start;
      current_profile->end = primary->end;
      current_profile->length = (primary->end - primary->start + 1);
      current_profile->strand = alignment->strand;

      // Store profile data if allowed by parameters -> memory reduction
      if ((gsl_stats_max(primary->profile, 1, (primary->end - primary->start + 1)) >= arguments->min_reads) &&  // Contig has more than r reads
          ((primary->end - primary->start + 1) >= arguments->min_len))                                          // Contig has, at most, M nucleotides
      {
        current_profile->valid = 1;
        current_profile->profile = (double*) malloc(sizeof(double) * (primary->end - primary->start + 1));
        for (i = 0; i < (primary->end - primary->start + 1); i++) {
          if (arguments->replicate_treat == REPLICATE_MEAN)
            current_profile->profile[i] = primary->profile[i] / ((double) arguments->number_replicates);
          else
            current_profile->profile[i] = primary->profile[i];
        }
      }
      else
        current_profile->valid = 0;

      // Increment index pointer
      (*index)++;

      // Restart primary alignment
      primary->start = alignment->start;
      primary->end = alignment->end;
      strncpy(primary->chromosome, alignment->chromosome, MAX_FEATURE);
      double *profile_realloc = (double*) realloc(primary->profile, sizeof(double) * (primary->end - primary->start + 1));
      if (profile_realloc == NULL)
        return(1);
      primary->profile = profile_realloc;
      if ((arguments->replicate_treat != REPLICATE_REPLICATE) || ((arguments->replicate_treat == REPLICATE_REPLICATE) && ((arguments->replicate_number - 1) == alignment->replicate)))
        for (i = 0; i < (primary->end - primary->start + 1); i++) primary->profile[i] = alignment->nreads;
      for (i = 0; i < arguments->number_replicates; i++) primary->nreads[i] = 0;
      primary->nreads[alignment->replicate] += alignment->nreads;
    }
  }
  
  return(0);
}


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
 *  Function to copy values from one alignment_struct to another
 *
 * @arg alignment_struct* destiny
 *   Pointer to the alignment_struct that is updated
 *
 * @arg alignment_struct* source
 *   Pointer to the alignment_struct that will be copied
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
  args_p_struct arguments;                              // Struct for handling command line parameters
  FILE *profiles_file, *contigs_file;                   // Output file descriptors
  samfile_t* replicate_file[MAX_REPLICATES];            // Array of BAM file descriptors
  alignment_struct current_alignments[MAX_REPLICATES];  // Array of alignment_struct struct
  int result;                                           // Result of any operation
  int results[MAX_REPLICATES];                          // Replicate-specific results
  char* error_message;                                  // Error message to display in case of abnormal termination
  profile_struct* profiles;                             // Array of profile_struct structs
  contig_struct contig_fwd, contig_rev;                 // Struct for building profiles
  int i, j, stop, index, chridx, blkidx, maxstart;      // Multi-purpose indexes and checkpoint variables
  char curr_chrom[MAX_FEATURE];                         // Current chromosome
  int curr_len;                                         // Current length

  // Initialize options with default values. Parse command line.
  // Exit if command is not well-formed.
  arguments.min_len = MIN_LEN;
  arguments.max_len = MAX_LEN;
  arguments.spacing = SPACING;
  arguments.min_reads = MIN_READS;
  arguments.trimming = TRIMMING;
  arguments.replicate_treat = REPLICATE_POOL;
  arguments.replicate_number = REPLICATE_NUMBER;
  arguments.number_replicates = MAX_REPLICATES;
  arguments.idr_method = IDR_COMMON;
  arguments.idr_cutoff = CUTOFF;
  arguments.read_length = MAX_READ_LENGTH;
  if (parse_command_line_p(argc, argv, &error_message, &arguments) < 0) {
    fprintf(stderr, "%s\n", error_message);
    if ((strcmp(error_message, PROFILES_HELP_MSG) == 0) || (strcmp(error_message, VERSION_MSG) == 0))
      return(0);
    fprintf(stderr, "%s\n", ERR_PROFILES_HELP_MSG);
    return(1);
  }

  // Open output files for writing results.
  // Exit if output files do not exist or are not readable.
  char *profiles_file_name = malloc((MAX_PATH + strlen(PROFILES_SUFFIX) + 2) * sizeof(char));
  char *contigs_file_name = malloc((MAX_PATH + strlen(CONTIGS_SUFFIX) + 2) * sizeof(char));
  strncpy(profiles_file_name, arguments.output_f_path, MAX_PATH);
  strcat(profiles_file_name, PATH_SEPARATOR);
  strcat(profiles_file_name, PROFILES_SUFFIX);
  strncpy(contigs_file_name, arguments.output_f_path, MAX_PATH);
  strcat(contigs_file_name, PATH_SEPARATOR);
  strcat(contigs_file_name, CONTIGS_SUFFIX);
  profiles_file = fopen(profiles_file_name, "w");
  contigs_file = fopen(contigs_file_name, "w");
  if (!profiles_file || !contigs_file) {
    fprintf(stderr, "%s\n", ERR_OUTPUT_F_NOT_WRITABLE);
    return (1);
  }
  free(profiles_file_name);
  free(contigs_file_name);

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
  // BIG TODO -> replicates can have different number of chromosomes
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

  // Allocate memory for profiles
  profiles = (profile_struct*) malloc(MAX_CONTIGS * sizeof(profile_struct));

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
      results[i] = next_alignment(replicate_file[i], &current_alignments[i], i);
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

    alignments = (alignment_struct*) malloc (MAX_BLOCK * arguments.read_length * 2 * arguments.number_replicates * sizeof(alignment_struct));//TODO

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
          results[i] = next_alignment(replicate_file[i], &current_alignments[i], i);
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
        if (parse_alignment(&arguments, algn, &profiles[index], &contig_fwd, &index) < 0) {
          fprintf(stderr, "%s\n", ERR_REALLOC_FAILED);
          return(1);
        }
      }
      else {
        if (parse_alignment(&arguments, algn, &profiles[index], &contig_rev, &index) < 0) {
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
  profiles[index].nreads = (int*) malloc(sizeof(int) * arguments.number_replicates);
  for (i = 0; i < arguments.number_replicates; i++) profiles[index].nreads[i] = contig_fwd.nreads[i];
  strncpy(profiles[index].chromosome, contig_fwd.chromosome, MAX_FEATURE);
  profiles[index].start = contig_fwd.start;
  profiles[index].end = contig_fwd.end;
  profiles[index].length = (contig_fwd.end - contig_fwd.start + 1);
  profiles[index].strand = FWD_STRAND;

  // Flush the profile contents in the forward contig struct if allowed by parameters
  if ((gsl_stats_max(contig_fwd.profile, 1, (contig_fwd.end - contig_fwd.start + 1)) >= arguments.min_reads) &&  // Contig has more than r reads
      ((contig_fwd.end - contig_fwd.start + 1) >= arguments.min_len))                                            // Contig has, at most, M nucleotides
  {
    profiles[index].valid = 1;
    profiles[index].profile = (double*) malloc(sizeof(double) * (contig_fwd.end - contig_fwd.start + 1));
    for (i = 0; i < (contig_fwd.end - contig_fwd.start + 1); i++) {
      if (arguments.replicate_treat == REPLICATE_MEAN)
        profiles[index].profile[i] = contig_fwd.profile[i] / ((double) arguments.number_replicates);
      else
        profiles[index].profile[i] = contig_fwd.profile[i];
    }
  }
  else
    profiles[index].valid = 0;
  index++;

  // Flush the contents in the reverse contig struct
  profiles[index].nreads = (int*) malloc(sizeof(int) * arguments.number_replicates);
  for (i = 0; i < arguments.number_replicates; i++) profiles[index].nreads[i] = contig_rev.nreads[i];
  strncpy(profiles[index].chromosome, contig_rev.chromosome, MAX_FEATURE);
  profiles[index].start = contig_rev.start;
  profiles[index].end = contig_rev.end;
  profiles[index].length = (contig_rev.end - contig_rev.start + 1);
  profiles[index].strand = REV_STRAND;

  // Flush the profile contents in the reverse contig struct if allowed by parameters
  if ((gsl_stats_max(contig_rev.profile, 1, (contig_rev.end - contig_rev.start + 1)) >= arguments.min_reads) &&  // Contig has more than r reads
      ((contig_rev.end - contig_rev.start + 1) >= arguments.min_len))                                            // Contig has, at most, M nucleotides
  {
    profiles[index].valid = 1;
    profiles[index].profile = (double*) malloc(sizeof(double) * (contig_rev.end - contig_rev.start + 1));
    for (i = 0; i < (contig_rev.end - contig_rev.start + 1); i++) {
      if (arguments.replicate_treat == REPLICATE_MEAN)
        profiles[index].profile[i] = contig_rev.profile[i] / ((double) arguments.number_replicates);
      else
        profiles[index].profile[i] = contig_rev.profile[i];
    }
  }
  else
    profiles[index].valid = 0;
  index++;

  // Parse all the contigs and calculate irreproducibility scores
  fprintf(stderr, "[LOG] CALCULATING IRREPRODUCIBILITY SCORES\n");
  if (arguments.idr_method == IDR_SERE)
    calculate_sere_scores(profiles, index, arguments.number_replicates, arguments.idr_cutoff);
  else if (arguments.idr_method == IDR_COMMON)
    calculate_common_scores(profiles, index, arguments.number_replicates);
  else if (arguments.idr_method == IDR_IDR)
    calculate_npidr_scores(profiles, index, arguments.number_replicates, arguments.idr_cutoff);

  // Print contigs
  fprintf(stderr, "[LOG] FILTERING AND PRINTING PROFILES\n");
  for(i = 0; i < index; i++) {
    if (profiles[i].strand == FWD_STRAND)
      fprintf(contigs_file, "%s\t%d\t%d\t+", profiles[i].chromosome, profiles[i].start, profiles[i].end);
    else
      fprintf(contigs_file, "%s\t%d\t%d\t-", profiles[i].chromosome, profiles[i].start, profiles[i].end);

    for (j = 0; j < arguments.number_replicates; j++)
      fprintf(contigs_file, "\t%d", profiles[i].nreads[j]);

    fprintf(contigs_file, "\t%f\n", profiles[i].idr_score);
  }

  // Print profiles
  for (i = 0; i < index; i++) {
    int startfp, endtp;

    // Trim 5' end of the profile
    if (profiles[i].valid) {
      j = 0; stop = 0;
      while (!stop) {
        if (profiles[i].profile[j] > arguments.trimming) {
          stop++;
          startfp = j;
        }
        else if (j == profiles[i].end - profiles[i].start) stop++;
        else j++;
      }

      // Trim 3' end of the profile
      j = profiles[i].end - profiles[i].start; stop = 0;
      while (!stop) {
        if (profiles[i].profile[j] > arguments.trimming) {
          stop++;
          endtp = j;
        }
        else if (j == 0) stop++;
        else j--;
      }

      profiles[i].length = endtp - startfp + 1;
      profiles[i].valid = ((profiles[i].length >= arguments.min_len) && (profiles[i].length <= arguments.max_len));
    }

    // Print profile if allowed by parameters
    if (profiles[i].valid) {
      if (profiles[i].strand == FWD_STRAND)
        fprintf(profiles_file, "%s:%d-%d:+", profiles[i].chromosome, profiles[i].start + startfp, profiles[i].start + endtp);
      else
        fprintf(profiles_file, "%s:%d-%d:-", profiles[i].chromosome, profiles[i].start + startfp, profiles[i].start + endtp);

      if (profiles[i].strand == FWD_STRAND)
        for (j = startfp; j <= endtp; j++) fprintf(profiles_file, "\t%f", profiles[i].profile[j]);
      else
        for (j = endtp; j >= startfp; j--) fprintf(profiles_file, "\t%f", profiles[i].profile[j]);

      fprintf(profiles_file, "\n");
    }
  }

  // Close file descriptors and free
  free(profiles);
  for(i = 0; i < arguments.number_replicates; i++)
    samclose(replicate_file[i]);
  fclose(profiles_file);
  fclose(contigs_file);
  return(0);
}
/* END OF FILE */
