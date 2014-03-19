#include <sdrnaprofiles.h>


/*
 * parse_command_line:
 *   Parses the command line
 *
 * @arg int argc
 *   Number of arguments in the command line
 * @arg char** argv
 *   Array containing the command line arguments
 * @arg error_message
 *   Pointer to a char array where to store the error message
 * @arg args_struct* arguments
 *   Pointer to a the argument handler
 *
 * @return -1 if an error occurred. 0 otherwise.
 */
int parse_command_line(int argc, char** argv, char** error_message, args_struct* arguments)
{
  opterr = 0;
  char carg;
  int terminate = 0;

  while((carg = getopt(argc, argv, "x:m:M:s:r:t:")) != -1) {
    switch (carg) {
      case 'x':
        arguments->cutoff = (double) atof(optarg);
        if (arguments->cutoff <= (double) 0 || arguments->cutoff >= (double) 1) {
          terminate--;
          *error_message = ERR_INVALID_x_VALUE;
        }
        break;
      case 'm':
        arguments->min_len = atoi (optarg);
        if (arguments->min_len < 5) {
          terminate--;
          *error_message = ERR_INVALID_m_VALUE;
        }
        break;
      case 'M':
        arguments->max_len = atoi(optarg);
        if (arguments->max_len < 6) {
          terminate--;
          *error_message = ERR_INVALID_M_VALUE;
        }
        break;
      case 's':
        arguments->spacing = atoi(optarg);
        if (arguments->spacing < 0) {
          terminate--;
          *error_message = ERR_INVALID_s_VALUE;
        }
        break;
      case 'r':
        arguments->min_reads = (double) atof(optarg);
        if (arguments->min_reads < (double) 0) {
          terminate--;
          *error_message = ERR_INVALID_r_VALUE;
        }
        break;
      case 't':
        arguments->trimming = atoi(optarg);;
        if (arguments->trimming < 0) {
          terminate--;
          *error_message = ERR_INVALID_t_VALUE;
        }
        break;
      case '?':
        terminate--;
        *error_message = ERR_INVALID_ARGUMENT;
    }
  }
  if (!terminate && arguments->max_len <= arguments->min_len) {
    terminate--;
    *error_message = ERR_INVALID_LENGTH;
  }
  else if (!terminate && (argc - optind) != (MAX_REPLICATES + 1)) {
    terminate--;
    *error_message = ERR_INVALID_NUMBER_ARGUMENTS;
  }
  else if (!terminate) {
    int i, j = 0;
    for (i = optind; i < argc - 1; i++)
      strncpy(arguments->replicate_f_path[j++], argv[i], MAX_PATH);
    strncpy(arguments->output_f_path, argv[argc - 1], MAX_PATH);
  }
  return(terminate);  
}


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
int next_alignment(samfile_t *bam_file, alignment_struct *alignment)
{
  int r;
  bam1_t *bam_alignment = bam_init1();
  
  if ((r = samread(bam_file, bam_alignment)) >= 0) {
    int32_t length = bam_alignment->core.l_qseq;
    int32_t pos = bam_alignment->core.pos + 1;
    int32_t end = bam_alignment->core.pos;
    int32_t flag = bam_alignment->core.flag;
    char *chr = bam_file->header->target_name[bam_alignment->core.tid];
    uint32_t *cigar = bam1_cigar(bam_alignment);
    uint32_t clipping_start = 0, clipping_end = 0;
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
      alignment->chromosome = chr;
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
int parse_alignment(const args_struct *arguments, const alignment_struct *alignment, profile_struct *current_profile, contig_struct *primary, int *index, FILE *output)
{
  // FIRST alignment
  if (primary->start == 0) {
    primary->start = alignment->start;
    primary->end = alignment->end;
    strncpy(primary->chromosome, alignment->chromosome, MAX_FEATURE);
    primary->profile = (double*) malloc(sizeof(double) * (primary->end - primary->start + 1));
    int i;
    for (i = 0; i < (primary->end - primary->start + 1); i++) primary->profile[i] = 1;
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
    if (alignment->end <= primary->end && alignment->end > primary->start) {
      int i;
      for (i = alignment->start - primary->start; i < (alignment->end - primary->start + 1); i++) primary->profile[i]++;
    }

    // Cases:
    //       current contig    chr A  |----------------|
    //
    //       alignments        chr A               |--------|
    //                         chr A  |---------------------|
    //                         chr A                   |----|
    //                         chr A                    spacing |-------|
    else if (alignment->start <= primary->end + arguments->spacing && alignment->end > primary->start) {
      double *profile_realloc = (double*) realloc(primary->profile, sizeof(double) * (alignment->end - primary->start + 1));
      if (profile_realloc == NULL)
        return(-1);
      primary->profile = profile_realloc;
      int i;
      for (i = (primary->end - primary->start + 1); i < (alignment->end - primary->start + 1); i++) primary->profile[i] = 0;
      for (i = (alignment->start - primary->start); i < (alignment->end - primary->start + 1); i++) primary->profile[i]++;
      primary->end = alignment->end; 
    }

    // Cases:
    //       current contig    chr A             |----------------|
    //
    //       alignments        chr A                                > spacing   |------|
    //                         chr B  |--------|
    else {
      int startfp, endtp, stop, i;

      // Trim 5'-end
      i = 0; stop = 0;
      while (!stop) {
        if (primary->profile[i] > arguments->trimming) {
          stop++;
          startfp = i;
        }
        else if (i == primary->end - primary->start) stop++;
        else i++;
      }

      // Trim 3'-end
      i = primary->end - primary->start; stop = 0;
      while (!stop) {
        if (primary->profile[i] > arguments->trimming) {
          stop++;
          endtp = i;
        }
        else if (i == 0) stop++;
        else i--;
      }

      // Load profile if allowed by parameters
      if ((gsl_stats_max(primary->profile, 1, (primary->end - primary->start + 1)) >= arguments->min_reads) &&  // Contig has more than r reads
          ((endtp - startfp + 1) >= arguments->min_len)                                       &&  // Contig has, at least, m nucleotides
          ((endtp - startfp + 1) <= arguments->max_len))                                          // Contig has, at most, M nucleotides
      {
        current_profile->profile = (double*) malloc(sizeof(double) * (endtp - startfp + 1));
        strncpy(current_profile->chromosome, primary->chromosome, MAX_FEATURE);
        current_profile->start = primary->start + startfp;
        current_profile->end = primary->start + endtp;
        current_profile->length = (endtp - startfp + 1);
        current_profile->strand = alignment->strand;
        if (alignment->strand == FWD_STRAND) {
          fprintf(output, "%s:%d-%d:+", primary->chromosome, primary->start + startfp, primary->start + endtp);
          for (i = 0; i < (endtp - startfp + 1); i++) {
            current_profile->profile[i] = primary->profile[i + startfp];
            fprintf(output, "\t%f", primary->profile[i + startfp]);
          }
          fprintf(output, "\n");
        }
        else {
          fprintf(output, "%s:%d-%d:-", primary->chromosome, primary->start + startfp, primary->start + endtp);
          int j = 0;
          for (i = endtp; i >= startfp; i--) {
            current_profile->profile[j] = primary->profile[i];
            j++;
            fprintf(output, "\t%f", primary->profile[i]);
          }
          fprintf(output, "\n");
        } 
        (*index)++;
      }
      primary->start = alignment->start;
      primary->end = alignment->end;
      strncpy(primary->chromosome, alignment->chromosome, MAX_FEATURE);
      double *profile_realloc = (double*) realloc(primary->profile, sizeof(double) * (primary->end - primary->start + 1));
      if (profile_realloc == NULL)
        return(1);
      primary->profile = profile_realloc;
      for (i = 0; i < (primary->end - primary->start + 1); i++) primary->profile[i] = 1;
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
 * Application entry point
 */
int main(int argc, char** argv)
{
  // Define and declare variables
  args_struct arguments;                                // Struct for handling command line parameters
  FILE *profiles_file, *crosscor_file;                  // Output file descriptors
  samfile_t* replicate_file[MAX_REPLICATES];            // Array of BAM file descriptors
  alignment_struct alignment;                           // Structure containing useful information for one alignment
  int index;                                            // Index for filling coverage vector
  int result;                                           // Result of any operation
  char* error_message;                                  // Error message to display in case of abnormal termination
  profile_struct profiles[MAX_CONTIGS];                 // Array of profile_struct structs
  contig_struct contig_fwd, contig_rev;                 // Struct for building profiles
  int i, j, stop;                                       // Multi-purpose indexes and checkpoint variables
  int startfp, endtp;                                   // Auxiliar variables for trimming

  // Initialize options with default values. Parse command line.
  // Exit if command is not well-formed.
  arguments.cutoff = CUTOFF;
  arguments.min_len = MIN_LEN;
  arguments.max_len = MAX_LEN;
  arguments.spacing = SPACING;
  arguments.min_reads = MIN_READS;
  if (parse_command_line(argc, argv, &error_message, &arguments) < 0) {
    MSG_PRINT(1, error_message, ERR_MSG);
    DEBUG_PRINT(error_message);
    return(1);
  }
  
  // Open output files for writing results.
  // Exit if output files do not exist or are not readable.
  char *profiles_file_name = malloc((MAX_PATH + strlen(PROFILES_SUFFIX) + 2) * sizeof(char));
  char *crosscor_file_name = malloc((MAX_PATH + strlen(CROSSCOR_SUFFIX) + 2) * sizeof(char));
  strncpy(profiles_file_name, arguments.output_f_path, MAX_PATH);
  strcat(profiles_file_name, PATH_SEPARATOR);
  strcat(profiles_file_name, PROFILES_SUFFIX);
  strncpy(crosscor_file_name, arguments.output_f_path, MAX_PATH);
  strcat(crosscor_file_name, PATH_SEPARATOR);
  strcat(crosscor_file_name, CROSSCOR_SUFFIX);
  profiles_file = fopen(profiles_file_name, "w");
  crosscor_file = fopen(crosscor_file_name, "w");
  if (!profiles_file || !crosscor_file) {
    MSG_PRINT(0, ERR_OUTPUT_F_NOT_WRITABLE, "");
    DEBUG_PRINT(ERR_OUTPUT_F_NOT_WRITABLE);
    return (1);
  }
  free(profiles_file_name);
  free(crosscor_file_name);

  // Open replicate file for reading and read BAM header
  // Exit if replicate BAM files do not exist, are not readable or do not have header
  for (index = 0; index < MAX_REPLICATES; index++) {
    if ((replicate_file[index] = samopen(arguments.replicate_f_path[index], "rb", 0)) == 0) {
      MSG_PRINT(0, ERR_BAM_F_NOT_READABLE, "");
      DEBUG_PRINT(ERR_BAM_F_NOT_READABLE);
      return(1);
    }
    if (replicate_file[index]->header == 0) {
      MSG_PRINT(0, ERR_BAM_F_NOT_HEADER, "");
      DEBUG_PRINT(ERR_BAM_F_NOT_HEADER);
      return(1);
    }
  }
  
  // Read BAM file and build sRNA profiles
  // Exit if BAM files are truncated or ill-formed, or if not enough memory
  fprintf(stderr, "Generating contigs... ");
  contig_fwd.start = 0; contig_fwd.end = 0; contig_rev.start = 0; contig_rev.end = 0;
  index = 0;
  while ((result = next_alignment(replicate_file[0], &alignment)) > -1) {
    if (alignment.valid) {
      if (alignment.strand == FWD_STRAND) {
        if (parse_alignment(&arguments, &alignment, &profiles[index], &contig_fwd, &index, profiles_file) < 0) {
          MSG_PRINT(0, ERR_REALLOC_FAILED, "");
          DEBUG_PRINT(ERR_REALLOC_FAILED);
          return(1);
        }
      }
      else {
        if (parse_alignment(&arguments, &alignment, &profiles[index], &contig_rev, &index, profiles_file) < 0) {
          MSG_PRINT(0, ERR_REALLOC_FAILED, "");
          DEBUG_PRINT(ERR_REALLOC_FAILED);
          return(1);
        }
      }
    }
  }
  if (result < -1) {
    MSG_PRINT(0, ERR_BAM_F_TRUNCATED, "");
    DEBUG_PRINT(ERR_BAM_F_TRUNCATED);
    return(1);
  }

  // Trim 5'-end of the forward contig struct
  i = 0; stop = 0;
  while (!stop) {
    if (contig_fwd.profile[i] > arguments.trimming) {
      stop++;
      startfp = i;
    }
    else if (i == contig_fwd.end - contig_fwd.start) stop++;
    else i++;
  }

  // Trim 3'-end of the forward contig struct
  i = contig_fwd.end - contig_fwd.start; stop = 0;
  while (!stop) {
    if (contig_fwd.profile[i] > arguments.trimming) {
      stop++;
      endtp = i;
    }
    else if (i == 0) stop++;
    else i--;
  }

  // Flush the contents in the forward contig struct
  if ((gsl_stats_max(contig_fwd.profile, 1, (contig_fwd.end - contig_fwd.start + 1)) >= arguments.min_reads) &&  // Contig has more than r reads
      ((endtp - startfp + 1) >= arguments.min_len)                                                           &&  // Contig has, at least, m nucleotides
      ((endtp - startfp + 1) <= arguments.max_len))                                                              // Contig has, at most, M nucleotides
  {
    profiles[index].profile = (double*) malloc(sizeof(double) * (endtp - startfp + 1));
    strncpy(profiles[index].chromosome, contig_fwd.chromosome, MAX_FEATURE);
    profiles[index].start = contig_fwd.start + startfp;
    profiles[index].end = contig_fwd.start + endtp;
    profiles[index].length = (endtp - startfp + 1);
    profiles[index].strand = alignment.strand;
    fprintf(profiles_file, "%s:%d-%d:+", contig_fwd.chromosome, contig_fwd.start + startfp, contig_fwd.start + endtp);
    for (i = 0; i < (endtp - startfp + 1); i++) {
      profiles[index].profile[i] = contig_fwd.profile[i + startfp];
      fprintf(profiles_file, "\t%f", contig_fwd.profile[i + startfp]);
    }
    fprintf(profiles_file, "\n");
    index++;
  }

  // Trim 5'-end of the reverse contig struct
  i = 0; stop = 0;
  while (!stop) {
    if (contig_rev.profile[i] > arguments.trimming) {
      stop++;
      startfp = i;
    }
    else if (i == contig_rev.end - contig_rev.start) stop++;
    else i++;
  }

  // Trim 3'-end of the reverse contig struct
  i = contig_rev.end - contig_rev.start; stop = 0;
  while (!stop) {
    if (contig_rev.profile[i] > arguments.trimming) {
      stop++;
      endtp = i;
    }
    else if (i == 0) stop++;
    else i--;
  }

  // Flush the contents in the reverse contig struct
  if ((gsl_stats_max(contig_rev.profile, 1, (contig_rev.end - contig_rev.start + 1)) >= arguments.min_reads) &&  // Contig has more than r reads
      ((endtp - startfp + 1) >= arguments.min_len)                                                           &&  // Contig has, at least, m nucleotides
      ((endtp - startfp + 1) <= arguments.max_len))                                                              // Contig has, at most, M nucleotides
  {
    profiles[index].profile = (double*) malloc(sizeof(double) * (endtp - startfp + 1));
    strncpy(profiles[index].chromosome, contig_rev.chromosome, MAX_FEATURE);
    profiles[index].start = contig_rev.start + startfp;
    profiles[index].end = contig_rev.start + endtp;
    profiles[index].length = (endtp - startfp + 1);
    profiles[index].strand = alignment.strand;
    fprintf(profiles_file, "%s:%d-%d:-", contig_rev.chromosome, contig_rev.start + startfp, contig_rev.start + endtp);
    j = 0;
    for (i = endtp; i >= startfp; i--) {
      profiles[index].profile[j] = contig_rev.profile[i];
      j++;
      fprintf(profiles_file, "\t%f", contig_rev.profile[i]);
    }
    fprintf(profiles_file, "\n");
    index++;
  }

  // Close file descriptors
  samclose(replicate_file[0]);
  fclose(profiles_file);
  fprintf(stderr, "        [OK] => %d contigs\n", index);

  // Calculate cross-correlations for pairwise signals
  fprintf(stderr, "Calculating X-correlations... ");
  for (i = 0; i < index; i++) {
    for (j = 0; j < index; j++) {

      char istrand, jstrand;
      if (profiles[i].strand == FWD_STRAND)
        istrand = '+';
      else
        istrand = '-';
      if (profiles[j].strand == FWD_STRAND)
        jstrand = '+';
      else
        jstrand = '-';

      if (j == i) {
        fprintf(crosscor_file, "%s:%d-%d:%c\t%s:%d-%d:%c\t1.0\n", profiles[i].chromosome, profiles[i].start, profiles[i].end, istrand, profiles[j].chromosome, profiles[j].start, profiles[j].end, jstrand);
      }
      else if (j > i) {
        double* corr;
        if (profiles[i].length >= profiles[j].length)
          corr = (double*) malloc(profiles[i].length * sizeof(double));
        else
          corr = (double*) malloc(profiles[j].length * sizeof(double));
        int max_index = nxcorr(corr, profiles[i].profile, profiles[i].length, profiles[j].profile, profiles[j].length);
        fprintf(crosscor_file, "%s:%d-%d:%c\t%s:%d-%d:%c\t%f\n", profiles[i].chromosome, profiles[i].start, profiles[i].end, istrand, profiles[j].chromosome, profiles[j].start, profiles[j].end, jstrand, corr[max_index]);
        free(corr);
      }
      else if (j < i) {
        double* corr;
        if (profiles[j].length >= profiles[i].length)
          corr = (double*) malloc(profiles[j].length * sizeof(double));
        else
          corr = (double*) malloc(profiles[i].length * sizeof(double));
        int max_index = nxcorr(corr, profiles[j].profile, profiles[j].length, profiles[i].profile, profiles[i].length);
        fprintf(crosscor_file, "%s:%d-%d:%c\t%s:%d-%d:%c\t%f\n", profiles[i].chromosome, profiles[i].start, profiles[i].end, istrand, profiles[j].chromosome, profiles[j].start, profiles[j].end, jstrand, corr[max_index]);
        free(corr);
      }
    }
  }
  fprintf(stderr, "[OK]\n");

  // Close output file and exit
  fclose(crosscor_file);
  DEBUG_PRINT(SUCCESS);
  return(0);
}
/* END OF FILE */
