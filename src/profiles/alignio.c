#include <profiles/alignio.h>


/*
 * next_alignment
 *
 * @see include/profiles/alignio.h
 */
int next_alignment(samfile_t *bam_file, alignment_struct *alignment, int replica, args_p_struct *arguments)
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
    // - it has the minimum length specified in the arguments
    if (!(spliced)            &&
        !(flag & BAM_FPAIRED) &&
        !(flag & BAM_FUNMAP)  &&
        !(flag & BAM_FQCFAIL) &&
        !(flag & BAM_FDUP)    &&
        ((end - pos + 1) >= arguments->min_read_len))
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
 * parse_alignment
 *
 * @see include/profiles/alignio.h
 */
int parse_alignment(const args_p_struct *arguments, const alignment_struct *alignment, contig_struct *primary, FILE* tmprofiles_file)
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
    else if (strcmp(alignment->chromosome, primary->chromosome) == 0 && alignment->start <= (primary->end + arguments->spacing + 1) && alignment->end > primary->start) {
      double *profile_realloc = (double*) realloc(primary->profile, sizeof(double) * (alignment->end - primary->start + 1));
      if (profile_realloc == NULL) {
        free(primary->profile);
        return(-1);
      }
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
      fprintf(tmprofiles_file, "%s", primary->chromosome);
      fprintf(tmprofiles_file, "\t%d", primary->start);
      fprintf(tmprofiles_file, "\t%d", primary->end);
      fprintf(tmprofiles_file, "\t%d", alignment->strand);
      fprintf(tmprofiles_file, "\t%d", arguments->number_replicates);
      for (i = 0; i < arguments->number_replicates; i++) fprintf(tmprofiles_file, "\t%d", primary->nreads[i]);

      // Store profile data if allowed by parameters -> memory reduction
      if ((gsl_stats_max(primary->profile, 1, (primary->end - primary->start + 1)) >= arguments->min_reads) &&  // Contig has more than r reads
          ((primary->end - primary->start + 1) >= arguments->min_len))                                          // Contig has, at most, M nucleotides
      {
        fprintf(tmprofiles_file, "\t%d", primary->end - primary->start + 1);
        for (i = 0; i < (primary->end - primary->start + 1); i++) {
          int idx = i;
          if (alignment->strand == REV_STRAND)
            idx = primary->end - primary->start - i;
          if (arguments->replicate_treat == REPLICATE_MEAN)
            fprintf(tmprofiles_file, "\t%f", primary->profile[idx] / ((double) arguments->number_replicates));
          else
            fprintf(tmprofiles_file, "\t%f", primary->profile[idx]);
        }
        fprintf(tmprofiles_file, "\n");
      }
      else
        fprintf(tmprofiles_file, "\t%d\n", 0);

      // Restart primary alignment
      primary->start = alignment->start;
      primary->end = alignment->end;
      strncpy(primary->chromosome, alignment->chromosome, MAX_FEATURE);
      double *profile_realloc = (double*) realloc(primary->profile, sizeof(double) * (primary->end - primary->start + 1));
      if (profile_realloc == NULL) {
        free(primary->profile);
        return(-1);
      }
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
 * next_tmprofile
 *
 * @see include/profiles/alignio.h
 */
int next_tmprofile(FILE* fp, profile_struct* profile)
{
  char *line = NULL;
  char *token;
  size_t len = 0;
  ssize_t read;
  int i, nreads, length;

  profile->free = 0;
  profile->valid = 0;

  read = getline(&line, &len, fp);

  if (read < 0) {
    free(line);
    return(0);
  }

  if ((token = strtok(line, "\t")) == NULL) {
    free(line);
    return(-1);
  }
  strcpy(profile->chromosome, token);

  if ((token = strtok(NULL, "\t")) == NULL) {
    free(line);
    return(-1);
  }
  profile->start = atoi(token);

  if ((token = strtok(NULL, "\t")) == NULL) {
    free(line);
    return(-1);
  }
  profile->end = atoi(token);
  profile->length = profile->end - profile->start + 1;
  profile->olength = profile->length;

  if ((token = strtok(NULL, "\t")) == NULL) {
    free(line);
    return(-1);
  }
  profile->strand = atoi(token);

  if ((token = strtok(NULL, "\t")) == NULL) {
    free(line);
    return(-1);
  }
  nreads = atoi(token);

  profile->nreads = (int*) malloc(nreads * sizeof(int));
  for (i = 0; i < nreads; i++) {
    if ((token = strtok(NULL, "\t")) != NULL)
      profile->nreads[i] = atoi(token);
    else {
      free(line);
      return(-1);
    }
  }

  if ((token = strtok(NULL, "\t")) == NULL) {
    free(line);
    return(-1);
  }
  length = atoi(token);

  if (length == 0) {
    free(line);
    profile->valid = 0;
    return(1);
  }

  profile->profile = (double*) malloc((profile->end - profile->start + 1) * sizeof(double));
  profile->free = 1;
  for (i = 0; i < (profile->end - profile->start + 1); i++) {
    if ((token = strtok(NULL, "\t")) != NULL)
      profile->profile[i] = (double) atof(token);
    else {
      free(line);
      return(-1);
    }
  }

  profile->valid = 1;

  free(line);
  return(1);
}
