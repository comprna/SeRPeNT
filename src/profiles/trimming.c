#include <profiles/trimming.h>

/*
 * trim
 *
 * @see src/include/profiles/trimming.h
 */
void trim(profile_struct* profile, args_p_struct* arguments)
{
  int i, j, startfp, endtp, stop;
  double max_height, threshold;

  // Calculate maximum height
  max_height = gsl_stats_max(profile->profile, 1, profile->length);
  threshold = max_height * arguments->trim_threshold;

  // Trim 5' end of the profile
  i = 0; stop = 0; startfp = i;
  while (!stop) {
    if ((profile->profile[i] >= arguments->trim_max) ||
        ((profile->profile[i] > arguments->trim_min) && (profile->profile[i] > threshold))) {
      stop++;
      startfp = i;
    }
    else if (i == profile->end - profile->start) stop++;
    else i++;
  }

  // Trim 3' end of the profile
  i = profile->end - profile->start; stop = 0; endtp = i; j = 0;
  while (!stop) {
    if ((profile->profile[i] >= arguments->trim_max) ||
        ((profile->profile[i] > arguments->trim_min) && (profile->profile[i] > threshold))) {
      stop++;
      endtp = i;
    }
    else if (i == 0) stop++;
    else {
      i--;
      j++;
    }
  }

  profile->tstart = startfp;
  profile->tend = endtp;
  profile->length = profile->tend - profile->tstart + 1;

  if (profile->strand == FWD_STRAND) {
    profile->end = profile->start + endtp;
    profile->start = profile->start + startfp;
  }
  else {
    profile->start = profile->start + j;
    profile->end = profile->end - startfp;
  }

  profile->valid = ((profile->length >= arguments->min_len) && (profile->length <= arguments->max_len) && (max_height >= arguments->min_reads));
}
