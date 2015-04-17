#include <profiles/trimming.h>

/*
 * trim
 *
 * @see src/include/profiles/trimming.h
 */
void trim(profile_struct* profile, args_p_struct* arguments, double** fp, double** tp)
{
  int i, startfp, endtp, stop;
  double max_height, threshold;

  // Calculate maximum height
  max_height = gsl_stats_max(profile->profile, 1, profile->length);
  threshold = max_height * arguments->trim_threshold;

  // Trim 5' end of the profile
  i = 0; stop = 0; startfp = i;
  while (!stop) {
    if ((profile->profile[i] >= arguments->trim_max) ||
        ((profile->profile[i] > arguments->trim_min) && (profile->profile[i] > arguments->trim_threshold))) {
      stop++;
      startfp = i;
    }
    else if (i == profile->end - profile->start) stop++;
    else i++;
  }

  // Trim 3' end of the profile
  i = profile->end - profile->start; stop = 0; endtp = i;
  while (!stop) {
    if ((profile->profile[i] >= arguments->trim_max) ||
        ((profile->profile[i] > arguments->trim_min) && (profile->profile[i] > arguments->trim_threshold))) {
      stop++;
      endtp = i;
    }
    else if (i == 0) stop++;
    else i--;
  }

  profile->length = endtp - startfp + 1;
  profile->start = profile->start + startfp;
  profile->end = profile->start + endtp;
  profile->valid = ((profile->length >= arguments->min_len) && (profile->length <= arguments->max_len));
}
