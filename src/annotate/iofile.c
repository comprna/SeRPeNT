#include <annotate/iofile.h>

/*
 * wcl
 * 
 * @see src/include/annotate/iofile.h
 */
int wcl (FILE *fp)
{
  int c;        
  int newline_count = 0;

  while ((c=fgetc(fp)) != EOF)
    if ( c == '\n' )
      newline_count++;

  return(newline_count++);
}

/*
 * next_profile
 *
 * @see src/include/annotate/iofile.h
 */
int next_profile(FILE* fp, profile_struct* profile)
{
  char *line = NULL;
  char *token, *feature, *cline;
  size_t len = 0;
  ssize_t read;
  int i;

  read = getline(&line, &len, fp);

  if (read < 0)
    return(0);

  cline = (char*) malloc((strlen(line) + 1) * sizeof(char));
  strcpy(cline, line);

  if ((token = strtok(line, "\t")) == NULL) {
    free(cline);
    return(-1);
  }

  if ((feature = strtok(token, ":")) == NULL)  {
    free(cline);
    return(-1);
  }
  strcpy(profile->chromosome, feature);

  if ((feature = strtok(NULL, "-")) == NULL)  {
    free(cline);
    return(-1);
  }
  profile->start = atoi(feature);

  if ((feature = strtok(NULL, ":")) == NULL)  {
    free(cline);
    return(-1);
  }
  profile->end = atoi(feature);
  profile->length = profile->end - profile->start + 1;

  if ((feature = strtok(NULL, ":")) == NULL)  {
    free(cline);
    return(-1);
  }

  if (strcmp(feature, "+") == 0)
    profile->strand = FWD_STRAND;
  else if (strcmp(feature, "-") == 0)
    profile->strand = REV_STRAND;
  else {
    free(cline);
    return(-1);
  }

  profile->additional = 0;

  profile->anscore = INT_MIN;

  profile->profile = (double*) malloc((profile->end - profile->start + 1) * sizeof(double));

  if ((token = strtok(cline, "\t")) == NULL) {
    free(cline);
    return(-1);
  }

  for (i = 0; i < (profile->end - profile->start + 1); i++) {
    if ((token = strtok(NULL, "\t")) != NULL)
      profile->profile[i] = (double) atof(token);
    else {
      free(cline);
      return(-1);
    }
  }

  profile->max_height = gsl_stats_max(profile->profile, 1, profile->length);
  profile->mean = gsl_stats_mean(profile->profile, 1, profile->length);
  profile->variance = gsl_stats_variance(profile->profile, 1, profile->length);

  free(cline);
  return(1);
}

/*
 * next_additional_profile
 *
 * @see src/include/annotate/iofile.h
 */
int next_additional_profile(FILE* fp, profile_struct* profile, char* species)
{
  char *line = NULL;
  char *token, *feature, *cline;
  size_t len = 0;
  ssize_t read;
  int i;

  read = getline(&line, &len, fp);

  if (read < 0)
    return(0);

  cline = (char*) malloc((strlen(line) + 1) * sizeof(char));
  strcpy(cline, line);

  if ((token = strtok(line, "\t")) == NULL) {
    free(cline);
    return(-1);
  }

  if ((feature = strtok(token, ":")) == NULL)  {
    free(cline);
    return(-1);
  }
  strcpy(profile->chromosome, feature);

  if ((feature = strtok(NULL, "-")) == NULL)  {
    free(cline);
    return(-1);
  }
  profile->start = atoi(feature);

  if ((feature = strtok(NULL, ":")) == NULL)  {
    free(cline);
    return(-1);
  }
  profile->end = atoi(feature);
  profile->length = profile->end - profile->start + 1;

  if ((feature = strtok(NULL, ":")) == NULL)  {
    free(cline);
    return(-1);
  }

  if (strcmp(feature, "+") == 0)
    profile->strand = FWD_STRAND;
  else if (strcmp(feature, "-") == 0)
    profile->strand = REV_STRAND;
  else {
    free(cline);
    return(-1);
  }

  profile->additional = 1;
  strcpy(profile->species, species);

  profile->anscore = INT_MIN;

  profile->profile = (double*) malloc((profile->end - profile->start + 1) * sizeof(double));

  if ((token = strtok(cline, "\t")) == NULL) {
    free(cline);
    return(-1);
  }

  for (i = 0; i < (profile->end - profile->start + 1); i++) {
    if ((token = strtok(NULL, "\t")) != NULL)
      profile->profile[i] = (double) atof(token);
    else {
      free(cline);
      return(-1);
    }
  }

  profile->max_height = gsl_stats_max(profile->profile, 1, profile->length);
  profile->mean = gsl_stats_mean(profile->profile, 1, profile->length);
  profile->variance = gsl_stats_variance(profile->profile, 1, profile->length);

  free(cline);
  return(1);
}


/*
 * next_feature
 *
 * @see include/annotate/iofile.h
 */
int next_feature(FILE* bedf, feature_struct* feature)
{
  char *line = NULL;
  char *token, *index;
  size_t len = 0;
  ssize_t read;

  read = getline(&line, &len, bedf);

  if (read < 0)
    return(0);

  if ((token = strtok(line, "\t")) == NULL) 
    return(-1);
  strncpy(feature->chromosome, token, MAX_FEATURE);

  if ((token = strtok(NULL, "\t")) == NULL)
    return(-1);
  feature->start = atoi(token);

  if ((token = strtok(NULL, "\t")) == NULL)
    return(-1);
  feature->end = atoi(token);

  if ((token = strtok(NULL, "\t")) == NULL)
    return(-1);
  strncpy(feature->name, token, MAX_FEATURE);

  if ((token = strtok(NULL, "\t")) == NULL)
    return(-1);

  if ((token = strtok(NULL, "\t")) == NULL)
    return(-1);

  // Check if bed file is BED6 format and remove end of line character
  index = token;
  while(index[1]) ++index;
  if (*index == '\n') *index = '\0';

  if (strcmp(token, "+") == 0) feature->strand = FWD_STRAND;
  else feature->strand = REV_STRAND;

  return(1);
}
