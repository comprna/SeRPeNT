#include <annotate/paramclust.h>

/*
 * parse_command_line_c
 *
 * @see src/include/annotate/paramclust.h
 */
int parse_command_line_c(int argc, char** argv, char** error_message, args_a_struct* arguments)
{
  opterr = 0;
  char carg;
  int terminate = 0;

  while(((carg = getopt(argc, argv, "hva:x:c:o:")) != -1) && (terminate >= 0)) {
    switch (carg) {
      case 'h':
        terminate--;
        *error_message = ANNOTATE_HELP_MSG;
        break;
      case 'v':
        terminate--;
        *error_message = VERSION_MSG;
        break;
      case 'a':
        terminate = parse_annotation_parameters(optarg, error_message, arguments);
        break;
      case 'x':
        terminate = parse_additional_profiles_parameters(optarg, error_message, arguments);
        break;
      case 'c':
        terminate = parse_threshold_parameters(optarg, error_message, arguments);
        break;
      case 'o':
        terminate = parse_overlapping_parameters(optarg, error_message, arguments);
        break;
      case '?':
        terminate--;
        *error_message = ERR_INVALID_ARGUMENT;
    }
  }

  if (!terminate && (argc - optind) != 2) {
    terminate--;
    *error_message = ERR_INVALID_NUMBER_ARGUMENTS;
  }
  else if (!terminate) {
    strncpy(arguments->profiles_f_path, argv[argc -2], MAX_PATH);
    strncpy(arguments->output_f_path, argv[argc - 1], MAX_PATH);
  }
  return(terminate);
}

/*
 * parse_annotation_parameters
 *
 * @see src/include/annotate/paramclust.h
 */
int parse_annotation_parameters(char* option, char** error_message, args_a_struct* arguments)
{
  char* token;

  if ((token = strtok(option, ":")) == NULL) {
    strncpy(arguments->annotation_f_path[arguments->annotation++], option, MAX_PATH);
    return(0);
  }
  strncpy(arguments->annotation_f_path[arguments->annotation++], token, MAX_PATH);

  while((token = strtok(NULL, ":")) != NULL)
    strncpy(arguments->annotation_f_path[arguments->annotation++], token, MAX_PATH);

  return(0);
}

/*
 * parse_additional_profiles_parameters
 *
 * @see src/include/annotate/paramclust.h
 */
int parse_additional_profiles_parameters(char* option, char** error_message, args_a_struct* arguments)
{
  char* token;
  
  if ((token = strtok(option, ":")) == NULL) {
    *error_message = ERR_INVALID_x_VALUE;
    return(-1);
  }
  strncpy(arguments->species, token, MAX_FEATURE);

  if ((token = strtok(NULL, ":")) == NULL) {
    *error_message = ERR_INVALID_x_VALUE;
    return(-1);
  }
  strncpy(arguments->additional_profiles_f_path, token, MAX_PATH);

  arguments->additional_profiles = 1;

  return(0);
}

/*
 * parse_cutoff_parameters
 *
 * @see src/include/annotate/paramclust.h
 */
int parse_threshold_parameters(char* option, char** error_message, args_a_struct* arguments)
{
  arguments->cluster_cutoff = atof(option);
  if (arguments->cluster_cutoff < 0) {
    *error_message = ERR_INVALID_c_VALUE;
    return(-1);
  }
  return(0);
}

/*
 * parse_overlapping_parameters
 *
 * @see include/annotate/paramclust.h
 */
int parse_overlapping_parameters(char* option, char** error_message, args_a_struct* arguments)
{
  char* token;

  if (strchr(option, ':') == NULL) {
    *error_message = ERR_INVALID_o_VALUE;
    return(-1);
  }

  if ((token = strtok(option, ":")) != NULL) {
    arguments->overlap_ftop = atof(token);
    if (arguments->overlap_ftop < 0 || arguments->overlap_ftop > 1) {
      *error_message = ERR_INVALID_ovpftop_VALUE;
      return(-1);
    }
  }
  else {
    *error_message = ERR_INVALID_o_VALUE;
    return(-1);
  }

  if ((token = strtok(NULL, ":")) != NULL) {
    arguments->overlap_ptof = atof(token);
    if (arguments->overlap_ptof < 0 || arguments->overlap_ptof > 1) {
      *error_message = ERR_INVALID_ovpptof_VALUE;
      return(-1);
    }
  }
  else {
    *error_message = ERR_INVALID_o_VALUE;
    return(-1);
  }

   if (((token = strtok(NULL, ":")) != NULL)) {
    *error_message = ERR_INVALID_o_VALUE;
    return(-1);
  }

  return(0);
}
