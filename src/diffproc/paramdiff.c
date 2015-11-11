#include <diffproc/paramdiff.h>

/*
 * parse_command_line_d
 *
 * @see include/diffproc/paramdiff.h
 */
int parse_command_line_d(int argc, char** argv, char** error_message, args_d_struct* arguments)
{
  opterr = 0;
  char carg;
  int terminate = 0;

  while(((carg = getopt(argc, argv, "hvg:")) != -1) && (terminate >= 0)) {
    switch (carg) {
      case 'h':
        terminate--;
        *error_message = DIFFPROC_HELP_MSG;
        break;
      case 'v':
        terminate--;
        *error_message = VERSION_MSG;
        break;
      case 'g':
        terminate = parse_filter_output_parameters(optarg, error_message, arguments);
        break;
      case '?':
        terminate--;
        *error_message = ERR_INVALID_ARGUMENT;
    }
  }
  if (!terminate && (argc - optind) != 5) {
    terminate--;
    *error_message = ERR_INVALID_NUMBER_ARGUMENTS;
  }
  else if (!terminate) {
    strncpy(arguments->profiles_a_f_path, argv[argc - 5], MAX_PATH);
    strncpy(arguments->clusters_a_f_path, argv[argc - 4], MAX_PATH);
    strncpy(arguments->profiles_b_f_path, argv[argc - 3], MAX_PATH);
    strncpy(arguments->clusters_b_f_path, argv[argc - 2], MAX_PATH);
    strncpy(arguments->output_f_path, argv[argc - 1], MAX_PATH);
  }

  return(terminate);
}


/*
 * parse_filter_output_parameters
 *
 * @see include/diffproc/paramdiff.h
 */
int parse_filter_output_parameters(char* option, char** error_message, args_d_struct* arguments)
{
  char* token;

  if (strchr(option, ':') == NULL) {
    *error_message = ERR_INVALID_g_VALUE;
    return(-1);
  }

  if ((token = strtok(option, ":")) != NULL) {
    arguments->pvalue = atof(token);
    if (arguments->pvalue < 0 || arguments->pvalue > 1) {
      *error_message = ERR_INVALID_pvalue_VALUE;
      return(-1);
    }
  }
  else {
    *error_message = ERR_INVALID_g_VALUE;
    return(-1);
  }

  if ((token = strtok(NULL, ":")) != NULL) {
    arguments->overlap = atof(token);
    if (arguments->overlap < 0 || arguments->overlap > 1) {
      *error_message = ERR_INVALID_cluster_overlapping_VALUE;
      return(-1);
    }
  }
  else {
    *error_message = ERR_INVALID_g_VALUE;
    return(-1);
  }

   if (((token = strtok(NULL, ":")) != NULL)) {
    *error_message = ERR_INVALID_g_VALUE;
    return(-1);
  }

  return(0);
}
