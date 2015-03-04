#include <profiles/paramprof.h>

/*
 * parse_command_line
 *
 * @see include/profiles/paramprof.h
 */
int parse_command_line_p(int argc, char** argv, char** error_message, args_p_struct* arguments)
{
  opterr = 0;
  char carg;
  int terminate = 0;

  while(((carg = getopt(argc, argv, "hvf:p:r:i:")) != -1) && (terminate >= 0)) {
    switch (carg) {
      case 'h':
        terminate--;
        *error_message = PROFILES_HELP_MSG;
        break;
      case 'v':
        terminate--;
        *error_message = VERSION_MSG;
        break;
      case 'f':
        terminate = parse_filter_parameters(optarg, error_message, arguments);
        break;
      case 'i':
        terminate = parse_irreproducibility_parameters(optarg, error_message, arguments);
        break;
      case 'r':
        terminate = parse_replicates_parameters(optarg, error_message, arguments);
        break;
      case 'p':
        terminate = parse_profiles_parameters(optarg, error_message, arguments);
        break;
      case '?':
        terminate--;
        *error_message = ERR_INVALID_ARGUMENT;
    }
  }
  if (!terminate && (argc - optind) < 2) {
    terminate--;
    *error_message = ERR_INVALID_NUMBER_ARGUMENTS;
  }
  else if (!terminate && (argc - optind) > (MAX_REPLICATES + 1)) {
    terminate--;
    *error_message = ERR_INVALID_NUMBER_REPLICATES;
  }
  else if (!terminate) {
    int i, j = 0;
    for (i = optind; i < argc - 1; i++)
      strncpy(arguments->replicate_f_path[j++], argv[i], MAX_PATH);
    strncpy(arguments->output_f_path, argv[argc - 1], MAX_PATH);
    arguments->number_replicates = argc - optind - 1;
  }

  return(terminate);
}


/*
 * parse_filter_parameters
 *
 * @see include/profiles/paramprof.h
 */
int parse_filter_parameters(char* option, char** error_message, args_p_struct* arguments)
{
  arguments->min_read_len = atoi(option);
  if (arguments->min_read_len < 0) {
    *error_message = ERR_INVALID_minreadlen_VALUE;
    return(-1);
  }

  return(0);
}


/*
 * parse_irreproducibility_parameters
 *
 * @see include/profiles/paramprof.h
 */
int parse_irreproducibility_parameters(char* option, char** error_message, args_p_struct* arguments)
{
  char* token;

  if (strcmp(option, IDR_COMMON_STR) == 0)
    arguments->idr_method = IDR_COMMON;
  else if (strcmp(option, IDR_NONE_STR) == 0)
    arguments->idr_method = IDR_NONE;
  else {
    // Read method
    if ((token = strtok(option, ":")) != NULL) {
      if (strcmp(token, IDR_SERE_STR) == 0)
        arguments->idr_method = IDR_SERE;
      else if (strcmp(token, IDR_IDR_STR) == 0)
        arguments->idr_method = IDR_IDR;
      else {
        *error_message = ERR_INVALID_i_VALUE;
        return(-1);
      }
    }
    else {
      *error_message = ERR_INVALID_i_VALUE;
      return(-1);
    }

    // Read cutoff
    if ((token = strtok(NULL, ":")) != NULL) {
      arguments->idr_cutoff = atof(token);
      if (arguments->trimming < 0) {
        *error_message = ERR_INVALID_i_VALUE;
        return(-1);
      }
    }
    else {
      *error_message = ERR_INVALID_i_VALUE;
      return(-1);
    }

    // No more args
    if ((token = strtok(NULL, ":")) != NULL) {
      *error_message = ERR_INVALID_i_VALUE;
      return(-1);
    }
  }
  return(0);
}


/*
 * parse_replicates_parameters
 *
 * @see include/profiles/paramprof.h
 */
int parse_replicates_parameters(char* option, char** error_message, args_p_struct* arguments)
{
  char* token;

  if (strcmp(option, REPLICATE_POOL_STR) == 0)
    arguments->replicate_treat = REPLICATE_POOL;
  else if (strcmp(option, REPLICATE_MEAN_STR) == 0)
    arguments->replicate_treat = REPLICATE_MEAN;
  else {
    // Read treatment
    if ((token = strtok(option, ":")) != NULL) {
      if (strcmp(token, REPLICATE_REPLICATE_STR) == 0)
        arguments->replicate_treat = REPLICATE_REPLICATE;
      else {
        *error_message = ERR_INVALID_r_VALUE;
        return(-1);
      }
    }
    else {
      *error_message = ERR_INVALID_r_VALUE;
      return(-1);
    }

    // Read number
    if ((token = strtok(NULL, ":")) != NULL) {
      arguments->replicate_number = atoi(token);
      if ((arguments->replicate_number <= 0) ||
          (arguments->replicate_number > arguments->number_replicates))
      {
        *error_message = ERR_INVALID_repnumber_VALUE;
        return(-1);
      }
    }
    else {
      *error_message = ERR_INVALID_r_VALUE;
      return(-1);
    }

    // No more args
    if ((token = strtok(NULL, ":")) != NULL) {
      *error_message = ERR_INVALID_r_VALUE;
      return(-1);
    }
  }
  return(0);
}


/*
 * parse_profiles_parameters
 *
 * @see include/profiles/paramprof.h
 */
int parse_profiles_parameters(char* option, char** error_message, args_p_struct* arguments)
{
  char* token;

  if ((token = strtok(option, ":")) != NULL) {
    arguments->min_len = atoi(token);
    if (arguments->min_len < 5) {
      *error_message = ERR_INVALID_minlen_VALUE;
      return(-1);
    }
  }
  else {
    *error_message = ERR_INVALID_p_VALUE;
    return(-1);
  }

  if ((token = strtok(NULL, ":")) != NULL) {
    arguments->max_len = atoi(token);
    if (arguments->max_len < arguments->min_len) {
      *error_message = ERR_INVALID_maxlen_VALUE;
      return(-1);
    }
  }
  else {
    *error_message = ERR_INVALID_p_VALUE;
    return(-1);
  }

  if ((token = strtok(NULL, ":")) != NULL) {
      arguments->spacing = atoi(token);
      if (arguments->spacing < 0) {
        *error_message = ERR_INVALID_spacing_VALUE;
        return(-1);
      }
  }
  else {
    *error_message = ERR_INVALID_p_VALUE;
    return(-1);
  }

  if ((token = strtok(NULL, ":")) != NULL) {
    arguments->min_reads = atoi(token);
    if (arguments->min_reads <= 0) {
      *error_message = ERR_INVALID_minheight_VALUE;
      return(-1);
    }
  }
  else {
    *error_message = ERR_INVALID_p_VALUE;
    return(-1);
  }

  if ((token = strtok(NULL, ":")) != NULL) {
    arguments->trimming = atoi(token);
    if (arguments->trimming < 0) {
      *error_message = ERR_INVALID_trimming_VALUE;
      return(-1);
    }
  }
  else {
    *error_message = ERR_INVALID_p_VALUE;
    return(-1);
  }

  if (((token = strtok(NULL, ":")) != NULL)) {
    *error_message = ERR_INVALID_p_VALUE;
    return(-1);
  }

  return(0);
}
