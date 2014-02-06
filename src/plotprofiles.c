#include <plotprofiles.h>

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
int parse_command_line(int argc, char **argv, char **error_message, args_struct *arguments)
{
  opterr = 0;
  char carg;
  int terminate = 0;

  while((carg = getopt(argc, argv, "o:")) != -1) {
    switch (carg) {
      case 'o':
        arguments->output_mode = OUTPUT_PNG;
        strncpy(arguments->output_png_path, optarg, MAX_PATH);
        break;
      case '?':
        terminate--;
        *error_message = ERR_INVALID_ARGUMENT;
    }
  }

  if (!terminate && (((argc - optind) < 2)) || ((argc - optind) > MAX_NUMBER_FEATURES + 1)) {
    terminate--;
    *error_message = ERR_INVALID_NUMBER_ARGUMENTS;
  }
  else if (!terminate) {
    int i, j = 0;
    strncpy(arguments->input_profile_path, argv[optind], MAX_PATH);
    for (i = optind + 1; i < argc; i++)
      strncpy(arguments->identifiers[j++], argv[i], MAX_FEATURE);
  }

  return(terminate);
}

/*
 * Application entry point
 */
int main(int argc, char** argv)
{
  // Define and declare variables
  args_struct arguments;            // Struct for handling command line parameters
  FILE *input_file, *output_file;   // Input and output (if needed) file descriptors
  char *error_message;              // Error message to display in case of abnormal termination
  

  // Initialize options with default values. Parse command line.
  // Exit if command is not well-formed.
  arguments.output_mode = OUTPUT_SCREEN;
  if (parse_command_line(argc, argv, &error_message, &arguments) < 0) {
    MSG_PRINT(1, error_message, ERR_MSG_PLOTPROFILES);
    DEBUG_PRINT(error_message);
    return(1);
  }

  // Open input file
  // Exit if file does not exist or is not readable
  input_file = fopen(arguments.input_profile_path, "r");
  if (!input_file) {
    MSG_PRINT(0, ERR_INPUT_F_NOT_READABLE, "");
    DEBUG_PRINT(ERR_INPUT_F_NOT_READABLE);
    return (1);
  }

  // Parse input file and extract features to plot
  size_t linesize = 0;
  char *linebuf = 0;
  ssize_t linelength = 0;
  while ((linelength = getline(&linebuf, &linesize, input_file)) > 0) {
    fprintf(stderr, "%s", linebuf);
  }

  // Execute gnuplot on generated files

  // Close file descriptors and return
  return(0);
}
