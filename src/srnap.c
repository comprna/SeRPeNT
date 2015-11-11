#include <srnap.h>

/*
 * Application entry point
 */
int main(int argc, char** argv)
{
  // Make sure the user at least entered a subcommand
  if (argc < 2) {
    fprintf(stderr, "\n%s\n", ERR_INVALID_NUMBER_ARGUMENTS);
    fprintf(stderr, "\n%s\n", ERR_HELP_MSG);
    return(-1);
  }

  // Parse user subcommand and run module if needed
  else if (strcmp(argv[1], PROFILES_SUBCOMMAND) == 0)
    return profiles_sc(argc - 1, argv + 1);
  else if (strcmp(argv[1], ANNOTATE_SUBCOMMAND) == 0)
    return annotate_sc(argc - 1, argv + 1);
  else if (strcmp(argv[1], DIFFPROC_SUBCOMMAND) == 0)
    return diffproc_sc(argc - 1, argv + 1);
  else if ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0)) {
    fprintf(stderr, "%s\n", GENERAL_HELP_MSG);
    return(0);
  }
  else if ((strcmp(argv[1], "-v") == 0) || (strcmp(argv[1], "--version") == 0)) {
    fprintf(stderr, "%s\n", VERSION_MSG);
    return(0);
  }
  else {
    fprintf(stderr, "\n%s\n", ERR_INVALID_SUBCOMMAND);
    fprintf(stderr, "\n%s\n", ERR_HELP_MSG);
    return(-1);
  }
}
