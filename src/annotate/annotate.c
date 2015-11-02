#include <annotate/annotate.h>

/*
 * Application entry point
 */
int annotate_sc(int argc,  char **argv)
{
  // Define and declare variables
  args_a_struct arguments;                               // Struct for handling command line parameters
  FILE *profiles_file;                                   // Profiles file descriptor
  FILE *additional_profiles_file;                        // Annotation file descriptor and additional profiles file descriptor
  FILE *correlations_file;                               // Correlations file descriptor
  FILE *xcorr_file, *clusters_file, *annotation_o_file;  // Output file descriptors
  int result;                                            // Result of any operation
  char* error_message;                                   // Error message to display in case of abnormal termination
  int nprofiles;                                         // Total number of profiles
  int naprofiles;                                        // Total number of additional profiles
  int nclusters;                                         // Total number of clusters
  double** xcorr;                                        // 2-dimensional matrix containing correlations between profiles
  double** xspeciescorr;                                 // 2-dimensional matrix containing correlations between intra-species profiles
  int i, j, index;                                       // Multi-purpose indexes
  profile_struct_annotation* profiles;                   // Array of profiles
  profile_struct_annotation* additional_profiles;        // Array of additional profiles
  map_struct map;                                        // Profile map
  char categories[2][6] = {"NOVEL\0", "KNOWN\0"};        // Array for printing category
  char strands[2][2] = {"+\0", "-\0"};                   // Array for printing strand
  double cutoff;                                         // Branching calculated cutoff value

  // Initialize options with default values. Parse command line.
  // Exit if command is not well-formed.
  arguments.annotation = ANNOTATION_CONDITION;
  arguments.additional_profiles = ADDITIONAL_P_CONDITION;
  arguments.cluster_cutoff = CLUSTER_CUTOFF;
  arguments.overlap_ftop = OVERLAP_FTOP;
  arguments.overlap_ptof = OVERLAP_PTOF;
  arguments.correlations = CORRELATIONS_CONDITION;
  if (parse_command_line_c(argc, argv, &error_message, &arguments) < 0) {
    fprintf(stderr, "%s\n", error_message);
    if ((strcmp(error_message, ANNOTATE_HELP_MSG) == 0) || (strcmp(error_message, VERSION_MSG) == 0))
      return(0);
    fprintf(stderr, "%s\n", ERR_ANNOTATE_HELP_MSG);
    return(1);
  }

  // Open output files for writing results.
  // Exit if output files do not exist or are not readable.
  char *clusters_file_name = malloc((MAX_PATH + strlen(CLUSTERS_SUFFIX) + 2) * sizeof(char));
  strncpy(clusters_file_name, arguments.output_f_path, MAX_PATH);
  strcat(clusters_file_name, PATH_SEPARATOR);
  strcat(clusters_file_name, CLUSTERS_SUFFIX);
  clusters_file = fopen(clusters_file_name, "w");
  if (!clusters_file) {
    fprintf(stderr, "%s\n", ERR_OUTPUT_F_NOT_WRITABLE);
    return (1);
  }
  free(clusters_file_name);

  // Open annotation output file if annotation input file is provided
  if (arguments.annotation) {
    char *annotation_file_output_name = malloc((MAX_PATH + strlen(ANNOTATION_O_SUFFIX) + 2) * sizeof(char));
    strncpy(annotation_file_output_name, arguments.output_f_path, MAX_PATH);
    strcat(annotation_file_output_name, PATH_SEPARATOR);
    strcat(annotation_file_output_name, ANNOTATION_O_SUFFIX);
    annotation_o_file = fopen(annotation_file_output_name, "w");
    if (!annotation_o_file) {
      fprintf(stderr, "%s\n", ERR_OUTPUT_F_NOT_WRITABLE);
      return (1);
    }
    free(annotation_file_output_name);
  }

  // Open profiles file for reading and check number of lines.
  // Exit if profiles file do not exist or is not readable
  profiles_file = fopen(arguments.profiles_f_path, "r");
  if (!profiles_file) {
    fprintf(stderr, "%s - %s\n", ERR_PROFILE_F_NOT_READABLE, arguments.profiles_f_path);
    return(1);
  }
  nprofiles = wcl(profiles_file);
  fclose(profiles_file);

  // Check if additional profile file is provided.
  // Open additional profiles file for reading and check number of lines.
  // Exit if additional profiles file do not exist or is not readable
  if (arguments.additional_profiles) {
    additional_profiles_file = fopen(arguments.additional_profiles_f_path, "r");
    if (!additional_profiles_file) {
      fprintf(stderr, "%s - %s\n", ERR_PROFILE_F_NOT_READABLE, arguments.additional_profiles_f_path);
      return(1);
    }
    naprofiles = wcl(additional_profiles_file);
    fclose(additional_profiles_file);
  }

  // Allocate memory for profiles
  profiles = (profile_struct_annotation*) malloc(nprofiles * sizeof(profile_struct_annotation));
  map_init(&map);

  // Open profiles file for reading and load them into memory
  // Map the profiles if annotation is provided and no additional profiles file is provided
  // Exit if profiles file does not exist, is not readable or is ill-formatted
  fprintf(stderr, "[LOG] LOADING PROFILES\n");
  index = 0;
  profiles_file = fopen(arguments.profiles_f_path, "r");
  if (!profiles_file) {
    fprintf(stderr, "%s - %s\n", ERR_PROFILE_F_NOT_READABLE, arguments.profiles_f_path);
    return(1);
  }
  while((result = next_profile(profiles_file, &profiles[index++]) > 0)) {
    if (arguments.annotation && !arguments.additional_profiles)
      map_add_profile(&map, &profiles[index - 1]);
  }
  if (result < 0) {
    fprintf(stderr, "%s - %s\n", ERR_PROFILE_F_NOT_READABLE, arguments.profiles_f_path);
    return(1);
  }

  // Check if additional profile file is provided.
  // Allocate memory for additional profiles
  // Open additional profiles file for reading and load them into memory.
  // Exit if profiles file does not exist, is not readable or is ill-formatted.
  if (arguments.additional_profiles) {
    additional_profiles = (profile_struct_annotation*) malloc(naprofiles * sizeof(profile_struct_annotation));
    index = 0;
    fprintf(stderr, "[LOG] LOADING ADDITIONAL PROFILES FROM %s\n", arguments.species);
    additional_profiles_file = fopen(arguments.additional_profiles_f_path, "r");
    if (!additional_profiles_file) {
      fprintf(stderr, "%s - %s\n", ERR_PROFILE_F_NOT_READABLE, arguments.additional_profiles_f_path);
      return(1);
    }
    while((result = next_additional_profile(additional_profiles_file, &additional_profiles[index++], arguments.species) > 0)) {
      if (arguments.annotation)
        map_add_profile(&map, &additional_profiles[index - 1]);
    }
    if (result < 0) {
      fprintf(stderr, "%s - %s\n", ERR_PROFILE_F_NOT_READABLE, arguments.additional_profiles_f_path);
      return(1);
    }
  }

  // Check if annotation files are provided
  // Read annotation files and annotate profiles if no additional profiles provided
  // Read annotation files and annotate additional profiles if additional profiles provided
  if (arguments.annotation) {
    feature_struct feature;

    fprintf(stderr, "[LOG] LOADING ANNOTATIONS\n");
    for (i = 0; i < arguments.annotation; i++) {
      FILE* annotation_i_file;

      fprintf(stderr, "[LOG]   Reading %s\n", arguments.annotation_f_path[i]);
      annotation_i_file = fopen(arguments.annotation_f_path[i], "r");
      if (!annotation_i_file) {
        fprintf(stderr, "%s - %s\n", ERR_ANNOTATION_F_NOT_READABLE, arguments.annotation_f_path[i]);
        return(1);
      }
      while((result = next_feature(annotation_i_file, &feature) > 0))
        map_annotate(&arguments, &map, feature.chromosome, feature.start, feature.end, feature.strand, feature.name);
      if (result < 0) {
        fprintf(stderr, "%s - %s\n", ERR_ANNOTATION_F_NOT_READABLE, arguments.annotation_f_path[i]);
        return(1);
      }
      fclose(annotation_i_file);
    }
  }

  // Destroy map and free memory
  map_destroy(&map);

  // Allocate memory for correlation
  xcorr = (double**) malloc(nprofiles * sizeof(double*));
  for (i = 0; i < nprofiles; i++) {
    xcorr[i] = (double*) malloc(nprofiles * sizeof(double));
    profiles[i].anscore = 0;
  }

  // Allocate memory for cross-species correlations if additional profiles provided
  if (arguments.additional_profiles) {
    xspeciescorr = (double**) malloc(nprofiles * sizeof(double*));
    for (i = 0; i < nprofiles; i++)
      xspeciescorr[i] = (double*) malloc(naprofiles * sizeof(double));
  }

  // Read correlations file and store data
  if (arguments.correlations) {
    double score;

    fprintf(stderr, "[LOG] LOADING DISTANCE SCORES\n");
    correlations_file = fopen(arguments.correlations_f_path, "r");
    if (!correlations_file) {
      fprintf(stderr, "%s - %s\n", ERR_CORRELATIONS_F_NOT_READABLE, arguments.correlations_f_path);
      return(1);
    }
    i = 0; j = i + 1;
    while((result = next_correlation(correlations_file, &score) > 0)) {
      xcorr[i][j] = score;
      xcorr[j][i] = score;
      j++;
      if (j == nprofiles) {
        xcorr[i][i] = 0;
        i++;
        j = i + 1;
      }
    }
    if (result < 0) {
      fprintf(stderr, "%s - %s\n", ERR_CORRELATIONS_F_NOT_READABLE, arguments.correlations_f_path);
      return(1);
    }

    xcorr[nprofiles - 1][nprofiles - 1] = 0.0f;

    fclose(correlations_file);
  }

  // Calculate xcorrelations
  else {
    fprintf(stderr, "[LOG] CALCULATING DISTANCE SCORES\n");
    for (i = 0; i < (nprofiles - 1); i++) {
      xcorr[i][i] = (double) 0.0f;
      for (j = i + 1; j < nprofiles; j++) {
        double corr = xdtw(&profiles[i], &profiles[j]);
        if (corr < 0)
          corr = 0;
        xcorr[i][j] = 1 - corr;
        xcorr[j][i] = 1 - corr;
      }
    }
    xcorr[nprofiles - 1][nprofiles - 1] = 0.0f;

    // Print xcorrelations
    char *xcorr_file_name = malloc((MAX_PATH + strlen(CROSSCOR_SUFFIX) + 2) * sizeof(char));
    strncpy(xcorr_file_name, arguments.output_f_path, MAX_PATH);
    strcat(xcorr_file_name, PATH_SEPARATOR);
    strcat(xcorr_file_name, CROSSCOR_SUFFIX);
    xcorr_file = fopen(xcorr_file_name, "w");

    if (!xcorr_file) {
      fprintf(stderr, "%s\n", ERR_OUTPUT_F_NOT_WRITABLE);
      return (1);
    }
    free(xcorr_file_name);

    for (i = 0; i < (nprofiles - 1); i++) {
      for (j = i + 1; j < nprofiles; j++) {
        if (profiles[i].strand == FWD_STRAND)
          fprintf(xcorr_file, "%s:%d-%d:+\t", profiles[i].chromosome, profiles[i].start, profiles[i].end);
        else
          fprintf(xcorr_file, "%s:%d-%d:-\t", profiles[i].chromosome, profiles[i].start, profiles[i].end);
        if (profiles[j].strand == FWD_STRAND)
          fprintf(xcorr_file, "%s:%d-%d:+\t", profiles[j].chromosome, profiles[j].start, profiles[j].end);
        else
          fprintf(xcorr_file, "%s:%d-%d:-\t", profiles[j].chromosome, profiles[j].start, profiles[j].end);
        fprintf(xcorr_file, "%f\n", xcorr[i][j]);
      }
    }

    fclose(xcorr_file);
  }

  // Clustering by dpClust
  fprintf(stderr, "[LOG] PERFORMING DP-CLUSTERING\n");
  nclusters = dclust(xcorr, nprofiles, profiles, -1, 0);

  // Calculate cross-species correlations if additional profiles provided
  if (arguments.additional_profiles) {
    fprintf(stderr, "[LOG] CALCULATING DISTANCE SCORES BETWEEN SPECIES\n");
    for (i = 0; i < nprofiles; i++) {
      for (j = 0; j < naprofiles; j++) {
        double corr = xdtw(&profiles[i], &additional_profiles[j]);
        if (corr < 0)
          corr = 0;
        xspeciescorr[i][j] = 1 - corr;
      }
    }
  }

  // Annotate unknown profiles if no additional profiles provided
  if (arguments.annotation && !arguments.additional_profiles) {
    cluster_annotate(nclusters, nprofiles, profiles);
  }

  // Annotate unknown profiles if additional profiles provided
  if (arguments.annotation && arguments.additional_profiles) {
    fprintf(stderr, "[LOG] ANNOTATING UNKNOWN PROFILES\n");

    annotation_struct** matrix = (annotation_struct**) malloc(sizeof(annotation_struct*) * nprofiles);

    for (i = 0; i < nprofiles; i++) {
      matrix[i] = (annotation_struct*) malloc(sizeof(annotation_struct) * naprofiles);
      for (j = 0; j < naprofiles; j++) {
        matrix[i][j].score = xspeciescorr[i][j];
        matrix[i][j].index_i = i;
        matrix[i][j].index_j = j;
      }
    }

    // TODO
    if (arguments.cluster_cutoff < 0)
      cutoff = 0.01;
    else
      cutoff = arguments.cluster_cutoff;

    xspeciescorr_annotate(matrix, nprofiles, profiles, naprofiles, additional_profiles, cutoff);

    for (i = 0; i < nprofiles; i++) free(matrix[i]);
    free(matrix);
  }

  // Print annotated profiles in BED format
  if (arguments.annotation) {
    for (i = 0; i < nprofiles; i++) {
      profile_struct_annotation p = profiles[i];
      fprintf(annotation_o_file, "%s\t%d\t%d\t%s\t%f\t%s\t%s\t%d\n", p.chromosome, p.start, p.end, p.annotation, p.anscore, strands[p.strand], categories[p.category], p.cluster);
    }
  }

  // Close descriptors, free structures and exit
  for (i = 0; i < nprofiles; i++)
    free(profiles[i].profile);
  free(profiles);
  for (i = 0; i < nprofiles; i++)
    free(xcorr[i]);
  free(xcorr);
  fclose(clusters_file);
  fclose(profiles_file);
  if (arguments.annotation)
    fclose(annotation_o_file);
  if (arguments.additional_profiles) {
    for (i = 0; i < naprofiles; i++)
      free(additional_profiles[i].profile);
    free(additional_profiles);
    for (i = 0; i < nprofiles; i++)
      free(xspeciescorr[i]);
    free(xspeciescorr);
    fclose(additional_profiles_file);
  }
  return(0);
}
