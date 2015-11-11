#include <diffproc/diffproc.h>

/*
 * double comparison function for qsort
 */
int cmpds(const void *x, const void *y)
{
  double xx = *(double*)x, yy = *(double*)y;
  if (xx < yy) return -1;
  if (xx > yy) return  1;
  return 0;
}


/*
 * add_profile
 *   Inserts a profile into the array of profile structs
 */
void insert_profile(profile_struct_diffproc** condition, int* condition_n, profile_struct_diffproc profile, feature_struct_diffproc feature)
{
  int i, j, k;
  
  i = feature.cluster - 1;
  j = condition_n[feature.cluster - 1];

  condition[i][j].profile = (double*) malloc(profile.length * sizeof(double));
  for (k = 0; k < profile.length; k++) condition[i][j].profile[k] = profile.profile[k];
  strncpy(condition[i][j].chromosome, profile.chromosome, MAX_FEATURE);
  condition[i][j].start = profile.start;
  condition[i][j].end = profile.end;
  condition[i][j].length = profile.length;
  condition[i][j].strand = profile.strand;
  strncpy(condition[i][j].annotation, feature.name, MAX_FEATURE);
  for (k = 0; k < MAX_PROFILE_LENGTH; k++) condition[i][j].noise[k] = profile.noise[k];
  condition[i][j].cluster = feature.cluster;
  condition[i][j].position = j;
  condition[i][j].differential = 0;
  condition[i][j].partner = NULL;

  condition_n[feature.cluster - 1]++;
}


/*
 * assess_dp_c
 *   Assesses if a couple of clusters are differentially processed
 *
 * When using Mann-Whitney U test, certain sample sizes are too small to reject the null hypothesis under a certain alpha level.
 * In this case, means are compared and considered as differentially processed if 2fold.
 *
 * @arg double* a
 *   Sample of intra-cluster distances for cluster in condition A
 * @arg int n
 *   Number of observations in sample a
 * @arg double* b
 *   Sample of intra-cluster distances for cluster in condition B
 * @arg int m
 *   Number of observations in sample b
 * @arg double* x
 *   Sample of inter-cluster distances
 * @arg int y
 *   Number of observations in sample x
 * @arg double pvalue
 *   Alpha level
 *
 * @return
 *   1 if clusters are differentially processed under a pvalue level. 0 otherwise.
 *   -1 if differential proccessing is not assessable.
 */
int assess_dp_c(double* a, int n, double* b, int m, double* x, int y, double pval)
{
  double mean_a, mean_b, mean_x;

  // Differential processing is assessable by Mann-Whitney U test
  if (mannwhitney_a(n, y, pval) && mannwhitney_a(m, y, pval)) {
    if (y > 20)
      return (mannwhitney_i(a, n, x, y, pval) && mannwhitney_i(b, m, x, y, pval));
    else {
      if (n > 20 && m > 20)
        return (mannwhitney_i(a, n, x, y, pval) && mannwhitney_i(b, m, x, y, pval));
      if (n > 20 && m <= 20)
        return (mannwhitney_i(a, n, x, y, pval) && mannwhitney_d(b, m, x, y, pval));
      if (n <= 20 && m > 20)
        return (mannwhitney_d(a, n, x, y, pval) && mannwhitney_i(b, m, x, y, pval));
      else
        return (mannwhitney_d(a, n, x, y, pval) && mannwhitney_d(b, m, x, y, pval));
    }
  }

  // Differential processing is not assessable by Mann-Whitney U test.
  // Only one profile in cluster
  if (n == 0 || m == 0)
    return -1;

  // Differential processing is not assessable by Mann-Whitney U test
  //TODO => At least 2 profiles per cluster
  // Calculate means and check fold
  mean_a = gsl_stats_mean(a, 1, n);
  mean_b = gsl_stats_mean(b, 1, m);
  mean_x = gsl_stats_mean(x, 1, y);
  if ((mean_x >= 2.0 * mean_a) && (mean_x >= 2.0 * mean_b))
    return 1;

  return(0);
}


int assess_dp_p(double* a, int n, double* ba, int mn, double* b, int m, double* ab, int nm, double pval)
{
  double mean_a, mean_b, mean_ab, mean_ba;
  int mwab, mwba;

  // Differential processing is not assessable by Mann-Whitney U test.
  // Only one profile in cluster
  if (n == 0 || m == 0 || nm == 0 || mn == 0)
    return -1;

  // Differential processing is assessable by Mann-Whitney U test
  // between a and ba
  if (mannwhitney_a(n, mn, pval)) {
    if (n > 20 || mn > 20)
      mwba = mannwhitney_i(a, n, ba, mn, pval);
    else
      mwba = mannwhitney_d(a, n, ba, mn, pval);
  } 
  else {
    mean_a = gsl_stats_mean(a, 1, n);
    mean_ba = gsl_stats_mean(ba, 1, mn);
    mwba = (mean_ba >= 2.0 * mean_a);//TODO
  }

  // Differential processing is assessable by Mann-Whitney U test
  // between b and ab
  if (mannwhitney_i(b, m, ab, nm, pval)) {
    if (m > 20 || nm > 20)
      mwab = mannwhitney_i(b, m, ab, nm, pval);
    else
      mwab =  mannwhitney_d(b, m, ab, nm, pval);
  }
  else {
    mean_b = gsl_stats_mean(b, 1, m);
    mean_ab = gsl_stats_mean(ab, 1, nm);
    mwab = (mean_ab >= 2.0 * mean_b); //TODO
  }

  if (mwab > 0 && mwba > 0)
    return 1;

  return 0;
}


/*
 * Application entry point
 */
int diffproc_sc(int argc, char **argv)
{
  //Declare and define variables
  args_d_struct arguments;                 // Struct for handling command line parameters
  char* error_message;                     // Error message to display in case of abnormal termination
  FILE *clusters_a_file, *clusters_b_file; // Clusters file descriptors for conditions A and B
  FILE *profiles_a_file, *profiles_b_file; // Profile file descriptors for conditions A and B
  FILE *profiles_file, *clusters_file;     // Output file descriptors
  int nclusters_a, nclusters_b;            // Total number of clusters in conditions A and B
  int nprofiles_a, nprofiles_b;            // Total number of profiles in conditions A and B
  profile_struct_diffproc** cond_a;        // Array of profile structs for condition A. Each position in the array is a cluster number (0-based)
  profile_struct_diffproc** cond_b;        // Array of profile structs for condition B. Each position in the array is a cluster number (0-based)
  int result;                              // Result of any operation
  feature_struct_diffproc feature;         // Feature from the cluster file
  profile_struct_diffproc profile;         // Profile from the profile file
  int *cond_a_n, *cond_b_n;                // Array with numbers of profiles per cluster for conditions A and B
  int i, j, k;                             // General purpose variables
  double** intra_a;                        // Intracluster distances for condition A
  double** intra_b;                        // Intracluster distances for condition B
  char strands[2][2] = {"+\0", "-\0"};     // Array for printing strand

  // Initialize options with default values
  arguments.pvalue = P_VALUE;
  arguments.overlap = DF_OVERLAP;

  // Parse command line
  // Exit if command is not well-formed
  if (parse_command_line_d(argc, argv, &error_message, &arguments) < 0) {
    fprintf(stderr, "%s\n", error_message);
    if ((strcmp(error_message, DIFFPROC_HELP_MSG) == 0) || (strcmp(error_message, VERSION_MSG) == 0))
      return(0);
    fprintf(stderr, "%s\n", ERR_DIFFPROC_HELP_MSG);
    return(1);
  }

  // Open clusters file from condition A and check number of clusters
  // Exit if clusters file does not exist or is not readable
  clusters_a_file = fopen(arguments.clusters_a_f_path, "r");
  if (!clusters_a_file) {
    fprintf(stderr, "%s - %s\n", ERR_CLUSTER_F_NOT_READABLE, arguments.clusters_a_f_path);
    return(1);
  }
  nclusters_a = find_clusters(clusters_a_file);
  fclose(clusters_a_file);

  // Open clusters file from condition B and check number of clusters
  // Exit if clusters file does not exist or is not readable
  clusters_b_file = fopen(arguments.clusters_b_f_path, "r");
  if (!clusters_b_file) {
    fprintf(stderr, "%s - %s\n", ERR_CLUSTER_F_NOT_READABLE, arguments.clusters_b_f_path);
    return(1);
  }
  nclusters_b = find_clusters(clusters_b_file);
  fclose(clusters_b_file);

  // Allocate memory for arrays with number of profiles per cluster
  cond_a_n = (int*) malloc(nclusters_a * sizeof(int));
  cond_b_n = (int*) malloc(nclusters_b * sizeof(int));
  for (i = 0; i < nclusters_a; i++) cond_a_n[i] = 0;
  for (i = 0; i < nclusters_b; i++) cond_b_n[i] = 0;

  // Open clusters file from condition A and check number of profiles per cluster
  fprintf(stderr, "[LOG] LOADING CLUSTERS FOR CONDITION A\n");
  clusters_a_file = fopen(arguments.clusters_a_f_path, "r");
  allocate_clusters(clusters_a_file, cond_a_n);
  fclose(clusters_a_file);
  fprintf(stderr, "[LOG]   %d clusters loaded\n", nclusters_a);

  // Open clusters file from condition B and check number of profiles per cluster
  fprintf(stderr, "[LOG] LOADING CLUSTERS FOR CONDITION B\n");
  clusters_b_file = fopen(arguments.clusters_b_f_path, "r");
  allocate_clusters(clusters_b_file, cond_b_n);
  fclose(clusters_b_file);
  fprintf(stderr, "[LOG]   %d clusters loaded\n", nclusters_b);

  // Allocate profile arrays
  cond_a = (profile_struct_diffproc**) malloc(nclusters_a * sizeof(profile_struct_diffproc*));
  cond_b = (profile_struct_diffproc**) malloc(nclusters_b * sizeof(profile_struct_diffproc*));
  for (i = 0; i < nclusters_a; i++) cond_a[i] = (profile_struct_diffproc*) malloc(cond_a_n[i] * sizeof(profile_struct_diffproc));
  for (i = 0; i < nclusters_b; i++) cond_b[i] = (profile_struct_diffproc*) malloc(cond_b_n[i] * sizeof(profile_struct_diffproc));
  for (i = 0; i < nclusters_a; i++) cond_a_n[i] = 0;
  for (i = 0; i < nclusters_b; i++) cond_b_n[i] = 0;

  // Simultaneously open clusters and profile files from condition A. Read and store data.
  fprintf(stderr, "[LOG] LOADING PROFILES FOR CONDITION A\n");
  nprofiles_a = 0;
  clusters_a_file = fopen(arguments.clusters_a_f_path, "r");
  profiles_a_file = fopen(arguments.profiles_a_f_path, "r");
  result = 1;
  while(result > 0) {
    int r1 = next_diffproc_feature(clusters_a_file, &feature);
    int r2 = next_diffproc_profile(profiles_a_file, &profile);
    if (r1 > 0 && r2 > 0) {
      insert_profile(cond_a, cond_a_n, profile, feature);
      free(profile.profile);
      nprofiles_a++;
    }
    else
      result = MIN(r1, r2);
  }
  fclose(clusters_a_file);
  fclose(profiles_a_file);
  fprintf(stderr, "[LOG]   %d profiles loaded\n", nprofiles_a);

  // Simultaneously open clusters and profile files from condition B. Read and store data.
  fprintf(stderr, "[LOG] LOADING PROFILES FOR CONDITION B\n");
  nprofiles_b = 0;
  clusters_b_file = fopen(arguments.clusters_b_f_path, "r");
  profiles_b_file = fopen(arguments.profiles_b_f_path, "r");
  result = 1;
  while(result > 0) {
    int r1 = next_diffproc_feature(clusters_b_file, &feature);
    int r2 = next_diffproc_profile(profiles_b_file, &profile);
    if (r1 > 0 && r2 > 0) {
      insert_profile(cond_b, cond_b_n, profile, feature);
      free(profile.profile);
      nprofiles_b++;
    }
    else
      result = MIN(r1, r2);
  }
  fclose(clusters_b_file);
  fclose(profiles_b_file);
  fprintf(stderr, "[LOG]   %d profiles loaded\n", nprofiles_b);

  // Calculate intracluster distances for condition A
  fprintf(stderr, "[LOG] CALCULATING INTRA CLUSTER DISTANCES FOR CONDITION A\n");
  intra_a = (double**) malloc(nclusters_a * sizeof(double*));
  for (i = 0; i < nclusters_a; i++) {
    int nc = cond_a_n[i];
    int td = ((nc - 1) * nc) / 2;
    intra_a[i] = (double*) malloc (td * sizeof(double));
    int idx = 0;
    for (j = 0; j < (nc - 1); j++) {
      profile_struct_annotation pa;
      pa.profile = cond_a[i][j].profile;
      int l;
      for (l = 0; l < MAX_PROFILE_LENGTH; l++) pa.noise[l] = cond_a[i][j].noise[l];
      pa.length = cond_a[i][j].length;
      for (k = j + 1; k < nc; k++) {
        profile_struct_annotation pb;
        pb.profile = cond_a[i][k].profile;
        for (l = 0; l < MAX_PROFILE_LENGTH; l++) pb.noise[l] = cond_a[i][k].noise[l];
        pb.length = cond_a[i][k].length;
        double xcr = xdtw(&pa, &pb);
        if (xcr < 0) xcr = 0;
        intra_a[i][idx] = 1 - xcr;
        idx++;
      }
    }
    qsort(intra_a[i], td, sizeof(double), cmpds);
  }

  // Calculate intracluster distances for condition B
  fprintf(stderr, "[LOG] CALCULATING INTRA CLUSTER DISTANCES FOR CONDITION B\n");
  intra_b = (double**) malloc(nclusters_b * sizeof(double*));
  for (i = 0; i < nclusters_b; i++) {
    int nc = cond_b_n[i];
    int td = ((nc - 1) * nc) / 2;
    intra_b[i] = (double*) malloc (td * sizeof(double));
    int idx = 0;
    for (j = 0; j < (nc - 1); j++) {
      profile_struct_annotation pa;
      pa.profile = cond_b[i][j].profile;
      int l; 
      for (l = 0; l < MAX_PROFILE_LENGTH; l++) pa.noise[l] = cond_b[i][j].noise[l];
      pa.length = cond_b[i][j].length;
      for (k = j + 1; k < nc; k++) {
        profile_struct_annotation pb;
        pb.profile = cond_b[i][k].profile;
        for (l = 0; l < MAX_PROFILE_LENGTH; l++) pb.noise[l] = cond_b[i][k].noise[l];
        pb.length = cond_b[i][k].length;
        double xcr = xdtw(&pa, &pb);
        if (xcr < 0) xcr = 0;
        intra_b[i][idx] = 1 - xcr;
        idx++;
      }
    }
    qsort(intra_b[i], td, sizeof(double), cmpds);
  }

  // Open output file for differentially processed profiles 
  char *profile_output_name = malloc((MAX_PATH + strlen(DIFFPROC_PROFILE_O_SUFFIX) + 2) * sizeof(char));
  strncpy(profile_output_name, arguments.output_f_path, MAX_PATH);
  strcat(profile_output_name, PATH_SEPARATOR);
  strcat(profile_output_name, DIFFPROC_PROFILE_O_SUFFIX);
  profiles_file = fopen(profile_output_name, "w");
  if (!profiles_file) {
    fprintf(stderr, "%s\n", ERR_OUTPUT_F_NOT_WRITABLE);
    return (1);
  }
  free(profile_output_name);

  // Open output file for differentially processed clusters
  char *cluster_output_name = malloc((MAX_PATH + strlen(DIFFPROC_CLUSTER_O_SUFFIX) + 2) * sizeof(char));
  strncpy(cluster_output_name, arguments.output_f_path, MAX_PATH);
  strcat(cluster_output_name, PATH_SEPARATOR);
  strcat(cluster_output_name, DIFFPROC_CLUSTER_O_SUFFIX);
  clusters_file = fopen(cluster_output_name, "w");
  if (!clusters_file) {
    fprintf(stderr, "%s\n", ERR_OUTPUT_F_NOT_WRITABLE);
    return (1);
  }
  free(cluster_output_name);

  // Do cluster pairwise comparison
  // For each pair of clusters, calculate distributions and overlapping
  // For each pair of clusters with overlapping >= arguments.overlapping, assess differential processing
  fprintf(stderr, "[LOG] ASSESSING DIFFERENTIALLY PROCESSED CLUSTERS BETWEEN CONDITIONS\n");
  fprintf(stderr, "      pval < %f\n", arguments.pvalue);
  fprintf(stderr, "      overlap > %f\n", arguments.overlap);
  for (i = 0; i < nclusters_a; i++) {
    for (j = 0; j < nclusters_b; j++) {
      int idxi;
      int idxj;

      // Check overlap
      int elements = 0;
      double ovp;
      for (idxi = 0; idxi < cond_a_n[i]; idxi++) {
        profile_struct_diffproc pda = cond_a[i][idxi];
        for (idxj = 0; idxj < cond_b_n[j]; idxj++) {
          profile_struct_diffproc pdb = cond_b[j][idxj];
          if ((strcmp(pda.chromosome, pdb.chromosome) == 0) && (pda.strand == pdb.strand) &&
              (pdb.start <= pda.end) && (pda.start <= pdb.end)) {
            cond_a[i][idxi].partner = &cond_b[j][idxj];
            cond_b[j][idxj].partner = &cond_a[i][idxi];
            elements++;
          }
        }
      }
      ovp = ((double) (elements * 2)) / ((double) (cond_a_n[i] + cond_b_n[j]));

      // Calculate inter-cluster distances
      if (ovp >= arguments.overlap) {
        double* inter = (double*) malloc((cond_a_n[i] * cond_b_n[j]) * sizeof(double));
        int idx = 0;

        for (idxi = 0; idxi < cond_a_n[i]; idxi++) {
          profile_struct_diffproc pda = cond_a[i][idxi];
          profile_struct_annotation pa;
          pa.profile = pda.profile;
          int l; for (l = 0; l < MAX_PROFILE_LENGTH; l++) pa.noise[l] = pda.noise[l];
          pa.length = pda.length;

          for (idxj = 0; idxj < cond_b_n[j]; idxj++) {
            profile_struct_diffproc pdb = cond_b[j][idxj];
            profile_struct_annotation pb;
            pb.profile = pdb.profile;
            int l; for (l = 0; l < MAX_PROFILE_LENGTH; l++) pb.noise[l] = pdb.noise[l];
            pb.length = pdb.length;
            double xcr = xdtw(&pa, &pb);
            if (xcr < 0) xcr = 0;
            inter[idx] = 1 - xcr;
            idx++;
          }
        }
        qsort(inter, cond_a_n[i] * cond_b_n[j], sizeof(double), cmpds);
        int tda = ((cond_a_n[i] - 1) * cond_a_n[i]) / 2;
        int tdb = ((cond_b_n[j] - 1) * cond_b_n[j]) / 2;
        int adpc = assess_dp_c(intra_a[i], tda, intra_b[j], tdb, inter, cond_a_n[i] * cond_b_n[j], arguments.pvalue);
        if (adpc > 0) {
          /* DEBUGGING PURPOSES
          int r;
          for (r = 0; r < tda; r++) fprintf(stdout, "%dV%d\tA\t%f\n", i + 1, j + 1, intra_a[i][r]);
          for (r = 0; r < tdb; r++) fprintf(stdout, "%dV%d\tB\t%f\n", i + 1, j + 1, intra_b[j][r]);         
          for (r = 0; r < cond_a_n[i] * cond_b_n[j]; r++) fprintf(stdout, "%dV%d\tAB\t%f\n", i + 1, j + 1, inter[r]);
             END DEBUGGING PURPOSES */
          double mean_a = gsl_stats_mean(intra_a[i], 1, tda);
          double mean_b = gsl_stats_mean(intra_b[j], 1, tdb);
          double mean_x = gsl_stats_mean(inter, 1, cond_a_n[i] * cond_b_n[j]);
          fprintf(clusters_file, "%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\n", i + 1, j + 1, cond_a_n[i], cond_b_n[j], ovp, mean_x / mean_a, mean_x / mean_b);
        }
        else if (adpc < 0)
          fprintf(clusters_file, "%d\t%d\t%d\t%d\t%.2f\tNA\tNA\n", i + 1, j + 1, cond_a_n[i], cond_b_n[j], ovp);
        free(inter);
      }
    }
  }

  // Assess differentially processed profiles
  // For each profile in condition A that has a partner, calculate differential processing
  fprintf(stderr, "[LOG] ASSESSING DIFFERENTIALLY PROCESSED PROFILES BETWEEN CONDITIONS\n");
  fprintf(stderr, "      pval < %f\n", arguments.pvalue);
  for (i = 0; i < nclusters_a; i++) {
    int idxi;
    for (idxi = 0; idxi < cond_a_n[i]; idxi++) {

      // Profile has no partner
      if (cond_a[i][idxi].partner == NULL)
        fprintf(profiles_file, "%s:%d-%d:%s\t%s\tNA\tNA\tNA\tNA\n", cond_a[i][idxi].chromosome, cond_a[i][idxi].start, cond_a[i][idxi].end, strands[cond_a[i][idxi].strand], cond_a[i][idxi].annotation);

      // Profile has partner
      else {
        int j = cond_a[i][idxi].partner->cluster - 1;
        int idxj = cond_a[i][idxi].partner->position;

        int tda = ((cond_a_n[i] - 1) * cond_a_n[i]) / 2;
        int tdb = ((cond_b_n[j] - 1) * cond_b_n[j]) / 2;
        int tdab = cond_b_n[j];
        int tdba = cond_a_n[i];

        profile_struct_diffproc pda = cond_a[i][idxi];
        profile_struct_diffproc pdb = cond_b[j][idxj];

        profile_struct_annotation pa;
        pa.profile = pda.profile; 
        for (k = 0; k < MAX_PROFILE_LENGTH; k++) pa.noise[k] = pda.noise[k];
        pa.length = pda.length;

        profile_struct_annotation pb;
        pb.profile = pdb.profile;
        for (k = 0; k < MAX_PROFILE_LENGTH; k++) pb.noise[k] = pdb.noise[k];
        pb.length = pdb.length;

        // Calculate profile in A against cluster in B
        double* interab = (double*) malloc(tdab * sizeof(double));
        int idxjj;
        for(idxjj = 0; idxjj < tdab; idxjj++) {
          profile_struct_annotation pbb;
          pbb.profile = cond_b[j][idxjj].profile;
          for (k = 0; k < MAX_PROFILE_LENGTH; k++) pbb.noise[k] = cond_b[j][idxjj].noise[k];
          pbb.length = cond_b[j][idxjj].length;
          double xcr = xdtw(&pa, &pbb);
          if (xcr < 0) xcr = 0;
          interab[idxjj] = 1 - xcr;
        }
        qsort(interab, tdab, sizeof(double), cmpds);

        // Calculate profile in B against cluster in A
        double* interba = (double*) malloc(tdba * sizeof(double));
        int idxii;
        for(idxii = 0; idxii < tdba; idxii++) {
          profile_struct_annotation paa;  
          paa.profile = cond_a[i][idxii].profile;
          for (k = 0; k < MAX_PROFILE_LENGTH; k++) paa.noise[k] = cond_a[i][idxii].noise[k];
          paa.length = cond_a[i][idxii].length;
          double xcr = xdtw(&paa, &pb);
          if (xcr < 0) xcr = 0;
          interba[idxii] = 1 - xcr;
        }
        qsort(interba, tdba, sizeof(double), cmpds);

        // Assess differential processing of profile in A against cluster in B
        // Assess differential processing of profile in B against cluster in A
        int adpc = assess_dp_p(intra_a[i], tda, interba, tdba, intra_b[j], tdb, interab, tdab, arguments.pvalue);
        if (adpc > 0) {
          /* DEBUGGING PURPOSES
          int r;
          char ft[MAX_FEATURE];
          sprintf(ft, "%s:%d-%d:%s", pda.chromosome, pda.start, pda.end, strands[pda.strand]);
          for (r = 0; r < tda; r++) fprintf(stdout, "%s\tA\t%f\n", ft, intra_a[i][r]);
          for (r = 0; r < tdb; r++) fprintf(stdout, "%s\tB\t%f\n", ft, intra_b[j][r]);         
          for (r = 0; r < tdba; r++) fprintf(stdout, "%s\tBA\t%f\n", ft, interba[r]);
          for (r = 0; r < tdab; r++) fprintf(stdout, "%s\tAB\t%f\n", ft, interab[r]);
             END DEBUGGING PURPOSES */
          double mean_a = gsl_stats_mean(intra_a[i], 1, tda);
          double mean_b = gsl_stats_mean(intra_b[j], 1, tdb);
          double mean_ab = gsl_stats_mean(interab, 1, tdab);
          double mean_ba = gsl_stats_mean(interba, 1, tdba);
          fprintf(profiles_file, "%s:%d-%d:%s\t%s\t", pda.chromosome, pda.start, pda.end, strands[pda.strand], pda.annotation);
          fprintf(profiles_file, "%s:%d-%d:%s\t%s\t", pdb.chromosome, pdb.start, pdb.end, strands[pdb.strand], pdb.annotation);
          fprintf(profiles_file, "%.2f\t%.2f\n", mean_ab / mean_b, mean_ba / mean_a);
        }
        else if (adpc < 0) {
          fprintf(profiles_file, "%s:%d-%d:%s\t%s\t", pda.chromosome, pda.start, pda.end, strands[pda.strand], pda.annotation);
          fprintf(profiles_file, "%s:%d-%d:%s\t%s\t", pdb.chromosome, pdb.start, pdb.end, strands[pdb.strand], pdb.annotation);
          fprintf(profiles_file, "NA\tNA\n");
        }

        free(interab);
        free(interba);        
      }
    }
  }

  // Print each profile in condition B that does not have a partner
  for (i = 0; i < nclusters_b; i++) {
    int idxi;
    for (idxi = 0; idxi < cond_b_n[i]; idxi++) {
      profile_struct_diffproc pda = cond_b[i][idxi];
      if (pda.partner == NULL)
        fprintf(profiles_file, "NA\tNA\t%s:%d-%d:%s\t%s\tNA\tNA\n", pda.chromosome, pda.start, pda.end, strands[pda.strand], pda.annotation);
    }
  }

  // Close descriptors
  fclose(clusters_file);
  fclose(profiles_file);

  // Free pointer and exit
  for (i = 0; i < nclusters_a; i++) {
    for (j = 0; j < cond_a_n[i]; j++)
      free(cond_a[i][j].profile);
    free(cond_a[i]);
    free(intra_a[i]);
  }
  for (i = 0; i < nclusters_b; i++) {
    for (j = 0; j < cond_b_n[i]; j++)
      free(cond_b[i][j].profile);
    free(cond_b[i]);
    free(intra_b[i]);
  }
  free(cond_a);
  free(cond_b);
  free(cond_a_n);
  free(cond_b_n);
  free(intra_a);
  free(intra_b);
  return(0);
}
