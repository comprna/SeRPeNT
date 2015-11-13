#include <annotate/annotation.h>

/*
 * compare_annotation
 *   Compare two variables of type annotation_struct
 *
 * @arg const void *a
 *   Pointer to the first value to be compared
 * @arg const void *b
 *   Pointer to the second value to be compared
 *
 * @return -1 if a < b. 0 if a = b. 1 if a > b.
 */
int compare_annotation (const void *a, const void *b)
{
  annotation_struct aa = *(annotation_struct*) a;
  annotation_struct bb = *(annotation_struct*) b;

  if (aa.score < bb.score) return -1;
  if (aa.score > bb.score) return  1;
  return 0;
}


/*
 * xcorr_annotate
 *
 * @see include/annotate/annotation.c
 */
void xcorr_annotate(annotation_struct** xcorr, int nprofiles, profile_struct_annotation* profiles)
{
  int i, sindex, stotal;
  annotation_struct top_profiles[3];

  for (i = 0; i < nprofiles; i++) {
    if (strcmp(profiles[i].annotation, "unknown") == 0) {
      sindex = 0;
      stotal = 0;

      annotation_struct* scores = xcorr[i];
      qsort(scores, nprofiles, sizeof(annotation_struct), compare_annotation);

      while(stotal < 3) {
        if (scores[sindex].index_i != scores[sindex].index_j) {
          top_profiles[stotal].score = scores[sindex].score;
          top_profiles[stotal].index_i = scores[sindex].index_i;
          top_profiles[stotal].index_j = scores[sindex].index_j;
          stotal++;
        }
        sindex++;
      }

      if (((strcmp(profiles[top_profiles[0].index_j].annotation, profiles[top_profiles[1].index_j].annotation) == 0)  ||
           (strcmp(profiles[top_profiles[0].index_j].annotation, profiles[top_profiles[2].index_j].annotation) == 0)) &&
          (strcmp(profiles[top_profiles[0].index_j].annotation, "unknown") != 0)) {
        strncpy(profiles[i].annotation, profiles[top_profiles[0].index_j].annotation, MAX_FEATURE);
        profiles[i].anscore = top_profiles[0].score;
      }
      else if ((strcmp(profiles[top_profiles[1].index_j].annotation, profiles[top_profiles[2].index_j].annotation) == 0) &&
               (strcmp(profiles[top_profiles[1].index_j].annotation, "unknown") != 0)) {
        strncpy(profiles[i].annotation, profiles[top_profiles[1].index_j].annotation, MAX_FEATURE);
        profiles[i].anscore = top_profiles[1].score;
      }
      else if (strcmp(profiles[top_profiles[0].index_j].annotation, "unknown") != 0) {
        strncpy(profiles[i].annotation, profiles[top_profiles[0].index_j].annotation, MAX_FEATURE);
        profiles[i].anscore = top_profiles[0].score;
      }
    }
  }
}


/*
 * cluster_annotate
 *
 * @see include/annotate/annotation.c
 */
void cluster_annotate(int nclusters, int nprofiles, profile_struct_annotation* profiles)
{
  int i, j, max_profiles_per_cluster, nclasses;
  int* profiles_per_cluster;
  int** profiles_index;
  StrMap *sm;
  char buffer[MAX_FEATURE];

  // Initialize data structures
  max_profiles_per_cluster = nprofiles - (nclusters - 1);
  profiles_per_cluster = (int*) malloc(nclusters * sizeof(int));
  profiles_index = (int**) malloc(nclusters * sizeof(int*));
  for (i = 0; i < nclusters; i++) {
    profiles_per_cluster[i] = 0;
    profiles_index[i] = (int*) malloc(max_profiles_per_cluster * sizeof(int));
    for (j = 0; j < max_profiles_per_cluster; j++) profiles_index[i][j] = -1;
  }
  sm = sm_new(nprofiles);
  
  // Fill data structures
  nclasses = 0;
  for (i = 0; i < nprofiles; i++) {
    int cindex = profiles[i].cluster - 1;
    j = 0;

    profiles_per_cluster[cindex]++;
    while (profiles_index[cindex][j] >= 0) j++;
    profiles_index[cindex][j] = i;

    if (sm_get(sm, profiles[i].annotation, buffer, MAX_FEATURE) == 0) {
      sprintf(buffer, "%d", nclasses);
      sm_put(sm, profiles[i].annotation, buffer);
      nclasses++;
    }
  }

  // Annotate unknown profiles in a cluster with the majority class
  for (i = 0; i < nclusters; i++) {
    int max = 0;
    char class[MAX_FEATURE] = "unknown";
    int* profiles_per_class = (int*) malloc(sizeof(int) * nclasses);

    for (j = 0; j < nclasses; j++) profiles_per_class[j] = 0;

    for (j = 0; j < profiles_per_cluster[i]; j++) {
      int pindex = profiles_index[i][j];
if (strcmp(profiles[pindex].annotation, "unknown") != 0) {
      if (sm_get(sm, profiles[pindex].annotation, buffer, MAX_FEATURE) != 0) profiles_per_class[atoi(buffer)]++;
      if (max < profiles_per_class[atoi(buffer)]) {
        max = profiles_per_class[atoi(buffer)];
        strncpy(class, profiles[pindex].annotation, MAX_FEATURE);
      }
}
    }

    for (j = 0; j < profiles_per_cluster[i]; j++) {
if (strcmp(class, "unknown") != 0) {
      int pindex = profiles_index[i][j];
      if (strcmp(profiles[pindex].annotation, "unknown") == 0) {
        //if (profiles[pindex].halo == 0)
          strncpy(profiles[pindex].annotation, class, MAX_FEATURE);
      }
}
    }

    free(profiles_per_class);
  }

  // Free data structures
  for (i = 0; i < nclusters; i++) free(profiles_index[i]);
  free(profiles_index);
  free(profiles_per_cluster);
  sm_delete(sm);
}


/*
 * xcorrspecies_annotate
 *
 * @see include/annotate/annotation.h
 */
void xspeciescorr_annotate(annotation_struct** xcorr, int nprofiles, profile_struct_annotation* profiles, int naprofiles, profile_struct_annotation* aprofiles, double cutoff)
{
  int i, sindex, stotal;
  annotation_struct top_profiles[3];

  for (i = 0; i < nprofiles; i++) {
    sindex = 0;
    stotal = 0;

    annotation_struct* scores = xcorr[i];
    qsort(scores, naprofiles, sizeof(annotation_struct), compare_annotation);

    while((stotal < 3) && (sindex < naprofiles)) {
      if ((scores[sindex].score <= cutoff) && (strcmp(aprofiles[scores[sindex].index_j].annotation, "unknown") != 0)) {
        top_profiles[stotal].score = scores[sindex].score;
        top_profiles[stotal].index_i = scores[sindex].index_i;
        top_profiles[stotal].index_j = scores[sindex].index_j;
        stotal++;
      }
      sindex++;
    }

    if ((stotal < 3) && (stotal > 0)) {
      strncpy(profiles[i].annotation, aprofiles[top_profiles[0].index_j].annotation, MAX_FEATURE);
      profiles[i].anscore = top_profiles[0].score;
    }

    else if (stotal > 0) {
      if ((strcmp(aprofiles[top_profiles[0].index_j].annotation, aprofiles[top_profiles[1].index_j].annotation) == 0)  ||
          (strcmp(aprofiles[top_profiles[0].index_j].annotation, aprofiles[top_profiles[2].index_j].annotation) == 0)) { 
        strncpy(profiles[i].annotation, aprofiles[top_profiles[0].index_j].annotation, MAX_FEATURE);
        profiles[i].anscore = top_profiles[0].score;
      }
      else if (strcmp(aprofiles[top_profiles[1].index_j].annotation, aprofiles[top_profiles[2].index_j].annotation) == 0) {
        strncpy(profiles[i].annotation, aprofiles[top_profiles[1].index_j].annotation, MAX_FEATURE);
        profiles[i].anscore = top_profiles[1].score;
      }
      else {
        strncpy(profiles[i].annotation, aprofiles[top_profiles[0].index_j].annotation, MAX_FEATURE);
        profiles[i].anscore = top_profiles[0].score;
      }
    }
  }
}
