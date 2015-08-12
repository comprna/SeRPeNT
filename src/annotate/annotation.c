#include <annotate/annotation.h>

/*
 *  * compare_annotation
 *   *   Compare two variables of type annotation_struct
 *    *
 *     * @arg const void *a
 *      *   Pointer to the first value to be compared
 *       * @arg const void *b
 *        *   Pointer to the second value to be compared
 *         *
 *          * @return -1 if a < b. 0 if a = b. 1 if a > b.
 *           */
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

/*fprintf(stderr, "%s:%d-%d:%d = [", profiles[i].chromosome, profiles[i].start, profiles[i].end, profiles[i].strand);
fprintf(stderr, "%s:%d-%d:%d, ", profiles[top_profiles[0].index_j].chromosome, profiles[top_profiles[0].index_j].start, profiles[top_profiles[0].index_j].end, profiles[top_profiles[0].index_j].strand);
fprintf(stderr, "%s:%d-%d:%d, ", profiles[top_profiles[1].index_j].chromosome, profiles[top_profiles[1].index_j].start, profiles[top_profiles[1].index_j].end, profiles[top_profiles[1].index_j].strand);
fprintf(stderr, "%s:%d-%d:%d]\n", profiles[top_profiles[0].index_j].chromosome, profiles[top_profiles[2].index_j].start, profiles[top_profiles[2].index_j].end, profiles[top_profiles[2].index_j].strand);*/

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
