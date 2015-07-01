#include <annotate/hierarchical.h>

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
 * convert_tree
 *   Recursive function for hc_cluster
 */
void convert_tree(hcnode_struct* tree, int index)
{
  if (tree[index].left < 0) {
    int newindex = tree[index].left * (-1) - 1;
    tree[newindex].parent = (index + 1) * (-1);
    convert_tree(tree, newindex);
    tree[index].lleafs += tree[newindex].lleafs;
    tree[index].lleafs += tree[newindex].rleafs;
  }
  else
    tree[index].lleafs++;

  if (tree[index].right < 0) {
    int newindex = tree[index].right * (-1) - 1;
    tree[newindex].parent = (index + 1) * (-1);
    convert_tree(tree, newindex);
    tree[index].rleafs += tree[newindex].lleafs;
    tree[index].rleafs += tree[newindex].rleafs;
  }
  else
    tree[index].rleafs++;
}

/*
 * print_tree
 *   Recursive function for hc_print
 */
void print_tree(FILE* fp, hcnode_struct* hc, int index, profile_struct_annotation* profiles)
{
  fprintf(fp, "(");

  if (hc[index].left < 0) {
    int new_index = hc[index].left * (-1) - 1;
    print_tree(fp, hc, new_index, profiles);
  }
  else {
    profile_struct_annotation p = profiles[hc[index].left];
    if (p.strand == FWD_STRAND) fprintf(fp, "%s_%d-%d_+", p.chromosome, p.start, p.end);
    else fprintf(fp, "%s_%d-%d_-", p.chromosome, p.start, p.end);
  }

  fprintf(fp, ":%f,", hc[index].distance);

  if (hc[index].right < 0) {
    int new_index = hc[index].right * (-1) - 1;
    print_tree(fp, hc, new_index, profiles);
  }
  else {
    profile_struct_annotation p = profiles[hc[index].right];
    if (p.strand == FWD_STRAND) fprintf(fp, "%s_%d-%d_+", p.chromosome, p.start, p.end);
    else fprintf(fp, "%s_%d-%d_-", p.chromosome, p.start, p.end);
  }

  fprintf(fp, ":%f)", hc[index].distance);
}

/*
 * annotate_tree_r
 *   Recursive function for annotate_tree
 */
void annotate_tree_r(hcnode_struct* hc, int root, int start, int end, profile_struct_annotation* profiles, int* leafs)
{
  // Mark node as visited
  hc[root].visited = 1;

  // Left leaf
  if (root[hc].left >= 0)
    leafs[start] = hc[root].left;
  else
    annotate_tree_r(hc, hc[root].left * (-1) - 1, start, start + hc[root].lleafs, profiles, leafs);

  // Right leaf
  if (root[hc].right >= 0)
    leafs[end - 1] = hc[root].right;
  else
    annotate_tree_r(hc, hc[root].right * (-1) - 1, end - hc[root].rleafs, end, profiles, leafs);
}

/*
 * annotate_tree
 *   Auxiliar function for hc_annotate
 */
void annotate_tree(hcnode_struct* hc, int nprofiles, profile_struct_annotation* profiles, double** correlations, int root, int* cluster)
{
  int* leafs;
  int nleafs, i, add;
  
  // Initialize 
  nleafs = hc[root].lleafs + hc[root].rleafs;
  leafs = (int*) malloc(sizeof(int) * nleafs);
  add = 0;

  // Find clustered profiles
  annotate_tree_r(hc, root, 0, nleafs, profiles, leafs);

  // For each unknown leaf, sort remaining leafs by distance
  for (i = 0; i < nleafs; i++) {
    if ((strcmp(profiles[leafs[i]].annotation, "unknown") == 0) &&
        (strcmp(profiles[leafs[i]].tmp_annotation, "unknown") == 0)) {
      sprintf(profiles[leafs[i]].tmp_annotation, "cluster_%d", *cluster);
      profiles[leafs[i]].anscore = hc[root].distance;
      add++;
    }

    profiles[leafs[i]].cluster = *cluster;
  }

  // Increment cluster pointer
  if (add) (*cluster)++;

  // Free structures
  free(leafs);
}

/*
 * hc_cluster
 * 
 * @see include/annotate/hierarchical.h
 */
hcnode_struct* hc_cluster(double** correlation, int nprofiles)
{
  hcnode_struct* hc;
  int i;

  // Perform hierarchical clustering
  Node* tree = treecluster(nprofiles, nprofiles, NULL, NULL, NULL, 0, 'e', 'm', correlation);

  // Deep copy of nodes
  hc = (hcnode_struct*) malloc((nprofiles - 1) * sizeof(hcnode_struct));
  for (i = 0; i < (nprofiles - 1); i++) {
    hc[i].left = tree[i].left;
    hc[i].right = tree[i].right;
    hc[i].distance = tree[i].distance;
    hc[i].parent = nprofiles * (-1);
    hc[i].lleafs = 0;
    hc[i].rleafs = 0;
    hc[i].visited = 0;
  }
  free(tree);

  // Assign parents to nodes
  convert_tree(hc, nprofiles - 2);

  return(hc);
}

/*
 * hc_print
 *
 * @see include/annotate/hierarchical.h
 */
void hc_print(FILE* fp, hcnode_struct* hc, int nprofiles, profile_struct_annotation* profiles)
{
  print_tree(fp, hc, nprofiles - 2, profiles);
  fprintf(fp, ";");
}

/*
 * hc_annotate
 *
 * @see include/annotate/hierarchical.h
 */
void hc_annotate(hcnode_struct* hc, int nprofiles, profile_struct_annotation* profiles, double** correlations, double cutoff)
{
  int i, cluster;

  // Annotate each node in the tree ascendently by distance
  cluster = 1;
  for (i = 0; i < (nprofiles - 1); i++) {
    if (!hc[i].visited) {
      int root = i;
      int parent = hc[root].parent * (-1) - 1;
      while ((root < (nprofiles - 2)) && (hc[root].distance <= cutoff) && (hc[parent].distance <= cutoff)) {
        root = parent;
        parent = hc[root].parent * (-1) - 1;
      }
      if (hc[root].distance <= cutoff)
        annotate_tree(hc, nprofiles, profiles, correlations, root, &cluster);
    }
  }

  // Assign labels
  for (i = 0; i < nprofiles; i++) {
    if (profiles[i].cluster < 0)
      profiles[i].cluster = cluster;

    if ((strcmp(profiles[i].annotation, "unknown") == 0) &&
        (strcmp(profiles[i].tmp_annotation, "unknown") == 0)) {
      sprintf(profiles[i].annotation, "cluster_%d", cluster++);
    }
    else if ((strcmp(profiles[i].annotation, "unknown") == 0) &&
             (strcmp(profiles[i].tmp_annotation, "unknown") != 0))
      strncpy(profiles[i].annotation, profiles[i].tmp_annotation, MAX_FEATURE);
  }
}

/*
 * xcorr_annotate
 *
 * @see include/annotate/hierarchical.c
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
        strncpy(profiles[i].tmp_annotation, profiles[top_profiles[0].index_j].annotation, MAX_FEATURE);
        profiles[i].anscore = top_profiles[0].score;
      }
      else if ((strcmp(profiles[top_profiles[1].index_j].annotation, profiles[top_profiles[2].index_j].annotation) == 0) &&
               (strcmp(profiles[top_profiles[1].index_j].annotation, "unknown") != 0)) {
        strncpy(profiles[i].tmp_annotation, profiles[top_profiles[1].index_j].annotation, MAX_FEATURE);
        profiles[i].anscore = top_profiles[1].score;
      }
      else if (strcmp(profiles[top_profiles[0].index_j].annotation, "unknown") != 0) {
        strncpy(profiles[i].tmp_annotation, profiles[top_profiles[0].index_j].annotation, MAX_FEATURE);
        profiles[i].anscore = top_profiles[0].score;
      }
    }
  }
}
