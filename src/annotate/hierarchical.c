#include <annotate/hierarchical.h>

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
void print_tree(FILE* fp, hcnode_struct* hc, int index, profile_struct* profiles)
{
  fprintf(fp, "(");

  if (hc[index].left < 0) {
    int new_index = hc[index].left * (-1) - 1;
    print_tree(fp, hc, new_index, profiles);
  }
  else {
    profile_struct p = profiles[hc[index].left];
    if (p.strand == FWD_STRAND) fprintf(fp, "%s_%d-%d_+", p.chromosome, p.start, p.end);
    else fprintf(fp, "%s_%d-%d_-", p.chromosome, p.start, p.end);
  }

  fprintf(fp, ":%f,", hc[index].distance);

  if (hc[index].right < 0) {
    int new_index = hc[index].right * (-1) - 1;
    print_tree(fp, hc, new_index, profiles);
  }
  else {
    profile_struct p = profiles[hc[index].right];
    if (p.strand == FWD_STRAND) fprintf(fp, "%s_%d-%d_+", p.chromosome, p.start, p.end);
    else fprintf(fp, "%s_%d-%d_-", p.chromosome, p.start, p.end);
  }

  fprintf(fp, ":%f)", hc[index].distance);
}

/*
 * cut_tree
 *   Auxiliar function for hc_annotate
 */
int cut_tree(hcnode_struct* hc, int nprofiles, int index, double cutoff)
{
  int parent_index;
  int new_index = index;

  // Return -1 if node is alone
  if (hc[new_index].distance > cutoff)
    return -1;

  // Check parent distance
  parent_index = hc[new_index].parent * (-1) - 1;
  while ((new_index < (nprofiles - 2)) && (hc[parent_index].distance <= cutoff)) {
    new_index = parent_index;
    parent_index = hc[new_index].parent * (-1) - 1;
  }

  return new_index;
}

/*
 * annotate_tree_r
 *   Recursive function for annotate_tree
 */
void annotate_tree_r(hcnode_struct* hc, int root, int start, int end, profile_struct* profiles, double** contributions, int* leafs, char** labels)
{
  int i, j;

  // Fill matrix
  for (i = start; i < (end - 1); i++) {
    for (j = (i + 1); j < end; j++) {
      contributions[i][j] = (1 - hc[root].distance);
      contributions[j][i] = (1 - hc[root].distance);
    }
  }

  // Left leaf
  if (root[hc].left >= 0) {
    leafs[start] = hc[root].left;
    i = 0;
    while((strcmp(labels[i], "NULLPOINTER") != 0) && (strcmp(labels[i], profiles[leafs[start]].annotation) != 0)) i++;
    if (strcmp(labels[i], "NULLPOINTER") == 0) strncpy(labels[i], profiles[leafs[start]].annotation, MAX_FEATURE);
  }
  else
    annotate_tree_r(hc, hc[root].left * (-1) - 1, start, start + hc[root].lleafs, profiles, contributions, leafs, labels);

  // Right leaf
  if (root[hc].right >= 0) {
    leafs[end - 1] = hc[root].right;
    i = 0;
    while((strcmp(labels[i], "NULLPOINTER") != 0) && (strcmp(labels[i], profiles[leafs[end - 1]].annotation) != 0)) i++;
    if (strcmp(labels[i], "NULLPOINTER") == 0) strncpy(labels[i], profiles[leafs[end - 1]].annotation, MAX_FEATURE);
  }
  else
    annotate_tree_r(hc, hc[root].right * (-1) - 1, end - hc[root].rleafs, end, profiles, contributions, leafs, labels);
}

/*
 * annotate_tree
 *   Auxiliar function for hc_annotate
 */
void annotate_tree(hcnode_struct* hc, int nprofiles, profile_struct* profiles, int index, int root, double cutoff)
{
  double** contributions;
  int* leafs;
  char** labels;
  char** annotations;
  int i, j, nleafs;
  double max_score;
  char current_annotation[MAX_FEATURE];

  // Initialize 
  nleafs = hc[root].lleafs + hc[root].rleafs;
  leafs = (int*) malloc(sizeof(int) * nleafs);
  labels = (char**) malloc(sizeof(char*) * nleafs);
  annotations = (char**) malloc(sizeof(char*) * nleafs);
  contributions = (double**) malloc(sizeof(double*) * nleafs);
  for (i = 0; i < nleafs; i++) {
    contributions[i] = (double*) malloc(sizeof(double) * nleafs);
    contributions[i][i] = 0;
    labels[i] = (char*) malloc(sizeof(char) * MAX_FEATURE);
    strncpy(labels[i], "NULLPOINTER", MAX_FEATURE);
    annotations[i] = (char*) malloc(sizeof(char) * MAX_FEATURE);
  }

  // Calculate contributions
  annotate_tree_r(hc, root, 0, nleafs, profiles, contributions, leafs, labels);

  // Annotate profiles
  for (i = 0; i < nleafs; i++) {
    if (profiles[leafs[i]].annotation == NULL || strcmp(profiles[leafs[i]].annotation, "\0") == 0) {
      strncpy(current_annotation, labels[0], MAX_FEATURE);
      max_score = 0;
      for (j = 0; j < nleafs; j++) {
        if (j != i && strcmp(profiles[leafs[j]].annotation, labels[0]) == 0)
          max_score += contributions[i][j];
      }
      max_score = max_score / (nleafs - 1);

      int z = 1;
      while(z < nleafs && strcmp(labels[z], "NULLPOINTER") != 0) {
        double score = 0;
        for (j = 0; j < nleafs; j++) {
          if (j != i && strcmp(profiles[leafs[j]].annotation, labels[z]) == 0)
            score += contributions[i][j];
        }
        score = score / (nleafs - 1);
        if (score > max_score) {
          strncpy(current_annotation, labels[z], MAX_FEATURE);
          max_score = score;
        }
        z++;
      }

      strncpy(annotations[i], current_annotation, MAX_FEATURE);
      profiles[leafs[i]].anscore = max_score;
    }
  }

  // Copy annotations
  for (i = 0; i < nleafs; i++) {
    if (profiles[leafs[i]].annotation == NULL || strcmp(profiles[leafs[i]].annotation, "\0") == 0) {
      if (strcmp(annotations[i], "\0") == 0)
        strncpy(profiles[leafs[i]].annotation, "cluster", MAX_FEATURE);
      else
        strncpy(profiles[leafs[i]].annotation, annotations[i], MAX_FEATURE);
    }
  }

  // Free structures
  for(i = 0; i < nleafs; i++) {
    free(contributions[i]);
    free(labels[i]);
    free(annotations[i]);
  }
  free(leafs);
  free(labels);
  free(annotations);
  free(contributions);
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
void hc_print(FILE* fp, hcnode_struct* hc, int nprofiles, profile_struct* profiles)
{
  print_tree(fp, hc, nprofiles - 2, profiles);
  fprintf(fp, ";");
}

/*
 * hc_annotate
 *
 * @see include/annotate/hierarchical.h
 */
void hc_annotate(hcnode_struct* hc, int nprofiles, profile_struct* profiles, double cutoff)
{
  int i;

  for (i = 0; i < (nprofiles - 1); i++) {
    if ((hc[i].right >= 0 && (profiles[hc[i].right].annotation == NULL || strcmp(profiles[hc[i].right].annotation, "") == 0)) ||
        (hc[i].left >= 0 && (profiles[hc[i].left].annotation == NULL || strcmp(profiles[hc[i].left].annotation, "") == 0)))
    {
      int node = cut_tree(hc, nprofiles, i, cutoff);
      if (node >= 0)
        annotate_tree(hc, nprofiles, profiles, i, node, cutoff);
      else
        if ((hc[i].right >= 0) && ((profiles[hc[i].right].annotation == NULL) || (strcmp(profiles[hc[i].right].annotation, "") == 0)))
          strncpy(profiles[hc[i].right].annotation, "cluster", MAX_FEATURE);
        if ((hc[i].left >= 0) && ((profiles[hc[i].left].annotation == NULL) || (strcmp(profiles[hc[i].left].annotation, "") == 0)))
          strncpy(profiles[hc[i].left].annotation, "cluster", MAX_FEATURE);
    }
  }
}
