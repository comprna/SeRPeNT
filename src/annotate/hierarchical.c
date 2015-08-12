#include <annotate/hierarchical.h>

/*
 * convert_tree
 *   Recursive function for hc_cluster
 *
 * @arg hcnode_struct* tree
 *   A pointer to the hierarchical clustering solution
 * @arg index
 *   Index of the node being visited in the hierarchical clustering solution
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
 *
 * @arg FILE* fp
 *   Pointer to the output file descriptor
 * @arg hcnode_struct* hc
 *   A pointer to the hierarchical clustering solution
 * @arg index
 *   Index of the node being visited in the hierarchical clustering solution
 * @arg profile_struct_annotation* profiles
 *   An array of profiles
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
 * branch_tree
 *   Recursive function for hc_branch
 *
 * @arg hcnode_struct* hc
 *   A pointer to the hierarchical clustering solution
 * @arg index
 *   Index of the node being visited in the hierarchical clustering solution
 * @arg profile_struct_annotation* profiles
 *   An the array of profiles
 * @arg int cluster
 *   Cluster number
 */
void branch_tree(hcnode_struct* hc, int index, profile_struct_annotation* profiles, int cluster)
{
  // Mark node as visited
  hc[index].visited = 1;

  // Left leaf
  if (hc[index].left < 0)
    branch_tree(hc, hc[index].left * (-1) - 1, profiles, cluster);
  else
    profiles[hc[index].left].cluster = cluster;

  // Right leaf
  if (hc[index].right < 0)
    branch_tree(hc, hc[index].right * (-1) - 1, profiles, cluster);
  else
    profiles[hc[index].right].cluster = cluster;
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
 * hc_branch
 *
 * @see include/annotate/hierarchical.h
 */
int hc_branch(hcnode_struct* hc, int nprofiles, profile_struct_annotation* profiles, double cutoff)
{
  int i, cluster;

  // Assign cluster number to each branch
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
        branch_tree(hc, root, profiles, cluster++);
    }
  }

  // Assign cluster number to each individual non-visited leaf
  for (i = 0; i < nprofiles; i++)
    if (profiles[i].cluster < 0) profiles[i].cluster = cluster++;

  // Return number of clusters
  return (cluster - 1);
}


/*
 * hc_eval
 *
 * @see include/annotate/hierarchical.h
 */
double hc_eval(int nprofiles, profile_struct_annotation* profiles)
{
  int i, j, nclusters, nannotated, nclasses;
  double score, h_k, h_c, h_conk, h_konc, homogeneity, completeness;
  int *clusters, *classes;
  char** names;
  int** contingency;

  // Initialize variables
  clusters = (int*) malloc(sizeof(int) * nprofiles);
  classes = (int*) malloc(sizeof(int) * nprofiles);
  names = (char**) malloc(sizeof(char*) * nprofiles);
  contingency = (int**) malloc(sizeof(int*) * nprofiles);
  for (i = 0; i < nprofiles; i++) {
    clusters[i] = 0;
    classes[i] = 0;
    contingency[i] = (int*) malloc(sizeof(int) * nprofiles);
    for (j = 0; j < nprofiles; j++) contingency[i][j] = 0;
  }

  // Calculate number of classes and number of clusters
  nclusters = 0;
  nannotated = 0;
  nclasses = 0;
  for (i = 0; i < nprofiles; i++) {
    profile_struct_annotation p = profiles[i];
    if (strcmp(p.annotation, "unknown") != 0) {
      nannotated++;
      clusters[p.cluster - 1]++;
      j = 0;
      while ((j < nclasses) && (strcmp(names[j], p.annotation) != 0)) j++;
      if (j == nclasses) {
        names[j] = (char*) malloc(sizeof(char) * MAX_FEATURE);
        strncpy(names[j], p.annotation, MAX_FEATURE);
        nclasses++;
      }
      classes[j]++;
      contingency[p.cluster - 1][j]++;
    }
    if (nclusters < p.cluster) nclusters = p.cluster;
  }

  // Calculate H(K)
  h_k = 0;
  for (i = 0; i < nclusters; i++) {
    if (clusters[i] > 0) {
      double term = ((double) clusters[i]) / ((double) nannotated);
      double logterm = log2f(term);
      h_k -= term * logterm;
    }
  }

  // Calculate H(C)
  h_c = 0;
  for (i = 0; i < nclasses; i++) {
    if (classes[i] > 0) {
      double term = ((double) classes[i]) / ((double) nannotated);
      double logterm = log2f(term);
      h_c -= term * logterm;
    }
  }

  // Calculate H(C|K) and H(K|C)
  h_conk = 0;
  h_konc = 0;
  for (i = 0; i < nclusters; i++) {
    for (j = 0; j < nclasses; j++) {
      if (contingency[i][j] > 0) {
        double term_a = ((double) contingency[i][j]) / ((double) nannotated);
        double term_b = ((double) contingency[i][j]) / ((double) clusters[i]);
        double logterm = log2f(term_b);
        h_conk += term_a * logterm;
        term_b = ((double) contingency[i][j]) / ((double) classes[j]);
        logterm = log2f(term_b);
        h_konc += term_a * logterm;
      }
    }
  }
  h_conk *= -1;
  h_konc *= -1;

  // Calculate V-measure
  homogeneity = 1;
  completeness = 1;
  if (h_c != 0) homogeneity = 1 - (h_conk / h_c);
  if (h_k != 0) completeness = 1 - (h_konc / h_k);
  score = (2 * homogeneity * completeness) / (homogeneity + completeness);

  // Free arrays
  free(clusters);
  free(classes);
  for (i = 0; i < nclasses; i++) free(names[i]);
  for (i = 0; i < nprofiles; i++) free(contingency[i]);
  free(names);
  free(contingency);

  return (score);
}
