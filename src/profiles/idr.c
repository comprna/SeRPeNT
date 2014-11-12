#include <profiles/idr.h>

/*
 * update_hash
 *   Auxiliar function for updating hash table
 */
void update_hash(struct llist_struct **elems, int index, int field, int n_replicas, int n_contigs)
{
  int idx = index % (n_replicas * n_contigs);

  // No index in the hash table
  if (elems[idx] == NULL) {
    elems[idx] = (struct llist_struct*) malloc(sizeof(struct llist_struct));
    elems[idx]->index = index;
    elems[idx]->absolute = 1;
    elems[idx]->conditional = 0;
    elems[idx]->next = NULL;
  }

  // Index in the hash table
  else {
    struct llist_struct* pointer = elems[idx];
    while((pointer->index != index) && (pointer->next != NULL))
      pointer = pointer->next;

    // Index in the collision list
    if (pointer->index == index) {
      if (field == ABSOLUTE) pointer->absolute++;
      if (field == CONDITIONAL) pointer->conditional++;
    }

    // Index not in the collision list
    else if (pointer->next == NULL) {
      pointer->next = (struct llist_struct*) malloc(sizeof(struct llist_struct));
      pointer->next->index = index;
      pointer->next->absolute = 1;
      pointer->next->conditional = 0;
      pointer->next->next = NULL;
    }
  }
}

/*
 * get_hash
 *   Auxiliar function for retrieving value from hash table
 */
struct llist_struct* get_hash(struct llist_struct **elems, int index, int n_replicas, int n_contigs)
{
  int idx = index % (n_replicas * n_contigs);
  struct llist_struct* pointer = elems[idx];

  while((pointer != NULL) && (pointer->index != index)) pointer = pointer->next;

  return(pointer);
}

/*
 * destroy_hash
 *   Auxiliar function for destroying elements in a hash table
 */
void destroy_hash(struct llist_struct *element)
{
  if (element == NULL)
    return;

  if (element->next != NULL)
    destroy_hash(element->next);

  free(element);
}

/*
 * calculate_sere_scores
 *
 * @see include/profiles/idr.h
 */
void calculate_sere_scores(profile_struct* profiles, int n_contigs, int n_replicates, double cutoff)
{
  double total_reads_r[MAX_REPLICATES];
  double total_reads = 0;
  int i, j;

  for (i = 0; i < n_replicates; i++) total_reads_r[i] = 0;
  for (i = 0; i < n_contigs; i++) for(j = 0; j < n_replicates; j++) total_reads_r[j] += (double) profiles[i].nreads[j];
  for (i = 0; i < n_replicates; i++) total_reads += total_reads_r[i];

  for (i = 0; i < n_contigs; i++) {
    double sisqr = 0;
    double sumyij = 0;
    double ei = 0;

    for (j = 0; j < n_replicates; j++) ei += (double) profiles[i].nreads[j];

    for (j = 0; j < n_replicates; j++) {
      double yij = (double) profiles[i].nreads[j];
      double estyij = ei * (total_reads_r[j] / total_reads);
      sumyij += pow(yij - estyij, 2) / estyij;
    }

    sisqr = (1 / (n_replicates - 1)) * sumyij;

    profiles[i].idr_score = sqrt(sisqr);

    if (profiles[i].idr_score > cutoff) profiles[i].valid = 0;
  }
}


/*
 * calculate_common_scores
 *
 * @see include/profiles/idr.h
 */
void calculate_common_scores(profile_struct* profiles, int n_contigs, int n_replicates)
{
  int i, j;

  for (i = 0; i < n_contigs; i++) {
    int nreps = 0;
    for (j = 0; j < n_replicates; j++) if (profiles[i].nreads[j] > 0) nreps++;
      if (nreps < n_replicates)
        profiles[i].valid = 0;
  }
}


/*
 * calculate_npidr_scores
 *
 * @see include/profiles/idr.h
 */
void calculate_npidr_scores(profile_struct* profiles, int n_contigs, int n_replicates, double cutoff)
{
  int i, bin_size = 1, pool_type = 2;
  struct llist_struct** elems = (struct llist_struct**) malloc(sizeof(struct llist_struct*) * n_replicates * n_contigs);
  long* result = (long*) malloc(sizeof(long) * n_contigs);

  // Calculate
  for(i = 0; i < n_contigs; i++) {
    int n_zeroes = 0;
    long res = 0;
    int j = 0, k = 0;

    for (k = 0; k < n_replicates; k++) {
      profiles[i].nreads[k] = profiles[i].nreads[k] / bin_size;

      if (profiles[i].nreads[k] == 0)
        n_zeroes++;
      else
        j = k;

      if (pool_type == 1) res += profiles[i].nreads[k];
      else if (profiles[i].nreads[k] > res) res = profiles[i].nreads[k];

      update_hash(elems, profiles[i].nreads[k], ABSOLUTE, n_replicates, n_contigs);
    }

    if(n_zeroes == (n_replicates - 1))
      update_hash(elems, profiles[i].nreads[j], CONDITIONAL, n_replicates, n_contigs);

    result[i] = res;
  }

  // Assign npIDR scores
  for (i = 0; i < n_contigs; i++) {
    struct llist_struct* element = get_hash(elems, result[i], n_replicates, n_contigs);

    if (element == NULL)
      profiles[i].idr_score = 0;
    else if (element->absolute > 0)
      profiles[i].idr_score = ((double) element->conditional) / ((double) element->absolute);
    else
      profiles[i].idr_score = 0;
  }

  // Free structures
  for (i = 0; i < (n_replicates * n_contigs); i++) destroy_hash(elems[i]);
  free(elems);
  free(result);
}
