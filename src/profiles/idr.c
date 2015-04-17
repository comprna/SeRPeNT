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
 * calculate_sere_partials
 * 
 * @see include/profiles/idr.h
 */
void calculate_sere_score(profile_struct* profile, int** reads_per_contig, int n_contigs, int n_replicates, double cutoff)
{
  long total_reads = 0;
  int i, j;
  double sisqr = 0;
  double sumyij = 0;
  long ei = 0;
  long* reads_per_replicate;

  reads_per_replicate = (long*) malloc(n_replicates * sizeof(long));

  for (i = 0; i < n_replicates; i++) reads_per_replicate[i] = 0;
  for (i = 0; i < n_contigs; i++) for(j = 0; j < n_replicates; j++) reads_per_replicate[j] += reads_per_contig[i][j];
  for (i = 0; i < n_replicates; i++) total_reads += reads_per_replicate[i];

  for (i = 0; i < n_replicates; i++) ei += profile->nreads[i];

  for (i = 0; i < n_replicates; i++) {
    double yij = (double) profile->nreads[i];
    double estyij = (double) ei * ((double) reads_per_replicate[i] / (double) total_reads);
    sumyij += pow(yij - estyij, 2) / estyij;
  }

  sisqr = ((double) 1 / (double) (n_replicates - 1)) * sumyij;

  profile->idr_score = sqrt(sisqr);

  if (profile->idr_score > cutoff) profile->valid = 0;

  free(reads_per_replicate);
}


/*
 * calculate_common_scores
 *
 * @see include/profiles/idr.h
 */
void calculate_common_score(profile_struct* profile, int n_replicates)
{
  int i;
  int nreps = 0;

  for (i = 0; i < n_replicates; i++)
    if (profile->nreads[i] > 0)
      nreps++;

  if (nreps < n_replicates) {
    profile->valid = 0;
    profile->idr_score = 0;
  }
  else
    profile->idr_score = 1;
}


/*
 * calculate_npidr_scores
 *
 * @see include/profiles/idr.h
 */
void calculate_npidr_score(profile_struct* profile, int** reads_per_contig, int index, int n_contigs, int n_replicates, double cutoff)
{
  int i;
  struct llist_struct** elems = (struct llist_struct**) malloc(sizeof(struct llist_struct*) * n_replicates * n_contigs);
  long* result = (long*) malloc(sizeof(long) * n_contigs);

  // Initialize
  for (i = 0; i < n_replicates * n_contigs; i++)
    elems[i] = NULL;

  // Calculate
  for(i = 0; i < n_contigs; i++) {
    int n_zeroes = 0;
    long res = 0;
    int j = 0, k = 0;

    for (k = 0; k < n_replicates; k++) {
      if (reads_per_contig[i][k] == 0)
        n_zeroes++;
      else
        j = k;

      if (reads_per_contig[i][k] > res) res = reads_per_contig[i][k];

      update_hash(elems, reads_per_contig[i][k], ABSOLUTE, n_replicates, n_contigs);
    }

    if(n_zeroes == (n_replicates - 1))
      update_hash(elems, reads_per_contig[i][j], CONDITIONAL, n_replicates, n_contigs);

    result[i] = res;
  }

  // Assign npIDR scores
  struct llist_struct* element = get_hash(elems, result[index], n_replicates, n_contigs);

  if (element == NULL)
    profile->idr_score = 0;
  else if (element->absolute > 0)
    profile->idr_score = ((double) element->conditional) / ((double) element->absolute);
  else
    profile->idr_score = 0;

  if (profile->idr_score > cutoff) profile->valid = 0;

  // Free structures
  for (i = 0; i < (n_replicates * n_contigs); i++) destroy_hash(elems[i]);
  free(elems);
  free(result);
}
