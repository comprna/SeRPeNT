#ifndef STRUCTS_H
#define STRUCTS_H

#include <core/constants.h>

/*
 * Struct for handling profiles command line arguments
 */
typedef struct {
  char output_f_path[MAX_PATH];
  char replicate_f_path[MAX_REPLICATES][MAX_PATH];
  int min_len;
  int max_len;
  int spacing;
  int min_reads;
  int replicate_treat;
  int replicate_number;
  int number_replicates;
  int idr_method;
  double idr_cutoff;
  int read_length;
  int min_read_len;
  double trim_threshold;
  int trim_min;
  int trim_max;
} args_p_struct;

/*
 * Struct for handling clusters command line arguments
 */
typedef struct {
  char output_f_path[MAX_PATH];
  char profiles_f_path[MAX_PATH];
  int annotation;
  char annotation_f_path[MAX_ANNOTATIONS][MAX_PATH];
  char species[MAX_FEATURE];
  int additional_profiles;
  char additional_profiles_f_path[MAX_PATH];
  double cluster_cutoff;
  double overlap_ftop;
  double overlap_ptof;
  int correlations;
  char correlations_f_path[MAX_PATH];
} args_a_struct;

/*
 * Struct for handling sRNA profiles
 */
typedef struct {
  double *profile;
  char chromosome[MAX_FEATURE];
  int start;
  int end;
  int length;
  int32_t strand;
  int* nreads;
  int32_t valid;
  int32_t free;
  double idr_score;
  int tstart;
  int tend;
} profile_struct;

typedef struct {
  double *profile;
  char chromosome[MAX_FEATURE];
  int start;
  int end;
  int length;
  int32_t strand;
  int additional;
  char species[MAX_FEATURE];
  char annotation[MAX_FEATURE];
  char tmp_annotation[MAX_FEATURE];
  int cluster;
  double anscore;
  double max_height;
  double mean;
  double variance;
  double noise[MAX_PROFILE_LENGTH];
  int32_t category;
} profile_struct_annotation;

/*
 * Struct for handling BAM alighment
 */
typedef struct {
  char chromosome[MAX_FEATURE];
  int32_t start;
  int32_t end;
  int32_t strand;
  int32_t valid;
  int replicate;
  int nreads;
} alignment_struct;

/*
 * Struct for handling contigs
 */
typedef struct {
  int start;
  int end;
  double *profile;
  char chromosome[MAX_FEATURE];
  int* nreads;
} contig_struct;

/*
 * Structure for handling binary heap of alignments
 */
typedef struct {
  alignment_struct** alignments;
  int* heap;
  int heap_size;
} heap_struct;

/*
 * Typedef for wc utility
 */
typedef unsigned long count_t;

/*
 * Structure for handling nodes of interval trees
 */
struct itnode_struct {
  int low;
  int high;
  profile_struct_annotation* profile;
  struct itnode_struct* left;
  struct itnode_struct* right;
  int height;
};

/*
 * Typedef for interval tree node
 */
typedef struct itnode_struct itnode_struct;

/*
 * Structure for handling elements of profiles maps
 */
typedef struct {
  char identifier[MAX_FEATURE];
  itnode_struct* root;
} map_element_struct;

/*
 * Structure for handling profile maps
 */
typedef struct {
  map_element_struct* elements;
  int size;
} map_struct;

/*
 * Structure for handling features in a BED file
 */
typedef struct {
  char chromosome[MAX_FEATURE];
  int start;
  int end;
  char name[MAX_FEATURE];
  int strand;
} feature_struct;

/*
 * Structure for handling linked list for collisions in hash table
 * Used in npIDR method
 */
struct llist_struct {
  long index;
  int conditional;
  int absolute;
  struct llist_struct* next;
};

/*
 * Hierarchical cluster node
 */
typedef struct {
  int left;
  int right;
  int parent;
  double distance;
  int lleafs;
  int rleafs;
  int visited;
} hcnode_struct;

/*
 * Struct for fast SERE IDR calculation
 */
typedef struct {
  long* reads_per_replicate;
  long total_reads;
} sere_struct;

/*
 * Struct for fast SERE IDR calculation
 */
typedef struct {
  struct llist_struct** elems;
  long* result;
  long nelems;
} npidr_struct;

/*
 * Struct for profile annotation
 */
typedef struct {
  double score;
  int index_i;
  int index_j;
} annotation_struct;
#endif
