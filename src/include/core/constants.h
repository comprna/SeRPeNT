#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <utils/version.h>
#include <utils/help.h>
#include <utils/error.h>

/*
 * Default value for minlen parameter
 */
#define MIN_LEN 15

/*
 * Default value for maxlen parameter
 */
#define MAX_LEN 30

/*
 * Default value for spacing parameter
 */
#define SPACING 20

/*
 * Default value for minheight parameter
 */
#define MIN_READS 10.0f

/*
 * Default value for trimming parameter
 */
#define TRIMMING 1

/*
 * Default value for IDR cutoff
 */
#define CUTOFF 2.0f

/*
 * Default value for Replicate number
 */
#define REPLICATE_NUMBER 1

/*
 * Replicate treatment options
 */
#define REPLICATE_POOL_STR "pool" // default value
#define REPLICATE_MEAN_STR "mean"
#define REPLICATE_REPLICATE_STR "replicate"
#define REPLICATE_POOL 0
#define REPLICATE_MEAN 1
#define REPLICATE_REPLICATE 2

/*
 * IDR method options
 */
#define IDR_COMMON_STR "common" // default value
#define IDR_NONE_STR "none"
#define IDR_SERE_STR "sere"
#define IDR_IDR_STR "idr"
#define IDR_COMMON 0
#define IDR_NONE 1
#define IDR_SERE 2
#define IDR_IDR 3

/*
 * Maximum length of a contig
 */
#define MAX_CONTIG_LENGTH 200

/*
 * Maximum number of contigs
 */
#define MAX_CONTIGS 10000000

/*
 * Maximum size, in chars, of an error message
 */
#define MAX_ERR_MSG 100

/*
 * Maximum size, in chars, of a given path
 */
#define MAX_PATH 500

/*
 * Maximum size, in chars, of the 4th field (name) in the query BED file
 */
#define MAX_FEATURE 50

/*
 * Maximum number of replicates
 */
#define MAX_REPLICATES 10 

/*
 * Maximum number of alignments per heap
 */
#define MAX_ALIGN_HEAP 50000000

/*
 * Maximum read length
 */
#define MAX_READ_LENGTH 200

/*
 * Maximum number of block base pairs
 */
#define MAX_BLOCK 10000

/*
 * Alignment strand
 */
#define FWD_STRAND 0 // forward/watson
#define REV_STRAND 1 // reverse/crick

/*
 * Alignment validity
 */
#define VALID_ALIGNMENT 1
#define INVALID_ALIGNMENT 0

/*
 * Constants for npIDR method
 */
#define ABSOLUTE 0
#define CONDITIONAL 1 

/*
 * Path separator
 */
#ifdef __unix__
  #define PATH_SEPARATOR "/"
#else
  #define PATH_SEPARATOR "\\"
#endif

/*
 * Output file suffixes
 */
#define PROFILES_SUFFIX "profiles.dat"
#define CONTIGS_SUFFIX "contigs.dat"
#define CROSSCOR_SUFFIX "crosscor.dat"
#define CLUSTERS_SUFFIX "clusters.neWick"
#define ANNOTATION_O_SUFFIX "annotation.bed" 

/*
 * Maximum N limit for gaussian white noise generation
 */
#define MAX_GNOISE_N 20

/*
 * Cluster cutoff default value
 */
#define CLUSTER_CUTOFF 0.01f

/*
 * Condition for existence of annotation file
 */
#define ANNOTATION_CONDITION 0

/*
 * Condition for existence of additional profiles file
 */
#define ADDITIONAL_P_CONDITION 0

/*
 * Maximum number of annotation files
 */
#define MAX_ANNOTATIONS 10

#endif
