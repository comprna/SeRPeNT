#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <gsl/gsl_statistics_double.h>
#include <samtools/sam.h>
#include <debug.h>
#include <error.h>
#include <xcorr.h>

/*
 * Maximum length of a contig
 */
#define MAX_CONTIG_LENGTH 200
/*
 * Maximum number of contigs
 */
#define MAX_CONTIGS 100000 
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
#define MAX_REPLICATES 1 
/*
 * Maximum size, in base pairs, of a chromosome
 */
#define MAX_CHROMOSOME 250000000
/*
 * Default values for -c option
 */
#define CUTOFF 0.99f
/*
 * Default value for -m option
 */
#define MIN_LEN 15
/*
 * Default value for -M option
 */
#define MAX_LEN 30
/*
 * Default value for -s option
 */
#define SPACING 20
/*
 * Default value for -r option
 */
#define MIN_READS 10.0f
/*
 * Alignment in forward/watson strand
 */
#define FWD_STRAND 0
/*
 * Alignment in reverse/crick strand
 */
#define REV_STRAND 1
/*
 * Alignment is valid
 */
#define VALID_ALIGNMENT 1
/*
 * Alignment is not valid
 */
#define INVALID_ALIGNMENT 0
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
#define CROSSCOR_SUFFIX "crosscor.dat"

/*
 * Struct for handling command line arguments
 */
typedef struct {
  char output_f_path[MAX_PATH];
  char replicate_f_path[MAX_REPLICATES][MAX_PATH];
  int min_len;
  int max_len;
  double cutoff;
  int spacing;
  double min_reads;
} args_struct;

/*
 * Struct for handling sRNA profiles
 */
typedef struct {
  double *profile;
  char chromosome[MAX_FEATURE];
  int start;
  int end;
  int length;
  char strand;
} profile_struct;

/*
 * Struct for handling BAM alighment
 */
typedef struct {
  char *chromosome;
  int32_t start;
  int32_t end;
  char strand;
  char valid;
} alignment_struct;

/*
 * Struct for handling contigs
 */
typedef struct {
  int start;
  int end;
  double *profile;
  char chromosome[MAX_FEATURE];
} contig_struct;
