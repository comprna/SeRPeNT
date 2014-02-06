#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <debug.h>
#include <gnuplot_i.h>
#include <error.h>

/*
 * Output modes
 */
#define OUTPUT_SCREEN 1
#define OUTPUT_PNG 2
/*
 * Maximum path length
 */
#define MAX_PATH 500
/*
 * Maximum identifier length
 */
#define MAX_FEATURE 50
/*
 * Maximum number of identifiers
 */
#define MAX_NUMBER_FEATURES 5

/*
 * Struct for handling command line arguments
 */
typedef struct {
  int output_mode;
  char output_png_path[MAX_PATH];
  char input_profile_path[MAX_PATH];
  char identifiers[MAX_NUMBER_FEATURES][MAX_FEATURE];
} args_struct;
