#include <core/structs.h>

/*
 * parse_command_line
 *   Parses the command line
 *
 * @arg int argc
 *   Number of arguments in the command line
 * @arg char** argv
 *   Array containing the command line arguments
 * @arg error_message
 *   Pointer to a char array where to store the error message
 * @arg args_a_struct* arguments
 *   Pointer to the argument handler
 *
 * @return -1 if an error occurred. 0 otherwise.
 */
int parse_command_line_c(int argc, char** argv, char** error_message, args_a_struct* arguments);

/*
 * parse_annotation_parameters
 *   Parses the string defining the annotation file option
 *
 * @arg char* option
 *   String defining the replicates treatment options
 * @arg char** error_message
 *   Pointer to a char array where to store the error message
 * @args args_a_struct* arguments
 *   Pointer to the argument handler
 *
 * @return -1 if an error occurred. 0 otherwise.
 */
int parse_annotation_parameters(char* option, char** error_message, args_a_struct* arguments);

/*
 * parse_additional_profiles_parameters
 *   Parses the string defining the additional profiles file option
 *
 * @arg char* option
 *   String defining the replicates treatment options
 * @arg char** error_message
 *   Pointer to a char array where to store the error message
 * @args args_a_struct* arguments
 *   Pointer to the argument handler
 *
 * @return -1 if an error occurred. 0 otherwise.
 */
int parse_additional_profiles_parameters(char* option, char** error_message, args_a_struct* arguments);

/*
 * parse_cutoff_parameters
 *   Parses the string defining the cutoff option
 *
 * @arg char* option
 *   String defining the replicates treatment options
 * @arg char** error_message
 *   Pointer to a char array where to store the error message
 * @args args_a_struct* arguments
 *   Pointer to the argument handler
 *
 * @return -1 if an error occurred. 0 otherwise.
 */
int parse_threshold_parameters(char* option, char** error_message, args_a_struct* arguments);
