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
 * @arg args_p_struct* arguments
 *   Pointer to the argument handler
 *
 * @return -1 if an error occurred. 0 otherwise.
 */
int parse_command_line_d(int argc, char** argv, char** error_message, args_d_struct* arguments);

/*
 * parse_filter_output_parameters
 *   Parses the string defining the output filtering options
 *
 * @arg char* option
 *   String defining the output filtering options
 * @arg char** error_message
 *   Pointer to a char array where to store the error message
 * @args args_p_struct* arguments
 *   Pointer to the argument handler
 *
 * @return -1 if an error occurred. 0 otherwise.
 */
int parse_filter_output_parameters(char* option, char** error_message, args_d_struct* arguments);
