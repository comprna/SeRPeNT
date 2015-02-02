#include <core/structs.h>

/*
 * parse_contigs_parameters
 *   Parses the string defining the contig discovery options
 *
 * @arg char* option
 *   String defining the contig discovery options
 * @arg char** error_message
 *   Pointer to a char array where to store the error message
 * @args args_p_struct* arguments
 *   Pointer to the argument handler
 *
 * @return -1 if an error occurred. 0 otherwise.
 */
int parse_contigs_parameters(char* option, char** error_message, args_p_struct* arguments);

/*
 * parse_replicates_parameters
 *   Parses the string defining the replicates treatment options
 *
 * @arg char* option
 *   String defining the replicates treatment options
 * @arg char** error_message
 *   Pointer to a char array where to store the error message
 * @args args_p_struct* arguments
 *   Pointer to the argument handler
 *
 * @return -1 if an error occurred. 0 otherwise.
 */
int parse_replicates_parameters(char* option, char** error_message, args_p_struct* arguments);

/*
 * parse_irreproducibility_parameters
 *   Parses the string defining the irreproducibility check options
 *
 * @arg char* option
 *   String defining the irreproducibility check options
 * @arg char** error_message
 *   Pointer to a char array where to store the error message
 * @args args_p_struct* arguments
 *   Pointer to the argument handler
 *
 * @return -1 if an error occurred. 0 otherwise.
 */
int parse_irreproducibility_parameters(char* option, char** error_message, args_p_struct* arguments);

/*
 * parse_command_line
 *   Parses the command line
 *
 * @arg int argc
 *   Number of arguments in the command line
 * @arg char** argv
    Array containing the command line arguments
 * @arg error_message
 *   Pointer to a char array where to store the error message
 * @arg args_p_struct* arguments
 *   Pointer to the argument handler
 *
 * @return -1 if an error occurred. 0 otherwise.
 */
int parse_command_line_p(int argc, char** argv, char** error_message, args_p_struct* arguments);