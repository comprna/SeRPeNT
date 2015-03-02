#ifndef ERROR_H
#define ERROR_H

/*
 * ERROR : Bad syntax message for plotprofiles
 */
#define ERR_MSG_PLOTPROFILES "Tool      : plotprofiles\n\
Version   : v1.0\n\
Summary   : A tool to draw sRNA profiles discovered using the sdrnaprofiles tool\n\n\
Usage     : plotprofiles [OPTIONS] <profiles.dat> <profile_name_1> [<profile_name_2> ... <profile_name_n>]\n\n\
Output    : plotprofiles will open a X11 window displaying the profile or will generate a png file\n\n\
Options   :\n\
           -o    PNG file.\n\
                 The plot will be printed in the specified PNG file.\n\
                 - Default behaviour is to open the plot in a system window"
/*
 * ERROR : Bad syntax
 */
#define ERR_HELP_MSG "Please type <srnap -h> or <srnap --help> for help"
/*
 * ERROR : profiles - bad syntax
 */
#define ERR_PROFILES_HELP_MSG "Please type <srnap profiles -h> or <srnap profiles --help> for help"
/*
 * ERROR : clusters - bad syntax
 */
#define ERR_ANNOTATE_HELP_MSG "Please type <srnap annotate -h> or <srnap annotate --help> for help"
/*
 * ERROR : BED Format
 */
#define ERR_BED "BED format file:\n\
    chromosome  start_position  end_position  name  score  strand   ..."
/*
 * ERROR : Invalid argument for -c option
 */
#define ERR_INVALID_c_VALUE "Invalid argument for option -c"
/*
 * ERROR : Invalid argument for -r option
 */
#define ERR_INVALID_r_VALUE "Invalid argument for option -r"
/*
 * ERROR : Invalid argument for -i option
 */
#define ERR_INVALID_i_VALUE "Invalid argument for option -i"
/*
 * ERROR : Invalid argument for -t option
 */
#define ERR_INVALID_t_VALUE "Invalid argument for option -t"
/*
 * ERROR : Invalid argument for -a option
 */
#define ERR_INVALID_a_VALUE "Invalid argument for option -a"
/*
 * ERROR : Invalid argument for -p option
 */
#define ERR_INVALID_p_VALUE "Invalid argument for option -p"
/*
 * ERROR : Invalid minimum contig length
 */
#define ERR_INVALID_m_VALUE "Minimum contig length must be an integer number greater than 5"
/*
 * ERROR : Invalid maximum contig length
 */
#define ERR_INVALID_M_VALUE "Maximum contig length must be an integer number greater then the minimum length"
/*
 * ERROR : Invalid spacing value
 */
#define ERR_INVALID_s_VALUE "Spacing must be an integer number greater or equal to 0"
/*
 * ERROR : Invalid minimum height value
 */
#define ERR_INVALID_h_VALUE "Minimum height must be an integer number greater than 0"
/*
 * ERROR : Invalid argument for -t option
 */
#define ERR_INVALID_tr_VALUE "Trimming must be an integer number greater than 0"
/*
 * ERROR : Unrecognized option or argument missing
 */
#define ERR_INVALID_ARGUMENT "Unrecognized option or argument missing"
/*
 * ERROR : Minimum length is higher or equal than maximum length
 */
#define ERR_INVALID_LENGTH "Minimum length is higher or equal than maximum length"
/*
 * ERROR : Invalid number of arguments
 */
#define ERR_INVALID_NUMBER_ARGUMENTS "Invalid number of arguments"
/*
 * ERROR : Invalid number of replicates
 */
#define ERR_INVALID_NUMBER_REPLICATES "Invalid number of arguments. Only up to 10 replicates allowed"
/*
 * ERROR : Unrecognized subcommand
 */
#define ERR_INVALID_SUBCOMMAND "Subcommand not recognized"
/*
 * ERROR : Invalid replicate number
 */
#define ERR_INVALID_REPLICATE_NUMBER "Replicate number for option -r is greater than the number of actual replicates"
/*
 * ERROR : Cannot read input file
 */
#define ERR_INPUT_F_NOT_READABLE "Profiles file does not exist or is not readable"
/*
 * ERROR : Query file is not formatted as a BED file
 */
#define ERR_QUERY_F_NOT_FORMATTED "Query file is not a BED file"
/*
 * ERROR : Cannot read BAM files
 */
#define ERR_BAM_F_NOT_READABLE "BAM file do not exist or is not readable"
/*
 * ERROR : Cannot read BAM header file
 */
#define ERR_BAM_F_NOT_HEADER "Failed to read the header from BAM file"
/*
 * ERROR : BAM truncated
 */
#define ERR_BAM_F_TRUNCATED "BAM file is truncated or ill-formed"
/*
 * ERROR : Cannot read replicate files
 */
#define ERR_REPLICATE_F_NOT_READABLE "Replicate files do not exist or are not readable"
/*
 * ERROR : Cannot write in output file
 */
#define ERR_OUTPUT_F_NOT_WRITABLE "Output path is not valid or output file is not writable"
/*
 * ERROR : Cannot re-allocate memory
 */
#define ERR_REALLOC_FAILED "Not enough memory or memory is corrupted"
/*
 * ERROR : Cannot read profile file
 */
#define ERR_PROFILE_F_NOT_READABLE "Profile file does not exist or is not readable"
/*
 * ERROR : Cannot read annotation file
 */
#define ERR_ANNOTATION_F_NOT_READABLE "Annotation file does not exist or is not readable"
/*
 * ERROR : Invalid minimum read length value
 */
#define ERR_INVALID_MIN_READ_LEN_VALUE "Minimum read length must be a number greater or equal to 0"
/*
 * ERROR : Invalid option for -f value
 */
#define ERR_INVALID_f_VALUE "Invalid argument for option -f"

#endif
