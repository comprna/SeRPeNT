#ifndef ERROR_H
#define ERROR_H

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
 * ERROR : Unrecognized subcommand
 */
#define ERR_INVALID_SUBCOMMAND "Subcommand not recognized"
/*
 * ERROR : Unrecognized option or argument missing
 */
#define ERR_INVALID_ARGUMENT "Unrecognized option or argument missing"
/*
 * ERROR : Invalid number of arguments
 */
#define ERR_INVALID_NUMBER_ARGUMENTS "Invalid number of arguments"

/*
 * ERROR : Invalid number of replicates
 */
#define ERR_INVALID_NUMBER_REPLICATES "Invalid number of arguments. Only up to 10 replicates allowed"

/*
 * ERROR : Invalid option for -f value
 */
#define ERR_INVALID_f_VALUE "Invalid argument for option -f"
/*
 * ERROR : Invalid minimum read length <minreadlen> value
 */
#define ERR_INVALID_minreadlen_VALUE "Minimum read length <minreadlen> must be a number greater or equal to 0"

/*
 * ERROR : Invalid argument for -r option
 */
#define ERR_INVALID_r_VALUE "Invalid argument for option -r"
/*
 * ERROR : Invalid replicate number
 */
#define ERR_INVALID_repnumber_VALUE "Replicate number <repnumber> must be an integer number between 1 and the number of actual replicates"

/*
 * ERROR : Invalid argument for -i option
 */
#define ERR_INVALID_i_VALUE "Invalid argument for option -i"

/*
 * ERROR : Invalid argument for -p option
 */
#define ERR_INVALID_p_VALUE "Invalid argument for option -p"
/*
 * ERROR : Invalid minimum contig length
 */
#define ERR_INVALID_minlen_VALUE "Minimum contig length <minlen> must be an integer number greater than 5"
/*
 * ERROR : Invalid maximum contig length
 */
#define ERR_INVALID_maxlen_VALUE "Maximum contig length <maxlen> must be an integer number greater than the minimum length"
/*
 * ERROR : Invalid spacing value
 */
#define ERR_INVALID_spacing_VALUE "Maximum distance between contigs <spacing> must be an integer number greater or equal to 0"
/*
 * ERROR : Invalid minimum height value
 */
#define ERR_INVALID_minheight_VALUE "Minimum height <minheight> must be an integer number greater than 0"

/*
 * ERROR : Invalid argument for -t option
 */
#define ERR_INVALID_t_VALUE "Invalid argument for option -t"
/*
 * ERROR : Invalid trimming threshold
 */
#define ERR_INVALID_trimthreshold_VALUE "Trimming threshold <trim_threshold> must be a number between 0 and 1"
/*
 * ERROR : Invalid trimming minimum value
 */
#define ERR_INVALID_trimin_VALUE "Trimming minimum height <trim_min> must be an integer number greater or equal to 0"
/*
 * ERROR : Invalid maximum trimming value
 */
#define ERR_INVALID_trimax_VALUE "Trimming maximum height <trim_max> must be an integer number greater or equal than <trim_min>"

/*
 * ERROR : Invalid argument for -a option
 */
#define ERR_INVALID_a_VALUE "Invalid argument for option -a"

/*
 * ERROR : Invalid argument for -x option
 */
#define ERR_INVALID_x_VALUE "Invalid argument for option -x"

/*
 * ERROR : Invalid argument for -c option
 */
#define ERR_INVALID_c_VALUE "Invalid argument for option -c"

/*
 * ERROR : BED Format
 */
#define ERR_BED "BED format file:\n\
    chromosome  start_position  end_position  name  score  strand   ..."
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

#endif
