/*
 * SUCCESS : Program normally terminated
 */
#define SUCCESS "Success"

/*
 * ERROR : Bad syntax error message
 */
#define ERR_MSG "Tool      : sdRNAprofiles\n\
Version   : v1.0\n\
Summary   : A tool for the de-novo discovery of small RNAs and derived RNA molecules, from small RNA-Seq data\n\n\
Usage     : sdrnaprofiles [OPTIONS] replicate_1.bam ... replicate_n.bam output_file\n\n\
Options   :\n\
           -x     Minimum cross-correlation\n\
                  Threshold for the cross-correlation value of the same contig between different replicates.\n\
                  Must be between [0,1]\n\
                  - Default is 0.99\n\n\
           -m     Minimum length\n\
                  Minimum length required to report a contig.\n\
                  Must be > 5.\n\
                  - Default is 15\n\n\
           -M     Maximum length\n\
                  Maximum length required to report a contig.\n\
                  Must be > 6.\n\
                  If -m option is present, maximum length must be greater than minimum length.\n\
                  - Default is 30\n\
           -s     Spacing\n\
                  Maximum space required between two reads to be considered part of the same contig.\n\
                  Must be 0 or a positive number.\n\
                  - Default is 20\n\
           -t     Trimming threshold\n\
                  Trim nucleotides in the 5'- and 3'-ends of the contig that have -t reads or less.\n\
                  Must be 0 or a positive number.\n\
                  - Default is 2"
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
 * ERROR : BED Format
 */
#define ERR_BED "BED format file:\n\
    chromosome  start_position  end_position  name  score  strand   ..."
/*
 * ERROR : Invalid argument for -x option
 */
#define ERR_INVALID_x_VALUE "Invalid argument for option -x"
/*
 * ERROR : Invalid argument for -m option
 */
#define ERR_INVALID_m_VALUE "Invalid argument for option -m"
/*
 * ERROR : Invalid argument for -M option
 */
#define ERR_INVALID_M_VALUE "Invalid argument for option -M"
/*
 * ERROR : Invalid argument for -s option
 */
#define ERR_INVALID_s_VALUE "Invalid argument for option -s"
/*
 * ERROR : Invalid argument for -r option
 */
#define ERR_INVALID_r_VALUE "Invalid argument for option -r"
/*
 * ERROR : Invalid argument for -t option
 */
#define ERR_INVALID_t_VALUE "Invalid argument for option -t"
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
