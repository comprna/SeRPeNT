#ifndef HELP_H
#define HELP_H

#define GENERAL_HELP_MSG "srnap - A suite of tools for ncRNA profiling, classification and analysis\n\n\
The srnap subcommands include:\n\n\
[ Profiling tools ]\n\
   profiles    ncRNA discovery and profiling from small RNA-Seq data\n\
   annotate    ncRNA clustering, classification and annotation from profile data\n\
   diffproc    ncRNA differential processing from profile and clustering data between two conditions\n\n\
[ General help ]\n\
   -h          Print this help menu\n\
   -v          What version of srnap are you using?"

#define PROFILES_HELP_MSG "Tool      : profiles\n\n\
Summary   : ncRNA discovery and profiling from small RNA-Seq data\n\n\
Usage     : srnap profiles [OPTIONS] replicate_1.bam ... replicate_n.bam output_folder\n\n\
Options   :\n\
           -f   Read filtering\n\
                Format is <minreadlen>, where:\n\
                  - <minreadlen> is the minimum required length for the reads. Reads shorter than <readlen> are discarded. If 0, no reads are discarded. Must be >= 0.\n\
                [ Default is 0 ]\n\n\
           -i   Irreproducibility control for contigs\n\
                Format is <method:cutoff> | <common> | <none>, where:\n\
                  - <method> is the irreproducibility control method. Options are:\n\
                    - sere : Single-parameter quality control. (Schulze et al. BMC Genomics 2012)\n\
                    - idr  : Non-parametric irreproducibility discovery rate. (Dobin et al. Bioinformatics 2013)\n\
                  - <cutoff> is the cutoff value. Contigs that have an irreproducibility score higher than <cutoff> are not reported. Must be > 0.\n\
                  - <common> : Contigs that do not overlap in all the replicates are not reported.\n\
                  - <none>   : All contigs are reported.\n\
                [ Default is common ]\n\n\
           -r   Replicates treatment\n\
                Format is <pool> | <mean> | <replicate:repnumber>, where:\n\
                  - <pool> : Profiles are built by pooling the reads of all the replicates.\n\
                  - <mean> : Profiles are built by averaging the reads of all the replicates.\n\
                  - <replicate:repnumber> : Profiles are built by using only the reads of the <repnumber> replicate.\n\
                [ Default is pool ]\n\n\
           -t   Trimming\n\
                Format is <trim_threshold:trim_min:trim_max>\n\
                  - <trim_percentage> is the trimming threshold. Nucleotides in the ends of the profile having less than <trim_precentage> percent of reads compared to the\n\
                                      maximum height will be trimmed.\n\
                  - <trim_min> is the trimming minimum height. All nucleotides in both ends of the profile having less than <trim_min> reads will be trimmed.\n\
                  - <trim_max> is the trimming maximum height. No nucleotides in both ends of the profile having more than <trim_max> reads will be trimmed.\n\
                [ Default is 0.1:2:10 ]\n\n\
           -p   Profile definition\n\
                Format is <minlen:maxlen:spacing:minheight:trimming>, where:\n\
                  - <minlen> is the minimum length of the profile after trimming. Profiles shorter than <minlen> are not reported. Must be > 5.\n\
                  - <maxlen> is the maximum length of the profile after trimming. Profiles longer than <maxlen> are not reported. Must be >= minlen.\n\
                  - <spacing> is the maximum distance between profiles. Profiles separated by <spacing> or less bp are merged into one single profile. Must be >= 0.\n\
                  - <minheight> is the minimum number of piled-up reads. Profiles that have less than <minheight> piled-up reads are not reported. Must be > 0.\n\
                [ Default is 16:200:20:50 ]\n\n\
Output :\n\
           output_folder/profiles.dat : List of ncRNA profiles with per-base heights\n\
           output_folder/contigs.dat  : List of unfiltered contigs\n\n\
Examples :\n\
           srnap profiles -f 20 -i sere:2 -r pool -t 0.1:5:20 -p 20:200:39:100 replicate1.bam replicate2.bam output_dir"

#define ANNOTATE_HELP_MSG "Tool      : annotate\n\n\
Summary   : ncRNA clustering, classification and annotation from profile data\n\n\
Usage     : srnap annotate [OPTIONS] profiles_file.dat output_folder\n\n\
Options   :\n\
            -a   Annotation file\n\
                 Format is <annotation_file>, where:\n\
                   - <annotation_file> is a BED file with annotated features\n\
                 [ No default value ]\n\n\
            -d   Additional profiles file\n\
                 Format is <species:profiles_file>, where:\n\
                   - <species> is the name of the species where the profiles comes from. e.g. hsap.\n\
                   - <profiles_file> is a profiles file\n\
                 [ No default value ]\n\n\
            -o   Overlapping parameters\n\
                 Format is <feature_to_profile:profile_to_feature>, where:\n\
                   - <feature_to_profile> is the percentage of nucleotides from the feature that overlap the profile\n\
                   - <profile_to_feature> is the percentage of nucleotides from the profile that overlap the feature\n\
                 [ Default is 0.9:0.5 ]\n\n\
            -c   Cutoff for branching\n\
                 Format is <cutoff>, where:\n\
                   - <cutoff> is the distance threshold for branching the hierarchical clustering solution\n\
                 If no -c option is specified, srnap calculates the optimal cutoff\n\
                 [ No default value ]\n\n\
            -x   Distance file\n\
                 Format is <distance_file>, where:\n\
                   - <distance_file> is the file with pairwise distances between profiles\n\
                 When -x option is specified, distances are not calculated and directly taken from the provided file\n\
                 [ No default value ]\n\n\
Output    :\n\
            output_folder/crosscorr.dat    : List of distances between pairs of profiles (only if no distance file is provided)\n\
            output_folder/clusters.neWick  : List of clustered profiles\n\
            output_folder/annotation.bed   : List of annotated features in BED file (only if annotation file is provided)\n\n\
Examples  :\n\
            srnap annotate -a hsap_micrornas.bed profiles.dat output_dir\n\
            srnap annotate -a hsap_micrornas.bed -x crosscor.dat profiles.dat output_dir"

#define DIFFPROC_HELP_MSG "Tool      : diffproc\n\n\
Summary   : ncRNA differential processing from profile and clustering data between two conditions\n\n\
Usage     : srnap diffproc [OPTIONS] profiles_file_1.dat clustering_file_1.bed profile_file_2.dat clustering_file_2.bed output_folder\n\n\
Options   :\n\
            -g   p-value and overlap\n\
                 Format is <pvalue:overlap> where:\n\
                   - <pvalue> is the p-value threshold for filtering differentially processed profiles\n\
                   - <overlap> is the overlap threshold for filtering differentially processed clusters\n\
                 [ Default is 0.01:0.5 ]\n\n\
Output    :\n\
            output_folder/diffprofiles.dat : List of differentially processed profiles\n\
            output_folder/diffclusters.dat : List of differentially processed clusters\n\n\
Examples  :\n\
            srnap diffproc wild_type/profiles.dat wild_type/annotation.bed treated/profiles.dat treated/annotation.bed output_dir\n\
            srnap diffproc -g 0.01:0.2 wild_type/profiles.dat wild_type/annotation.bed treated/profiles.dat treated/annotation.bed output_dir"
#endif
