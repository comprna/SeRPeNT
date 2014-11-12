#ifndef HELP_H
#define HELP_H

#define GENERAL_HELP_MSG "srnap - A suite of tools for ncRNA profiling and classification\n\n\
The srnap subcommands include:\n\n\
[ Profiling tools ]\n\
   profiles    ncRNA discovery and profiling from small RNA-Seq data\n\
   annotate    ncRNA clustering, classification and annotation from profile data\n\n\
[ General help ]\n\
   -h, --help     Print this help menu\n\
   -v, --version  What version of srnap are you using?"

#define PROFILES_HELP_MSG "Tool      : profiles\n\n\
Summary   : ncRNA discovery and profiling from small RNA-Seq data\n\n\
Usage     : srnap profiles [OPTIONS] replicate_1.bam ... replicate_n.bam output_folder\n\n\
Options   :\n\
           -c, --contigs     Profile definition\n\
                             Format is <minlen:maxlen:spacing:minheight:trimming>, where:\n\
                               - <minlen> is the minimum length of the contig. Contigs shorter than <minlen> are not reported. Must be > 5.\n\
                               - <maxlen> is the maximum length of the contig. Contigs longer than <maxlen> are not reported. Must be > minlen.\n\
                               - <spacing> is the maximum distance between contigs. Contigs separated by <spacing> or less bp are considered the same contig. Must be >= 0.\n\
                               - <minheight> is the minimum number of piled-up reads. Contigs that have less than <minheight> piled-up reads are not reported. Must be > 0.\n\
                               - <trimming> is the number of bases to trim. Bases in the 5' and 3' sites of the contig that have less than <trimming> piled-up reads are trimmed. Must be >= 0.\n\
                             [ Default is 16:200:20:50:5 ]\n\n\
           -r, --replicates  Replicates treatment\n\
                             Format is <pool> | <mean> | <replicate:number>, where:\n\
                               - <pool> : Profiles are built by pooling the reads of all the replicates.\n\
                               - <mean> : Profiles are built by averaging the reads of all the replicates.\n\
                               - <replicate:number> : Profiles are built by using only the reads of one replicate.\n\
                             [ Default is pool ]\n\n\
           -i, --ic          Irreproducibility control\n\
                             Format is <method:cutoff> | <common> | <none>, where:\n\
                               - <method> is the irreproducibility control method. Options are:\n\
                                 - sere : Single-parameter quality control. (Schulze et al. BMC Genomics 2012)\n\
                                 - idr  : Irreproducibility discovery rate. (Li et al. The Annals of Applied Statistics 2011)\n\
                               - <cutoff> is the cutoff value. Contigs that have an irreproducibility score higher than <cutoff> are not reported.\n\
                               - <common> : Contigs that do not overlap in all the replicates are not reported.\n\
                               - <none>   : All contigs are reported.\n\
                             [ Default is common ]\n\n\
Output    :\n\
           output_folder/profiles.dat : List of ncRNA profiles with per-base heights\n\
           output_folder/contigs.dat  : List of unfiltered contigs\n\n\
Example   :\n\
           srnap profiles -c 20:200:20:50:5 -r pool -i sere:2"

#define ANNOTATE_HELP_MSG "Tool      : annotate\n\n\
Summary   : ncRNA clustering, classification and annotation from profile data\n\n\
Usage     : srnap annotate [OPTIONS] profiles_file.dat output_folder\n\n\
Options   : \n\
            -a, --afile  Annotation file\n\
                         Format is <annotation_file>, where:\n\
                           - <annotation_file> is a BED file with annotated features\n\
                         [ No default value ]\n\n\
            -p, --pfile  Additional profiles file\n\
                         Format is <species:profiles_file>, where:\n\
                           - <species> is the name of the species where the profiles comes from. e.g. hsap.\n\
                           - <profiles_file> is a profiles file\n\
                         [ No default value ]\n\n\
            -c, --cutoff  Cutoff for label propagation\n\
                          Format is <cutoff>, where:\n\
                            - <cutoff> is the cutoff of the hierarchical clustering distance\n\
                          [ Default is 0.01 ]\n\n\
Output    :\n\
           output_folder/crosscorr.dat    : List of cross correlation coefficients between pairs of profiles\n\
           output_folder/clusters.neWick  : List of clustered profiles\n\
           output_folder/annotation.bed   : List of annotated features in BED file (only if annotation file is provided)\n\n\
Example   :\n\
           srnap annotate -c 0.05 -a hsap_micrornas.bed profiles.dat output_dir\n\
           srnap annotate -a hsap_micornas.bed -p hsap:hsap_profiles.dat profiles.dat output_dir"
#endif
