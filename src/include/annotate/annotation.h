#include <core/structs.h>
#include <annotate/strmap.h>

void xcorr_annotate(annotation_struct** xcorr, int nprofiles, profile_struct_annotation* profiles);
void cluster_annotate(int nclusters, int nprofiles, profile_struct_annotation* profiles);
void xspeciescorr_annotate(annotation_struct** xcorr, int nprofiles, profile_struct_annotation* profiles, int naprofiles, profile_struct_annotation* aprofiles, double cutoff);
