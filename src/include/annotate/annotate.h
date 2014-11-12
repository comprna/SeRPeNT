#include <annotate/xcorr.h>
#include <annotate/iofile.h>
#include <annotate/hierarchical.h>
#include <annotate/paramclust.h>
#include <annotate/profilemap.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

int annotate_sc(int argc, char** argv);
