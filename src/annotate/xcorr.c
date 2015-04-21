#include <annotate/xcorr.h>

/*
 * nxcorr
 *
 * @see src/include/annotate/xcorr.h
 */
double nxcorr(profile_struct_annotation* p1, profile_struct_annotation* p2)
{
  int index, lag, i, j, length;
  double rxx = 0, ryy = 0, rxy = 0, rnm, noise, nrxy;
  double* corr;

  srand(time(NULL));
  index = 0;
  length = abs(p1->length - p2->length) + 1;
  corr = (double*) malloc(length * sizeof(double));

  // profile1.length > profile2.length
  if (p1->length > p2->length) {
    for (i = 0; i < p1->length; i++)
      rxx += p1->profile[i] * p1->profile[i];

    for (lag = 0; lag < length; lag++) {
      ryy = 0; rxy = 0; j = 0; i = 0;

      while(i < lag) {
        noise = p2->noise[rand() % MAX_PROFILE_LENGTH];
        ryy += noise * noise;
        rxy += p1->profile[i] * noise;
        i++;
      }
      while(j < p2->length) {
        ryy += p2->profile[j] * p2->profile[j];
        rxy += p1->profile[i] * p2->profile[j];
        i++; j++;
      }
      while(i < p1->length) {
        noise = p2->noise[rand() % MAX_PROFILE_LENGTH];
        ryy += noise * noise;
        rxy += p1->profile[i] * noise;
        i++;
      }

      rnm = sqrt(rxx * ryy);
      corr[lag] = rxy / rnm;
      if (corr[lag] > corr[index]) index = lag;
    }
  }

  // profile1.length < profile2.length
  else if (p1->length < p2->length) {
    for (i = 0; i < p2->length; i++)
      ryy += p2->profile[i] * p2->profile[i];

    for (lag = 0; lag < length; lag++) {
      rxx = 0; rxy = 0; j = 0; i = 0;

      while(i < lag) {
        noise = p1->noise[rand() % MAX_PROFILE_LENGTH];
        rxx += noise * noise;
        rxy += noise * p2->profile[i];
        i++;
      }
      while(j < p1->length) {
        rxx += p1->profile[j] * p1->profile[j];
        rxy += p1->profile[j] * p2->profile[i];
        i++; j++;
      }
      while(i < p2->length) {
        noise = p1->noise[rand() % MAX_PROFILE_LENGTH];
        rxx += noise * noise;
        rxy += noise * p2->profile[i];
        i++;
      }

      rnm = sqrt(rxx * ryy);
      corr[lag] = rxy / rnm;
      if (corr[lag] > corr[index]) index = lag;
    }
  }

  // profile1.length == profile2.length
  else {
    i = 0; j = 0;
    while(i < p1->length) {
      rxx += p1->profile[i] * p1->profile[i];
      ryy += p2->profile[i] * p2->profile[i];
      rxy += p1->profile[i] * p2->profile[i];
      i++;
    }
    rnm = sqrt(rxx * ryy);
    corr[0] = rxy / rnm;
    index = 0;
  }

  nrxy = corr[index];
  free(corr);

  return(nrxy);
}
