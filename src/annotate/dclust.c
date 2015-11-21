#include <annotate/dclust.h>

int cmpd(const void *x, const void *y)
{
  double xx = *(double*)x, yy = *(double*)y;
  if (xx < yy) return -1;
  if (xx > yy) return  1;
  return 0;
}

int cmpi(const void *x, const void *y)
{
  int xx = *(int*)x, yy = *(int*)y;
  if (xx < yy) return -1;
  if (xx > yy) return  1;
  return 0;
}

int cmpsm (const void *x, const void *y)
{
  rho_struct xx = *(rho_struct*) x;
  rho_struct yy = *(rho_struct*) y;

  if (xx.value < yy.value) return 1;
  if (xx.value > yy.value) return -1;
  return 0;
}


/*
 * dcoptimize
 *
 * @see include/annotate/dclust.h
 */
double dcoptimize(double** dist, int n, double* max)
{
  double dc, z, h, hmin, lower, upper, sigma;
  double* allds;
  int i, j, counter;
  double* potentials;

  // Initialize variables
  potentials = (double*) malloc(n * sizeof(double));
  allds = (double*) malloc(sizeof(double) * (n * (n - 1) / 2));
  hmin = DBL_MAX;
  dc = 0.0;

  // Convert dissimilarity matrix into a vector
  // Sort dissimilarity vector ascendently
  counter = 0;
  for (i = 0; i < n; i++) {
    for(j = i + 1; j < n; j++) {
      allds[counter] = dist[i][j];
      counter++;
    }
  }
  qsort(allds, n * (n - 1) / 2, sizeof(double), cmpd);
  *max = allds[counter - 1];

  // Calculate lower and upper boundaries for sigma (impact factor)
  //   lower : minimum distance > 0
  //   upper : distance at 10 percentile
  counter = 0;
  while(allds[counter] == 0) counter++;
  lower = allds[counter];
  upper = gsl_stats_quantile_from_sorted_data(allds, 1, (n * (n - 1) / 2), 0.10);

  // Calculate entropy for values of sigma between lower and upper in increments of 0.005
  for (sigma = lower; sigma <= upper; sigma += 0.005f) {

    // Calculate potential for each point
    for(i = 0; i < n; i++) {
      potentials[i] = 0.0f;
      for(j = 0; j < n; j++) {
        double sq = (dist[i][j] / sigma) * (dist[i][j] / sigma);
        double expsq = exp(-sq);
        if (j != i) potentials[i] += expsq;
      }
    }

    // Calculate Z (normalization factor)
    z = 0;
    for (i = 0; i < n; i++)
      z += potentials[i];

    // Calculate H (entropy)
    h = 0;
    for (i = 0; i < n; i++) {
      double a = potentials[i] / z;
      double b = log(a);
      h += a*b;
    }
    h = h * (-1);

    // Store dc and hmin
    if (h < hmin) {
      hmin = h;
      dc = sigma;
    }
  }

  free(potentials);
  free(allds);

  return((3/sqrt(2)) * dc);
}

/*
 * dclust
 *
 * @see include/annotate/dclust.h
 */
int dclust(double** dist, int n, profile_struct_annotation* profiles, double cutoff, int gaussian)
{
  int *cl, *halo;
  double *rho, *delta, *sortrho, *sortdelta, *bord_rho;
  int i, j, nclust;
  double dc, maxd, rhothreshold, deltathreshold;
  rho_struct *strho;

  // Initialize structures
  cl = (int*) malloc(sizeof(int) * n);
  halo = (int*) malloc(sizeof(int) * n); 
  rho = (double*) malloc(sizeof(double) * n);
  delta = (double*) malloc(sizeof(double) * n);
  sortrho = (double*) malloc(sizeof(double) * n);
  sortdelta = (double*) malloc(sizeof(double) * n);
  strho = (rho_struct*) malloc(sizeof(rho_struct) * n);
  for (i = 0; i < n; i++) rho[i] = 0.0;

  // Calculate optimal dc
  // Calculate maximum distance
  dc = dcoptimize(dist, n, &maxd);
  if (cutoff > 0)
    dc = cutoff;
  fprintf(stderr, "        Distance cutoff is %f\n", dc);

  // Calculate RHO[i] per point using gaussian kernel
  //   RHO[i] = sum {exp(-(dist(i,j)/dc)^2)}
  if (gaussian) {
    for(i = 0; i < n; i++) {
      for(j = i + 1; j < n; j++) {
        double sq = -(dist[i][j] / (double)dc) * (dist[i][j] / (double)dc);
        double expsq = exp(sq);
        rho[i] += expsq;
        rho[j] += expsq;
      }
    }
  }

  // Calculate RHO per point
  //   RHO[i] = number of points j that satisfy dist(i,j) < dc
  if (!gaussian) {
    for(i = 0; i < n; i++) {
      for(j = i + 1; j < n; j++) {
        if(dist[i][j] < dc) {
          rho[i]++;
          rho[j]++;
        }
      }
    }
  }

  // Calculate DELTA per point
  //   DELTA[i] = minimum {dist(i,j) if RHO[j] > RHO[i]}
  for (i = 0; i < n; i++) {
    delta[i] = maxd;
    for (j = 0; j < n; j++) {
      if ((j != i) && (rho[i] < rho[j]) && (delta[i] > dist[i][j]))
        delta[i] = dist[i][j];
    }
  }

  // Calculation of RHO and DELTA threshold
  //   deltathreshold = 3rd quartile
  //   rhothreshold = 1st quartile 
  for (i = 0; i < n; i++) {
    sortrho[i] = rho[i];
    sortdelta[i] = delta[i];
  }
  qsort(sortrho, n, sizeof(double), cmpd);
  qsort(sortdelta, n, sizeof(double), cmpd);
  rhothreshold = gsl_stats_quantile_from_sorted_data(sortrho, 1, n, 0.25);
  deltathreshold = gsl_stats_quantile_from_sorted_data(sortdelta, 1, n, 0.75);
  fprintf(stderr, "        Rho cutoff is %f\n", rhothreshold);
  fprintf(stderr, "        Delta cutoff is %f\n", deltathreshold);

  // Identification of cluster centers
  //   i is a center <=> RHO[i] > rhothreshold and DELTA[i] > deltathreshold
  nclust = 0;
  for (i = 0; i < n; i++) {
    profiles[i].center = 0;
    cl[i] = 0;
    if (rho[i] > rhothreshold && delta[i] > deltathreshold) {
      cl[i] = ++nclust;
      profiles[i].center = 1;
    }
  }

  // Order points according RHO values in ORDRHO descendently
  for (i = 0; i < n; i++) {
    strho[i].value = rho[i];
    strho[i].index = i;
  }
  qsort(strho, n, sizeof(rho_struct), cmpsm);

  // Assign clusters
  // Point is assigned to the same cluster as its nearest neighbor of higher density
  for (i = 0; i < n; i++) {
    int idxi = strho[i].index;
    if (cl[idxi] == 0) {
      double mindist = maxd;
      int minidx = n;
      for (j = 0; j < i; j++) {
        int idxj = strho[j].index;
        if ((rho[idxi] < rho[idxj]) && (dist[idxi][idxj] < mindist)) {
          cl[idxi] = cl[idxj];
          mindist = dist[idxi][idxj];
          minidx = idxj;
        }
        else if ((rho[idxi] < rho[idxj]) && (dist[idxi][idxj] == mindist) && (idxj < minidx)) {
          cl[idxi] = cl[idxj];
          mindist = dist[idxi][idxj];
          minidx = idxj;
        }
      }
    }
  }

  // Find border densities per cluster
  bord_rho = (double*) malloc((nclust + 1) * sizeof(double));
  for (i = 1; i <= nclust; i++) {
    bord_rho[i] = 0;
    for (j = 0; j < n; j++) {
      if (cl[j] == i) {
        int k;
        for (k = 0; k < n; k++) {
          if ((cl[k] != cl[j]) && (dist[j][k] <= dc) && (rho[j] > bord_rho[i]))
            bord_rho[i] = rho[j];
        }
      }
    }
  }

  // Generate Halo
  for (i = 0; i < n; i++) {
    if (rho[i] >= bord_rho[cl[i]])
      halo[i] = 0;
    else
      halo[i] = 1;
  }

  // Asign cluster number and halo to profiles
  for(i = 0; i < n; i++) {
    profiles[i].cluster = cl[i];
    profiles[i].halo = halo[i];
  }

  // Free structures and exit
  free(cl);
  free(halo);
  free(rho);
  free(sortrho);
  free(delta);
  free(sortdelta);
  free(strho);
  free(bord_rho);

  return nclust;
}


int dclustr_f(double** dist, int n, double dc, int gaussian, profile_struct_annotation** profiles, int ncluster)
{
  double *rho;
  int i, j, nclust, grhoidx;
  double maxrho;

  // Initialize structures
  rho = (double*) malloc(sizeof(double) * n);
  for (i = 0; i < n; i++) rho[i] = 0.0;

  // Calculate RHO[i] per point using gaussian kernel
  //   RHO[i] = sum {exp(-(dist(i,j)/dc)^2)}
  if (gaussian) {
    for(i = 0; i < n; i++) {
      for(j = i + 1; j < n; j++) {
        double sq = -(dist[i][j] / (double)dc) * (dist[i][j] / (double)dc);
        double expsq = exp(sq);
        rho[i] += expsq;
        rho[j] += expsq;
      }
    }
  }

  // Calculate RHO per point
  //   RHO[i] = number of points j that satisfy dist(i,j) < dc 
  if (!gaussian) {
    for(i = 0; i < n; i++) {
      for(j = i + 1; j < n; j++) {
        if(dist[i][j] < dc) {
          rho[i]++;
          rho[j]++;
        }
      }
    }
  }

  // Find point with greater RHO and assign cluster
  maxrho = 0;
  grhoidx = -1;
  for (i = 0; i < n; i++) {
    if (maxrho < rho[i]) {
      grhoidx = i;
      maxrho = rho[i];
    }
  }

  // Assign same cluster to profiles that are at a distance <= dc
  nclust = 0;
  for (i = 0; i < n; i++) {
    if (dist[grhoidx][i] <= dc) {
      nclust++;
      profiles[i]->cluster = ncluster;
    }
  }

  free(rho);
  return(nclust); 
}


/*
 * dclustr
 *
 * @see include/annotate/dclust.h
 */
int dclustr(double** dist, int n, profile_struct_annotation* profiles, double cf, int gaussian)
{
  double dc, maxd;//, cutoff;
  int i, j;//, counter;
  int ncluster, nvisited, stop;
  //double* allds;

  /*
  // Calculate distance cutoff
  if (cutoff < 0)
    dc = dcoptimize(dist, n, &maxd);
  else
    dc = cutoff;
  */

  /*
  allds = (double*) malloc((n * (n - 1) / 2) * sizeof(double));
  counter = 0;
  for (i = 0; i < (n - 1); i++) {
    for(j = i + 1; j < n; j++) {
      allds[counter] = dist[i][j];
      counter++;
    }
  }
  qsort(allds, n * (n - 1) / 2, sizeof(double), cmpd);
  cutoff = gsl_stats_quantile_from_sorted_data(allds, 1, (n * (n - 1)) / 2, 0.02);
  free(allds);
  */

  // Perform dclustr_f till an empty cluster is found
  nvisited = 0;
  stop = 0;
  ncluster = 1;
  while (!stop) {
    // Prepare
    profile_struct_annotation** prf = (profile_struct_annotation**) malloc((n - nvisited) * sizeof(profile_struct_annotation*));
    double** dm = (double**) malloc ((n - nvisited) * sizeof(double*));
    for (i = 0; i < (n - nvisited); i++) dm[i] = (double*) malloc((n - nvisited) * sizeof(double));
    int idxi = 0;
    for (i = 0; i < n; i++) {
      if (profiles[i].cluster < 0) {
        prf[idxi] = &profiles[i];
        int idxj = 0;
        for (j = 0; j < n; j++) {
          if (profiles[j].cluster < 0) {
            dm[idxi][idxj] = dist[i][j];
            idxj++;
          }
        }
        idxi++;
      }
    }

    // Cluster
    dc = dcoptimize(dm, (n - nvisited), &maxd);
    int nv = 0;
    if (dc <= cf)//cutoff)
      nv = dclustr_f(dm, (n - nvisited), dc, gaussian, prf, ncluster);
    else//if (nv == 1) //else
      stop++;

    // Finish
    for (i = 0; i < (n - nvisited); i++) free(dm[i]);
    nvisited += nv;
    ncluster++;
    free(dm);
    free(prf);
  }

  // Assign remaining profiles to clusters
  for (i = 0; i < n; i++) {
    if (profiles[i].cluster < 0) profiles[i].cluster = ncluster++;
  }

  // Free and return
  return (ncluster - 1);
}
