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
  if (xx.value > yy.value) return 0;
  return 0;
}


/*
 * dclust
 *
 * @see include/annotate/dclust.h
 */
int dclust(double** dist, int n, int k, profile_struct_annotation* profiles)
{
  int *cl, *icl, *halo, *cds, *ordrho, *nneigh, *sites;
  double *allds, *rho, *delta, *sortrho, *rhotdelta, *bord_rho;
  int i, j, counter, pos, nclust, idx;
  double dc, maxd, maxdelta, maxrho, rhotdeltamin, rho_aver;
  rho_struct *strho;

  // Initialize structures
  cl = (int*) malloc(sizeof(int) * n);
  icl = (int*) malloc(sizeof(int) * n);
  halo = (int*) malloc(sizeof(int) * n); 
  cds = (int*) malloc(sizeof(int) * k);
  allds = (double*) malloc(sizeof(double) * (n * (n - 1) / 2));
  rho = (double*) malloc(sizeof(double) * n);
  delta = (double*) malloc(sizeof(double) * n);
  sortrho = (double*) malloc(sizeof(double) * n);
  ordrho = (int*) malloc(sizeof(int) * n);
  strho = (rho_struct*) malloc(sizeof(rho_struct) * n);
  nneigh = (int*) malloc(sizeof(int) * n);
  rhotdelta = (double*) malloc(sizeof(double) * n);
  sites =  (int*) malloc(sizeof(int) * n);
  bord_rho = (double*) malloc(sizeof(double) * n);
  for (i = 0; i < n; i++) rho[i] = 0.0f;

  // Convert dissimilarity matrix into a vector
  counter = 0;
  for (i = 0; i < n; i++) {
    for(j = i + 1; j < n; j++) {
      allds[counter] = dist[i][j];
      counter++;
    }
  }

  // STEP
  qsort(allds, n * (n - 1) / 2, sizeof(double), cmpd);
  pos = (int) round(((double) counter) * 2.0f / 100.0f) - 1;
  dc = allds[pos];

  // STEP
  for(i = 0; i < n; i++) {
    for(j = i + 1; j < n; j++) {
      double sq = -(dist[i][j] / (double)dc) * (dist[i][j] / (double)dc);
      double expsq = exp(sq);
      rho[i] += expsq;
      rho[j] += expsq;
    }
  }

  /* METHOD B
  for(i = 0; i < n; i++) {
    for(j = i + 1; j < n; j++) {
      if(dist[i][j] < dc) {
        rho[i]++;
        rho[j]++;
      }
    }
  } END METHOD B */

  // STEP
  for (i = 0; i < n; i++){
    strho[i].value = rho[i];
    sortrho[i] = rho[i];
    strho[i].index = i;
    ordrho[i] = i;
  }
  qsort(strho, n, sizeof(rho_struct), cmpsm);
  for (i = 0; i < n; i++) ordrho[i] = strho[i].index;

  // STEP
  maxd = allds[counter-1];
  delta[ordrho[0]] = -1;
  nneigh[ordrho[0]] = 0;
  for(i = 1; i < n; i++) {
    delta[ordrho[i]] = maxd;
    for(j = 0; j < i; j++) {
      if(dist[ordrho[i]][ordrho[j]] < delta[ordrho[i]]) {
        delta[ordrho[i]] = dist[ordrho[i]][ordrho[j]];
        nneigh[ordrho[i]] = ordrho[j];
      }
    }
  }

  delta[ordrho[0]] = gsl_stats_max (delta, 1, n); // delta[ordrho[0]] = *max_element(delta.begin(),delta.end());
  maxdelta = gsl_stats_max (delta, 1, n); // maxdelta = *max_element(delta.begin(),delta.end());
  maxrho = gsl_stats_max (rho,  1, n); // maxrho = *max_element(rho.begin(),rho.end());

  // STEP
  for(i = 0; i < n; i++) {
    rhotdelta[i] = sqrt(pow((delta[i]/maxdelta), 2) + pow((rho[i] / maxrho), 2));
    strho[i].value = rhotdelta[i];
    sites[i] = i;
    strho[i].index = i;
  }
  qsort(strho, n, sizeof(rho_struct), cmpsm); // sort(sites.begin(),sites.end(),mtd(rhotdelta));
  for (i = 0; i < n; i++) sites[i] = strho[i].index;

  /* CLUSTER NUMBER CALCULATION NOT NEEDED */
  // If k = 0, a way to calculate number of clusters
  if (k == 0) {
    double sum = 0;
    for (idx = 0; idx < n; idx++) sum += rhotdelta[idx]; //accumulate(rhotdelta.begin(), rhotdelta.end(), 0.0);
    double mean = sum / n; //rhotdelta.size();

    // Returns the result of accumulating init with the inner products of the pairs formed by the elements of two ranges starting at first1 and first2
    double sq_sum = 0; //inner_prodiuct(rhotdelta.begin(), rhotdelta.end(), rhotdelta.begin(), 0.0);
    for (idx = 0; idx < n; idx++) sq_sum += rhotdelta[idx] * rhotdelta[idx];
    
    double stdev = (double) sqrt(sq_sum / n - mean * mean); //(double) sqrt(sq_sum / rhotdelta.size() - mean * mean);
    double stdevby10 = stdev / 10.0;

    rhotdeltamin = gsl_stats_max(rhotdelta, 1, n) - 0.000001;;//*max_element(rhotdelta.begin(),rhotdelta.end()) - 0.000001;

    double rhotdeltamininit = rhotdeltamin;
    int kk = 1;
    int exit = 0; //bool exit = false;

    while(!exit) {
      int pkk = kk;
      rhotdeltamin -= stdevby10;
      while(rhotdelta[sites[kk]] > rhotdeltamin) kk++;

      if(kk-pkk > 5*pkk) { // Exit if numclusters increases 5 times the previous amount
        rhotdeltamin += stdev;
        exit = 1;//exit = true;
      }
      else {
        if((kk-pkk == 0) && (k > 1)) // Exit if no new clusters and already more than 1
          exit = 1;//exit = true;
      }

      if(rhotdeltamin < rhotdeltamininit - 2 * stdev) // Exit if we are getting close to the mean
        exit = 1;//exit = true;
    }
  }
  else 
    rhotdeltamin = rhotdelta[sites[k]];

  // STEP
  nclust = 0;
  for(i = 0; i < n; i++)
    cl[i] = -1;  // push_back : dds a new element at the end of the vector, after its current last element

  idx = 0;
  for(i = 0; i < n; i++) {
    if(rhotdelta[i] > rhotdeltamin) {
      cl[i] = nclust;
      icl[idx++] = i;
      nclust = nclust + 1; 
    }
  }

  // Assignation
  for(i = 0; i < n; i++)
    if(cl[ordrho[i]] == -1)
      cl[ordrho[i]] = cl[nneigh[ordrho[i]]];//cl[ordrho[i]] = cl[nneigh[ordrho[i]]];

  // Halo
  idx = 0;
  for(i = 0; i < n; i++)
    halo[idx++] = cl[i];

  rho_aver = 0.0;
  if (nclust > 1) {
    for(i = 0;i < nclust; i++)
      bord_rho[i] = 0.;

    for(i = 0; i < n-1; i++) {
      for(j = i+1; j < n; j++) {
        if ((cl[i] != cl[j]) && (dist[i][j] <= dc)) {
          rho_aver=(rho[i]+rho[j]) / 2.;
          if (rho_aver > bord_rho[cl[i]])
            bord_rho[cl[i]]=rho_aver;
          if (rho_aver > bord_rho[cl[j]])
            bord_rho[cl[j]]=rho_aver;
        }
      }
    }

    for(i = 0; i < n; i++) {
      if (rho[i] < bord_rho[cl[i]]) {
        halo[i]=-1;
      }
    }
  }

  /* PRINT 
  for(i = 0; i < nclust; i++) {
    int nc = 0;
    int nh = 0;
    for(j = 0; j < n; j++) {
      if (cl[j] == i)
        nc = nc + 1;
      if (halo[j]==i)
        nh = nh + 1;
    }
    fprintf(stderr, "CLUSTER: %d CENTER: %d ELEMENTS: %d CORE: %d HALO: %d \n", i, icl[i], nc, nh, nc - nh);
  }
  END PRINT */

  for(i = 0; i < n; i++)
    //fprintf(stderr, "%i %i %i\n", i, cl[i], halo[i]);
    profiles[i].cluster = cl[i] + 1;

  // Free structures and exit
  free(cl);
  free(icl);
  free(halo);
  free(cds);
  free(allds);
  free(rho);
  free(delta);
  free(sortrho);
  free(ordrho);
  free(strho);
  free(nneigh);
  free(rhotdelta);
  free(sites);
  free(bord_rho);

  return nclust;
}
