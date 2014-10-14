//  Author: Jinglin Zhao, 2014
//  Institute: University of Goettingen
//  @CM restrict keyword added
//  @CM progress bar format changed to match other code
//dist array problem solved
//include RSS

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#if __STDC_VERSION__ < 199901L
#define restrict
#endif
/*
 *  Natural Splines
 *  ---------------
 *  Here the end-conditions are determined by setting the second
 *  derivative of the spline at the end-points to equal to zero.
 *
 *  There are n-2 unknowns (y[i]'' at x[2], ..., x[n-1]) and n-2
 *  equations to determine them.  Either Choleski or Gaussian
 *  elimination could be used.
 */

/*-------------------------------------------------------------------------*/

static void natural_spline(
  int     n,
  double * restrict x,
  double * restrict y,
  double * restrict b,
  double * restrict c,
  double * restrict d
){
    int nm1, i;
    double t;

    x--; y--; b--; c--; d--;

    if(n < 2) {
  errno = EDOM;
  return;
    }

    if(n < 3) {
  t = (y[2] - y[1]);
  b[1] = t / (x[2] - x[1]);
  b[2] = b[1];
  c[1] = c[2] = d[1] = d[2] = 0.0;
  return;
    }

    nm1 = n - 1;

    /* Set up the tridiagonal system */
    /* b = diagonal, d = offdiagonal, c = right hand side */

    d[1] = x[2] - x[1];
    c[2] = (y[2] - y[1]) / d[1];
    for( i = 2 ; i < n ; i++) {
  d[i]  = x[i+1] - x[i];
  b[i]  = 2.0 * (d[i-1] + d[i]);
  c[i+1]  = (y[i+1] - y[i]) / d[i];
  c[i]  = c[i+1] - c[i];
    }

    /* Gaussian elimination */

    for(i = 3; i < n; i++) {
  t = d[i-1] / b[i-1];
  b[i] = b[i] - t * d[i-1];
  c[i] = c[i] - t * c[i-1];
    }

    /* Backward substitution */

    c[nm1] = c[nm1] / b[nm1];
    for(i = n-2; i > 1; i--)
  c[i] = (c[i] - d[i] * c[i+1]) / b[i];

    /* End conditions */

    c[1] = c[n] = 0.0;

    /* Get cubic coefficients */

    b[1] = (y[2] - y[1]) / d[1] - d[i] * c[2];
    c[1] = 0.0;
    d[1] = c[2] / d[1];
    b[n] = (y[n] - y[nm1]) / d[nm1] + d[nm1] * c[nm1];
    for(i = 2; i < n; i++) {
  b[i] = (y[i+1] - y[i]) / d[i] - d[i] * (c[i+1] + 2.0 * c[i]);
  d[i] = (c[i+1] - c[i]) / d[i];
  c[i] = 3.0 * c[i];
    }
    c[n] = 0.0;
    d[n] = 0.0;

    return;
}

/*-------------------------------------------------------------------------*/

void spline_eval(
  int    * restrict nu, //evaluation: number of points to be evaluated
  double * restrict u,  //evaluation: input array
  double * restrict v,  //evaluation: output array
  int    * restrict n,  //interpolation: number of interporating points
  double * restrict x,  //interpolation: input array
  double * restrict y,  //interpolation: output array
  double * restrict b,  //interpolation coefficient
  double * restrict c,  //interpolation coefficient
  double * restrict d   //interpolation coefficient
){

/* Evaluate  v[l] := spline(u[l], ...),       l = 1,..,nu, i.e. 0:(nu-1)
 * Nodes x[i], coef (y[i]; b[i],c[i],d[i]); i = 1,..,n , i.e. 0:(*n-1)
 */
    const int n_1 = *n - 1;
    int i, j, k, l;
    double ul, dx;

    for(l = 0; l < *nu; l++)
      v[l] = u[l];

    i = 0;
    for(l = 0; l < *nu; l++) {
       ul = v[l];
       if(ul < x[i] || (i < n_1 && x[i+1] < ul)) {
           /* reset i  such that  x[i] <= ul <= x[i+1] : */
           i = 0;
           j = *n;
           do {
              k = (i+j)/2;
              if(ul < x[k]) j = k;
              else i = k;
           }
           while(j > i+1);
       }
       dx = ul - x[i];
       v[l] = y[i] + dx*(b[i] + dx*(c[i] + dx*d[i]));
    }
}

/*-------------------------------------------------------------------------*/

void random_wave(double * restrict wav_org, double * restrict dist, int N, int flag, double * restrict y, int NS){

/* wav_org  - array of input wavelegnth
   dist     - array of photon number distribution
   N        - size of wav_org and dist
   flag     - interpolation degree : 0 (for line list) or 3 (for continious spectra)
   y        - output: randomly sampled wavelength
   NS       - size of y (number of sampling points)
*/

    clock_t start = clock();
    srand48((long)time(NULL));
    int n, N_plus = N + 1;          //number of extended data

    if (flag == 0){               // no interpolation, step-like CDF, emission-line spectra
  double *CDF  = malloc(sizeof(double) * N_plus);
  CDF[0] = 0;
  for (n = 0; n < N; n++)
      CDF[n+1] = CDF[n] + dist[n];
  double total = CDF[N];
  for (n = 0; n <= N; n++)
      CDF[n] = CDF[n] / total;        //normalization

  unsigned long long ns;
  double random;
  int i, j, k;
  //double *num_den = malloc(sizeof(double) * N);
  for (ns = 0; ns < NS; ns ++) {
      random = drand48();
      i = 0;
      j = N;
      do{                 // search algorithm (bi-section) to find position of random in CDF
    k = (i+j)/2;
    if (random < CDF[k])
        j = k;
    else i = k;
      } while(j > i+1);             //until random number falls into (CDF[i], CDF[i+1])
      y[ns] = wav_org[i];
      //num_den[i] ++;
  }

  /*---data analysis---*/
/*
  FILE *num_den_file;
      num_den_file = fopen("./num_den_cum.dat", "w+");
  FILE *normalized_PDF_file;
      normalized_PDF_file = fopen("./normalized_PDF.dat", "w+");
  printf("Writing statistics to file\n");
  for (n = 0; n < N; n++) {
      fprintf(num_den_file, "%lf\t%le\n", wav_org[n], (double) num_den[n] / NS);  //normalization
      fprintf(normalized_PDF_file, "%lf\t%le\n", wav_org[n], dist[n] / total);
  }
  fclose(num_den_file);
  fclose(normalized_PDF_file);
*/
    }

  if (flag == 3){             //continious spectra, data interpolation

      /*---CDF; Data extension---*/

      double *wav = malloc(sizeof(double) * N_plus);
      double *CDF = malloc(sizeof(double) * N_plus);

      CDF[0] = 0;
      CDF[1] = dist[0];
      wav[0] = 1.5 * wav_org[0] - 0.5 * wav_org[1];
      for (n = 1; n < N; n++) {
    CDF[n+1] = CDF[n] + dist[n];
    wav[n] = (wav_org[n-1] + wav_org[n]) / 2;
      }
      wav[N] = 1.5 * wav_org[N-1] - 0.5 * wav_org[N-2];
      double total = CDF[N];

      unsigned long long ns;
      double *b, *c, *d, *x;
      b = malloc(sizeof(double) * N_plus);
      c = malloc(sizeof(double) * N_plus);
      d = malloc(sizeof(double) * N_plus);
      natural_spline(N_plus, CDF, wav, b, c, d);    //wav as a function of CDF
      int j;

      if (NS <= pow(10, 7)){
  x = malloc(sizeof(double) * NS);
  for (ns = 0; ns < NS; ns ++)
      x[ns] = total * drand48();
  spline_eval(&NS, x, y, &N_plus, CDF, wav, b, c, d); //returns y to be the sampling wavelength
  free(x);
      }

      else {              //(NS > pow(10, 7))
    int new_NS = pow(10, 7);
    int progress, times, TIMES = (int)(NS / new_NS), tail = NS - new_NS * TIMES;
    double ratio;
    x = malloc(sizeof(double) * new_NS);
    //printf("-----------------\nsampling progress\n-----------------\n");
    for (times = 0; times < TIMES; times ++){
      for (ns = 0; ns < new_NS; ns ++)
        x[ns] = total * drand48();

      spline_eval(&new_NS, x, y, &N_plus, CDF, wav, b, c, d);
      y += new_NS;

      // Show the load bar.
      ratio = (double)(times + 1) / TIMES;
      printf("    %3d%% [", (int)(ratio*100) );
      for (progress = 0; progress < (int)(ratio * 50) ; progress++)
        printf("=");
      for (progress = (int)(ratio * 50); progress < 50; progress++)
        printf(" ");
      printf("]\r");
      fflush(stdout);
    }
    free(x);
    double *x2      = malloc(sizeof(double) * tail);
    for (ns = 0; ns < tail; ns ++)
      x2[ns] = total * drand48();
    spline_eval(&tail, x2, y, &N_plus, CDF, wav, b, c, d);
    printf("\n");
    free(x2);
    y -= new_NS * TIMES;
    }

/*
    FILE *wavelegnth_sam;
  wavelegnth_sam = fopen("./wavelegnth_sam.dat", "w+");
    for (ns = 0; ns < NS; ns++)
  fprintf(wavelegnth_sam, "%lf\n", y[ns]);
    fclose(wavelegnth_sam);
*/

    /*---deallocate array memory---*/

    free(b);
    free(c);
    free(d);
    free(wav);
    free(CDF);

  }

    clock_t ends = clock();
    //printf("      Running time = %.2f s for %.2le samples\n",(double) (ends - start) / CLOCKS_PER_SEC, (double)NS);

}

/*-------------------------------------------------------------------------*/

int main(int argc, char *argv[]){
  /* calling
    gcc sampling_cum_clean.c && ./a.out 0

    SYNTAX: ./a.out  [intpol_deg]
    ./a.out 0
    intpol_deg : 0 or 3 interpolation degree
  */

    long long NS = 1*pow(10, 8);      //number of sampling points
    int intpol_deg;
    double *y = malloc(sizeof(double) * NS);
    int n, N;

    /*---interpolation mode---*/

    FILE *file;
    if (argc == 2) {
      intpol_deg = (int) atoi(argv[1]);
      if (intpol_deg == 3)        //spectrum interpolation
  file = fopen("./wav_int_sim.dat", "r");
      else if (intpol_deg == 0)
  file = fopen("./ThAr.dat", "r");
      else
  printf("Please input the correct interpolation degree\n");
    }
    else if (argc == 3){
      file = fopen(argv[1], "r");
      intpol_deg = (int) atoi(argv[2]);
    }
    else printf("Please follow the correct input format: <filename>(optional) <interpolation degree>\n");

    /*---number of data points---*/

    if (intpol_deg == 3)        //spectrum interpolation
      N = 500001;
    else if (intpol_deg == 0)       //no interpolation
      N = 619;
    else printf("N is not designated\n");

    /*---read in data---*/

    double *wav_org = malloc(sizeof(double) * N);
    double *dist_org  = malloc(sizeof(double) * N);
    for(n = 0; n < N; n++)
  fscanf(file, "%lf%lf", &wav_org[n], &dist_org[n]);
    random_wave(wav_org, dist_org, N, intpol_deg, y, NS);

    return 0;
}
