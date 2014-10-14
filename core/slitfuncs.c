#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

#if __STDC_VERSION__ < 199901L
#define restrict
#endif


// progress bar, inspired by Jinglin Zhao
// k:iteration
// n:total number
inline void progress_bar(int k, int n) {
  int res = 20;     // progress bar resolution
  int window = 50;  // progress bar width
  int counter;
  int progress;
  double ratio;

  if (n < res) return;
  if ( (k != n) && (k + 1) % (n/res) != 0) return;
  ratio = (double)(k + 1) / n;
  counter = ratio * window;
  printf("    [");
  for (progress=0; progress<window; ++progress) {
    if (progress < counter) {
      printf("=");
    }
    else if (progress == counter) {
      printf(">");
    }
    else {
      printf(" ");
    }
  }
  printf("\b] %6.2f %%\r", (ratio*100.0) );
  fflush(stdout);
}


inline double cradians(double angle) {
    return angle * M_PI / 180.0;
}

void slit_uniform_psf(
    const int n,
    const double seeing,
    const double mu_x,
    const double mu_y,
    const double tau_0,
    const double slit_width,
    const double slit_height,
    double* restrict x_out,
    double* restrict y_out) {
  int cond;
  int i;
  double x;
  double y;
  const double tau = cradians(tau_0);
  const double r2 = 0.25*seeing*seeing;
  const double sw = slit_width * 0.5;
  const double sh = slit_height * 0.5;
  long seed = time(NULL);
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_taus);   /* global rng generator */
  gsl_rng_set(rng, seed);                       /* seed the rng */

  for (i=0; i<n; ++i) {
    // CIRCLE + SLIT FILTER
    cond = 0;
    while (cond != 1) {
      x = gsl_rng_uniform(rng) * seeing - seeing*0.5;
      y = gsl_rng_uniform(rng) * seeing - seeing*0.5;
      if (x*x + y*y < r2) {
        x += mu_x;
        y += mu_y;
        if (x < sw && x > -sw && y < sh && y > -sh)
          cond = 1;
      }
    }
    // ROTATE
    x_out[i] = x * cos(tau) + y * sin(tau);
    y_out[i] = -x * sin(tau) + y * cos(tau);
    progress_bar(i, n);
  }
	printf("\n");
  gsl_rng_free(rng);
}


void slit_gaussian_psf(
    const int n,
    const double mu_x,
    const double mu_y,
    const double sig_x,
    const double sig_y,
    const double tau_0,
    const double slit_width,
    const double slit_height,
    double* restrict x_out,
    double* restrict y_out) {
  int cond;
  int i;
  double x;
  double y;
  const double tau = cradians(tau_0);
  const double sw = slit_width * 0.5;
  const double sh = slit_height * 0.5;
  long seed = time(NULL);
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_taus);   /* global rng generator */
  gsl_rng_set(rng, seed);                       /* seed the rng */

  for (i=0; i<n; ++i) {
    // SLIT FILTER
    cond = 0;
    while (cond == 0) {
      x = gsl_ran_gaussian_ziggurat(rng, sig_x) + mu_x;
      y = gsl_ran_gaussian_ziggurat(rng, sig_y) + mu_y;
      if (x < sw && x > -sw && y < sh && y > -sh)
        cond += 1;
    }
    // ROTATE
    x_out[i] = x * cos(tau) + y * sin(tau);
    y_out[i] = -x * sin(tau) + y * cos(tau);
    progress_bar(i, n);
  }
	printf("\n");
  gsl_rng_free(rng);
}