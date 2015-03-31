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
inline void progress_bar(int k, int n) {
  int res = 20;     // progress bar resolution
  int window = 50;  // progress bar width
  int counter;
  int progress;
  double ratio;

  if (n < res) return;
  if (k != n && (k + 1) % (n/res) != 0) return;
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


void raytrace_interp_bin(
    const int NXPIX,
    const int NYPIX,
    const double DPIX,
    const double XDL_0,                       /* left detector x-offset */
    const double XDM_0,                       /* middle detector y-offset */
    const double XDR_0,                       /* right detector x-offset */
    const double YDL_0,                       /* left detector y-offset */
    const double YDM_0,                       /* middle detector y-offset */
    const double YDR_0,                       /* right detector y-offset */
    const double TAU_DL,                      /* left detector rotation */
    const double TAU_DM,                      /* middle detector rotation */
    const double TAU_DR,                      /* right detector rotation */
    const double slit_ratio,                  /* slit_height / slit_width */
    const unsigned long n,                    /* size of slit, waves */
    const unsigned int cn,                    /* interpolation size */
    const double* restrict cwl,               /* interpolation wavelengths */
    const double* restrict cxb,
    const double* restrict cxm,
    const double* restrict cxt,
    const double* restrict cyb,
    const double* restrict cym,
    const double* restrict cyt,
    const double* restrict cphi,
    const double* restrict wl,                /* wavelengths, size n */
    const double* restrict slitx,             /* slit location */
    const double* restrict slity,             /* slit location */
    unsigned long* restrict out) {
  const double YPIX_1_2 = (double)NYPIX/2.0;
  const double X_OFF_L = (double)NXPIX/6.0;
  const double X_OFF_M = X_OFF_L + (double)NYPIX;
  const double X_OFF_R = X_OFF_M + (double)NYPIX;
  const int NXPIX_1_3 = NXPIX / 3;
  const int NXPIX_2_3 = 2 * NXPIX_1_3;
  const double TAUDL = cradians(TAU_DL);
  const double TAUDM = cradians(TAU_DM);
  const double TAUDR = cradians(TAU_DR);
  long ix;                                  /* x coordinate [pix] */
  long iy;                                  /* y coordinate [pix] */
  unsigned long i;                            /* wavelength iterator */
  double xd, yd;
  double sx, sy, sw, sh;
  double xb, xm, xt, yb, ym, yt, phi;
  double xdp, ydp;

  // CUBIC SPLINES
  gsl_interp_accel* acc1 = gsl_interp_accel_alloc();
  gsl_interp_accel* acc2 = gsl_interp_accel_alloc();
  gsl_interp_accel* acc3 = gsl_interp_accel_alloc();
  gsl_interp_accel* acc4 = gsl_interp_accel_alloc();
  gsl_interp_accel* acc5 = gsl_interp_accel_alloc();
  gsl_interp_accel* acc6 = gsl_interp_accel_alloc();
  gsl_interp_accel* acc7 = gsl_interp_accel_alloc();
  gsl_spline* spline1 = gsl_spline_alloc(gsl_interp_akima, cn);
  gsl_spline* spline2 = gsl_spline_alloc(gsl_interp_akima, cn);
  gsl_spline* spline3 = gsl_spline_alloc(gsl_interp_akima, cn);
  gsl_spline* spline4 = gsl_spline_alloc(gsl_interp_akima, cn);
  gsl_spline* spline5 = gsl_spline_alloc(gsl_interp_akima, cn);
  gsl_spline* spline6 = gsl_spline_alloc(gsl_interp_akima, cn);
  gsl_spline* spline7 = gsl_spline_alloc(gsl_interp_akima, cn);
  gsl_spline_init(spline1, cwl, cxb, cn);
  gsl_spline_init(spline2, cwl, cxm, cn);
  gsl_spline_init(spline3, cwl, cxt, cn);
  gsl_spline_init(spline4, cwl, cyb, cn);
  gsl_spline_init(spline5, cwl, cym, cn);
  gsl_spline_init(spline6, cwl, cyt, cn);
  gsl_spline_init(spline7, cwl, cphi, cn);

  for (i=0; i<n; ++i) {
    phi = gsl_spline_eval(spline7, wl[i], acc7);
    sx = slitx[i]*cos(phi) - slity[i]*sin(phi);
    sy = slitx[i]*sin(phi) + slity[i]*cos(phi);
    xb = gsl_spline_eval(spline1, wl[i], acc1);
    xm = gsl_spline_eval(spline2, wl[i], acc2);
    xt = gsl_spline_eval(spline3, wl[i], acc3);
    yb = gsl_spline_eval(spline4, wl[i], acc4);
    ym = gsl_spline_eval(spline5, wl[i], acc5);
    yt = gsl_spline_eval(spline6, wl[i], acc6);
    sh = sqrt((xt-xb)*(xt-xb) + (yt-yb)*(yt-yb));
    sw = sh / slit_ratio;
    xd = xm + sx * sw;
    yd = ym + sy * sh;

    /* LEFT DETECTOR */
    xdp =  xd-XDL_0 + (xd-XDL_0)*cos(TAUDL) + (yd-YDL_0)*sin(TAUDL);
    ydp =  yd-YDL_0 - (xd-XDL_0)*sin(TAUDL) + (yd-YDL_0)*cos(TAUDL);
    xdp /= (DPIX*2.0);  // NOT QUITE SURE WHY WE NEED THIS
    ydp /= (DPIX*2.0);
    ix = (int)floor(xdp + X_OFF_L);
    if (ix >= 0 && ix < NXPIX_1_3) {
      iy = (int)floor(ydp + YPIX_1_2);
      if (iy >= 0 && iy < NYPIX) {
        out[ix+NXPIX*iy] += 1;
      }
    }
    /* MIDDLE DETECTOR */
    xdp = xd-XDM_0 + (xd-XDM_0)*cos(TAUDM) + (yd-YDM_0)*sin(TAUDM);
    ydp = yd-YDM_0 - (xd-XDM_0)*sin(TAUDM) + (yd-YDM_0)*cos(TAUDM);
    xdp /= (DPIX*2.0);
    ydp /= (DPIX*2.0);
    ix = (int)floor(xdp + X_OFF_M);
    if (ix >= NXPIX_1_3 && ix < NXPIX_2_3) {
      iy = (int)floor(ydp + YPIX_1_2);
      if (iy >= 0 && iy < NYPIX) {
        out[ix+NXPIX*iy] += 1;
      }
    }
    /* RIGHT DETECTOR */
    xdp = xd-XDR_0 + (xd-XDR_0)*cos(TAUDR) + (yd-YDR_0)*sin(TAUDR);
    ydp = yd-YDR_0 - (xd-XDR_0)*sin(TAUDR) + (yd-YDR_0)*cos(TAUDR);
    xdp /= (DPIX*2.0);
    ydp /= (DPIX*2.0);
    ix = (int)floor(xdp + X_OFF_R);
    if (ix >= NXPIX_2_3 && ix < NXPIX) {
      iy = (int)floor(ydp + YPIX_1_2);
      if (iy >= 0 && iy < NYPIX) {
        out[ix+NXPIX*iy] += 1;
      }
    }
    progress_bar(i, n);
  }
  gsl_spline_free(spline1);
  gsl_spline_free(spline2);
  gsl_spline_free(spline3);
  gsl_spline_free(spline4);
  gsl_spline_free(spline5);
  gsl_spline_free(spline6);
  gsl_spline_free(spline7);
  gsl_interp_accel_free(acc1);
  gsl_interp_accel_free(acc2);
  gsl_interp_accel_free(acc3);
  gsl_interp_accel_free(acc4);
  gsl_interp_accel_free(acc5);
  gsl_interp_accel_free(acc6);
  gsl_interp_accel_free(acc7);
  printf("\n");
}


// Physical model is based on implementation for CARMENES
// (Marvin, 2012; master's thesis)
// http://www.astro.physik.uni-goettingen.de/~cmarvin/masters_thesis/mthesis_chris.pdf
//
// which is based on Ballester & Rosa (1997)
// http://adsabs.harvard.edu/abs/1997A%26AS..126..563B
//
// which draws heavily from Astronomical Optics, 2nd Edition by Schroeder (1987).
void raytrace_solve_general(
    int BLAZE_FLAG,
    int RETURN_MODE,
    unsigned long n,      /* number of slit/wavelength elements */
    unsigned short m,
    int NXPIX,
    int NYPIX,
    double F_COL,
    double F_COL_2,
    double ALPHA_E,
    double BLAZE_E,
    double GAMMA_E,
    double SIGMA_E,
    double ALPHA_G,
    double SIGMA_G,
    double F_CAM,
    double F_CAM_1,
    double DPIX,
    const double XDL_0,                       /* left detector x-offset */
    const double XDM_0,                       /* middle detector y-offset */
    const double XDR_0,                       /* right detector x-offset */
    const double YDL_0,                       /* left detector y-offset */
    const double YDM_0,                       /* middle detector y-offset */
    const double YDR_0,                       /* right detector y-offset */
    const double TAU_DL,                      /* left detector rotation */
    const double TAU_DM,                      /* middle detector rotation */
    const double TAU_DR,                      /* right detector rotation */
    double* restrict xslit,          /* slit location */
    double* restrict yslit,          /* slit location */
    double* restrict lamb,            /* wavelengths */
    double* restrict returnx,
    double* restrict returny,
    double* restrict outwaves,                        /* ccd array */
    unsigned long* restrict outcounts) {              /* counts array */

  /* pre-calculate constants */
  ALPHA_E = cradians(ALPHA_E);
  GAMMA_E = cradians(GAMMA_E);
  ALPHA_G = cradians(ALPHA_G);
  const double ORDER = (double)m;
  const double F_COL2 = F_COL * F_COL;
  const double MU_E0 = ALPHA_E - M_PI - cradians(-3.55);
  const double NU_E0 = GAMMA_E;
  const double MU_E1 = ALPHA_E + M_PI + cradians(-3.55);
  const double NU_E1 = -GAMMA_E;
  const double NU_G0 = cradians(0.0);//ALPHA_G;
  const double NU_G1 = ALPHA_G + M_PI;
  const double NU_G2 = ALPHA_G + M_PI;
  const double COS_MU_E0 = cos(MU_E0);
  const double SIN_MU_E0 = sin(MU_E0);
  const double COS_MU_E1 = cos(MU_E1);
  const double SIN_MU_E1 = sin(MU_E1);
  const double COS_NU_E0 = cos(NU_E0);
  const double SIN_NU_E0 = sin(NU_E0);
  const double COS_NU_E1 = cos(NU_E1);
  const double SIN_NU_E1 = sin(NU_E1);
  const double COS_NU_G0 = cos(NU_G0);
  const double SIN_NU_G0 = sin(NU_G0);
  const double COS_NU_G1 = cos(NU_G1);
  const double SIN_NU_G1 = sin(NU_G1);
  const double COS_NU_G2 = cos(NU_G2);
  const double SIN_NU_G2 = sin(NU_G2);
  const double YPIX_1_2 = (double)NYPIX/2.0;
  const double X_OFF_L = (double)NXPIX/6.0;
  const double X_OFF_M = X_OFF_L + (double)NYPIX;
  const double X_OFF_R = X_OFF_M + (double)NYPIX;
  const int NXPIX_1_3 = NXPIX / 3;
  const int NXPIX_2_3 = 2 * NXPIX_1_3;
  const double TAUDL = cradians(TAU_DL);
  const double TAUDM = cradians(TAU_DM);
  const double TAUDR = cradians(TAU_DR);

  /* local variables */
  double veclen;
  double li;
  double wi;
  double n_g;
  double beta;
  double blaze_eff;
  long ix; /* x coordinate [pix] */
  long iy; /* y coordinate [pix] */
  double x, x0, xi, xd, xt, y, y0, yd, yi, yt, z, z0, zt, fcc;
  double xdpld, xdpmd, xdprd;
  unsigned long i;   /* iterator */

  // RANDOM NUMBER GENERATOR INITIALIZATION
  double u;
  gsl_rng* r = gsl_rng_alloc(gsl_rng_taus); /* global rng generator */
  long seed = time(NULL);
  gsl_rng_set(r, seed);                     /* seed the rng */

  for (i=0; i<n; ++i) {
    xi = xslit[i];       /* x coord */
    yi = yslit[i];       /* y coord */
    li = lamb[i];        /* wavelength */

    /* LAUNCH RAYS FROM SLIT (normalized vectors) */
    veclen = sqrt( xi*xi + yi*yi + F_COL2 );
    x0 = xi / veclen;
    y0 = yi / veclen;
    z0 = F_COL / veclen;

    /* BLAZE FUNCTION */
    //xb = xx;

    /* :::::::::::::::::::: GRATING :::::::::::::::::::::::::::::::::*/
    /* INTO PLANE RF */
    //xt = x0;
    //yt = y0 * COS_NU_G0 - z0 * SIN_NU_G0;
    //zt = y0 * SIN_NU_G0 + z0 * COS_NU_G0;
    /* PLANE RELATION */
    //x = xt / n_g;
    //y = yt / n_g;
    /* NORMALIZATION AFTER PLANE RELATION (REFRACTION) */
    //z = (zt / fabs(zt)) * sqrt(1.0 - x*x - y*y);

    /* INTO GRISM RF and GRISM RELATION */
    xt = x0;
    yt = (-li / SIGMA_G) + (COS_NU_G1 * y0) - (SIN_NU_G1 * z0);
    zt = (SIN_NU_G1 * y0) + (COS_NU_G1 * z0);
    /* NORMALIZATION AFTER GRATING RELATION */
    zt = (zt / fabs(zt)) * sqrt(1.0 - xt*xt - yt*yt);

    /* OUT OF GRISM RF */
    x = xt;
    y = (COS_NU_G2 * yt) - (SIN_NU_G2 * zt);
    z = (SIN_NU_G2 * yt) + (COS_NU_G2 * zt);
    /* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

    /* CAM-COLL */
    fcc =  F_COL_2 / F_CAM_1 ;
    veclen = sqrt( x*x + y*y + fcc*fcc );
    x = x /veclen;
    y = y / veclen;
    z = fcc / veclen;

    /* :::::::::::::::::::: ECHELLE :::::::::::::::::::::::::::::::::::::*/
    /* INTO ECHELLE RF */
    xt = (ORDER * li / SIGMA_E) - (COS_MU_E0 * x) + (SIN_MU_E0*SIN_NU_E0 * y) + (SIN_MU_E0*COS_NU_E0 * z);
    yt = -(COS_NU_E0 * y) + (SIN_NU_E0 * z);
    zt = -(SIN_MU_E0 * x) - (COS_MU_E0*SIN_NU_E0 * y) - (COS_MU_E0*COS_NU_E0 * z);
    /* NORMALIZATION AFTER ECHELLE RELATION */
    zt = (zt / fabs(zt)) * sqrt(1.0 - yt*yt - xt*xt);

    /* OUT OF ECHELLE RF */
    x = (COS_MU_E1 * xt) - (SIN_MU_E1 * zt);
    y = -(SIN_MU_E1*SIN_NU_E1) * xt + (COS_NU_E1 * yt - COS_MU_E1*SIN_NU_E1 * zt);
    z = (SIN_MU_E1*COS_NU_E1) * xt + (SIN_NU_E1 * yt + COS_MU_E1*COS_NU_E1 * zt);

    /* BLAZE FUNCTION*/
    if (BLAZE_FLAG == 1) {
      beta = M_PI * cos(ALPHA_E) * SIGMA_E * (x0 + x) / li;
      blaze_eff = (sin(beta)*sin(beta)) / (beta*beta);
    }
    /* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

    /* PROJECTION ONTO DETECTOR */
    xd = x/z * F_CAM;
    //yd = (y/z * F_CAM  + yd_0);
    yd = -(y/sqrt(x*x + z*z) * F_CAM);

    switch (RETURN_MODE) {
      /* BIN PIXELS */
      case 0:
        xd /= DPIX;
        yd /= DPIX;
        ix = (int)floor(xd + X_OFF_L);
        if (ix >= 0 && ix < NXPIX_1_3) {
          iy = (int)floor(yd + YPIX_1_2);
          if (iy >= 0 && iy < NYPIX) {
            if (BLAZE_FLAG == 1) {
              u = gsl_rng_uniform(r);
              if (u <= blaze_eff) {
                outwaves[ix+NXPIX*iy] += wi;
                outcounts[ix+NXPIX*iy] += 1;
              }
            }
            else {
              outwaves[ix+NXPIX*iy] += wi;
              outcounts[ix+NXPIX*iy] += 1;
            }
          }
        }
        break;

      /* RETURN EXACT LOCATIONS [mm] */
      case 1:
        returnx[i] = xd;
        returny[i] = yd;
        break;

      /* RETURN LOCATIONS 0-CENTERED IN [pix] */
      case 2:
        returnx[i] = xd/DPIX ;
        returny[i] = yd/DPIX ;
        break;

      /* RETURN LOCATIONS [pix] */
      case 3:
        returnx[i] = xd/DPIX + (double)NXPIX/2.0;
        returny[i] = yd/DPIX + (double)NYPIX/2.0;
        break;

      default:
        break;
    }
  }
  gsl_rng_free(r);
}
