"""Simulator.py

Simulator class.
"""
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
import scipy.integrate
import scipy.ndimage.filters
from astropy.io import fits
import ctypes as ct
import logging
import os
import sys
import time
import cinterface as ci
from defaults import *
from noise import add_noise
import physics
import slit
import wavefuncs as wf

log = logging.getLogger(__name__)


class Simulator(object):
    """Object that initializes and encapsulates simulation parameters and data.

    Parameters
    ----------
    inst : object_like
        Instance of Instrument object that contains data parameters.

    Attributes
    ----------
    source_spectrum : tuple
        Length-2 tuple consisting of the wavelengths in nm and pdf (flux
        converted to counts).
    slitfunc : function
        Function to produce slit psf.
    modelfunc: function
        Raytracing function to input source_spectrum and slitfunc.
    outarr : array_like
        Simulated image.

    """

    def __init__(self, inst, plot=False):

        log.info("Initializing simulator...")
        for k,v in inst.__dict__.iteritems():
            setattr(self, k, v)

        # IMPORT SOURCE SPECTRUM
        self.source_spectrum = self.import_source_spectrum()
        self.initialize_spectrum()

        # SETUP SLIT ILLUMINATION
        self.slitfunc = self.init_slitfunc()

        # SETUP RAYTRACING MODEL
        self.modelfunc = self.init_raytrace()

        # INITIALIZE DATA
        self.outarr = np.empty(self.det_dims)
        self.outarr = np.require(self.outarr, requirements=ci.req_out,
            dtype=np.uint)
        self.nrays_tot = 0
        self.nrays_per_order = []
        self.mean_rays_per_pixel = 0
        self.med_rays_per_pixel = 0
        self.min_rays_per_pixel = 0
        self.max_rays_per_pixel = 0


    def run(self):
        """
        TODO nrays keyword to break up long simulations into chunks
        this can either be to avoid memory errors or a quasi-simulation
        of time (psf tracking, stability/tilts, etc.) by altering
        simulation parameters per run"""
        for i in xrange(self.nruns):
            log.info("Simulation run: %s", i+1)
            self.simulate()


    def add_noise(self):
        if self.noise:
            add_noise(self)


    def simulate(self, plot=False, **kwargs):
        """
        TODO Different run parameters can be passed through kwargs
        """

        waves, pdf = self.source_spectrum[0], self.source_spectrum[1]
        pdf_tot = scipy.integrate.simps(pdf, waves)
        t0 = time.time()
        d0 = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())

        # START SIMULATION
        log.info("Beginning simulation, %s", d0)
        for m in self.orders:
            # pre sample wavelengths
            mwaves, mpdf = wf.feed_spectrum(self, m, waves, pdf)
            if mwaves.size == 0 or mpdf.size == 0:
                log.warning("Order %s failed. Source spectrum does not have wavelengths in this range.", m)
                continue
            pdf_ratio = scipy.integrate.simps(mpdf, mwaves) / pdf_tot
            mnrays = int(pdf_ratio * self.nrays)
            # just in case we are given 0 rays to simulate
            while mnrays <= 0:
                mnrays = np.random.poisson(mnrays)
            waves_in = wf.sample_cdf(mwaves, mpdf, mnrays)
            self.nrays_per_order.append(mnrays)

            # pre sample slit
            slit_x, slit_y = self.slitfunc(mnrays)
            # normalize to unity
            slit_x /= self.slit_height
            slit_y /= self.slit_height

            # input through model
            assert slit_x.size == slit_y.size == waves_in.size
            self.modelfunc(m, waves_in, slit_x, slit_y)
        self.sim_time = time.time() - t0
        inds = np.where(self.outarr != 0)
        self.mean_rays_per_pixel = np.mean(self.outarr[inds])
        self.med_rays_per_pixel = np.median(self.outarr[inds])
        self.min_rays_per_pixel = np.min(self.outarr[inds])
        self.max_rays_per_pixel = np.max(self.outarr[inds])
        self.nrays_tot = np.sum(self.outarr[inds])


    # =========================[ initfuncs ]===================================


    def import_source_spectrum(self):
        if len(self.source) > 2:
            log.error("More than 2 input source files specified.")
            sys.exit(0)
        elif len(self.source) == 2:
            return self.import_two_files()
        elif len(self.source) == 1:
            if self.source[0].lower() in "p phoenix".split():
                return self.phx_model()
            elif self.source[0].lower() in "f flatfield".split():
                return self.flatfield()
            else:
                return self.import_one_file()
        elif len(self.source) == 0:
            return self.phx_model()


    def init_slitfunc(self):
        log.info("Initializing source PSF: %s", self.psf)
        if self.psf == "gaussian":
            return self.psf_gaussian

    def init_raytrace(self):
        log.info("Initializing raytracing model: %s", self.model)
        if self.model == "interp":
            return self.interp
        elif self.model == "solve":
            return self.solve

    # =========================[ model methods ]===============================

    def interp(self, m, waves, slit_x, slit_y):
        fn = os.path.join(codevparsednpy_path % (self.band, self.echang, m))
        log.info("Loading '%s'", fn)
        _m, wl, xb, xmid, xt, yb, ymid, yt, slitheight, phi = np.load(fn).T
        inds = np.argsort(wl)
        wl = wl[inds]
        xbot = xb[inds]
        xmid = xmid[inds]
        xtop = xt[inds]
        ybot = yb[inds]
        ymid = ymid[inds]
        ytop = yt[inds]
        phi = phi[inds]
        buff = 2.0 * np.median(np.diff(wl))  # interpolation limit buffer
        log.info("Extending interpolation boundaries by %s Ang", buff)
        # EXTEND BOUNDARIES TO AVOID GSL INTERPOLATION ERROR
        fxb = InterpolatedUnivariateSpline(wl, xbot)
        fxm = InterpolatedUnivariateSpline(wl, xmid)
        fxt = InterpolatedUnivariateSpline(wl, xtop)
        fyb = InterpolatedUnivariateSpline(wl, ybot)
        fym = InterpolatedUnivariateSpline(wl, ymid)
        fyt = InterpolatedUnivariateSpline(wl, ytop)
        fphi = InterpolatedUnivariateSpline(wl, phi)
        # APPEND NEW ENDPOINTS
        wmin = np.concatenate((wl, waves)).min() - buff
        wmax = np.concatenate((wl, waves)).max() + buff
        wl = np.insert(wl, 0, wmin)
        wl = np.append(wl, wmax)
        xbot = np.insert(xbot, 0, fxb(wmin))
        xbot = np.append(xbot, fxb(wmax))
        xmid = np.insert(xmid, 0, fxm(wmin))
        xmid = np.append(xmid, fxm(wmax))
        xtop = np.insert(xtop, 0, fxt(wmin))
        xtop = np.append(xtop, fxt(wmax))
        ybot = np.insert(ybot, 0, fyb(wmin))
        ybot = np.append(ybot, fyb(wmax))
        ymid = np.insert(ymid, 0, fym(wmin))
        ymid = np.append(ymid, fym(wmax))
        ytop = np.insert(ytop, 0, fyt(wmin))
        ytop = np.append(ytop, fyt(wmax))
        phi = np.insert(phi, 0, fphi(wmin))
        phi = np.append(phi, fphi(wmax))
        # SEND TO C FUNCTION
        nxpix = self.det_dims[1]
        nypix = self.det_dims[0]
        dpix = self.dpix
        slit_ratio = self.slit_ratio
        n = slit_x.size
        cn = wl.size
        xdl_0 = self.xdl_0
        xdm_0 = self.xdm_0
        xdr_0 = self.xdr_0
        ydl_0 = self.ydl_0
        ydm_0 = self.ydm_0
        ydr_0 = self.ydr_0
        tau_dl = self.tau_dl
        tau_dm = self.tau_dm
        tau_dr = self.tau_dr
        func = ci.raytrace.raytrace_interp_bin
        func.argtypes = [
            ct.c_int,               # nxpix
            ct.c_int,               # nypix
            ct.c_double,            # dpix
            ct.c_double,            # xdl_0
            ct.c_double,            # xlm_0
            ct.c_double,            # xdr_0
            ct.c_double,            # ydl_0
            ct.c_double,            # ydm_0
            ct.c_double,            # ydr_0
            ct.c_double,            # tau_dl
            ct.c_double,            # tau_dm
            ct.c_double,            # tau_dr
            ct.c_double,            # slit_ratio
            ct.c_ulong,             # nslit
            ct.c_uint,              # cn
            ci.array_1d_double,     # cwl
            ci.array_1d_double,     # cxb
            ci.array_1d_double,     # cxm
            ci.array_1d_double,     # cxt
            ci.array_1d_double,     # cyb
            ci.array_1d_double,     # cym
            ci.array_1d_double,     # cyt
            ci.array_1d_double,     # cphi
            ci.array_1d_double,     # waves
            ci.array_1d_double,     # slit_x
            ci.array_1d_double,     # slit_y
            ci.array_2d_uint]       # outarr
        func.restype = None
        log.info("Raytracing order %s...", m)
        func(nxpix, nypix, dpix, xdl_0, xdm_0, xdr_0,
            ydl_0, ydm_0, ydr_0, tau_dl, tau_dm, tau_dr, slit_ratio, n, cn, wl,
            xbot, xmid, xtop, ybot, ymid, ytop, phi, waves, slit_x, slit_y,
            self.outarr)


    def solve(self, m, waves, slit_x, slit_y):
        # SEND TO C FUNCTION
        blaze_flag = int(self.blaze)
        return_mode = 0
        nxpix = self.det_dims[1]
        nypix = self.det_dims[0]
        dpix = self.dpix
        n = slit_x.size
        f_col_1 = self.f_col_1
        f_col_2 = self.f_col_2
        alpha_ech = self.alpha_ech
        blaze_ech = self.blaze_ech
        gamma_ech = self.gamma_ech
        sigma_ech = 1.0 / self.sigma_ech_inv
        alpha_cd = self.alpha_cd
        sigma_cd = 1.0 / self.sigma_cd_inv
        f_cam = self.f_cam
        f_cam_1 = self.f_cam_1
        returnx = np.empty(1)
        returny = np.empty(1)
        returnwaves = np.zeros(self.det_dims)
        xdl_0 = self.xdl_0
        xdm_0 = self.xdm_0
        xdr_0 = self.xdr_0
        ydl_0 = self.ydl_0
        ydm_0 = self.ydm_0
        ydr_0 = self.ydr_0
        tau_dl = self.tau_dl
        tau_dm = self.tau_dm
        tau_dr = self.tau_dr
        func = ci.raytrace.raytrace_solve_general
        func.argtypes = [
            ct.c_int,               # blaze_flag
            ct.c_int,               # return_mode
            ct.c_ulong,             # n (slit and waves)
            ct.c_uint,              # m
            ct.c_int,               # nxpix
            ct.c_int,               # nypix
            ct.c_double,            # f_col_1
            ct.c_double,            # f_col_2
            ct.c_double,            # alpha_ech
            ct.c_double,            # blaze_ech
            ct.c_double,            # gamma_ech
            ct.c_double,            # sigma_ech
            ct.c_double,            # alpha_cd
            ct.c_double,            # sigma_cd
            ct.c_double,            # f_cam
            ct.c_double,            # f_cam_1
            ct.c_double,            # dpix
            ct.c_double,            # xdl_0
            ct.c_double,            # xlm_0
            ct.c_double,            # xdr_0
            ct.c_double,            # ydl_0
            ct.c_double,            # ydm_0
            ct.c_double,            # ydr_0
            ct.c_double,            # tau_dl
            ct.c_double,            # tau_dm
            ct.c_double,            # tau_dr
            ci.array_1d_double,     # slit_x
            ci.array_1d_double,     # slit_y
            ci.array_1d_double,     # waves
            ci.array_1d_double,     # returnx
            ci.array_1d_double,     # returny
            ci.array_2d_double,     # returnwaves
            ci.array_2d_uint]       # returncounts
        func.restype = None
        log.info("Raytracing order %s...", m)
        func(blaze_flag, return_mode, n, m, nxpix, nypix, f_col_1, f_col_2,
            alpha_ech, blaze_ech, gamma_ech, sigma_ech, alpha_cd, sigma_cd,
            f_cam, f_cam_1, dpix, xdl_0, xdm_0, xdr_0, ydl_0, ydm_0, ydr_0,
            tau_dl, tau_dm, tau_dr, slit_x, slit_y, waves, returnx, returny,
            returnwaves, self.outarr)


    # =========================[ slit methods ]================================


    def psf_gaussian(self, nrays):
        n = int(nrays)
        mux = self.mu_x_psf
        muy = self.mu_y_psf
        sigx = self.sig_x_psf
        sigy = self.sig_y_psf
        tau = self.tau_s0
        sw = self.slit_width
        sh = self.slit_height
        return slit.slit_gaussian_psf(n, mux, muy, sigx, sigy, tau, sw, sh)


    def psf_uniform(self, nrays):
        pass


    def psf_point(self, offset=None):
        pass


    # ======================[ spectrum mtehods ]===============================


    def initialize_spectrum(self):
        # GET WAVELENGTHS INTO RIGHT UNITS
        if self.factor:
            log.info("Multiplying input wavelengths by factor %s.", self.factor)
            wavelengths = self.source_spectrum[0] * self.factor
        # REDSHIFT SPECTRA
        wavelengths = physics.redshift(self.source_spectrum[0]*1.e-9, self.rv) * 1.e9
        flux = self.source_spectrum[1]
        # ADD TELLURIC LINES
        if self.telluric and self.FITS_SOURCE == "STAR":
            flux = wf.convolve_telluric_lines(self.telluric, wavelengths, flux)
        # CONVERT FLUX TO NUMBER DENSITY
        if self.FITS_SOURCE == "STAR":
            pdf = physics.energy2counts(wavelengths, flux)
        else:
            pdf = flux
        # TRUNCATE TO SPECTROGRAPH LIMIT WITH A SMALL BUFFER AS CUSHION
        self.source_spectrum = wf.truncate_spectrum(self, wavelengths, pdf)


    def import_one_file(self):
        desc = "Importing source spectrum: object simulation"
        log.info("Initializing source spectrum: %s", self.source)
        fn = self.source[0]
        self.FITS_TYPE = "STAR"
        self.FITS_SOURCE = "OBJ"
        self.FITS_CATG = "OBS"
        self.FITS_SRC = 'SPECTRUM'
        self.FITS_INFILE1 = fn
        self.FITS_INFILE2 = None
        if os.path.isfile(fn):
            try:
                wavelengths, flux = fits.getdata(fn)
            except ValueError:
                log.exception("Input file not correct dimensions.", exc_info=True)
                sys.exit(0)
            except IOError:
                try:
                    wavelengths, flux = np.loadtxt(fn, unpack=True)
                except ValueError:
                    log.exception("Please specify an input source of the correct form.", exc_info=True)
                    sys.exit(0)
            return wavelengths, flux
        else:
            log.exception("Please specify a real input filename.", exc_info=True)
            sys.exit(0)

    def import_two_files(self):
        desc = "Importing source spectrum: object simulation"
        log.info("Initializing source spectrum: ")
        for fn in self.source:
            log.info("%s", fn)
        waves_fn = self.source[0]
        flux_fn = self.source[1]
        self.FITS_TYPE = "STAR"
        self.FITS_SOURCE = "OBJ"
        self.FITS_CATG = "OBS"
        self.FITS_SRC = 'SPECTRUM'
        self.FITS_INFILE1 = waves_fn
        self.FITS_INFILE2 = flux_fn
        if os.path.isfile(waves_fn):
            try:
                wavelengths = fits.getdata(waves_fn)
            except IOError:
                try:
                    wavelengths = np.loadtxt(waves_fn)
                    assert wavelengths.ndim == 1
                except AssertionError:
                    log.exception("Please specify an input source of the correct form.", exc_info=True)
                    sys.exit(0)
            try:
                flux = fits.getdata(flux_fn)
            except IOError:
                try:
                    flux = np.loadtxt(flux_fn)
                    assert flux.ndim == 1
                except AssertionError:
                    log.exception("Please specify an input source of the correct form.", exc_info=True)
                    sys.exit(0)
            return wavelengths, flux
        else:
            log.exception("Please specify a real input filename.", exc_info=True)
            sys.exit(0)


    def phx_model(self, plot=False):
        desc = "Importing source spectrum: PHOENIX synspec, T_eff=3000K [M/H]=0.0 log(g)=5.0"
        log.info(desc)
        self.FITS_TYPE = "STAR"
        self.FITS_SOURCE = "OBJ"
        self.FITS_CATG = "OBS"
        self.FITS_SRC = 'SPECTRUM'
        self.FITS_INFILE1 = os.path.basename(phx_waves_fn)
        self.FITS_INFILE2 = os.path.basename(phx_flux_fn)
        wavelengths = fits.open(phx_waves_path)[0].data * 0.1 # Ang to nm
        flux = fits.open(phx_flux_path)[0].data
        if plot:
            import matplotlib.pyplot as plt
            plt.plot(wavelengths, flux)
            plt.show()
        return wavelengths, flux


    def flatfield(self, plot=False):
        desc = "Importing source spectrum: flatfield lamp"
        self.FITS_TYPE = "FLAT"
        self.FITS_SOURCE = "HAL"
        self.FITS_CATG = "CAL"
        self.FITS_SRC = 'FLATFIELD'
        self.FITS_INFILE1 = None
        self.FITS_INFILE2 = None
        log.info(desc)
        wmin, wmax = self.wmin, self.wmax
        dw = self.dlamb
        n = np.int(np.ceil( (wmax-wmin) / dw ))
        wavelengths = np.linspace(wmin, wmax, num=n)
        flux = np.ones(wavelengths.size)
        if plot:
            import matplotlib.pyplot as plt
            plt.plot(wavelengths, flux)
            plt.show()
        return wavelengths, flux

    def wavemap(self):
        desc = "2D wavelength mapping"
        arm.fiber_description = 'wavelength mapping'
        arm.wavelengths = wf.calculate_wavelengths(arm, mode='CCD', nwaves=arm.nw)
        arm.intensities = arm.wavelengths
        arm.wavemap = True
        arm.fib_obstype[i] = 'WAVE'
        arm.fib_src[i] = 'WAVE'
        return wavelengths, wavelengths

    def wavetrace():
        desc = "1D wavelength tracing"

    def solve(inst, settings):
        log.info("Solving.")
        sys.exit(0)

    def spreadout(self, kernel=None):
        if not kernel:
            kernel = np.array([[0, 0, 1, 0, 0],
                               [0, 2, 2, 2, 0],
                               [1, 2, 5, 2, 1],
                               [0, 2, 2, 2, 0],
                               [0, 0, 1, 0, 0]], dtype='int16')

        self.outarr = scipy.ndimage.filters.convolve(self.outarr, kernel)
