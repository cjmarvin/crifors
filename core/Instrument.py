"""
Instrument.py

The Instrument object is a data container of the instrument parameters.
"""

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
import ConfigParser
import logging
import sys
from defaults import *
import cinterface as ci

# config file helpers
_bands = {
    "Y" : "yband",
    "J" : "jband",
    "H" : "hband",
    "K" : "kband",
    "L" : "lband",
    "M" : "mband"
    }
_sections = [
    "crossdisp",
    "detector",
    "echelle",
    "instrument",
    "model",
    #"noise",
    "settings",
    "slit",
    "slitpsf",
    ]
_cl = {
    "bglight" : "--bglight",
    "blaze" : "--blaze",
    "config" : "--config",
    "dlamb" : "--dlamb",
    "echang" : "--echang",
    "factor" : "--factor",
    "model" : "--model",
    "noise" : "--noise",
    "nrays" : "--nrays",
    "nruns" : "--nruns",
    "outfn" : "--outfn",
    "plotslit" : "--plot-slit",
    "psf" : "--psf",
    "rv" : "--rv",
    "seeing" : "--seeing",
    "slit_width" : "--slit-width",
    "source" : "SOURCE",
    "telluric" : "--telluric",
    "spread" : "--spread",
    }
log = logging.getLogger(__name__)


def merge(dict_1, dict_2):
    """Merge two dictionaries.
    Values that evaluate to true take priority over false values.
    `dict_1` takes priority over `dict_2`.
    """
    return dict(
        (str(key), dict_1.get(key) or dict_2.get(key))
        for key in set(dict_2) | set(dict_1))


class Instrument(object):
    """
    Container of instrument parameters.

    The main intention of this is a global object to pass around to modules.

    """
    def __init__(self, args, cfg=None):

        # FIRST, LOAD DEFAULT CONFIG FILES, ACCORDING TO BAND
        if cfg is None:
            cfg = dinst
        self.band = args["BAND"]
        log.info("Initializing %s band.", args["BAND"])
        log.info('Loading default instrument parameters.')
        log.debug("  '%s'", cfg)
        with open(cfg) as fp:
            config = ConfigParser.SafeConfigParser()
            config.readfp(fp)

            # SET SPECIFIC ATTRIBUTES FOR SPECTRAL BAND
            for item in config.items(_bands[self.band]):
                setattr(self, item[0], eval(item[1]))

            # SET GENERAL ATTRIBUTES
            for section in _sections:
                for item in config.items(section):
                    try:
                        setattr(self, item[0], eval(item[1])) # eval if float
                    except NameError:
                        setattr(self, item[0], item[1])
        log.debug("Initializing default instrument parameters:")
        for k,v in sorted(self.__dict__.iteritems()):
            log.debug("   %s.%s = %s", self.__class__.__name__, k, v)

        # TODO OVERRIDE DEFAULTS WITH COMMAND LINE CONFIG FILE
        # TODO need error handling of attributes
        if args["--config"]:
            log.info("Overriding defaults with input configuration file.")
            log.info("  %s", args["--config"])

            with open(args["--config"]) as fp:
                config = ConfigParser.SafeConfigParser()
                config.readfp(fp)

                # SET SPECIFIC ATTRIBUTES FOR SPECTRAL BAND
                for item in config.items(_bands[self.band]):
                    setattr(self, item[0], eval(item[1]))

                # SET GENERAL ATTRIBUTES ONLY IF ALREADY PRESENT
                for section in _sections:
                    for item in config.items(section):
                        if item[0] in self.__dict__.keys():
                            try:
                                setattr(self, item[0], eval(item[1])) # eval if float
                            except NameError:
                                setattr(self, item[0], item[1])
            log.debug("Initializing default instrument parameters:")
            for k,v in sorted(self.__dict__.iteritems()):
                log.info("   %s.%s = %s", self.__class__.__name__, k, v)


        # TODO OVERRIDE WITH COMMAND LINE ARGUMENTS
        # TODO need error handling of attributes
        log.info("Initializing command line arguments.")
        for k,v in sorted(_cl.iteritems()):
            log.debug("  %s.%s = %s", self.__class__.__name__, k, args[v])
            try:
                setattr(self, k, eval(args[v]))
            except (TypeError, NameError, AttributeError) as e:
                setattr(self, k, args[v])

        # ADDITIONAL PARAMETERS AND CONVENIENCE CALCULATIONS
        if self.telluric:
            self.telluric = tell_species
        self.slit_ratio = self.slit_height / self.slit_width
        self.sig_x_psf = self.seeing / (2. * np.sqrt(2. * np.e))
        self.sig_y_psf = self.seeing / (2. * np.sqrt(2. * np.e))
        self.det_dims = (self.nypix, self.nxpix*self.ndet)
        # set detector points wrt origin
        self.xdl_0 = self.xdl - (self.nxpix*self.dpix)
        self.xdm_0 = self.xdm + 0.0 # operation ensures copy (not sure but jic)
        self.xdr_0 = self.xdr + (self.nxpix*self.dpix)
        self.ydl_0 = self.ydl + 0.0
        self.ydm_0 = self.ydm + 0.0
        self.ydr_0 = self.ydr + 0.0
        self.xdl_le = self.xdl_0 - self.nxpix*self.dpix*0.5
        self.xdl_re = self.xdl_0 + self.nxpix*self.dpix*0.5
        self.xdm_le = self.xdm_0 - self.nxpix*self.dpix*0.5
        self.xdm_re = self.xdm_0 + self.nxpix*self.dpix*0.5
        self.xdr_le = self.xdr_0 - self.nxpix*self.dpix*0.5
        self.xdr_re = self.xdr_0 + self.nxpix*self.dpix*0.5
        self.xdlm = -(self.nxpix*self.dpix+self.xdl+self.xdm)*0.5
        self.xdmr = (self.nxpix*self.dpix+self.xdm+self.xdr)*0.5

        # RUN INITIALIZATION METHODS/PROCEDURES
        self.set_spectral_orders()
        self.set_ccd_limits()

        log.info("Instrument initialization complete.")
        return None

# =============================================================================
# :::::::::::::::::::::::: INITIALIZATION METHODS :::::::::::::::::::::::::::::
# =============================================================================


    def set_spectral_orders(self):
        self.orders = np.arange(self.m_range_min, self.m_range_max+1, 1)


    def set_ccd_limits(self):
        """Find the wavelength limits of the detector."""
        log.info("Model = %s", self.model)
        if self.model == "interp":
            self.find_ccd_limits_interp()

        elif self.model == "solve":
            self.find_ccd_limits_falloff()
            #self.find_ccd_limits_interp()


    def find_ccd_limits_interp(self, write=False):
        # TODO SHOULD THIS BE A WAVEFUNC METHOD ?
        n = self.orders.size
        det_wl_lim = np.empty([self.ndet*2, n])

        # detector x positions in mm
        # NOTE does not take detector tilt into account
        x_det = np.array([self.xdl_le, self.xdl_re, self.xdm_le,
            self.xdm_re, self.xdr_le, self.xdr_re])
        log.debug("Detector edges = %s mm", x_det)
        for i in xrange(n):
            m = self.orders[i]
            fn = os.path.splitext(codevparsed_path % (self.band, self.echang, m))[0] + ".npy"
            log.debug("Loading '%s'...", fn)
            _m, wl, xbot, xmid, xtop, yb, ymid, yt, slitheight, phi = np.load(fn).T
            inds = np.argsort(wl)
            wl = wl[inds]
            xbot = xbot[inds]
            xmid = xmid[inds]
            xtop = xtop[inds]
            fxbot = InterpolatedUnivariateSpline(xbot, wl)
            fxmid = InterpolatedUnivariateSpline(xmid, wl)
            fxtop = InterpolatedUnivariateSpline(xtop, wl)
            wlpbot = fxbot(x_det)
            wlpmid = fxmid(x_det)
            wlptop = fxtop(x_det)
            wlp = wlpmid
            det_wl_lim[:, i] = wlp
            log.info("Wavelength limits for order %s found.", m)
        self.det_wl_lim = det_wl_lim
        self.wmin = self.det_wl_lim.min()
        self.wmax = self.det_wl_lim.max()
        if write:
            # WRITE TO ETC OUTPUT FILE
            pass
        pass


    def find_ccd_limits_falloff(self, write=False):

        import matplotlib.pyplot as plt
        blaze_flag = 0
        return_mode = 1   # return in mm

        ndims = (21, self.orders.size)
        solvex = np.empty(ndims)
        solvey = np.empty(ndims)
        interpx = np.empty(ndims)
        interpy = np.empty(ndims)
        for i,m in enumerate(self.orders):
            print i
            fn = os.path.splitext(codevparsed_path % (self.band, self.echang, m))[0] + ".npy"
            _m, wl, xbot, xmid, xtop, yb, ymid, yt, slitheight, phi = np.load(fn).T
            n = wl.size

            print 'assigning parameters'
            slit_x = np.zeros(n)
            slit_y = np.zeros(n)
            nxpix = self.det_dims[1]
            nypix = self.det_dims[0]
            dpix = self.dpix
            f_col_1 = self.f_col_1
            f_col_2 = self.f_col_2
            alpha_ech = self.alpha_ech
            blaze_ech = self.blaze_ech
            gamma_ech = self.gamma_ech
            sigma_ech = 1.e7 / self.sigma_ech_inv
            alpha_cd = self.alpha_cd
            sigma_cd = 1.e7 / self.sigma_cd_inv
            f_cam = self.f_cam
            f_cam_1 = self.f_cam_1
            xdl_0 = self.xdl_0
            xdm_0 = self.xdm_0
            xdr_0 = self.xdr_0
            ydl_0 = self.ydl_0
            ydm_0 = self.ydm_0
            ydr_0 = self.ydr_0
            tau_dl = self.tau_dl
            tau_dm = self.tau_dm
            tau_dr = self.tau_dr

            print 'creating return arrays'
            returnx = np.empty(n)
            returny = np.empty(n)
            returnwaves = np.empty((1, 1))            # dummy
            returncounts = np.empty((1, 1), dtype=np.uint16)  # dummy

            # require
            slit_x = np.require(slit_x, requirements=ci.req_in, dtype=np.float64)
            slit_y = np.require(slit_y, requirements=ci.req_in, dtype=np.float64)
            wl = np.require(wl, requirements=ci.req_in, dtype=np.float64)
            returnx = np.require(returnx, requirements=ci.req_out, dtype=np.float64)
            returny = np.require(returny, requirements=ci.req_out, dtype=np.float64)
            returnwaves = np.require(returnwaves, requirements=ci.req_out, dtype=np.float64)
            returncounts = np.require(returncounts, requirements=ci.req_out, dtype=np.uint)

            func = ci.solve
            log.info("Raytracing order %s...", m)
            func(blaze_flag, return_mode, n, m, nxpix, nypix, f_col_1, f_col_2,
                alpha_ech, blaze_ech, gamma_ech, sigma_ech, alpha_cd, sigma_cd,
                f_cam, f_cam_1, dpix, xdl_0, xdm_0, xdr_0, ydl_0, ydm_0, ydr_0,
                tau_dl, tau_dm, tau_dr, slit_x, slit_y, wl, returnx, returny,
                returnwaves, returncounts)
            #print returnx, returny
            print m, 'done', n
            solvex[:, i] = returnx
            solvey[:, i] = returny
            interpx[:, i] = xmid
            interpy[:, i] = ymid


        solvexoff = (np.max(solvex[10, :]) + np.min(solvex[10, :])) / 2.0
        solveyoff = (np.max(solvey[10, :]) + np.min(solvey[10, :])) / 2.0
        interpxoff = (np.max(interpx[10, :]) + np.min(interpx[10, :])) / 2.0
        interpyoff = (np.max(interpy[10, :]) + np.min(interpy[10, :])) / 2.0
        solvex -= solvexoff
        solvey -= solveyoff
        interpx -= interpxoff
        interpy -= interpyoff
        solvex /= self.dpix
        solvey /= self.dpix
        interpx /= self.dpix
        interpy /= self.dpix
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='equal')
        ax.scatter(solvex, solvey, c='r', marker='x', label='solve')
        ax.scatter(interpx, interpy, c='b', marker='o', label='interp')
        plt.legend()
        plt.show()
        exit()
