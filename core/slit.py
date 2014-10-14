"""
slit.py

Generates slit functions according to parameters.

Functions to populate the slit psf are written in c and accessed via the
cinterface.py module.

"""

import numpy as np
import ctypes as ct
import logging
import os
import sys
import time
import cinterface as ci
from defaults import *

log = logging.getLogger(__name__)


def preview(ax, x, y, circle=False):
    import matplotlib.pylab as plt
    ax.scatter(x, y, alpha=0.5, edgecolors="none")
    if circle:
        an = np.linspace(0, 2*np.pi, 100)
        ax.plot(circle*np.cos(an), circle*np.sin(an), color="r", linewidth=3, alpha=0.6)
    #ax.set_xlabel(r"$x$ [mm]")
    #ax.set_ylabel(r"$y$ [mm]")
    #ax.set_title(title)
    return ax


def slit_filter(params):
    """Passes source psf through slit.

    Arguments
    ---------
    params : object
        Data container of instrument and simulation parameters.

    Returns
    -------
    slit_x: array_like
        x-coordinates
    slit_y: array_like
        y-coordinates
    """
    pass


def point_source(yoffset=0.0, plot=False):
    x = 0.0
    y = 0.0 + yoffset

    # scale slit function
    #x = x * arm.FN
    #y = y * arm.FN

    ## apply rotation
    #xp = (x * np.cos(arm.TAU_0)) + (y * np.sin(arm.TAU_0))
    #yp = (-x * np.sin(arm.TAU_0)) + (y * np.cos(arm.TAU_0))

    # convert to array for looping compatibility
    slit = np.array((x, y))
    log.info("Point source initialized.")

    # preview slit
    if plot:
        log.info("Opening preview plot of point source at slit center (%.1f, %.1f)." % (x, y))
        import matplotlib.pylab as plt
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='equal')
        ax.scatter(slit[0], slit[1], s=20, edgecolor=None)
        plt.title("0D Point Source PSF")
        plt.show()
    #arm.slittyp = "0D point-source"
    return slit


def slit_uniform_psf(n, seeing, mu_x, mu_y, tau_0, slit_width, slit_height, plot=False):
    """Returns x- and y- coordinate arrays of a 2D random uniformly distributed
    circle.

    Parameters
    ----------
    n : int
        Size of coordinate arrays.
    seeing: double
        Seeing of source psf in arcseconds.
    mu_x : double
        Center of PSF in x-coords.
    mu_y : double
        Center of PSF in y-coords.
    tau_0 : double
        Rotation about z-axis (tilt).
    slit_width : double
        Width of slit in arcseconds.
    slit_height : double
        Height of slit in arcseconds.

    Returns
    -------
    slit_x : array_like
        Array of x-coordinates.
    slit_y : array_like
        Array of y-coordinates.

    """
    desc = "Source psf: uniform, mux=%.2f muy=%.2f seeing=%.2f arcsec" % (mu_x, mu_y, seeing)
    log.info(desc)
    # initialize output arrays to send to c function
    slit_x = np.empty(n, dtype=np.float64)
    slit_y = np.empty(n, dtype=np.float64)
    slit_x = np.require(slit_x, requirements=ci.req_out, dtype=np.float64)
    slit_y = np.require(slit_y, requirements=ci.req_out, dtype=np.float64)
    func = ci.slitc.slit_uniform_psf
    func.argtypes = [
        ct.c_int,             # n
        ct.c_double,          # seeing
        ct.c_double,          # mu_x
        ct.c_double,          # mu_y
        ct.c_double,          # tau_0
        ct.c_double,          # slit_width
        ct.c_double,          # slit_height
        ci.array_1d_double,   # slit_x
        ci.array_1d_double]   # slit_y
    func.restype = None
    log.info("Slit Rejection Sampling: %s rays...", n)
    func(n, seeing, mu_x, mu_y, tau_0, slit_width, slit_height, slit_x, slit_y)
    # preview slit
    if plot:
        log.info("Opening preview plot of 2D uniformly random psf.")
        import matplotlib.pylab as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)#, aspect='equal')
        ax.scatter(slit_x, slit_y, s=20, edgecolor=None)
        plt.title("0D Point Source PSF")
        plt.show()
    return slit_x, slit_y


def slit_gaussian_psf(n, mu_x, mu_y, sig_x, sig_y, tau_0, slit_width, slit_height, plot=False):
    """Returns x- and y- coordinate arrays of a 2D random gaussian
    distribution.

    Parameters
    ----------
    n : int
        Size of coordinate arrays.
    mu_x : double
        Center of Gaussian in x-coords.
    mu_y : double
        Center of Gaussian in y-coords.
    sig_x : double
        Width of Gaussian in x-coords.
    sig_y : double
        Width of Gaussian in y-coords.
    tau_0 : double
        Rotation about z-axis (tilt).
    slit_width : double
        Width of slit in arcseconds.
    slit_height : double
        Height of slit in arcseconds.

    Returns
    -------
    slit_x : array_like
        Array of x-coordinates.
    slit_y : array_like
        Array of y-coordinates.

    """
    desc = "Source psf: gaussian, mux=%.2f muy=%.2f sigx=%.2f sigy=%.2f" % (mu_x, mu_y, sig_x, sig_y)
    log.info(desc)
    slit_x = np.empty(n, dtype=np.float64)
    slit_y = np.empty(n, dtype=np.float64)
    slit_x = np.require(slit_x, requirements=ci.req_out)
    slit_y = np.require(slit_y, requirements=ci.req_out)
    func = ci.slitc.slit_gaussian_psf
    func.argtypes = [
        ct.c_int,             # n
        ct.c_double,          # mu_x
        ct.c_double,          # mu_y
        ct.c_double,          # sig_x
        ct.c_double,          # sig_y
        ct.c_double,          # tau_0
        ct.c_double,          # slit_width
        ct.c_double,          # slit_height
        ci.array_1d_double,   # slit_x
        ci.array_1d_double]   # slit_y
    func.restype = None
    log.info("Slit Rejection Sampling: %s rays...", n)
    func(n, mu_x, mu_y, sig_x, sig_y, tau_0, slit_width, slit_height, slit_x, slit_y)
    # preview slit
    if plot:
        log.info("Opening preview plot of 2D Gaussian psf.")
        import matplotlib.pylab as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)#, aspect='equal')
        ax.scatter(slit_x, slit_y, s=20, edgecolor=None)
        plt.title("2D Gaussian PSF")
        plt.show()
    return slit_x, slit_y


# ==============================[ TESTS ]======================================


def _point_test():
    point_source(plot=True)


def _uniform_test():
    import matplotlib.pyplot as plt
    n = 10000
    slit_x = np.empty(n, dtype=np.float64)
    slit_y = np.empty(n, dtype=np.float64)
    tau = 0.0
    mu_x = 0.0
    mu_y = 0.0
    seeing = 1.5
    slit_width = 0.2
    slit_height = 10.0
    slit_x, slit_y = slit_uniform_psf(n, seeing, mu_x, mu_y, tau, slit_width, slit_height)
    log.info("x range: [%s, %s]", slit_x.min(), slit_x.max())
    log.info("y range: [%s, %s]", slit_y.min(), slit_y.max())
    plt.scatter(slit_x, slit_y)
    plt.fill([-slit_width/2, slit_width/2, slit_width/2, -slit_width/2],
             [-slit_height/2, -slit_height/2, slit_height/2, slit_height/2],
             'r',
             alpha=0.10,
             edgecolor='k')
    plt.title("Random uniform distribution")
    plt.show()


def _gaussian_test():
    import matplotlib.pyplot as plt
    n = 10000
    mu_x = 0.0
    mu_y = 0.0
    #sig_x, sig_y = 1.5, 1.5
    tau = 0.0
    seeing = 1.5
    sigma = seeing / (2. * np.sqrt(2. * np.e))
    slit_width = 0.2
    slit_height = 10.0
    slit_x = np.empty(n, dtype=np.float64)
    slit_y = np.empty(n, dtype=np.float64)
    slit_x, slit_y = slit_gaussian_psf(n, mu_x, mu_y, sigma, sigma, tau, slit_width, slit_height)
    log.info("x range: [%s, %s]", slit_x.min(), slit_x.max())
    log.info("y range: [%s, %s]", slit_y.min(), slit_y.max())
    plt.scatter(slit_x, slit_y, alpha=0.8)
    plt.fill([-slit_width/2, slit_width/2, slit_width/2, -slit_width/2],
             [-slit_height/2, -slit_height/2, slit_height/2, slit_height/2],
             'r',
             alpha=0.10,
             edgecolor='k')
    plt.gca().set_aspect("equal")
    plt.title("Gaussian distribution")
    plt.xlim([-slit_height/2., slit_height/2])
    plt.show()


def main():
    _point_test()
    _uniform_test()
    _gaussian_test()


if __name__ == "__main__":
    pname = os.path.splitext(os.path.basename(__file__))[0]
    log_path = os.path.join(logs_dir, pname)
    if not os.path.isdir(log_path):
        os.makedirs(log_path)
    d0 = time.strftime("%Y-%m-%d_%H.%M.%S", time.gmtime())
    fmt = '%(levelname)8s %(asctime)s %(name)16s %(filename)16s:%(lineno)-5s: %(message)s'
    logging.basicConfig(level=logging.DEBUG,
                        format=fmt,
                        datefmt='%m-%d-%y %H:%M',
                        filename=os.path.join(log_path, "%s-%s.log" % (pname, d0)),
                        filemode="w")
    log = logging.getLogger(__name__)
    formatter = logging.Formatter('%(levelname)8s %(name)16s: %(message)s')
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)
    main()