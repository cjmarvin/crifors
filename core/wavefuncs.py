"""wavefuncs.py

Module containing wavelength functions.
"""
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from astropy.io import fits
import ctypes as ct
import logging
import os
import pprint
import time
import cinterface as ci
from defaults import tell_dir, logs_dir

log = logging.getLogger(__name__)


def feed_spectrum(inst, m, wavelengths, pdf):
    """
    calculates the wavelengths and density function to be fed into
    the simulator for a given order

    Wavelengths MUST be in nm.
    """
    log.info("Feeding spectrum into spectral order %s.", m)
    i = m - inst.m_range_min
    wmin = inst.det_wl_lim[:, i].min()
    wmax = inst.det_wl_lim[:, i].max()
    idx = (wavelengths >= wmin) & (wavelengths <= wmax)
    wavelengths = np.array(wavelengths[idx])
    pdf = np.array(pdf[idx])
    return wavelengths, pdf


def sample_cdf(wavelengths, pdf, n, flag=3):
    """Returns an n sized array of wavelengths sampled from the pdf."""
    log.info("Sampling %s wavelengths from pdf...", n)
    out = np.empty(n, dtype=np.float64)
    wavelengths = np.require(wavelengths, requirements=ci.req_in, dtype=np.float64)
    pdf = np.require(pdf, requirements=ci.req_in, dtype=np.float64)
    out = np.require(out, requirements=ci.req_out, dtype=np.float64)
    func = ci.cdf.random_wave
    func.argtypes = [
        ci.array_1d_double,     # waves
        ci.array_1d_double,     # pdf
        ct.c_int,               # nw
        ct.c_int,               # interpolation flag
        ci.array_1d_double,     # outwaves
        ct.c_int]               # outwaves size
    func.restype = None
    func(wavelengths, pdf, pdf.size, flag, out, n)
    return out


def wavegrid(self, wavelengths=None, intensities=None, factor=10.0,
    buff=1.0e-7, kind=3, assign=False):
    """
    Re-evaluates input spectra on a pre-determined wavelength grid.
    """
    if wavelengths is None:
        wavelengths = self.wavelengths
    if intensities is None:
        intensities = self.intensities
    wmin = self.wmin - buff
    wmax = self.wmax + buff

    # determine number of samples based on desired resolution
    dw = self.nw * 1.e-7 # angstrom -> mm
    n = np.int(np.ceil( (wmax-wmin) / dw ))

    wavelengths2 = np.linspace(wmin, wmax, num=n)
    print "Re-evaluating wavelength grid to %i linearly spaced samples\n  dw = %4.e Angstrom" % (n, self.nw)

    f = InterpolatedUnivariateSpline(wavelengths, intensities)
    intensities2 = f(wavelengths2)

    if assign:
        self.wavelengths = wavelengths2
        self.intensities = intensities2


def load_telluric_lines(line_list, plot=False):
    """
    Adds telluric lines to spectra.
    """
    #wmin = self.wmin
    #wmax = self.wmax

    # import all fluxes
    log.info("Loading telluric spectra.")
    fluxes = []
    species_fn = [os.path.join(tell_dir, "LBLRTM_%s_+0.0.npy" % x) for x in line_list]
    for species in species_fn:
        log.info("Loading '%s'.", species)
        data = np.load(species)
        fluxes.append(data[1])
    waves = data[0]   # wavelengths are same for all species
    log.info("Combining species together.")
    all_spec = np.prod(fluxes, axis=0)  # multiply all species together
    if plot:
        log.info("Opening plot...")
        import matplotlib.pyplot as plt
        plt.plot(waves, all_spec)
        plt.xlabel("Angstrom")
        plt.ylabel("Normalized flux")
        #plt.xlim([wmin, wmax])
        plt.show()
    return waves, all_spec


def convolve_telluric_lines(line_list, wavelengths, intensities):
    """
    Multiplies spectra by telluric lines.

    Parameters
    ----------
    line_list : list
        List of line species names.
    wavelengths : array_like
        Wavelength array in Angstroms.
    intensities : array_like
        Intensity array in arbitrary units.

    Returns
    -------
    intensities : array_like
        intensities convolved with telluric species.

    """
    log.info("Loading telluric species: %s", line_list)
    tellwaves, tellflux = load_telluric_lines(line_list)
    log.info("Convolving spectrum with telluric species.")
    ftell = InterpolatedUnivariateSpline(tellwaves, tellflux)
    tellflux = ftell(wavelengths)
    return intensities * tellflux


def truncate_spectrum(sim, wavelengths, pdf, buff=10.0):
    log.info("Truncating spectrum to spectrograph limits.")
    wmin = sim.det_wl_lim.min() - buff
    wmax = sim.det_wl_lim.max() + buff
    idx = (wavelengths >= wmin) & (wavelengths <= wmax)
    wavelengths = np.array(wavelengths[idx])
    pdf = np.array(pdf[idx])
    return wavelengths, pdf


def main():
    load_telluric_lines("CO2 H2O O3 N2O CH4 O2".split(), plot=True)

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
    console.setLevel(logging.INFO)
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)
    main()
