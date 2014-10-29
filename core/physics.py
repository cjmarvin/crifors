"""physics.py
"""
import numpy as np
import logging

log = logging.getLogger(__name__)

c = 299792458.0       # [m/s]
h = 6.62606957e10-34  # m^2 kg / s

def lam_blaze_ech(m, sigma, alpha, beta):
    """Blaze Wavelength of Grating"""
    return 2.0 * sigma * np.sin((beta+alpha) / 2.) * np.cos((beta-alpha) / 2.) / m

def fsr(lam_blaze, m):
    """The wavelength limits for each echelle order"""
    return lam_blaze / m

def redshift(wavelengths, rv=0.0):
    """rv [m/s]"""
    log.info("Shifting spectra by %s m/s.", rv)
    return wavelengths * (1.0 + rv/c)


def energy2counts(wavelengths, flux):
    """Convert energy per wavelength to counts per wavelength.

    Parameters
    ----------
    wavelengths : array_like
        Wavelengths in m.
    flux : array_like
        Energy per wavelength in some form of ergs/s.

    Returns
    -------
    array_like
        Counts per wavelength.

    """
    assert wavelengths.shape == flux.shape
    #waves = wavelengths * 1.e-10 # Angstrom -> m
    return flux * wavelengths / (h*c)

# MISCELLANEOUS FUNCTIONS
def find_nearest(array, target):
    """Finds nearest value in an array. Returns the index and array value"""
    i = np.argmin(np.abs(array - target))
    return i, array[i]

