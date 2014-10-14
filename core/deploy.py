"""plant.py

This function pickles the included FITS files into the NumPy NPY binary format,
allowing for much faster reading and loading.
"""
import numpy as np
from astropy.io import fits
import glob
import logging
import os
import pprint
import sys
import time
try:
    from defaults import logs_dir, tell_dir, tell_species
except ImportError:
    sys.path.append(os.getcwd())
    from defaults import logs_dir, tell_dir, tell_species

log = logging.getLogger(__name__)

def plant_telluric_spectra():
    """
    Loads fits.gz files and saves to an npy file.

    wavelengths : mm
    flux : normalized to unity

    """
    spec_fn = glob.glob(os.path.join(tell_dir, "LBLRTM*0.0.fits.gz"))
    for species in spec_fn:
        log.info("Loading '%s'...", species)
        data = fits.getdata(species)
        #data[0] *= 1.0e-7    # convert Ang to mm
        plant_fn = species.replace("fits.gz", "npy")
        log.info("Writing '%s' to disk.", plant_fn)
        np.save(plant_fn, data)


# ============================== MAIN ========================================

def main():

    process = {
        'telluric' : plant_telluric_spectra}

    for arg in sys.argv[1:]:
        if arg in process:
            process[arg]()

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
