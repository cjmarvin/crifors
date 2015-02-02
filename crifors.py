#!/usr/bin/env python
"""
CriForS v{0}
CRIRES+ Forward Simulator feeds an input source spectrum into the future ESO
instrument CRIRES+ (Oliva et al. 2014) and outputs a raw 2D image.

Usage:
    crifors.py BAND [SOURCE]... [-bnt] [--bglight=BGLIGHT] [--config=CONFIG]
               [--dlamb=DLAMB] [--ds9] [--echang=ECHANG] [--factor=FACTOR]
               [--model=MODEL] [--nrays=NRAYS] [--nruns=NRUNS] [--outfn=OUTFN]
               [--plot | --plot-psf | plot-simple] [--psf=PSF] [--rv=RV]
               [--seeing=SEEING] [--slit-width=SLIT] [--verbose=LEVEL]
               [--spread]
    crifors.py [-h] | [--help] | [--version]

Arguments:
    BAND               spectral band
    SOURCE...          input file/s or source spectra [Default: P]

General options:
    -h, --help         show this help message and exit
    --version          show version and exit

Simulation options:
    -b, --blaze        include blaze efficiency
    -n, --noise        include noise in the simulation
    -t, --telluric     include telluric lines
    -m, --model=MODEL  computation model [Default: interp]
    --bglight=BGLIGHT  background light file
    --dlamb=DLAMB      wavelength grid resolution in nm [Default: 1e-5]
    --echang=ECHANG    incident echelle angle [Default: 63.5]
    --factor=FACTOR    factor to multiply wavelengths to convert to nm
    --nrays=NRAYS      number of rays [Default: 1e7]
    --nruns=NRUNS      number of simulation runs [Default: 1]
    --psf=PSF          psf [Default: gaussian]
    --rv=RV            radial velocity shift in m/s [Default: 0.0]
    --seeing=SEEING    seeing in arcseconds [Default: 1.5]
    --slit-width=SLIT  width of slit in arcseconds [Default: 0.2]
    --spread           spread out each ray by convolving with a kernel

Other options:
    --config=CONFIG    simulation config file
    --verbose=LEVEL    verbosity level [Default: INFO]
    --outfn=OUTFN      specific output filename
    --plot             open raytrace plot after simulation and exit
    --plot-psf         preview slit psf function before simulation and exit
    --plot-simple      open simple raytrace plot and exit
    --ds9              open simulated image in SAO-DS9

Argument details:
BAND
    This is a required argument. Must be one of the following spectral
    bands: (Y, J, H, K, L, M).
SOURCE...
    The input source of the spectrograph.
    Input choices:
    (p, phoenix) : included PHOENIX synthetic spectrum, T_eff = 3000 K,
        [M/H] = 0.0, log(g) = 5.0
        This is the default SOURCE argument.
    (f, flatfield) : ideal flatfield spectrum. Assumes a flat, uniform
        distribution (no lamp emission as of yet).
    (w, wavemap) : wavemap. Instead of count values as pixels, the wavelength
        value is given.  NOTE not currently supported.
    (<path/to/spectrum.fits/txt>) : path and filename of an input spectrum.
        It assumes a 2 column fits or txt file, with wavelengths units in [nm]
        in the 1st column and flux in the 2nd column.  The filetype can be
        either a fits or txt extension.
    (<path/to/wavelengths.fits/txt> <path/to/flux.fits/txt>) : paths and
            filenames of an input spectrum divided into 2 files.
        This assumes 2 fits or txt files, with wavelengths units in [nm]
        as the first filepath and flux as the 2nd filepath.

Option details:
-b, --blaze
    This flag turns on the option to include blaze efficiency in the
    simulation. It should be noted that this is only a rough approximation.
-n, --noise
    This flag turns on noise in the simulation. By default, shot noise is
    inherent in the raytracing process, but this option is needed to include
    dark current and readout noise.
-t, --tell
    This flag adds telluric lines to the input spectrum.  Lines are calculated
    using LBLRTM code (see Husser & Ulbrich, 2014).
-m, --model=MODEL
    This can be of the choices (interp, solve).  As of now, 'solve' is not
    currently implemented.
--bglight=BGLIGHT
    This is the background light file. Not currently supported.
--config=CONFIG
    Input configuration file to change instrument or simulation parameters.
--dlamb=DLAMB
    The wavelength grid resolution in nm.  This typically does not need to be
    altered.
--ds9
    Calls SAO-DS9 from the shell before exiting.
--echang=ECHANG
    The incident echelle angle.
--factor=FACTOR
    If the wavelengths file is in units other than nm, then supply the factor
    to get it in nm.
--nrays=NRAYS
    Input number of rays to simulate. Note that due to rays not hitting the
    detector, and order overlap, that will not be the exact number of rays
    that actually hit the detector. It is only an estimate.
--nruns=NRUNS
    Number of simulation iterations. This can be used to increase NRAYS
    without sacrificing memory. Also can be used to simulate observing over
    time, (ie. tracking, seeing changes, etc.) but this is not yet implemented.
--outfn=OUTFN
    Output filename or path.
--plot
    Opens an interactive plot of the simulated image and input spectrum
    before finishing.
--plot-psf
    Opens an interactive plot of the slit psf function and exits.
--psf=PSF
    Point spread function of source. Currently only 'gaussian' is supported.
    (point, uniform, gaussian)
--rv=RV
    Radial velocity shift in m/s of source.
--seeing=SEEING
    Seeing in arcseconds.
--slit-width=SLIT
    Width of slit in arcseconds (0.2, 0.4)
--spread
    tbw

Examples:
>> python crifors.py J --telluric --noise --nrays=1e9
    This will simulate a J band image using the included PHOENIX spectrum.
    It will include telluric lines and noise, and will begin with an input
    number of rays of 1e9.

>> python crifors.py L ~/Desktop/spectrum.fits --plot
    This will create an L band image of the source spectrum in the user's
    Desktop, and will open an interactive plot of the image and the spectrum
    prior to finishing.

>> python crifors.py K --seeing=2.0 --factor=0.1 ~/documents/waves.txt \
~/documents/flux.fits
    This will decrease the seeing to 2 arcseconds, and will convert the
    wavelength file, given in Angstroms, to nm.

"""
__author__ = 'Christopher J. Marvin, Mathias Zechmeister'

import numpy as np
import astropy.io.fits as fits
import os
import subprocess
import sys
import time

from version import __version__
import core
from defaults import paths
import logging

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                              MAIN PROGRAM
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

def main():

    # PROGRAM START TIME + DATE STAMPS
    t0 = time.time()
    d0 = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())

    # PARSE COMMAND LINE ARGUMENTS
    _doc_fmt = [__version__]
    args = core.docopt(__doc__.format(*_doc_fmt), version=__version__)

    log.info("CRIRES+ Forward Simulator (CRIForS) v%s", __version__)
    log.info("Start time: %s", d0)

    # INITIALIZE INSTRUMENT
    instrument = core.Instrument(args)

    # INITIALIZE SIMULATION
    simulator = core.Simulator(instrument)

    # RUN SIMULATION
    simulator.run()

    # SPREAD OUT EACH RAY
    if args["--spread"]:
        simulator.spreadout()

    # ADD NOISE
    if args["--noise"]:
        core.add_noise(simulator)

    if args["--plot"]:
        log.info("Opening plot...")
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
        fig = plt.figure("CRIRES+ Simulation results")
        gs = gridspec.GridSpec(2, 1)
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])
        ax1.imshow(simulator.outarr, origin="lower", interpolation='nearest', cmap="hot")
        ax2.plot(simulator.source_spectrum[0], simulator.source_spectrum[1])
        ax2.set_xlabel("Wavelength (nm)")
        ax2.set_ylabel("PDF")
        plt.tight_layout()
        plt.show()
    t1 = time.time()
    d1 = time.strftime("%Y-%m-%d_%H.%M.%S", time.gmtime())
    log.info("End time: %s", d1)
    log.info("Run time: %.1f s", t1-t0)

    # WRITE TO FITS FILE
    core.write_to_fits(simulator)

    # SPAWN FITS OUTPUT IMAGE IN SAO-DS9
    if args["--ds9"]:
        if os.path.isfile(simulator.outpath + ".fits.gz"):
            shell_call = "ds9 %s.fits.gz" % simulator.outpath
        elif os.path.isfile(simulator.outpath + ".fits"):
            shell_call = "ds9 %s.fits" % simulator.outpath
        try:
            log.info("Executing '%s'", shell_call)
            subprocess.check_call(shell_call.split())
        except OSError, e:
            log.error(e, exc_info=False)
            log.info("Shell call failed. Just run the following line:", exc_info=False)
            log.info(shell_call, exc_info=False)

    sys.exit(0)

if __name__ == "__main__":
    pname = os.path.splitext(os.path.basename(__file__))[0]
    log_path = os.path.join(paths.logs_dir, pname)
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
