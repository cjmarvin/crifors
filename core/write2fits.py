"""write2fits.py

"""
import numpy as np
import astropy.io.fits as fits
import logging
import os
import subprocess
import time
import astrodate
from defaults import output_dir, output_fn
from version import __version__

log = logging.getLogger(__name__)


def recurs_makedirs(out_dir):
    """Makes a directory.

    If out_dir exists, creates out_dir+1 so that the original out_dir is not
    overwritten.
    """
    try:
        os.makedirs(out_dir)
    except OSError:
        if out_dir[-1].isdigit():
            out_dir = out_dir[:-1] + str(int(out_dir[-1]) + 1)
        else:
            out_dir += "1"
        return recurs_makedirs(out_dir)
    else:
        return out_dir


def output_path(out_fn):
    """Sets up output path."""
    if out_fn is None:
        out_fn = output_fn % (
            time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime()))
    return os.path.join(output_dir, out_fn)


def add_simulation_keywords(header, sim):
    log.info("Adding simulation keywords.")
    # FOR CONVENIENCE AND CLARITY, ASSIGN HEADER VARIABLES UP HERE
    DATE = time.strftime("%Y-%m-%d", time.localtime())
    FILENAME = os.path.basename(sim.outpath)
    DATE_OBS = time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime())
    MJD_OBS = astrodate.AstroDate().mjd
    INFILE1 = sim.FITS_INFILE1
    INFILE2 = sim.FITS_INFILE2
    # HIERARCH BUG IN ASTROPY
    a = len('HIERARCH ESO SIMU IN')
    c = len('Infile 1')
    limit = 73
    if INFILE1:
        b = len(INFILE1)
        tot = a+b+c
        if a+b+c > limit:
            newlength = limit-tot
            INFILE1 = INFILE1[:newlength]
    if INFILE2:
        b = len(INFILE2)
        tot = a+b+c
        if a+b+c > limit:
            newlength = limit-tot
            INFILE2 = INFILE2[:newlength]
    if sim.telluric:
        CH4 = 'CH4' in sim.telluric
        CO2 = 'CO2' in sim.telluric
        H2O = 'H2O' in sim.telluric
        N2O = 'N2O' in sim.telluric
        O2 = 'O2' in sim.telluric
        O3 = 'O3' in sim.telluric
    else:
        CH4 = False
        CO2 = False
        H2O = False
        N2O = False
        O2 = False
        O3 = False
    TIME = np.float("%.2f" % sim.sim_time)

    header['HIERARCH ESO SIMU VERSION'] = (__version__, 'Simulation version')
    header['HIERARCH ESO SIMU FILENAME'] = (FILENAME, 'Filename')
    header['HIERARCH ESO SIMU TYPE'] = (sim.FITS_TYPE, 'Type')
    header['HIERARCH ESO SIMU CATG'] = (sim.FITS_CATG, 'Category')
    header['HIERARCH ESO SIMU MJD'] = (MJD_OBS, 'MJD')
    header['HIERARCH ESO SIMU SOURCE'] = (sim.FITS_SOURCE, 'Source')
    header['HIERARCH ESO SIMU MODEL'] = (sim.model, 'Simulation model')
    header['HIERARCH ESO SIMU ECH ANG'] = (sim.echang, 'Echelle angle')
    header['HIERARCH ESO SIMU TOBS'] = (sim.tobs, '[s] Duration of observation')
    header['HIERARCH ESO SIMU INPUT SRC'] = (sim.FITS_SRC, 'Input source')
    header['HIERARCH ESO SIMU IN 1'] = (INFILE1, 'Infile1 (spectrum or wavelengths)')
    header['HIERARCH ESO SIMU IN 2'] = (INFILE2, 'Infile2 (flux)')
    header['HIERARCH ESO SIMU DW'] = (sim.dlamb, '[nm] Input wavelength grid sampling')
    header['HIERARCH ESO SIMU OMIN'] = (sim.m_range_min, 'Minimum echelle order')
    header['HIERARCH ESO SIMU OMAX'] = (sim.m_range_max, 'Maximum echelle order')
    header['HIERARCH ESO SIMU RV'] = (sim.rv, '[m/s] Input radial velocity')
    header['HIERARCH ESO SIMU BLAZE'] = (sim.blaze, 'Blaze function flag')
    header['HIERARCH ESO SIMU TELL CH4'] = (CH4, 'Telluric species included')
    header['HIERARCH ESO SIMU TELL CO2'] = (CO2, 'Telluric species included')
    header['HIERARCH ESO SIMU TELL H2O'] = (H2O, 'Telluric species included')
    header['HIERARCH ESO SIMU TELL N2O'] = (N2O, 'Telluric species included')
    header['HIERARCH ESO SIMU TELL O2'] = (O2, 'Telluric species included')
    header['HIERARCH ESO SIMU TELL O3'] = (O3, 'Telluric species included')
    header['HIERARCH ESO SIMU SLIT WIDTH'] = (sim.slit_width, '[arcsec] Slit width')
    header['HIERARCH ESO SIMU SRC PSF'] = (sim.psf, 'Source PSF')
    header['HIERARCH ESO SIMU SRC PSF MUX'] = (sim.mu_x_psf, 'Source PSF x-center')
    header['HIERARCH ESO SIMU SRC PSF MUY'] = (sim.mu_y_psf, 'Source PSF y-center')
    header['HIERARCH ESO SIMU SRC PSF SIGX'] = (sim.sig_x_psf, 'Source PSF x-seeing')
    header['HIERARCH ESO SIMU SRC PSF SIGY'] = (sim.sig_y_psf, 'Source PSF y-seeing')
    header['HIERARCH ESO SIMU IN NRAYS'] = (sim.nrays, 'Number of input rays at cmdline')
    for i,N in enumerate(sim.nrays_per_order):
        header['HIERARCH ESO SIMU NRAYS %s' % (i+sim.m_range_min)] = (N, 'Input number of rays for order %s' % (i+sim.m_range_min))
    header['HIERARCH ESO SIMU MEAN NRAYS'] = (sim.mean_rays_per_pixel, 'Mean rays per raw pixel')
    header['HIERARCH ESO SIMU MED NRAYS'] = (sim.med_rays_per_pixel, 'Median rays per raw pixel')
    header['HIERARCH ESO SIMU MIN NRAYS'] = (sim.min_rays_per_pixel, 'Min rays per raw pixel')
    header['HIERARCH ESO SIMU MAX NRAYS'] = (sim.max_rays_per_pixel, 'Max rays per raw pixel')
    header['HIERARCH ESO SIMU TOT NRAYS'] = (sim.nrays_tot, 'Total rays at detector')
    header['HIERARCH ESO SIMU DL ID'] = (sim.dl_id, 'Left detector ID')
    header['HIERARCH ESO SIMU DL XD0'] = (sim.xdl, 'Left detector x-offset')
    header['HIERARCH ESO SIMU DL YD0'] = (sim.ydl, 'Left detector y-offset')
    header['HIERARCH ESO SIMU DL TAU'] = (sim.tau_dl, 'Left detector rotation')
    header['HIERARCH ESO SIMU DL RON'] = (sim.dl_ron, '[e-] Left detector readout noise')
    header['HIERARCH ESO SIMU DL DARK'] = (sim.dl_dc, '[e-/s] Left detector dark current')
    header['HIERARCH ESO SIMU DL BIAS'] = (sim.dl_bias, '[e-] Left detector bias')
    header['HIERARCH ESO SIMU DM ID'] = (sim.dm_id, 'Middle detector ID')
    header['HIERARCH ESO SIMU DM XD0'] = (sim.xdm, 'Middle detector x-offset')
    header['HIERARCH ESO SIMU DM YD0'] = (sim.ydm, 'Middle detector y-offset')
    header['HIERARCH ESO SIMU DM TAU'] = (sim.tau_dm, 'Middle detector rotation')
    header['HIERARCH ESO SIMU DM RON'] = (sim.dm_ron, '[e-] Middle detector readout noise')
    header['HIERARCH ESO SIMU DM DARK'] = (sim.dm_dc, '[e-/s] Middle detector dark current')
    header['HIERARCH ESO SIMU DM BIAS'] = (sim.dm_bias, '[e-] Middle detector bias')
    header['HIERARCH ESO SIMU DR ID'] = (sim.dr_id, 'Right detector ID')
    header['HIERARCH ESO SIMU DR XD0'] = (sim.xdr, 'Right detector x-offset')
    header['HIERARCH ESO SIMU DR YD0'] = (sim.ydr, 'Right detector y-offset')
    header['HIERARCH ESO SIMU DR TAU'] = (sim.tau_dr, 'Right detector rotation')
    header['HIERARCH ESO SIMU DR RON'] = (sim.dr_ron, '[e-] Right detector readout noise')
    header['HIERARCH ESO SIMU DR DARK'] = (sim.dr_dc, '[e-/s] Right detector dark current')
    header['HIERARCH ESO SIMU DR BIAS'] = (sim.dr_bias, '[e-] Right detector bias')
    header['HIERARCH ESO SIMU INV GAIN'] = (sim.inv_gain, '[e-/DN] inverse gain')
    header['HIERARCH ESO SIMU SIM TIME'] = (TIME, '[s] simulation time')


def add_default_keywords(header):
    log.info("Adding default keywords.")
    header['HIERARCH ESO OBS TARG NAME'] = ('input star spectrum', 'OB target name')
    header['HIERARCH ESO DPR CATG'] = ('SCIENCE', 'Observation category')
    header['HIERARCH ESO DPR TECH'] = ('SPECTRUM')
    header['HIERARCH ESO DPR TYPE'] = ('OBJECT', 'Observation type')
    header['HIERARCH ESO INS MODE'] = ('SCIENCE', 'Instrument mode')
    header['HIERARCH ESO INS FILT1 ENC'] = ('2666', 'Absolute position [Enc]')
    header['HIERARCH ESO INS FILT1 NAME'] = ('Ks', 'Element name')
    header['HIERARCH ESO INS FILT1 NO'] = (5, 'Element number')
    header['HIERARCH ESO INS GRAT ENC'] = (71999, 'Absolute position [Enc]')
    header['HIERARCH ESO INS GRAT ORDER'] = (12, 'Grating order')
    header['HIERARCH ESO INS GRAT SETNO'] = ('', 'Fixed echelle settings')
    header['HIERARCH ESO INS CROSS ENC']  = ('', 'Angle for the cross-disperser')
    header['HIERARCH ESO INS CROSS ORDER'] = ('', 'Order for the cross-disperser') # PROBABLY OBSOLETE
    header['HIERARCH ESO INS CROSS BAND'] = ('', 'Band name for the cross-disperser')
    header['HIERARCH ESO INS SLIT1 ENC'] = (5564, 'Absolute position [Enc]')
    header['HIERARCH ESO INS SLIT1 POS'] = (0.201, 'Position [arcsec]')
    header['HIERARCH ESO INS SLIT1 WID'] = (0.210, 'Slit width [arcsec]')
    header['HIERARCH ESO INS SLIT2 ENC'] = (5010, 'Absolute position [Enc]')
    header['HIERARCH ESO INS SLIT2 POS'] = (0.570, 'Position [mm]')
    header['HIERARCH ESO INS SLIT2 DECKER'] = ('', 'Decker layout (position)')
    header['HIERARCH ESO INS POL POS'] = ('', 'Polarimeter selection')
    header['HIERARCH ESO INS POL ANGLE'] = ('', 'Polarimeter rotation angle (0 or 180)')
    header['HIERARCH ESO DET CON OPMODE'] = ('NORMAL', 'Operational Mode')
    header['HIERARCH ESO DET DIT'] = (10.0, 'Integration Time')
    header['HIERARCH ESO DET FRAM NO'] = (1, 'Frame number')
    header['HIERARCH ESO DET FRAM TYPE'] = ('INT', 'Frame type')
    header['HIERARCH ESO DET NDIT'] = (6, '# of Sub-Integrations')


def write_to_fits(sim, gzip=True):
    """Write simulation results to fits file.

    Write header and data to FITS HDU. Compress if necessary.
    """

    # CREATE OUTPUT PATH
    sim.outpath = output_path(sim.outfn)

    # create PrimaryHDU object to encapsulate data
    log.info("Creating HDU...")
    hdu = fits.PrimaryHDU(np.asarray(sim.outarr, dtype=np.uint16))

    # create ImageHDU objects for detector images
    hdu_dl = fits.ImageHDU(sim.outarr[:, :sim.nxpix])
    hdu_dm = fits.ImageHDU(sim.outarr[:, sim.nxpix:2*sim.nxpix])
    hdu_dr = fits.ImageHDU(sim.outarr[:, 2*sim.nxpix:3*sim.nxpix])
    hdulist = fits.HDUList([hdu, hdu_dl, hdu_dm, hdu_dr])
    header = hdu.header

    # DEFAULT KEYWORDS
    add_default_keywords(header)

    # SIMULATION KEYWORDS
    add_simulation_keywords(header, sim)

    # WRITE TO FILE
    try:
        os.makedirs(output_dir)
    except OSError:
        pass
    log.info("Writing HDU to '%s.fits'", sim.outpath)
    hdulist.writeto('%s.fits' % sim.outpath, clobber=True)

    # COMPRESS
    if gzip:
        log.info("Compressing to '%s.fits.gz'.", sim.outpath)
        shcall = "gzip %s.fits" % sim.outpath
        subprocess.check_call(shcall.split())
    return 0
