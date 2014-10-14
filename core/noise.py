"""noise.py

Functions to simulate noise in spectral observations.
"""
import numpy as np
import logging

log = logging.getLogger(__name__)

def det_bias(bias, det=None):
    if det is not None:
        det_str = " to %s detector:" % det
    else:
        det_str = ":"
    log.info("Adding bias%s %s [e-].", det_str, bias)
    return bias

def shot_noise(image):
    print "Adding shot noise"
    return np.random.poisson(lam=image)

def readout_noise(rn, shape, det=None):
    if det is not None:
        det_str = " to %s detector:" % det
    else:
        det_str = ":"
    log.info("Adding readout noise%s rms = %s [e-].", det_str, rn)
    return np.random.normal(scale=rn, size=shape)

def dark_current(sig_dn, t, shape, det=None):
    if det is not None:
        det_str = " to %s detector:" % det
    else:
        det_str = ":"
    log.info("Adding dark noise%s rms = %s [e-].", det_str, sig_dn)
    return np.random.normal(scale=sig_dn*t, size=shape)

def gain(arr, invg):
    log.info("Converting e- to DN, gain = %s [DN / e-].", 1./invg)
    return arr / invg

def add_noise(sim):
    det_left = sim.outarr[:, :sim.nxpix]
    det_mid = sim.outarr[:, sim.nxpix:2*sim.nxpix]
    det_right = sim.outarr[:, 2*sim.nxpix:3*sim.nxpix]
    shape = det_left.shape

    det_left += det_bias(sim.dl_bias, det="left")
    det_mid += det_bias(sim.dm_bias, det="middle")
    det_right += det_bias(sim.dr_bias, det="right")

    det_left += readout_noise(sim.dl_ron, shape, det="left")
    det_mid += readout_noise(sim.dm_ron, shape, det="middle")
    det_right += readout_noise(sim.dr_ron, shape, det="right")

    det_left += dark_current(sim.dl_dc, sim.tobs, shape, det="left")
    det_mid += dark_current(sim.dm_dc, sim.tobs, shape, det="middle")
    det_right += dark_current(sim.dr_dc, sim.tobs, shape, det="right")

    sim.outarr = gain(sim.outarr, sim.inv_gain)

    if sim.outarr.max() > np.iinfo(np.uint16).max:
        log.info("Clipping array values larger than %s.", np.iinfo(np.uint16).max)
        sim.outarr[sim.outarr > np.iinfo(np.uint16).max] = np.iinfo(np.uint16).max
    sim.outarr = np.asarray(sim.outarr, dtype=np.uint16)    
    log.info("Converting image array back to %s.", sim.outarr.dtype)
    