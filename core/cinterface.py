"""cinterface.py

Interface between C code and Python.
"""

import numpy as np
import ctypes as ct
import os

# load c libraries
_path = os.path.dirname(os.path.abspath(__file__))
cdf = np.ctypeslib.load_library('cdf', _path)
raytrace = np.ctypeslib.load_library('raytrace', _path)
slitc = np.ctypeslib.load_library('slitfuncs', _path)

# define numpy array types
array_1d_double = np.ctypeslib.ndpointer(
    dtype=np.double,
    ndim=1,
    flags='C_CONTIGUOUS')
array_1d_int = np.ctypeslib.ndpointer(
    dtype=np.int,
    ndim=1,
    flags='C_CONTIGUOUS')
array_2d_double = np.ctypeslib.ndpointer(
    dtype=np.double,
    ndim=2,
    flags='C_CONTIGUOUS')
array_2d_int = np.ctypeslib.ndpointer(
    dtype=np.int,
    ndim=2,
    flags='C_CONTIGUOUS')
array_2d_uint = np.ctypeslib.ndpointer(
    dtype=np.uint,
    ndim=2,
    flags='C_CONTIGUOUS')

# numpy array requirements
req_in = "A C O".split()
req_out = "A C O W".split()

# ORIGINAL COMPUTE
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# declare arguments and return types for compute
# _lib_rt.compute_c_general.argtypes = [
#     ct.c_int,    # ARM ID
#     ct.c_int,    # BLAZE_FLAG
#     ct.c_int,    # LOC_FLAG
#     ct.c_int,    # GAP_FLAG
#     ct.c_ulong,  # nw
#     ct.c_ulong,  # nslit
#     ct.c_uint,   # m
#     ct.c_double, # gap
#     ct.c_double, # gapoffset
#     ct.c_double, # xd_0
#     ct.c_double, # yd_0
#     array_1d_double, # n_sell
#     array_1d_double, # slitx
#     array_1d_double, # slity
#     array_1d_double, # lamb
#     array_1d_double, # weights
#     array_2d_double, # ccd
#     array_2d_uint,   # counts
#     array_2d_uint,   # m_list
#     array_1d_double, # returnx
#     array_1d_double  # returny
#     ]
# _lib_rt.compute_c_general.restype = None
#
# def compute(
#     arm_flag,
#     blaze_flag,
#     location_flag,
#     gap_flag,
#     nwaves,
#     ns,
#     m,
#     gap,
#     gapoffset,
#     xd_0,
#     yd_0,
#     n_sell,
#     slitx,
#     slity,
#     waves,
#     weights,
#     ccd,
#     counts,
#     m_list,
#     returnx,
#     returny):
#
#     _lib_rt.compute_c_general(
#         arm_flag,
#         blaze_flag,
#         location_flag,
#         gap_flag,
#         nwaves,
#         ns,
#         m,
#         gap,
#         gapoffset,
#         xd_0,
#         yd_0,
#         n_sell,
#         slitx,
#         slity,
#         waves,
#         weights,
#         ccd,
#         counts,
#         m_list,
#         returnx,
#         returny)
#     return 0
# # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#
# # COMPUTE GRID (compute c without general options)
# # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# # declare arguments and return types for compute
# _lib_rt.compute_c_grid.argtypes = [
#     ct.c_int,    # ARM ID
#     ct.c_int,    # BLAZE_FLAG
#     ct.c_int,    # GAP_FLAG
#     ct.c_ulong,  # nw
#     ct.c_ulong,  # nslit
#     ct.c_uint,   # m
#     ct.c_double, # gap
#     ct.c_double, # gapoffset
#     ct.c_double, # xd_0
#     ct.c_double, # yd_0
#     array_1d_double, # n_sell
#     array_1d_double, # slitx
#     array_1d_double, # slity
#     array_1d_double, # lamb
#     array_1d_double, # weights
#     array_2d_double, # ccd
#     array_2d_uint   # counts
#     ]
# _lib_rt.compute_c_grid.restype = None
#
# def compute_grid(
#     arm_flag,
#     blaze_flag,
#     gap_flag,
#     nwaves,
#     ns,
#     m,
#     gap,
#     gapoffset,
#     xd_0,
#     yd_0,
#     n_sell,
#     slitx,
#     slity,
#     waves,
#     weights,
#     ccd,
#     counts):
#
#     _lib_rt.compute_c_grid(
#         arm_flag,
#         blaze_flag,
#         gap_flag,
#         nwaves,
#         ns,
#         m,
#         gap,
#         gapoffset,
#         xd_0,
#         yd_0,
#         n_sell,
#         slitx,
#         slity,
#         waves,
#         weights,
#         ccd,
#         counts)
#     return 0
# # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#
# # COMPUTE GRID WITH SLIT PERTURBATIONS
# # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# # declare arguments and return types for compute
# _lib_rt.compute_c_grid_perturb.argtypes = [
#     ct.c_int,    # ARM ID
#     ct.c_int,    # BLAZE_FLAG
#     ct.c_int,    # GAP_FLAG
#     ct.c_ulong,  # nw
#     ct.c_ulong,  # nslit
#     ct.c_uint,   # m
#     ct.c_double, # gap
#     ct.c_double, # gapoffset
#     ct.c_double, # xd_0
#     ct.c_double, # yd_0
#     ct.c_double, # perturb
#     array_1d_double, # n_sell
#     array_1d_double, # slitx
#     array_1d_double, # slity
#     array_1d_double, # lamb
#     array_1d_double, # weights
#     array_2d_double, # ccd
#     array_2d_uint   # counts
#     ]
# _lib_rt.compute_c_grid_perturb.restype = None
#
# def compute_grid_perturb(
#     arm_flag,
#     blaze_flag,
#     gap_flag,
#     nwaves,
#     ns,
#     m,
#     gap,
#     gapoffset,
#     xd_0,
#     yd_0,
#     perturb,
#     n_sell,
#     slitx,
#     slity,
#     waves,
#     weights,
#     ccd,
#     counts):
#     """c - wrapper"""
#
#     _lib_rt.compute_c_grid_perturb(
#         arm_flag,
#         blaze_flag,
#         gap_flag,
#         nwaves,
#         ns,
#         m,
#         gap,
#         gapoffset,
#         xd_0,
#         yd_0,
#         perturb,
#         n_sell,
#         slitx,
#         slity,
#         waves,
#         weights,
#         ccd,
#         counts)
#     return 0
# # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#
# # COMPUTE FULL MONTE CARLO METHOD USING REJECTION SAMPLING FOR SED SAMPLING
# # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# # declare arguments and return types for compute
# _lib_rt.compute_c_mc_rejsamp.argtypes = [
#     ct.c_int,    # ARM ID
#     ct.c_int,    # BLAZE_FLAG
#     ct.c_int,    # GAP_FLAG
#     ct.c_int,    # SLIT_FLAG
#     ct.c_ulong,  # nw
#     ct.c_ulong,  # nphotons
#     ct.c_uint,   # m
#     ct.c_double, # gap
#     ct.c_double, # gapoffset
#     ct.c_double, # xd_0
#     ct.c_double, # yd_0
#     ct.c_double, # flux max
#     ct.c_double, # offset
#     ct.c_double, # slit_locx
#     ct.c_double, # slit_locy
#     ct.c_double, # slit_scalex
#     ct.c_double, # slit_scaley
#     array_1d_double, # lamb
#     array_1d_double, # weights
#     array_2d_uint,   # ccd
#     array_2d_double  # wavemap values
#     ]
# _lib_rt.compute_c_mc_rejsamp.restype = None
#
# def compute_mc_rejsamp(
#     arm_flag,
#     blaze_flag,
#     gap_flag,
#     slit_flag,
#     nwaves,
#     nphotons,
#     m,
#     gap,
#     gapoffset,
#     xd_0,
#     yd_0,
#     fluxmax,
#     offset,
#     slit_locx,
#     slit_locy,
#     slit_scalex,
#     slit_scaley,
#     waves,
#     weights,
#     ccd,
#     wavemap_values):
#     """c - wrapper"""
#
#     _lib_rt.compute_c_mc_rejsamp(
#         arm_flag,
#         blaze_flag,
#         gap_flag,
#         slit_flag,
#         nwaves,
#         nphotons,
#         m,
#         gap,
#         gapoffset,
#         xd_0,
#         yd_0,
#         fluxmax,
#         offset,
#         slit_locx,
#         slit_locy,
#         slit_scalex,
#         slit_scaley,
#         waves,
#         weights,
#         ccd,
#         wavemap_values)
#     return 0
# # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#
# # COMPUTE FULL MONTE CARLO METHOD USING CDF FOR SED SAMPLING
# # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# # sampling CDF
# _lib_rt.compute_c_mc_cdf.argtypes = [
#     ct.c_int,    # ARM ID
#     ct.c_int,    # BLAZE_FLAG
#     ct.c_int,    # GAP_FLAG
#     ct.c_int,    # SLIT_FLAG
#     ct.c_ulong,  # nphotons
#     ct.c_uint,   # m
#     ct.c_double, # gap
#     ct.c_double, # gapoffset
#     ct.c_double, # xd_0
#     ct.c_double, # yd_0
#     ct.c_double, # offset
#     ct.c_double, # slit_locx
#     ct.c_double, # slit_locy
#     ct.c_double, # slit_scalex
#     ct.c_double, # slit_scaley
#     array_1d_double, # lamb
#     array_2d_uint,   # ccd
#     array_2d_double  # wavelength values
#     ]
# _lib_rt.compute_c_mc_cdf.restype = None
#
# def compute_mc_cdf(
#     arm_flag,
#     blaze_flag,
#     gap_flag,
#     slit_flag,
#     nphotons,
#     m,
#     gap,
#     gapoffset,
#     xd_0,
#     yd_0,
#     offset,
#     slit_locx,
#     slit_locy,
#     slit_scalex,
#     slit_scaley,
#     waves,
#     ccd,
#     wavemap_values):
#     """c - wrapper"""
#
#     _lib_rt.compute_c_mc_cdf(
#         arm_flag,
#         blaze_flag,
#         gap_flag,
#         slit_flag,
#         nphotons,
#         m,
#         gap,
#         gapoffset,
#         xd_0,
#         yd_0,
#         offset,
#         slit_locx,
#         slit_locy,
#         slit_scalex,
#         slit_scaley,
#         waves,
#         ccd,
#         wavemap_values)
#     return 0
# # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#
# # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# # declare arguments and return types for wavelength_grid_test
# _lib_rt.wavelength_rejection_sampling.argtypes = [
#     ct.c_ulong,  # nw
#     ct.c_ulong,  # nwout
#     array_1d_double, # waves
#     array_1d_double, # weights
#     array_1d_double, # env
#     array_1d_double # outwaves
#     ]
# _lib_rt.wavelength_rejection_sampling.restype = None
#
# def wavelength_rejection_sampling(
#     nw,
#     nwout,
#     waves,
#     weights,
#     env,
#     outwaves
#     ):
#     """c - wrapper"""
#     waves = np.ascontiguousarray(waves, dtype=np.float64)
#     weights = np.ascontiguousarray(weights, dtype=np.float64)
#     env = np.ascontiguousarray(env, dtype=np.float64)
#     outwaves = np.ascontiguousarray(outwaves, dtype=np.float64)
#     _lib_rt.wavelength_rejection_sampling(
#         nw,
#         nwout,
#         waves,
#         weights,
#         env,
#         outwaves
#         )
#     return 0
# # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#
# # CDF SAMPLING FUNCTION
# # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# # sampling CDF
# _lib_cdf.random_wave.argtypes = [
#     array_1d_double, # waves
#     array_1d_double, # weights
#     ct.c_int,    # nw
#     ct.c_int,    # interpolation flag
#     array_1d_double, # outwaves
#     ct.c_int     # outwaves size
#     ]
# _lib_cdf.random_wave.restype = None
#
# def random_wave_cdf(
#     waves,
#     weights,
#     nw,
#     flag,
#     outwaves,
#     nwout
#     ):
#     """c - wrapper"""
#     waves = np.ascontiguousarray(waves, dtype=np.float64)
#     weights = np.ascontiguousarray(weights, dtype=np.float64)
#     outwaves = np.ascontiguousarray(outwaves, dtype=np.float64)
#     _lib_cdf.random_wave(
#         waves,
#         weights,
#         nw,
#         flag,
#         outwaves,
#         nwout
#         )
#     return 0
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


if __name__ == "__main__":
    pass