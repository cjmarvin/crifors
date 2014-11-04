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


def interp(*args):
    func = raytrace.raytrace_interp_bin
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
        array_1d_double,     # cwl
        array_1d_double,     # cxb
        array_1d_double,     # cxm
        array_1d_double,     # cxt
        array_1d_double,     # cyb
        array_1d_double,     # cym
        array_1d_double,     # cyt
        array_1d_double,     # cphi
        array_1d_double,     # waves
        array_1d_double,     # slit_x
        array_1d_double,     # slit_y
        array_2d_uint]       # outarr
    func.restype = None
    func(*args)

def solve(*args):
    func = raytrace.raytrace_solve_general
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
        array_1d_double,     # slit_x
        array_1d_double,     # slit_y
        array_1d_double,     # waves
        array_1d_double,     # returnx
        array_1d_double,     # returny
        array_2d_double,     # returnwaves
        array_2d_uint]       # returncounts
    func.restype = None
    func(*args)