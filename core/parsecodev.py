"""parsecodev.py

Parses the Code V text files. Option to plot or write to disk.
"""

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
import glob
import logging
import os
import pprint
import re
import sys
import time
try:
    from defaults import logs_dir, codev_dir, codevparsed_path
except ImportError:
    sys.path.append(os.getcwd())
    from defaults import logs_dir, codev_dir, codevparsed_path

codev_fn = "criresplus_%s_v6_echelle_angle_%s_order_%s.txt"
write_dir = os.path.join(codev_dir, "parsed")
write_path = os.path.join(write_dir, codev_fn)
writenpy_path = os.path.splitext(write_path)[0] + ".npy"
bands = "Y J H K L M".split()
angles = ["{0}".format(x) for x in np.arange(60.0, 70.5, 0.5)]
headerline = "ORDER WAVELENGTH_NM XBOT_PIX XMID_PIX XTOP_PIX YBOT_PIX YMID_PIX YTOP_PIX SLITHEIGHT_PIX PHI_RAD".split()
headerfmt = "#{:<7} {:<15} {:<18} {:<18} {:<18} {:<18} {:<18} {:<18} {:<18} {:<18}\n"
datafmt = "{:<8} {:<15} {:<18} {:<18} {:<18} {:<18} {:<18} {:<18} {:<18} {:<18}\n"
log = logging.getLogger(__name__)


def linterp(X, X1, X2, Q11, Q21):
    """Linear interpolation."""
    return (X2-X)/(X2-X1)*Q11 + (X-X1)/(X2-X1)*Q21


class Mode():
    """Testing dummy class"""
    def __init__(self, sw=0.116, sh=5.8):
        self.SLIT_WIDTH = sw
        self.SLIT_HEIGHT = sh
        self.TAU_0 = 0.0
        self.DGAP = 0.0 # [mm]
        self.DPIX = 0.018 # [mm/pix]
        self.NDET = 3
        self.NXPIX = 2048
        self.NYPIX = 2048
        self.gap = self.DGAP / self.DPIX # ccd gap [pix]

def y2x(y, sw=0.2):
    """Ratio of x to y is 1 to 50 for 0.2 arcsecond slit,
    1 to 25 for 0.4 arcsecond slit."""
    if sw == 0.2:
        return 0.02 * y
    elif sw == 0.4:
        return 0.04 * y
    else:
        raise ValueError("Incorrect slit width.")

def grep_text(fn):
    with open(fn) as f:
        text = f.readlines()
    wlmatch = "WL ="
    xbmatch = "Xpos bottom of slit"
    xtmatch = "Xpos top of slit"
    ybmatch = "Ypos bottom of slit"
    ytmatch = "Ypos top of slit"
    wlind = []
    xbind = []
    xtind = []
    ybind = []
    ytind = []
    for i,l in enumerate(text):
        if re.search(wlmatch, l):
            wlind.append(i)
        if re.search(xbmatch, l):
            xbind.append(i+1)
        if re.search(xtmatch, l):
            xtind.append(i+1)
        if re.search(ybmatch, l):
            ybind.append(i+1)
        if re.search(ytmatch, l):
            ytind.append(i+1)
    wl = np.array([text[i].split("=")[1].strip(" ").strip("\n").strip("\r") for i in wlind], np.float)
    xb = np.array([text[i].split("=")[1].strip(" ").strip("\n").strip("\r") for i in xbind], np.float)
    xt = np.array([text[i].split("=")[1].strip(" ").strip("\n").strip("\r") for i in xtind], np.float)
    yb = np.array([text[i].split("=")[1].strip(" ").strip("\n").strip("\r") for i in ybind], np.float)
    yt = np.array([text[i].split("=")[1].strip(" ").strip("\n").strip("\r") for i in ytind], np.float)
    return wl, xb, xt, yb, yt

def load_data_tables(fn):
    log.info("Loading %s", fn)
    data = np.load(fn).T
    return data

def make_data_tables(band, angle, plot=False, write=False):
    log.info("Parsing data table: Band: %s, ech angle: %s", band, angle)
    log.debug("Plot tables  = %s", plot)
    log.debug("Write tables = %s", write)

    # find all orders
    fns = sorted(glob.glob(os.path.join(codev_dir,
        "original",
        codev_fn % (band, angle, "*"))))

    # create dummy instance
    mode = Mode()

    # PREPARE PLOT
    if plot:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_aspect("equal")
        for i in xrange(mode.NDET):
            i = i - mode.NDET/2
            xd = mode.NXPIX/2 * mode.DPIX
            yd = mode.NYPIX/2 * mode.DPIX
            xoff = i * (mode.gap * mode.DPIX + mode.NXPIX * mode.DPIX)
            ax.fill([-xd-xoff, xd-xoff, xd-xoff, -xd-xoff],
                    [-yd, -yd, yd, yd],
                    'r',
                    alpha=0.10,
                    edgecolor='k')

    # LOOP OVER ORDERS
    n = len(fns)
    for fn in fns:
        order = fn.split("_")[-1].replace(".txt", "")
        # READ DATA FROM CODE V OUTPUT
        log.debug("Parsing %s", fn)
        wl, yb, yt, xb, xt = grep_text(fn) # reverses x and y
        log.debug("Interpolating data...")
        xmid = (xb + xt) / 2.
        ymid = (yb + yt) / 2.
        slitheight = np.sqrt((yt-yb)**2 + (xt-xb)**2)
        slitwidth = slitheight * 0.2
        # slit image at origin
        xt0 = xt - xmid
        xb0 = xb - xmid
        yt0 = yt - ymid
        yb0 = yb - ymid
        # get tilt angle with slit centered at origin
        phi = np.pi/2.0 - np.arctan2(yt0-yb0, xt0-xb0)
        data = zip([int(order)]*wl.size, wl, xb, xmid, xt, yb, ymid, yt, slitheight, phi)

        # WRITE DATA TO FILE
        if write:
            if not os.path.isdir(write_dir):
                log.debug("Creating directory: %s", write_dir)
                os.makedirs(write_dir)
            wfn = write_path % (band, angle, order)
            wnpyfn = writenpy_path % (band, angle, order)
            with open(wfn, "w") as f1:
                log.debug("Writing data to %s", wfn)
                f1.write(headerfmt.format(*headerline))
                logging.debug(headerfmt.format(*headerline))
                for dat in data:
                    f1.write(datafmt.format(*dat))
                    logging.debug(datafmt.format(*dat))
            # binary data for loading faster
            log.debug("Writing data to %s", wnpyfn)
            np.save(wnpyfn, data)
        # PLOT
        if plot:
            log.info("Plotting...")
            px = zip(xb, xt)
            py = zip(yb, yt)
            for x,y in zip(px, py):
                plt.plot(x, y, c="b")

            # scale slit height
            pyt = yt0 * slitheight / mode.SLIT_HEIGHT
            pyb = yb0 * slitheight / mode.SLIT_HEIGHT

            nxt = (xt0) * np.cos(phi) - (pyt) * np.sin(phi) + xmid
            nyt = (xt0) * np.sin(phi) + (pyt) * np.cos(phi) + ymid
            nxb = (xb0) * np.cos(phi) - (pyb) * np.sin(phi) + xmid
            nyb = (xb0) * np.sin(phi) + (pyb) * np.cos(phi) + ymid
            plt.scatter(nxt, nyt, c='g')
            plt.scatter(nxb, nyb, c='r')

            for i,pt in enumerate(wl):
                if wl[i] == wl[0] or wl[i] == np.median(wl) or wl[i] == wl[-1]:
                    plt.annotate(str(wl[i])+" nm", xy=((xt[i]+xb[i])/2,(yt[i]+yb[i])/2))
    if plot:
        log.info("Opening plot...")
        plt.title("CRIRES+ %s band echangle = %s" % (band, angle))
        plt.xlabel("mm")
        plt.ylabel("mm")
        plt.show()
        plt.close()



def get_codev_files(self, m):
    ang_min = 60.0
    ang_max = 70.0
    ang_step = 0.5
    ang_range = np.arange(ang_min, ang_max+ang_step, ang_step)
    if self.echang not in ang_range:
        idx = np.searchsorted(ang_range, self.echang)
        ang_lower = ang_range[idx-1]
        ang_upper = ang_range[idx]
        fn1 = os.path.splitext(codevparsed_path % (self.band, ang_lower, m))[0] + ".npy"
        fn2 = os.path.splitext(codevparsed_path % (self.band, ang_upper, m))[0] + ".npy"
        log.debug("Loading '%s'...", fn1)
        log.debug("Loading '%s'...", fn2)
        table1 = np.load(fn1).T
        table2 = np.load(fn2).T
        out = np.empty(table1.shape)
        inds1 = np.argsort(table1[:, 1])
        inds2 = np.argsort(table2[:, 1])
        for i,(c1, c2) in enumerate(zip(table1, table2)):
            out[i] = linterp(self.echang, ang_lower, ang_upper, c1, c2)
        return out
    else:
        fn = os.path.splitext(codevparsed_path % (self.band, self.echang, m))[0] + ".npy"
        return np.load(fn).T


# =============================== MAIN =======================================

def main():
    # MAKESHIFT COMMAND LINE
    xb = bands
    xa = angles
    if "-p" in sys.argv[1:]:
        plot = True
    else:
        plot = False
    if "-w" in sys.argv[1:]:
        write = True
    else:
        write = False
    for arg in sys.argv[1:]:
        if "--bands" in arg:
            xb = [x for x in arg.split("=")[-1]] # not separated
        if "--angs" in arg:
            xa = [x for x in arg.split("=")[-1].split(",")] # separated by ,
    for band in xb:
        for angle in xa:
            make_data_tables(band, angle, plot=plot, write=write)

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
