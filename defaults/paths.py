import os.path

# DIRECTORIES
core_dir = "core"
data_dir = "data"
defaults_dir = "defaults"
output_dir = "output"
logs_dir = "logs"
codev_dir = os.path.join(data_dir, "codevtrace")
codevdata_dir = os.path.join(codev_dir, "original")
codevparsed_dir = os.path.join(codev_dir, "parsed")
spectra_dir = os.path.join(data_dir, "spectra")
phx_dir = os.path.join(spectra_dir, "phoenix")
tell_dir = os.path.join(spectra_dir, "telluric")

# FILENAMES
codev_fn = "criresplus_%s_v6_echelle_angle_%s_order_%s.txt"
phx_waves_fn = "WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"
phx_flux_fn = "lte03000-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"
output_fn = "crifors.%s"


# FULL PATHS
codevdata_path = os.path.join(codev_dir, "original", codev_fn)
codevparsed_path = os.path.join(codev_dir, "parsed", codev_fn)
codevparsednpy_path = os.path.splitext(codevparsed_path)[0] + ".npy"
phx_waves_path = os.path.join(phx_dir, phx_waves_fn)
phx_flux_path = os.path.join(phx_dir, phx_flux_fn)


# SETTINGS
dsettings_fn = "settings.cfg"
dsettings = os.path.join(defaults_dir, dsettings_fn)
tell_species = "CH4 CO2 H2O N2O O2 O3".split()


# INST PARAMS
dinst_fn = "inst.cfg"
dinst = os.path.join(defaults_dir, dinst_fn)
