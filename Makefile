#
# 'make'	build 'functions.so' file
# 'make clean'	removes all .so files
#

# compiler
CC = gcc
GCCVERSIONGTEQ4 := $(shell expr `gcc -dumpversion | cut -f1` \> 4.2.1)

# compiler flags
CFLAGS = -Wall -O3 -fPIC -g -pipe

# if GCC version > 4.2.1
ifeq "$(GCCVERSIONGTEQ4)" "1"
    CFLAGS +=-march=native
endif

# link libraries
LIBS = -lgsl -lgslcblas

# telluric spectra
TELLDIR = data/spectra/telluric
TELLSRC = $(TELLDIR)/LBLRTM_CH4_+0.0.fits.gz \
	$(TELLDIR)/LBLRTM_CO2_+0.0.fits.gz \
	$(TELLDIR)/LBLRTM_H2O_+0.0.fits.gz \
	$(TELLDIR)/LBLRTM_N2O_+0.0.fits.gz \
	$(TELLDIR)/LBLRTM_O2_+0.0.fits.gz \
	$(TELLDIR)/LBLRTM_O3_+0.0.fits.gz
TELLOBJS = $(TELLDIR)/LBLRTM_CH4_+0.0.npy \
	$(TELLDIR)/LBLRTM_CO2_+0.0.npy \
	$(TELLDIR)/LBLRTM_H2O_+0.0.npy \
	$(TELLDIR)/LBLRTM_N2O_+0.0.npy \
	$(TELLDIR)/LBLRTM_O2_+0.0.npy \
	$(TELLDIR)/LBLRTM_O3_+0.0.npy

# interpolation files
CODEV = data/codevtrace
CODEVSRC = $(CODEV)/original
CODEVOBJ = $(CODEV)/parsed

# c source files
COREDIR = core
SRC1 := $(COREDIR)/raytrace.c
OBJ1 := $(COREDIR)/raytrace.so
SRC2 := $(COREDIR)/slitfuncs.c
OBJ2 := $(COREDIR)/slitfuncs.so
SRC3 := $(COREDIR)/cdf.c
OBJ3 := $(COREDIR)/cdf.so

# object file
OBJS = $(OBJ1) $(OBJ2) $(OBJ3) $(TELLOBJS) $(CODEVOBJ)

.PHONY: all
all: $(OBJS)
	@echo "Done."

$(OBJ1): $(SRC1)
	$(CC) -shared -o $(OBJ1) $(CFLAGS) $(LIBS) $(SRC1)
	@echo  "Compilation of '$(SRC1)' to '$(OBJ1)' successful."

$(OBJ2): $(SRC2)
	$(CC) -shared -o $(OBJ2) $(CFLAGS) $(LIBS) $(SRC2)
	@echo  "Compilation of '$(SRC2)' to '$(OBJ2)' successful."

$(OBJ3): $(SRC3)
	$(CC) -shared -o $(OBJ3) $(CFLAGS) $(LIBS) $(SRC3)
	@echo  "Compilation of '$(SRC3)' to '$(OBJ3)' successful."

$(TELLOBJS): $(TELLSRC)
	python core/deploy.py telluric

$(CODEVOBJ): $(CODEVSRC)
	python core/parsecodev.py -w

clean:
	$(RM) -r *~ *.so *.so.dSYM core/*~ core/*.so core/*.so.dSYM $(CODEV)/parsed logs/* data/spectra/telluric/*.npy
	@find . -name ".*_*" -o -name "*.kate-swp" -o -name "*.swp" -o -name "*.pyc" -ls -exec rm -rf {} \;
	@echo  "Directory cleaned."

cleanlogs:
	$(RM) -r logs/*
	@echo  "logs directory cleaned."
