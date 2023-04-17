.SUFFIXES: .F90 .f90 .c .o

# Makefile to compile programs computing the
# RMS deviation between observations and model
#
# For HLRN-Emmy I load the following modules to build this:
# openmpi/intel/4.1.4
# intel/2022.2
# hdf5-parallel/ompi/intel/1.12.1
# netcdf-parallel/ompi/intel.22/4.9.1

# Variable OBJECTS containing all object files
OBJECTS_SST = timer.o rmse_sst.o
OBJECTS_SST_BAL = timer.o rmse_sst_baltic.o
OBJECTS_CHL_BAL = timer.o rmse_chl_baltic.o
OBJECTS_CHL_NSEA = timer.o rmse_chl_nsea.o

# Name of the executable
EXE_SST     = rmse_sst
EXE_SST_BAL = rmse_sst_baltic
EXE_CHL_BAL = rmse_chl_baltic
EXE_CHL_NSEA = rmse_chl_nsea

# Fortran compiler and compiler flags
FC      = ifort
FCFLAGS = -qmkl -O2 -check bounds

# Linker and linker flags
LINK    = ifort
LDFLAGS = $(FCFLAGS)

# Defines for preprocessor
#CPP_DEFS = -DIBM
CPP_DEFS = -cpp

# NetCDF library and include paths

# HLRN-LISE (perhaps - need to check)
#NCLIB = /sw/dataformats/netcdf/intel.18/4.7.3/skl/lib -lnetcdff -fPIC -L/sw/dataformats/hdf5/intel.18/1.10.5/skl/lib -L/sw/dataformats/szip/intel.18/2.1.1/skl/lib -L/sw/dataformats/netcdf/intel.18/4.7.3/skl/lib -lnetcdf -lnetcdf -lm -ldl
#NCINC = /sw/dataformats/netcdf/intel.18/4.7.3/skl/include

#HLRN-EMMY
NCLIB = /sw/dataformats/netcdf-parallel/ompi/intel.22/4.9.1/skl -lnetcdff -lnetcdf -lm -ldl -L/sw/dataformats/hdf5-parallel/ompi/intel.22/1.12.1/skl/lib
NCINC = /sw/dataformats/netcdf-parallel/ompi/intel.22/4.9.1/skl/include

LINK_LIBS = #-L/usr/lib/lapack -llapack -L/usr/lib/libblas -lblas

#--------------------------------------------------
$(EXE_SST):	 	$(OBJECTS_SST)
		$(LINK) -g -o $(EXE_SST) \
		$(LDFLAGS) $(OBJECTS_SST) $(MODULES) -L$(NCLIB) -lnetcdff $(LINK_LIBS)

$(EXE_SST_BAL):	 	$(OBJECTS_SST_BAL)
		$(LINK) -g -o $(EXE_SST_BAL) \
		$(LDFLAGS) $(OBJECTS_SST_BAL) $(MODULES) -L$(NCLIB) -lnetcdff $(LINK_LIBS)

$(EXE_CHL_BAL):	 	$(OBJECTS_CHL_BAL)
		$(LINK) -g -o $(EXE_CHL_BAL) \
		$(LDFLAGS) $(OBJECTS_CHL_BAL) $(MODULES) -L$(NCLIB) -lnetcdff $(LINK_LIBS)

$(EXE_CHL_NSEA):	 	$(OBJECTS_CHL_NSEA)
		$(LINK) -g -o $(EXE_CHL_NSEA) \
		$(LDFLAGS) $(OBJECTS_CHL_NSEA) $(MODULES) -L$(NCLIB) -lnetcdff $(LINK_LIBS)

.F90.o:
		$(FC) $< \
		$(FCFLAGS) -c -I$(NCINC) -fpp

.f90.o:
		$(FC) -c $< \
		$(FCFLAGS)  

.c.o:
		cc -c -I$(LASPACKDIR) stf_solver.c

#--------------------------------------------------
new:		clean make

clean:
	rm -f *.o $(EXE) $(EXE_SST) $(EXE_SST_BAL) $(EXE_CHL_NSEA) $(EXE_CHL_BAL) *.mod
