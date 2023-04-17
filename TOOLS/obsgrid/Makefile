.SUFFIXES: .F90 .f90 .c .o

# Makefile to compile program generating model fields
# on the observation grid (for chl...)

# For HLRN-Emmy I load the following modules to build this:
# openmpi/intel/4.1.4
# intel/2022.2
# hdf5-parallel/ompi/intel/1.12.1
# netcdf-parallel/ompi/intel.22/4.9.1

# Variable OBJECTS containing all object files
OBJECTS_CHL_Bal_obsgrid = timer.o write_chl_obsgrid_baltic.o
OBJECTS_CHL_No_obsgrid = timer.o write_chl_obsgrid_nsea.o
OBJECTS_SST_obsgrid = timer.o write_sst_obsgrid.o
OBJECTS_CHL_Bal_obsgrid_ens = timer.o write_chl_obsgrid_baltic_ens.o


# Name of the executable
EXE_CHL_Bal_obsgrid = write_chl_obsgrid_baltic
EXE_CHL_No_obsgrid = write_chl_obsgrid_nsea
EXE_SST_obsgrid = write_sst_obsgrid
EXE_CHL_Bal_obsgrid_ens = write_chl_obsgrid_baltic_ens

# Fortran compiler and compiler flags
FC      = ifort
FCFLAGS = -qmkl -O2 -check bounds

# Linker and linker flags
LINK    = ifort
LDFLAGS = $(FCFLAGS)

# Defines for preprocessor
CPP_DEFS = -cpp

# NetCDF library and include paths

#HLRN-EMMY
NCLIB = /sw/dataformats/netcdf-parallel/ompi/intel.22/4.9.1/skl -lnetcdff -lnetcdf -lm -ldl -L/sw/dataformats/hdf5-parallel/ompi/intel.22/1.12.1/skl/lib
NCINC = /sw/dataformats/netcdf-parallel/ompi/intel.22/4.9.1/skl/include

LINK_LIBS = #-L/usr/lib/lapack -llapack -L/usr/lib/libblas -lblas

#--------------------------------------------------

$(EXE_CHL_Bal_obsgrid): $(OBJECTS_CHL_Bal_obsgrid)
		$(LINK) -g -o $(EXE_CHL_Bal_obsgrid) \
		$(LDFLAGS) $(OBJECTS_CHL_Bal_obsgrid) $(MODULES) -L$(NCLIB) -lnetcdff $(LINK_LIBS)

$(EXE_CHL_No_obsgrid): $(OBJECTS_CHL_No_obsgrid)
		$(LINK) -g -o $(EXE_CHL_No_obsgrid) \
		$(LDFLAGS) $(OBJECTS_CHL_No_obsgrid) $(MODULES) -L$(NCLIB) -lnetcdff $(LINK_LIBS)

$(EXE_SST_obsgrid): $(OBJECTS_SST_obsgrid)
		$(LINK) -g -o $(EXE_SST_obsgrid) \
		$(LDFLAGS) $(OBJECTS_SST_obsgrid) $(MODULES) -L$(NCLIB) -lnetcdff $(LINK_LIBS)

$(EXE_CHL_Bal_obsgrid_ens): $(OBJECTS_CHL_Bal_obsgrid_ens)
		$(LINK) -g -o $(EXE_CHL_Bal_obsgrid_ens) \
		$(LDFLAGS) $(OBJECTS_CHL_Bal_obsgrid_ens) $(MODULES) -L$(NCLIB) -lnetcdff $(LINK_LIBS)

all: $(EXE_CHL_Bal_obsgrid) $(EXE_CHL_No_obsgrid) $(EXE_SST_obsgrid)

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
 
all: write_chl_obsgrid_baltic write_chl_obsgrid_baltic_ens

clean:
	rm -f *.o *.mod $(EXE_CHL_Bal_obsgrid) $(EXE_CHL_No_obsgrid) $(EXE_SST_obsgrid) $(EXE_CHL_Bal_obsgrid_ens)
