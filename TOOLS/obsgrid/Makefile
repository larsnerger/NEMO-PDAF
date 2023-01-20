.SUFFIXES: .F90 .f90 .c .o

# Makefile to compile program generating model fields
# on the observation grid (for chl...)

# For HLRN I load the following modules to build this:
# intel/18.0.5
# netcdf/intel/4.7.3  
# impi/2018.4


# Variable OBJECTS containing all object files
OBJECTS_CHL_Bal_obsgrid = timer.o write_chl_obsgrid_baltic.o
OBJECTS_CHL_Bal_obsgrid_ens = timer.o write_chl_obsgrid_baltic_ens.o


# Name of the executable
EXE_CHL_Bal_obsgrid = write_chl_obsgrid_baltic
EXE_CHL_Bal_obsgrid_ens = write_chl_obsgrid_baltic_ens

# Fortran compiler and compiler flags
FC      = ifort
#FCFLAGS = -O3 -qsuffix=f=f90 -qarch=pwr4 -qtune=auto -q64 -qrealsize=8
FCFLAGS = -mkl -O2 -check bounds
#FCFLAGS = -qextname -O3 -qstrict -qsuffix=f=f90 -qarch=auto -qtune=auto -q64 -qrealsize=8

#FCFLAGS = -fast -xtypemap=real:64  -xtarget=ultra3cu -xarch=v9b 

# Linker and linker flags
LINK    = ifort
LDFLAGS = $(FCFLAGS)

# Defines for preprocessor
#CPP_DEFS = -DIBM
CPP_DEFS = -cpp

# NetCDF library and include paths

# HLRN
NCLIB = /sw/dataformats/netcdf/intel.18/4.7.3/skl/lib -lnetcdff -fPIC -L/sw/dataformats/hdf5/intel.18/1.10.5/skl/lib -L/sw/dataformats/szip/intel.18/2.1.1/skl/lib -L/sw/dataformats/netcdf/intel.18/4.7.3/skl/lib -lnetcdf -lnetcdf -lm -ldl
NCINC = /sw/dataformats/netcdf/intel.18/4.7.3/skl/include

# OLLIE - need to 'module load netcdf/4.6.2_intel'
#NCLIB = -L/global/AWIsoft/netcdf/4.6.2_intel_18.0.5/lib -lnetcdff -L/global/AWIsoft/hdf5/1.10.2_gnu/lib -Wl,-rpath=/global/AWIsoft/hdf5/1.10.2_gnu/lib -L/global/AWIsoft/netcdf/4.6.2_intel_18.0.5/lib -Wl,-rpath=/global/AWIsoft/netcdf/4.6.2_intel_18.0.5/lib  -lnetcdf -lnetcdf -lhdf5_hl -lhdf5 -ldl -lsz -lz -lcurl
#NCINC = /global/AWIsoft/netcdf/4.6.2_intel_18.0.5/include

LINK_LIBS = #-L/usr/lib/lapack -llapack -L/usr/lib/libblas -lblas
#LINK_LIBS = -lessl -L/iblade/soft/lapack/3.2.1 -llapack_pwr6

#--------------------------------------------------

write_chl_obsgrid_baltic:	 	$(OBJECTS_CHL_Bal_obsgrid)
		$(LINK) -g -o $(EXE_CHL_Bal_obsgrid) \
		$(LDFLAGS) $(OBJECTS_CHL_Bal_obsgrid) $(MODULES) -L$(NCLIB) -lnetcdff $(LINK_LIBS)

write_chl_obsgrid_baltic_ens:	 	$(OBJECTS_CHL_Bal_obsgrid_ens)
		$(LINK) -g -o $(EXE_CHL_Bal_obsgrid_ens) \
		$(LDFLAGS) $(OBJECTS_CHL_Bal_obsgrid_ens) $(MODULES) -L$(NCLIB) -lnetcdff $(LINK_LIBS)

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
	rm -f *.o $(EXE) $(EXE_CHL_Bal) *.mod
