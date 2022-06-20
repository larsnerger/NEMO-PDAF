#######################################################
# Generic Makefile for to build #
#######################################################

.SUFFIXES: .F90 .f .o 

######################################################

# Include machine-specific definitions
# For available include files see directory make.arch
FC = ftn
CC = CC

PDAF_INC =
NETCDFF_INC = -I$(NETCDF_DIR)/include -I$(HDF5_DIR)/include
CPP_DEFS =
MPI_INC =

INC=$(MPI_INC) $(PDAF_INC) $(NETCDFF_INC) $(CPP_DEFS)

FCFLAGS = -O2 -cpp -fallow-argument-mismatch -fdefault-real-8 -fcray-pointer -ffree-line-length-none
FFLAGS  = -O2 -cpp -fallow-argument-mismatch -fdefault-real-8 -fcray-pointer -ffree-line-length-none


NETCDFF_LIB = -L$(HDF5_DIR)/lib -L$(NETCDF_DIR)/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz
PDAF_LIB =
LAPACK_LIB =

LDLIBS = -lstdc++ -lmpichf90 $(NETCDFF_LIB) $(PDAF_LIB) $(LAPACK_LIB)


######################################################
# Define objects
######################################################

# Generic routines in PDAF
OBJ_F = mod_kind_pdaf.o mod_io_pdaf.o mod_statevector_pdaf.o gen_cov.o


######################################################

generate_covar: $(OBJ_F)
	@echo "++++++ Generate Filter library ++++++"
	$(FC) -o $@ $^ $(LDLIBS)


.F90.o :
	$(FC) $(FFLAGS) $(OPT) $(INC) $(CPP_DEFS) -c $*.F90

# Cleans
clean :
	rm -f *.o *.mod generate_covar
