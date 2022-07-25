# $Id: Makefile 746 2009-08-04 12:16:28Z lnerger $

#######################################################
# Generic Makefile for to build PDAF with dummy model #
# To choose the architecture set $PDAF_ARCH           #
#######################################################

######################################################

# User specifications
# 1. Set BASEDIR, the directory where the PDAF package resides
# 2. Set PDAF_ARCH to include compile definitions
#    (See directory BASEDIR/make.arch for files. PDAF_ARCH is filename without .h)

# # Root directory of PDAF package
BASEDIR = ../../PDAF_V2.0
# DEF_NCPAR = -DUSE_NCPAR -Dkey_PDAF_offline

# # Include machine-specific definitions
# # For available include files see directory make.arch
# # To choose a file, set PDAF_ARCH either here or by an
# # environment variable.
# include $(BASEDIR)/make.arch/$(PDAF_ARCH).h
# Include machine-specific definitions
# For available include files see directory make.arch
FC = mpif90
CC = mpif90
LD = mpif90

PDAF_INC = -I$(BASEDIR)/include
NETCDFF_INC = 
CPP_DEFS = -Dkey_PDAF_offline
MPI_INC =

INC= -I. $(MPI_INC) $(PDAF_INC) $(NETCDFF_INC) $(CPP_DEFS)

OPT  = -O3 -cpp -Wall -Wextra --std=f2008


NETCDFF_LIB = 
PDAF_LIB = -L$(BASEDIR)/lib -lpdaf-d
LAPACK_LIB = -llapack  -lblas -lm

LDLIBS = $(NETCDFF_LIB) $(PDAF_LIB) $(LAPACK_LIB)

# Name of executable
EXE = PDAF_offline

# End of user specifications
######################################################

.SUFFIXES: .F90 .o 

# Modules used for the model part
MODULES =  	mod_kind_pdaf.o\
		mod_parallel_pdaf.o \
		timer.o \
		mod_memcount_pdaf.o \
		mod_mppini_pdaf_offline.o

MOD_ERGOM = mod_aux_pdaf.o \
	    	mod_io_pdaf.o \
	    	diag_histogram.o \
	    	diag_ensstats.o

# Module required for assimilation
MOD_ASSIM = mod_assimilation_pdaf.o \
		mod_nemo_pdaf.o \
		mod_statevector_pdaf.o

# Model routines
OBJ_MODEL =	initialize.o

# Routines of observation handling (PDAF-OMI)
OBJ_USER_PDAFOMI = mod_obs_ssh_mgrid_pdafomi.o \
		callback_obs_pdafomi.o
# 		obs_prof_pdafomi.o \
		obs_sst_cmems_pdafomi.o \

# Interface to PDAF - model sided
OBJ_PDAF_INT =  init_parallel_pdaf.o \
		mod_util_pdaf.o \
		mod_init_pdaf.o \
		mod_diagnostics_pdaf.o

# Generic user-supplied routines
OBJ_USER_GEN  = collect_state_pdaf.o \
		init_ens_pdaf.o \
		prepoststep_ens_pdaf.o \


# User-supplied routines for localized analysis (LSEIK)
OBJ_USER_LOCAL = init_n_domains_pdaf.o \
		init_dim_l_pdaf.o \
		g2l_state_pdaf.o \
		l2g_state_pdaf.o

# Full list of user-supplied routines for online modes
OBJ_PDAF_USER = $(OBJ_USER_PDAFOMI) $(OBJ_USER_LOCAL)

######################################################

$(EXE) : $(MODULES) $(MOD_ASSIM) $(MOD_ERGOM) $(OBJ_MODEL) $(OBJ_PDAF_USER) \
	$(OBJ_PDAF_INT)  $(OBJ_USER_GEN) main_offline.o
	$(LD)  $(OPT_LNK) $(NC_LIB) -o $@ $^ $(LDLIBS)
	@echo "++++++ Done ++++++"

######################################################

libpdaf-d.a: 
	@echo "++++++ Generate Filter library ++++++"
	@cd $(BASEDIR)/src; make;

######################################################

.F90.o :
	$(FC) $(OPT) $(CPP_DEFS) $(INC) -I$(BASEDIR)/include -c $*.F90

.f90.o :
	$(FC) $(OPT) $(CPP_DEFS) $(INC) -c $*.f90

# For older compilers one might need to separate the
# preprocessing from the compilation as defined below:
#.F90.o :
#	$(CPP) -P -C $(MPI_INC) $(CPP_DEFS) $*.F90 $*.f90
#	$(FC) $(OPT) $(MPI_INC) -c $*.F90
#	@rm -f  $*.f90

cleanall : cleanpdaf clean

clean :
	rm -f *.o *.mod *.a $(EXE)

cleanpdaf:
	@echo "+++ Clean up PDAF directory"
	cd $(BASEDIR)/src; make clean
