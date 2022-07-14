# Generate covariance matrix

This program generates covariance matrix from a trajectory of model snapshots used to perform EOF decomposition by PDAF routine. The covariance matrix can then be used to generate initial ensemble perturbation for the model field. The program is specifically tailored for ocean model, NEMO.

Some details can be found at [PDAF documentatopm](https://pdaf.awi.de/trac/wiki/PDAF_eofcovar)

## Install
The `Makefile` is tailored for running on Archer2 program. 

1. Add the path to PDAF library for `PDAF_INC` and `PDAF_LIB` in `Makefile`.
2. `make clean && make`

On other platforms, it is necessary to change the path to dependent libraries.


## Input filenames
It supports two ways to provide filenames to the program. The filename can be given by directly specify `sfields(i)%file`. If it is preferred to provide a list of filenames, `sfields(i)%read_from_list` should be set to `.true.` and `sfields(i)%file` is then the filename of the list.

An example namelist for the read from list case is:
```
&state_vector
  sfields(1)%read_from_list = .true.
  sfields(1)%file='zos.dat'
/
```
The corresponding `zos.dat` file should be:
```
ORCA2_5d_00010101_00011231_grid_T.nc
ORCA2_5d_00020101_00021231_grid_T.nc
...
```

## Namelist:
The namelist is needed for constructing the state vector and PDAF_EOFcovar options.

- n_statevector:
  - nfields --- This specifies the number of state vectors
- state_vector: --- This namelist initialise an array of `type(state_field)` in `mod_statevector_pdaf.F90`
  - sfields(1)%variable =  'zos'
  - sfields(2)%variable =  'thetao'
  - sfields(3)%variable =  'so'
  - sfields(4)%variable =  'uo'
  - sfields(5)%variable =  'vo'

  - sfields(1)%read_from_list = .true.
  - sfields(1)%file='zos.dat'
  - ! sfields(2)%file='ORCA2_5d_00010101_00011231_grid_T.nc'
  - ! sfields(3)%file='ORCA2_5d_00010101_00011231_grid_T.nc'
  - ! sfields(4)%file='ORCA2_5d_00010101_00011231_grid_U.nc'
  - ! sfields(5)%file='ORCA2_5d_00010101_00011231_grid_V.nc'

- cov_options --- options for the resulting singular values and vectors
  - maxtimes = 73 ---maximum allowing time steps used to generate covariance matrix
  - do_mv = 1  --- Set to 1 to to perform multivariate normalization.
  - remove_mean = 0 --- Set to 1 to compute mean state and subtract it from states
  - limit = 1.0e-12 --- smallest singular values and its corresponding vectors in the covariance matrix
  - hwindow = 6 --- Half time window for running mean (running mean window of 12)