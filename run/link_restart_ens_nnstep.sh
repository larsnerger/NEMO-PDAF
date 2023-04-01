#!/bin/bash

# Script to link restart files unsing the time index nnstep into the drietory intialstate 
# to be used in a model restart.
# The scrip tshould be executed in the main directory of an experiment

np_nemo=186
NENS=30

wdir=`pwd`
restart_out='output/restarts'

nnstep=00012480

np=$(( $np_nemo - 1 ))
for((i=1;i<=$NENS;i++)); do
  ENSstr=`printf %03d $i`

  rm -f $wdir/$ENSstr/initialstate/restart*.nc

  for n in `seq -f "%04g" 0 $np`;do
    ln -s $wdir/$ENSstr/$restart_out/'NORDIC_'$nnstep'_restart_out_'$n'.nc'     $wdir/$ENSstr/'initialstate/restart_in_'$n'.nc'     
    ln -s $wdir/$ENSstr/$restart_out/'NORDIC_'$nnstep'_restart_ice_out_'$n'.nc'     $wdir/$ENSstr/'initialstate/restart_ice_in_'$n'.nc'     
    ln -s $wdir/$ENSstr/$restart_out/'NORDIC_'$nnstep'_restart_trc_out_'$n'.nc'     $wdir/$ENSstr/'initialstate/restart_trc_in_'$n'.nc'     
   done
done
