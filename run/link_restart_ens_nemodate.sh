#!/bin/bash

# Script to link restart files from an archive directory (restartdir) into the directory
# initialstate/ of an experiment.

np_nemo=1 #186
NENS=2 #30

wdir=`pwd`
restartdir='/scratch/projects/hbk00095/restarts/free_N30'
exp='DA-SST_noSalt'

date_nemo='y2015m02d01'

echo "Link ensemble restart files from $restartdir to $exp"

np=$(( $np_nemo - 1 ))
for((i=1;i<=$NENS;i++)); do
  ENSstr=`printf %03d $i`
  wdir=`pwd`/$exp/${ENSstr}

  echo "- process member $i"

#  rm -f $wdir/'initialstate/'restart*.nc

  for n in `seq -f "%04g" 0 $np`;do
    ln -s $restartdir/$ENSstr'/restart_in_'$n'_'$date_nemo'.nc'     $wdir/'initialstate/restart_in_'$n'.nc'
    ln -s $restartdir/$ENSstr'/restart_ice_in_'$n'_'$date_nemo'.nc'     $wdir/'initialstate/restart_ice_in_'$n'.nc'
    ln -s $restartdir/$ENSstr'/restart_trc_in_'$n'_'$date_nemo'.nc'     $wdir/'initialstate/restart_trc_in_'$n'.nc'
  done
done
