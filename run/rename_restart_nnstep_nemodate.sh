#!/bin/bash

# Script to rename restart files holding the step index to restart files hold the NEMO date string.
# This is used when a job fails and one has to do the renaming afterwards. If a simulation job is
# successful this renaming is done in the postprocessing part of the job.

np_nemo=186
NENS=30

wdir=`pwd`
restart_out='output/restarts'

exp='DA-SST_noSalt'
nnstep=00013440
date_nemo='y2015m01d15'

echo "rename restart files from $nnstep to date $date_nemo".

np=$(( $np_nemo - 1 ))
for((i=1;i<=$NENS;i++)); do
  ENSstr=`printf %03d $i`

  echo Member $i

#  mkdir $savedir/$ENSstr

  for n in `seq -f "%04g" 0 $np`;do
    mv $wdir/$exp/$ENSstr/$restart_out/'NORDIC_'$nnstep'_restart_out_'$n'.nc' $wdir/$exp/$ENSstr/$restart_out/'restart_in_'$n'_'$date_nemo'.nc' 
    mv $wdir/$exp/$ENSstr/$restart_out/'NORDIC_'$nnstep'_restart_ice_out_'$n'.nc' $wdir/$exp/$ENSstr/$restart_out/'restart_ice_in_'$n'_'$date_nemo'.nc' 
    mv $wdir/$exp/$ENSstr/$restart_out/'NORDIC_'$nnstep'_restart_trc_out_'$n'.nc' $wdir/$exp/$ENSstr/$restart_out/'restart_trc_in_'$n'_'$date_nemo'.nc' 
   done
done
