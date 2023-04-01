#!/bin/bash

# Sript to store restart files from an experiment directory into an archive directory

np_nemo=186
NENS=30

exp='ERGOM-DA_free'
savedir='/scratch/projects/hbk00095/restarts/free_N30'

timestr='y2015m02d01'

#-------------------

echo "Save restart files for date $timestr from $wdir/$exp to $savedir".
wdir=`pwd`

restart_out='output/restarts'

np=$(( $np_nemo - 1 ))

for((i=1;i<=$NENS;i++)); do
  ENSstr=`printf %03d $i`

  echo Member $i

  mkdir $savedir/$ENSstr

  for n in `seq -f "%04g" 0 $np`;do
    cp -p $wdir/$exp/$ENSstr/$restart_out/restart_in*$n'_'$timestr'.nc' $savedir/$ENSstr/
    cp -p $wdir/$exp/$ENSstr/$restart_out/restart_ice_in*$n'_'$timestr'.nc' $savedir/$ENSstr/
    cp -p $wdir/$exp/$ENSstr/$restart_out/restart_trc_in*$n'_'$timestr'.nc' $savedir/$ENSstr/
   done
done
