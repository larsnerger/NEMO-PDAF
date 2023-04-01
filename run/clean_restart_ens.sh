#!/bin/bash

# Script to delete unnecessary restart files in an experiment directory.
# There are tow modes:
# rsttype=0
#    This removes all restart files NORDIC_*.nc
# rsttype=1
#    This removes the restart files restart_*_timestr


np_nemo=186
NENS=30

#  (0) for NORDIC*_0???.nc, (1) for using timestr
rsttype=0

# Name of experiment
exp='dev3'

# Timestring
timestr='y2015m01d18'

#------------------

wdir=`pwd`
restart_out='output/restarts'

if [ $rsttype -eq 1 ]; then
  echo "Clean restart files for $timestr from $wdir/$exp"
else
  echo "Clean NORDIC files from $wdir/$exp"
fi

np=$(( $np_nemo - 1 ))

for((i=1;i<=$NENS;i++)); do
  ENSstr=`printf %03d $i`

  echo clean directory $wdir/$exp/$ENSstr/$restart_out/

  if [ $rsttype -eq 1 ]; then
    rm $wdir/$exp/$ENSstr/$restart_out/restart*_$timestr.nc
  else
    rm $wdir/$exp/$ENSstr/$restart_out/NORDIC*_????.nc
  fi
done
