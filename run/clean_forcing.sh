#!/bin/bash

# Script to remove time dependent forcing files from an experiment

np=186
NENS=30

wdir=`pwd`
exp='ERGOM-DA_free_new'
exp='tmp'

echo "Clean time-dependent forcing in ensemble directories in $wdir/$exp"

for((i=1;i<=$NENS;i++)); do
  ENSstr=`printf %03d $i`

  rm -f $wdir/$exp/$ENSstr/forcing/FORCE_y*
  rm -f $wdir/$exp/$ENSstr/forcing/ERGOM_SBC_y*
  rm -f $wdir/$exp/$ENSstr/forcing/ERGOM_OBC_y*
  rm -f $wdir/$exp/$ENSstr/forcing/ERGOM_CBC_y*
  rm -f $wdir/$exp/$ENSstr/forcing/EHYPE_y*
  rm -f $wdir/$exp/$ENSstr/forcing/bdy_uvh_y*
  rm -f $wdir/$exp/$ENSstr/forcing/bdy_ts_y*

done

