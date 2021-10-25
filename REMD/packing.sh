#!/bin/bash

#This script moves each nvt_$i.tpr into a folder having its corresponding $i as its
# folder name, being $i the replica number.
#In each new folder the nvt_$i.tpr file is renamed as nvt_.tpr 

for ii in `ls *.tpr`
do
  num="$(echo "$ii" | cut -d'.' -f1 | cut -d'_' -f2)"
  base="$(echo "$ii" | cut -d'_' -f1)"
  mkdir -p ${num}
  mv ${ii} ${num}/${base}_.tpr
done

exit
