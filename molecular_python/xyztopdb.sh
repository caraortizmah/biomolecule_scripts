#!/bin/bash

python="/usr/bin/python"

for ii in `ls *.xyz`
do 
    suff=".xyz"
    name=${ii/%$suff}
    ${python} convert_xyztopdb.py ${ii} ${name}.pdb
done
