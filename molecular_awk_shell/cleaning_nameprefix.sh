#!/bin/bash

# Rename a file by removing the initial common pattern on the
#  name file

for ii in `ls *.xyz`
do

    # "cutting" the pattern of the name of the file
    pref="format_"
    name=${ii/#$pref}
    
    # rename file
    mv $ii $name.xyz
done
