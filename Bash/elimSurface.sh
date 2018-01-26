#!/bin/bash

read -e -p "Select a .surfsoap File. Make sure that the file is in the right format (surfsoap format): " FILE

read -e -p "Pic Soap threashold. [0.0 - Infinity]: " thresh



    ./elimSurf $FILE ${thresh} > "$FILE"_${thresh}_"$grid".surfUniq


read -p "Select the .xyz file corresponding to the .surfsoap file." FILE2
cp $FILE2 "check_$FILE2"

    while read p ; do sed -i "$(($p+2))s/^\s*[A-Za-z]./Lr/" "check_$FILE2" ; done <"$FILE"_${thresh}_"$grid".surfUniq

cp $FILE2.check check.xyz

    echo "${FILE}_${thresh}_$grid.surfUniq\nand \n check_$FILE2 \n was produced. Please check check.xyz file to see which surface atoms are Unique."
