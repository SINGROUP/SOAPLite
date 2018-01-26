#!/bin/bash

read -e -p "Select an xyz File: " FILE

read -e -p "Pic a grid density.Usually the smallest is accurate enough. [1-3]" grid

read -e -p "Select Atom Types.ex. type \"Au Cu\" without the quotation marks: " atoms

nAtoms=$(echo "$atoms" | wc -w)

echo "$nAtoms"

if [ $nAtoms -eq 2 ];then
    ./allTwo $FILE $atoms 0.5 5 3 9 $grid > "$FILE"_twoAt_05_5_3_9_"$grid".soapAll 

    echo "${FILE}_twoAt_05_5_3_9_$grid.soapAll \n was produced!"

 elif [ $nAtoms -eq 1 ]
    then
    ./allOne $FILE 0.5 5 3 9 $grid > "$FILE"_oneAt_05_5_3_9_"$grid".soapAll 

    echo "${FILE}_oneAt_05_5_3_3_$grid.soapAll \n was produced!"

 else
     echo "Error... Type1"
fi

