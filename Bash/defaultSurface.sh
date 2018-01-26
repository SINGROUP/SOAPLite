#!/bin/bash

read -e -p "Select a .xyz File." FILE

read -e -p "Pic a grid density. Usually the smallest is accurate enough:[1-3] " grid

read -e -p "Select Atom Types: e.x. type \"Au Cu \" but without the quotations: " atoms

nAtoms=$(echo "$atoms" | wc -w)

echo "$nAtoms"

if [ $nAtoms -eq 2 ];then
    ./surfaceTwo $FILE $atoms 0.5 5 3 9 $grid 2.8  > "$FILE"_twoAt_05_5_3_9_"$grid"_28.surfsoap

    echo "${FILE}_twoAt_05_5_3_9_$grid.surfsoap \nwas produced!"

 elif [ $nAtoms -eq 1 ]
    then
    ./surfaceOne $FILE 0.5 5 3 9 $grid 2.8  > "$FILE"_oneAt_05_5_3_9_"$grid"_28.surfsoap

    echo "${FILE}_oneAt_05_5_3_9_$grid_28.surfsoap\nwas produced!"

 else
     echo "Error... Type1"
fi


