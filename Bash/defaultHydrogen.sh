#!/bin/bash

read -e -p "Select a .xyz File: " FILE
read -e -p "Select a .Hscan File. You need to do a hydrogen scan first: " FILE2
read -e -p "Pic a grid density. Usually the smallest is accurate enough:[1-6] " grid
read -e -p "Select Atom Types. e.x. type \"Au Cu\" without the quotations: " atoms

nAtoms=$(echo "$atoms" | wc -w)

echo "$nAtoms"

if [ $nAtoms -eq 2 ];then
    ./Hsoap $FILE $atoms $FILE2 0.5 5 3 9 $grid > "$FILE"_twoAt_05_5_3_9_"$grid".Hsoap


    echo "${FILE}_twoAt_05_5_3_9_$grid.Hsoap \nwas produced!"

 elif [ $nAtoms -eq 1 ]
    then
    ./HsoapOne $FILE $FILE2 0.5 5 3 9 $grid > "$FILE"_oneAt_05_5_3_9_"$grid".Hsoap

    echo "${FILE}_oneAt_05_5_3_9_$grid.Hsoap\nwas produced!"

 elif [ $nAtoms -eq 3 ]
    then
    ./HsoapThree $FILE $atoms $FILE2 0.5 5 3 9 $grid > "$FILE"_ThreeAt_05_5_3_9_"$grid".Hsoap

    echo "${FILE}_ThreeAt_05_5_3_9_$grid.Hsoap\nwas produced!"

 else
     echo "Error: Check Number of Atoms in Input."
fi


