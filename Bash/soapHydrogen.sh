#!/bin/bash

read -e -p "Select a .xyz File: " FILE
read -e -p "Select a .Hscan or .etaHx File. You need to do a hydrogen scan first: " FILE2
read -e -p "Pic a grid density. Usually the smallest is accurate enough:[1-6] " grid
read -e -p "Select Atom Types. e.x. type \"Au Cu\" without the quotations: " atoms
read -e -p "Select the reach of the chemical environment.[1.0-10.0]" rcut

nAtoms=$(echo "$atoms" | wc -w)

echo "$nAtoms"

if [ $nAtoms -eq 2 ];then
    ./Hsoap $FILE $atoms $FILE2 0.5 $rcut 3 9 $grid > "$FILE"_twoAt_05_"$rcut"_3_9_"$grid".Hsoap


    echo "${FILE}_twoAt_05_"$rcut"_3_9_$grid.Hsoap \nwas produced!"

 elif [ $nAtoms -eq 1 ]
    then
    ./HsoapOne $FILE $FILE2 0.5 $rcut 3 9 $grid > "$FILE"_oneAt_05_"rcut"_3_9_"$grid".Hsoap

    echo "${FILE}_oneAt_05_"$rcut"_3_9_$grid.Hsoap\nwas produced!"

 else
     echo "Error... Type1"
fi


