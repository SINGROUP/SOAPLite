#!/bin/bash

read -e -p "Select an xyz File: " FILE
read -e -p "Pic a bohr radius. This defines the shape of the basis functions. [0.1-1.0]: " r0
read -e -p "Pic a Soap cuttof in Angs.[1.0-10.0]: " rcut
l=9
n=3
read -e -p "Pic the grid density. Usually the smallest is accurate enough. [1-6]: " grid
read -e -p "Select Atom Types:ex. type \"Au Cu\" without the double quotations: " atoms

nAtoms=$(echo "$atoms" | wc -w)

echo "$nAtoms"

if [ $nAtoms -eq 2 ];then
    ./allTwo $FILE $atoms $r0 $rcut $n $l $grid > "$FILE"_twoAt_05_${rcut}_${n}_${l}_"$grid".soapAll 

    echo "${FILE}_twoAt_${r0}_${rcut}_${n}_${l}_$grid.soapAll  was produced!"

 elif [ $nAtoms -eq 1 ]
    then
    ./allOne $FILE $r0 $rcut $n $l $grid > "$FILE"_oneAt_${r0}_${rcut}_${n}_${l}_"$grid".soapAll 

    echo "${FILE}_oneAt_${r0}_${rcut}_${n}_${l}_$grid.soapAll  was produced!"

 else
     echo "Error... Type1"
fi

