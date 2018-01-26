#!/bin/bash

read -e -p "Select a .xyz File: " FILE

read -e -p "Select a .Hscan File. You need to do a hydrogen scan first: " FILE2

read -e -p "Pic a bohr radius. This defines the shape of the basis functions. [0.01 - 1.0]" r0

read -e -p "Pic a Soap cuttof. MAX is 10.0, for large cutoff, make sure to use a dense grid.[1.0 - 10.0] " rcut

l=9
n=3

read -e -p "Pic a grid density. [1-7]" grid

read -e -p "Select Atom Types. e.x. type \"Au Cu \" without the quotations. " atoms

nAtoms=$(echo "$atoms" | wc -w)

echo "$nAtoms"

if [ $nAtoms -eq 2 ];then
    ./Hsoap $FILE $atoms $FILE2 ${r0} ${rcut} ${n} ${l} $grid > "$FILE"_twoAt_${r0}_${rcut}_${n}_${l}_"$grid".Hsoap


    if [ "$?" = -1 ] ; then
            zenity --error \
              --text="Update canceled."
    fi 

    echo "${FILE}_twoAt_${r0}_${rcut}_${n}_${l}_$grid.Hsoap \nwas produced!"

 elif [ $nAtoms -eq 1 ]
    then
    ./HsoapOne $FILE $FILE2 ${r0} ${rcut} ${n} ${l} $grid > "$FILE"_oneAt_${r0}_${rcut}_${n}_${l}_"$grid".Hsoap

    echo "${FILE}_oneAt_${r0}_${rcut}_${n}_${l}_$grid.Hsoap\nwas produced!"

 else
     echo "Error... Type1"
fi


