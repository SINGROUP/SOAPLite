#!/bin/bash

read -e -p "Select a .xyz File: " FILE

read -e -p "Pic a bohr radius. This defines the shape of the basis functions.[0.1 - 1.0] " r0

read -e -p "Pic a Soap cuttof. For large cutoff, make sure to use a dense grid.[1.0 - 10.0] " rcut

l=9

n=3

read -e -p "Pic a buble size . \n This will determine the surface.  make sure to check if the surface was defined correctly
afterwards. \n If too large, it might ignore valleys, if too small, it will create a bubble inside and think that there
is a surface there.[2.0 - ]: " b

read -e -p "Pic a grid density .Usually the smallest is accurate enough: [1-7] " grid

read -e -p "Select Atom Types.e.x. type \"Au Cu\" but without the quotations: " atoms

nAtoms=$(echo "$atoms" | wc -w)

echo "$nAtoms"

if [ $nAtoms -eq 2 ];then
    ./surfaceTwo $FILE $atoms $r0 $rcut $n $l $grid $b  > "$FILE"_twoAt_${r0}_${rcut}_${n}_${l}_"$grid"_$b.surfsoap 


    echo "${FILE}_twoAt_${r0}_${rcut}_${n}_${l}_${grid}_$b.surfsoap \nwas produced!"

 elif [ $nAtoms -eq 1 ]
    then
    ./surfaceOne $FILE $r0 $rcut $n $l $grid $b  > "$FILE"_oneAt_${r0}_${rcut}_${n}_${l}_"$grid"_$b.surfsoap 

    echo "${FILE}_oneAt_${r0}_${rcut}_${n}_${l}_${grid}_$b.surfsoap\nwas produced!"

 else
     echo "Error... Type1"
fi


