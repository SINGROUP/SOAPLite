#!/bin/bash

read -e -p "Select a .xyz File: " FILE
read -e -p "Select distance from atom type1: " dist1 
read -e -p "Select distance from atom type2(ignore if only one atom): " dist2 
read -e -p "Pic a buble size . \n This will determine the surface.  make sure to check if the surface was defined correctly
afterwards. \n If too large, it might ignore valleys, if too small, it will create a bubble inside and think that there
is a surface there.[2.0 - ]: " b

read -e -p "Select Atom Types.e.x. type \"Au Cu\" but without the quotations: " atoms

nAtoms=$(echo "$atoms" | wc -w)

echo "$nAtoms"

if [ $nAtoms -eq 2 ];then
time    ./Eta1 $FILE $atoms $dist1 $dist2 $b  > "$FILE".eta1H 

    echo "${FILE}.eta1H \nwas produced!"

 elif [ $nAtoms -eq 1 ]
    then
    ./Eta1One $FILE $dist1 $b  > "$FILE"One.eta1H 

    echo "${FILE}One.eta1H \nwas produced!"

 else
     echo "Error... Type1"
fi


