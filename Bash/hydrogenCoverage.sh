#!/bin/bash

read -e -p "Select a .soap File: " FILE

read -e -p "Select a .Hscan or .etaxH File. You need to do a hydrogen scan first: " FILE2

read -e -p "select soap threashold [d < 0.1]" d0

read -e -p " Number of Hydrogen Atoms [2,3]" nAtoms

echo "$nAtoms"

cat $FILE > bufferFile_DoNotTouch.xyz
 sed 's/^\s*[0-9].*$/H &/' $FILE2 >> bufferFile_DoNotTouch 

if [ $nAtoms -eq 2 ];then
    ./HydrogenCoverage2 $FILE $FILE2 $d0 > "$FILE"_"$FILE2"_2_"$d".Hcov

    if [ "$?" = -1 ] ; then
            zenity --error \
              --text="Update canceled."
    fi 
    echo "${FILE}_${FILE2}"_2_"${d}".Hcov \nwas produced!

elif [ $nAtoms -eq 3 ]
    then
    ./HydrogenCoverage3 $FILE $FILE2 $d0 > "$FILE"_"$FILE2"_3_"$d".Hcov

    echo "${FILE}_${FILE2}"_3_"${d}".Hcov \nwas produced!


 else
     echo "Error... Type1"
fi


