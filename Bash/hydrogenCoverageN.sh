#!/bin/bash

read -e -p "Select a .soap File: " FILE

read -e -p "Select a .Hscan or.etaxH File. You need to do a hydrogen scan first: " FILE2

read -e -p "select soap threashold [d < 0.1]" d0

read -e -p " Number of Hydrogen Atoms [n > 1]" nAtoms

echo "$nAtoms"

cat $FILE > bufferFile_DoNotTouch.xyz

 sed 's/^\s*[0-9].*$/H &/' $FILE2 >> bufferFile_DoNotTouch 

    ./HydrogenCoverageN $FILE $FILE2 $d0 $nAtoms > "${FILE}"_"${FILE2}"_"${d0}"_"${nAtoms}".Hcov

    echo "${FILE}_${FILE2}"_"${d0}"_"${nAtoms}".Hcov \nwas produced!


