#!/bin/bash

read -e -p "Select a list of atoms such as a .surfUniq file." FILE

read -e -p "Select distance from surf atom. [1.0 - 5.0]" dist
read -e -p "Select the .xyz file corresponding to the .surfUniq file: " FILE2

./putH $FILE $FILE2 $dist > ${FILE}_$dist.Hput 


echo "${FILE}_$dist.Hput was produced! Check it by catting it onto the .xyz file."
