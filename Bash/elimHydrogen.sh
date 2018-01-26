#!/bin/bash

read -e -p "Select a .Hsoap File. Make sure that the file is in the right format (Hsoap format.)" FILE

read -e -p "Pic Soap threashold. [0.0 - Infinity ] " thresh


    ./elimH $FILE ${thresh}  > "$FILE"_${thresh}_"$grid".Huniq


    echo "${FILE}_${thresh}_$grid.Huniq  was produced. which  atoms are Unique."



