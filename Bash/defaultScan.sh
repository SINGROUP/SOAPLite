#!/bin/bash

read -e -p "Select a .xyz File: " FILE
echo $FILE
./scanH $FILE 150 0.5 2.8 1.2 2.2 > "$FILE"_150_05_28_12_22.Hscan 

    echo "${FILE}_150_05_27_12_22.Hscan was produced!"


