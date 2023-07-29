#!/bin/bash

for folder in highlighting_v2 production_v2; do
    cd $folder
    latexmk
    rm *.fdb_latexmk # work-around hack
    latexmk
    cd ../
done
