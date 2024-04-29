#!/bin/bash

top=${PWD}

trail=('t01' 't02' 't03' 't04')
for ri in ligand complex; do
    for s in concerted; do

        for t in ${trail[@]}; do
            cd ${ri}/${s}/${t}
            echo "${ri} ${s} ${t}"
            if [ ! -d results ]; then
                 mkdir results
            fi
            for l in $(cat mbar.dat); do
                cd ${l}
                python3 ${top}/scripts/stdti_step2dats.py ti002.prod.out
                cp *.dat ../results
                cd ..
            done
            cd ../../..
        done
    done
done


python3 scripts/MakeLaTeX_unified.py
cd latex
latex Main.tex
latex Main.tex

if [ -e Main.dvi ]; then
   dvips -j0 -Ppdf -G0 -tletter Main.dvi
fi
if [ -e Main.ps ]; then
  ps2pdf -dCompatibilityLevel=1.3 -dMAxSubsetPct=100 -dSubsetFonts=true -dEmbedAllFonts=true -sPAPERSIZE=letter -dEPSCrop Main.ps
fi

exit
