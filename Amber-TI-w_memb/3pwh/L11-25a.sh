#!/bin/bash

morph_dir=(L11-25a)
systems=(complex ligand)
stages=(concerted)
trials=('t01' 't02' 't03' 't04')


for mdir in ${morph_dir[@]};do
    for system in ${systems[@]};do
        for stage in ${stages[@]};do
            for trial in ${trials[@]};do
                cd ${mdir}/${system}/${stage}/${trial}
                ./2_run.sh ${step}
                cd ../../../../
            done
        done
    done
done
