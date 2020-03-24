#!/bin/bash

l=16
w=3
seed=0
energy=8

d="../data/L_${l}/W_${w}"
mkdir -p ${d}
savefile="${d}/epsilon_${energy}_${seed}.txt"

qsub -N $l"_"$w"_"$seed"_"$energy -v W=$w,L=$l,SEED=$seed,ENERGY=$energy,NUMENERGIES=$l,SAVEFILE=$savefile mera_run.pbs
