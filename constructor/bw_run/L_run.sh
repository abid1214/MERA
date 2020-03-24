#!/bin/bash

l=32
N=9
for seed in $(seq -f "%02g" 0 99)
do
    for energy in $(seq -f "%02g" 0 $N)
    do
        echo "seed = ${seed}, energy = ${energy}"
        qsub -N $l"_"$seed"_"$energy -v L=$l,SEED=$seed,ENERGY=$energy,NUMENERGIES=$N mera_run.pbs

    done
done
