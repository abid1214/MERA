#!/bin/bash
source /opt/intel/bin/compilervars.sh intel64

d="../data"

L=16
l=3
N=10

Ld="${d}/L_${L}_l_${l}"
mkdir -p ${Ld}

for dseed in $(seq -f "%02g" 0 99)
do
    echo $dseed
    for W in 10 9 8 7 6 5 4 3 2 0.0001
    do
        Wd="${Ld}/W_${W}"
        mkdir -p ${Wd}
        for n in $(seq -f "%02g" 0 $N)
        do
            fn="${Wd}/epsilon_${n}_${dseed}.txt"
            OMP_NUM_THREADS=1 ./run $W $L $l ${dseed} ${n} $N > $fn
        done
    done
done
