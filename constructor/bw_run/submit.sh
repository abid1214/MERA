#!/bin/bash
source /opt/intel/bin/compilervars.sh intel64

d="../data"

L=64

Ld="${d}/L_${L}"
mkdir -p ${Ld}

W=${1}
for dseed in $(seq -f "%02g" 0 99)
do
    Wd="${Ld}/W_${W}"
    mkdir -p ${Wd}
    for n in $(seq -f "%02g" 0 $L)
    do
        OMP_NUM_THREADS=1 ./run $W $L ${dseed} ${n} $L > "${Wd}/epsilon_${n}_${dseed}.txt"
    done
done
