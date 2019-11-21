#!/bin/bash
source /opt/intel/bin/compilervars.sh intel64

L=64
dseed=31

W=${1}
#for W in `seq 1 10`
#do
d="data/last_few_bits_L_${L}_energy_avg/W_${W}"
mkdir -p ${d}
echo $d

echo "W=${W}"
for seed in $(seq -f "%02g" 0 99)
do
    echo "   ds=${dseed}   s=${seed}"
    OMP_NUM_THREADS=1 ./run $W $L ${dseed} ${seed} > $d/W_${W}_${seed}.txt
done
#done
