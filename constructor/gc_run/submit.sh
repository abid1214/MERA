#!/bin/bash

W=$1
L=$2
l=$3
seed=$4
n=$5
N=$6

fname="W_${W}_L_${L}_l_${l}.txt"

./run $W $L $l $seed $n $N > $fname
