#!/bin/bash

d="../data"

L=$1
l=$2
N=10

input="rerun_L_${L}_l_${l}.txt"

qsub -N "L_${L}_l_${l}" -v d=$d,input=$input,L=$L,l=$l,N=$N resume.pbs
