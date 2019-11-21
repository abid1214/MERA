#!/bin/bash
source /opt/intel/bin/compilervars.sh intel64
./run $1 > W_${1}.txt
