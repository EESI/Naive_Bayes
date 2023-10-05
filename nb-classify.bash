#! /bin/bash

NTHREADS=15
EXEC_ROOT=/mnt/HA/groups/rosenGrp/zz374/INBC/code_original/train/NB/NB.run
GEN_SRC=$1
GEN_TEST=$2
FOLD=$3
KMERSIZE=$4
MEMLIM=48000

echo "Classifying on fold $FOLD."
$EXEC_ROOT classify $GEN_TEST/fold$FOLD -s $GEN_SRC/save_$FOLD -t $NTHREADS -k $KMERSIZE -m $MEMLIM
