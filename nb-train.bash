#! /bin/bash

NTHREADS=15
EXEC_ROOT=/mnt/HA/groups/rosenGrp/zz374/INBC/code_original/train/NB/NB.run
GEN_SRC=$1
FOLD=$2
N_FOLDS=$3
KMERSIZE=$4
MEMLIM=48000

for ((N=1; N<=$N_FOLDS; N++)); do
    if [ $N -ne $FOLD ] && [ -d $GEN_SRC/fold$N ]
    then
        echo "Training on fold $N; doing all except $FOLD"
        $EXEC_ROOT train $GEN_SRC/fold$N -s $GEN_SRC/save_$FOLD -t $NTHREADS -k $KMERSIZE -m $MEMLIM
    fi
done
