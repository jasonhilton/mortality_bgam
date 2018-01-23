#! /bin/bash

# hold array id in environment variable
arrayid=$(qsub PBS/run_stan.pbs -v CONFIG=$1,SEX=$2)

# submitting stacking when all jobs in the array have completed
qsub -W depend=afterokarray:$arrayid PBS/run_stacking.pbs -v CONFIG=$1,SEX=$2
