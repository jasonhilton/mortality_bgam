#! /bin/bash

# hold array id in environment variable
# over-ride wall time in run_stan.pbs
arrayid=$(qsub PBS/run_stan.pbs -v CONFIG=$1,SEX=two_sex -l walltime=36:00:00)

# submitting stacking when all jobs in the array have completed
qsub -W depend=afterokarray:$arrayid PBS/run_stacking.pbs -v CONFIG=$1,SEX=two_sex
