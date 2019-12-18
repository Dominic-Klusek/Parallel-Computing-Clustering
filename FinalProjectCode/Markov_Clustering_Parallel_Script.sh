#!/bin/bash

#BATCH --job-name slurm_mpi
#SBATCH --nodes 7
#SBATCH --ntasks 49
#SBATCH --ntasks-per-node 7
#SBATCH --mem 20000
#SBATCH --partition partedu

cd $SLURM_SUBMIT_DIR
echo 0 > cap
mpirun -np 49 ./Markov
