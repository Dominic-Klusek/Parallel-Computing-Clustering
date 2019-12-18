#!/bin/bash

#BATCH --job-name slurm_mpi
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 20000
#SBATCH --partition partedu

cd $SLURM_SUBMIT_DIR
echo 0 > cap
mpirun -np 1 ./Markov_Seq
