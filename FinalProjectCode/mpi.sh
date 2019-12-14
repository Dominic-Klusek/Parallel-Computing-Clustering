#!/bin/bash

#BATCH --job-name slurm_mpi
#SBATCH --nodes 2
#SBATCH --ntasks 10
#SBATCH --ntasks-per-node 5
#SBATCH --mem 20000

mpirun -np 4 -npernode 2 ./test
