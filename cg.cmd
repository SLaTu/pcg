#!/bin/bash
#SBATCH --job-name=test
#SBATCH --workdir=.
#SBATCH --output=out_%j.out
#SBATCH --error=out_%j.err
#SBATCH --cpus-per-task=48
#SBATCH --ntasks=1
#SBATCH --time=00:05:00
#SBATCH --qos=debug

export OMP_NUM_THREADS=48



./cg.out 120 0 1 10000 1E-08 1
./cg.out 120 0 1 10000 1E-08 2

./cg.out 120 0 1 10000 1E-08 3 3 0 1
./cg.out 120 0 1 10000 1E-08 3 3 50 1
./cg.out 120 0 1 10000 1E-08 3 3 100 1
./cg.out 120 0 1 10000 1E-08 3 3 200 1
./cg.out 120 0 1 10000 1E-08 3 3 400 1
