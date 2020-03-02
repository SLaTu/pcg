#!/bin/bash
#SBATCH --job-name=cg_unstr
#SBATCH --workdir=.
#SBATCH --output=out_%j.out
#SBATCH --error=out_%j.err
#SBATCH --cpus-per-task=48
#SBATCH --ntasks=1
#SBATCH --time=00:01:00
#SBATCH --qos=debug

export OMP_NUM_THREADS=48


perf stat ./cg.out test_3.txt 0 1 1000 1E-08 4
perf stat ./cg.out smooth.txt 0 1 1000 1E-08 4

