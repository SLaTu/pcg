#!/bin/bash
#SBATCH --job-name=cg_unstr
#SBATCH --workdir=.
#SBATCH --output=out_%j.out
#SBATCH --error=out_%j.err
#SBATCH --cpus-per-task=48
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --qos=debug

export OMP_NUM_THREADS=48

perf stat ./cg.out UNST398.txt 0 1 1000 1E-08 3 3 0 1
perf stat ./cg.out UNST398.txt 0 1 1000 1E-08 3 3 1 1
perf stat ./cg.out UNST398.txt 0 1 1000 1E-08 3 3 2 1
perf stat ./cg.out UNST398.txt 0 1 1000 1E-08 3 3 3 1
perf stat ./cg.out UNST398.txt 0 1 1000 1E-08 3 3 4 1
perf stat ./cg.out UNST398.txt 0 1 1000 1E-08 3 3 5 1
perf stat ./cg.out UNST398.txt 0 1 1000 1E-08 3 3 6 1



# perf stat ./cg.out UNST398.txt 0 1 1000 1E-08 3 3 0 2
# perf stat ./cg.out UNST398.txt 0 1 1000 1E-08 3 3 1 2
# perf stat ./cg.out UNST398.txt 0 1 1000 1E-08 3 3 2 2
# perf stat ./cg.out UNST398.txt 0 1 1000 1E-08 3 3 3 2
# perf stat ./cg.out UNST398.txt 0 1 1000 1E-08 3 3 4 2
# perf stat ./cg.out UNST398.txt 0 1 1000 1E-08 3 3 5 2
# perf stat ./cg.out UNST398.txt 0 1 1000 1E-08 3 3 6 2
# perf stat ./cg.out UNST398.txt 0 1 1000 1E-08 3 3 7 2
# perf stat ./cg.out UNST398.txt 0 1 1000 1E-08 3 3 8 2
# 
# perf stat ./cg.out UNST398.txt 0 1 1000 1E-08 3 3 0 2

