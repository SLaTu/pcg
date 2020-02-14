#!/bin/bash
#SBATCH --job-name=cg_unstr
#SBATCH --workdir=.
#SBATCH --output=out_%j.out
#SBATCH --error=out_%j.err
#SBATCH --cpus-per-task=48
#SBATCH --ntasks=1
#SBATCH --time=00:35:00
#SBATCH --qos=debug

export OMP_NUM_THREADS=48

perf stat ./cg.out UNST398.txt 0 50 10000 1E-10 3 3 0 1
perf stat ./cg.out UNST398.txt 0 50 10000 1E-10 3 3 50 1
perf stat ./cg.out UNST398.txt 0 50 10000 1E-10 3 3 100 1
perf stat ./cg.out UNST398.txt 0 50 10000 1E-10 3 3 200 1
perf stat ./cg.out UNST398.txt 0 50 10000 1E-10 3 3 0 2
perf stat ./cg.out UNST398.txt 0 50 10000 1E-10 3 3 50 2
perf stat ./cg.out UNST398.txt 0 50 10000 1E-10 3 3 100 2
perf stat ./cg.out UNST398.txt 0 50 10000 1E-10 3 3 200 2

