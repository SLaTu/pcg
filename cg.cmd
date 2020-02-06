#!/bin/bash
#SBATCH --job-name=cg_test
#SBATCH --workdir=.
#SBATCH --output=out_%j.out
#SBATCH --error=out_%j.err
#SBATCH --cpus-per-task=48
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --qos=debug
 
module load python
module load gcc
module load papi
module load mkl


make

export OMP_NUM_THREADS=48

# 0.091211, 0.118150, 1.000000, 0.090771, 0.095210, 0.090771, 0.239770, 0.738339, 0.160703, 0.129669,
# ./cg.out ALYA1.txt 1 100 10000 1E-08 2


./cg.out ALYA1.txt 1 100 10000 1E-08 3 3 0 1
./cg.out ALYA1.txt 1 100 10000 1E-08 3 3 0 2

./cg.out ALYA1.txt 1 100 10000 1E-08 3 3 100 1


