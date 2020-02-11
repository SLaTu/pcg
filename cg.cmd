#!/bin/bash
#SBATCH --job-name=cg_tols
#SBATCH --workdir=.
#SBATCH --output=out_%j.out
#SBATCH --error=out_%j.err
#SBATCH --cpus-per-task=48
#SBATCH --ntasks=1
#SBATCH --time=05:00:00

module load python
module load intel
# module load papi
module load mkl


make

export OMP_NUM_THREADS=48


perf stat ./cg.out M_500.txt 1 50 10000 1E-8 3 3 0 1
perf stat ./cg.out M_500.txt 1 50 10000 1E-8 3 3 50 1
perf stat ./cg.out M_500.txt 1 50 10000 1E-8 3 3 100 1
perf stat ./cg.out M_500.txt 1 50 10000 1E-8 3 3 200 1
perf stat ./cg.out M_500.txt 1 50 10000 1E-8 3 3 300 1


perf stat ./cg.out M_500.txt 1 50 10000 1E-8 3 3 0 2
perf stat ./cg.out M_500.txt 1 50 10000 1E-8 3 3 50 2
perf stat ./cg.out M_500.txt 1 50 10000 1E-8 3 3 100 2
perf stat ./cg.out M_500.txt 1 50 10000 1E-8 3 3 200 2
perf stat ./cg.out M_500.txt 1 50 10000 1E-8 3 3 300 2


perf stat ./cg.out M_500.txt 1 50 10000 1E-8 3 3 0 3
perf stat ./cg.out M_500.txt 1 50 10000 1E-8 3 3 50 3



perf stat ./cg.out M_500.txt 1 50 10000 1E-10 3 3 0 1
perf stat ./cg.out M_500.txt 1 50 10000 1E-10 3 3 50 1
perf stat ./cg.out M_500.txt 1 50 10000 1E-10 3 3 100 1
perf stat ./cg.out M_500.txt 1 50 10000 1E-10 3 3 200 1
perf stat ./cg.out M_500.txt 1 50 10000 1E-10 3 3 300 1


perf stat ./cg.out M_500.txt 1 50 10000 1E-10 3 3 0 2
perf stat ./cg.out M_500.txt 1 50 10000 1E-10 3 3 50 2
perf stat ./cg.out M_500.txt 1 50 10000 1E-10 3 3 100 2
perf stat ./cg.out M_500.txt 1 50 10000 1E-10 3 3 200 2
perf stat ./cg.out M_500.txt 1 50 10000 1E-10 3 3 300 2


perf stat ./cg.out M_500.txt 1 50 10000 1E-10 3 3 0 3
perf stat ./cg.out M_500.txt 1 50 10000 1E-10 3 3 50 3



perf stat ./cg.out M_500.txt 1 50 10000 1E-12 3 3 0 1
perf stat ./cg.out M_500.txt 1 50 10000 1E-12 3 3 50 1
perf stat ./cg.out M_500.txt 1 50 10000 1E-12 3 3 100 1
perf stat ./cg.out M_500.txt 1 50 10000 1E-12 3 3 200 1
perf stat ./cg.out M_500.txt 1 50 10000 1E-12 3 3 300 1


perf stat ./cg.out M_500.txt 1 50 10000 1E-12 3 3 0 2
perf stat ./cg.out M_500.txt 1 50 10000 1E-12 3 3 50 2
perf stat ./cg.out M_500.txt 1 50 10000 1E-12 3 3 100 2
perf stat ./cg.out M_500.txt 1 50 10000 1E-12 3 3 200 2
perf stat ./cg.out M_500.txt 1 50 10000 1E-12 3 3 300 2


perf stat ./cg.out M_500.txt 1 50 10000 1E-12 3 3 0 3
perf stat ./cg.out M_500.txt 1 50 10000 1E-12 3 3 50 3








perf stat ./cg.out M3d_60.txt 1 50 10000 1E-8 3 3 0 1
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-8 3 3 50 1
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-8 3 3 100 1
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-8 3 3 200 1
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-8 3 3 300 1


perf stat ./cg.out M3d_60.txt 1 50 10000 1E-8 3 3 0 2
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-8 3 3 50 2
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-8 3 3 100 2
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-8 3 3 200 2
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-8 3 3 300 2


perf stat ./cg.out M3d_60.txt 1 50 10000 1E-8 3 3 0 3
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-8 3 3 50 3



perf stat ./cg.out M3d_60.txt 1 50 10000 1E-10 3 3 0 1
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-10 3 3 50 1
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-10 3 3 100 1
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-10 3 3 200 1
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-10 3 3 300 1


perf stat ./cg.out M3d_60.txt 1 50 10000 1E-10 3 3 0 2
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-10 3 3 50 2
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-10 3 3 100 2
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-10 3 3 200 2
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-10 3 3 300 2


perf stat ./cg.out M3d_60.txt 1 50 10000 1E-10 3 3 0 3
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-10 3 3 50 3



perf stat ./cg.out M3d_60.txt 1 50 10000 1E-12 3 3 0 1
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-12 3 3 50 1
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-12 3 3 100 1
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-12 3 3 200 1
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-12 3 3 300 1


perf stat ./cg.out M3d_60.txt 1 50 10000 1E-12 3 3 0 2
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-12 3 3 50 2
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-12 3 3 100 2
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-12 3 3 200 2
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-12 3 3 300 2


perf stat ./cg.out M3d_60.txt 1 50 10000 1E-12 3 3 0 3
perf stat ./cg.out M3d_60.txt 1 50 10000 1E-12 3 3 50 3








