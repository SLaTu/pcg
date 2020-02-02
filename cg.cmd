#!/bin/bash
#SBATCH --job-name=cg_test
#SBATCH --workdir=.
#SBATCH --output=out_%j.out
#SBATCH --error=out_%j.err
#SBATCH --cpus-per-task=48
#SBATCH --ntasks=1
#SBATCH --time=00:20:00
#SBATCH --qos=debug

module load python
module load gcc
module load papi
module load mkl

export OMP_NUM_THREADS=48
make

# python3 run_cg.py FULL_thermomech_TC.mtx 20 1E-08 0 1
# python3 run_cg.py FULL_thermomech_TC.mtx 20 1E-08 50 1
# python3 run_cg.py FULL_thermomech_TC.mtx 20 1E-08 100 1
# python3 run_cg.py FULL_thermomech_TC.mtx 20 1E-08 200 1
# python3 run_cg.py FULL_thermomech_TC.mtx 20 1E-08 300 1
# python3 run_cg.py FULL_thermomech_TC.mtx 20 1E-08 400 1
# 
# 
# 
python3 run_cg.py FULL_FEM3D_421875_75.mtx 50 1E-08 0 1
python3 run_cg.py FULL_FEM3D_421875_75.mtx 50 1E-08 50 1
# python3 run_cg.py FULL_FEM3D_421875_75.mtx 50 1E-08 100 1
# python3 run_cg.py FULL_FEM3D_421875_75.mtx 50 1E-08 200 1
# python3 run_cg.py FULL_FEM3D_421875_75.mtx 50 1E-08 300 1
# python3 run_cg.py FULL_FEM3D_421875_75.mtx 50 1E-08 400 1
# 
# 
# 
# python3 run_cg.py FULL_FEM_1000000_1000.mtx 50 1E-08 0 1
# python3 run_cg.py FULL_FEM_1000000_1000.mtx 50 1E-08 50 1
# python3 run_cg.py FULL_FEM_1000000_1000.mtx 50 1E-08 100 1
# python3 run_cg.py FULL_FEM_1000000_1000.mtx 50 1E-08 200 1
# python3 run_cg.py FULL_FEM_1000000_1000.mtx 50 1E-08 300 1
# python3 run_cg.py FULL_FEM_1000000_1000.mtx 50 1E-08 400 1
# 
# 
# 
# python3 run_cg.py FULL_FEM3D_1000000_100.mtx 50 1E-08 0 1
# python3 run_cg.py FULL_FEM3D_1000000_100.mtx 50 1E-08 50 1
# python3 run_cg.py FULL_FEM3D_1000000_100.mtx 50 1E-08 100 1
# python3 run_cg.py FULL_FEM3D_1000000_100.mtx 50 1E-08 200 1
# python3 run_cg.py FULL_FEM3D_1000000_100.mtx 50 1E-08 300 1
# python3 run_cg.py FULL_FEM3D_1000000_100.mtx 50 1E-08 400 1
# 
# 
# 
# python3 run_cg.py FULL_FEM_250000_500.mtx 50 1E-08 0 1
# python3 run_cg.py FULL_FEM_250000_500.mtx 50 1E-08 50 1
# python3 run_cg.py FULL_FEM_250000_500.mtx 50 1E-08 100 1
# python3 run_cg.py FULL_FEM_250000_500.mtx 50 1E-08 200 1
# python3 run_cg.py FULL_FEM_250000_500.mtx 50 1E-08 300 1
# python3 run_cg.py FULL_FEM_250000_500.mtx 50 1E-08 400 1
