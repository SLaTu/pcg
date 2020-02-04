#!/bin/bash
#SBATCH --job-name=cg_test
#SBATCH --workdir=.
#SBATCH --output=out_%j.out
#SBATCH --error=out_%j.err
#SBATCH --cpus-per-task=48
#SBATCH --ntasks=1
#SBATCH --time=02:00:00

module load python
module load gcc
module load papi
module load mkl

export OMP_NUM_THREADS=48
make

# 
# python3 run_cg.py FULL_FEM_250000_500.mtx 50 1E-08 0 1
# python3 run_cg.py FULL_FEM_250000_500.mtx 50 1E-08 50 1
# python3 run_cg.py FULL_FEM_250000_500.mtx 50 1E-08 100 1
# python3 run_cg.py FULL_FEM_250000_500.mtx 50 1E-08 200 1
# python3 run_cg.py FULL_FEM_250000_500.mtx 50 1E-08 300 1
# python3 run_cg.py FULL_FEM_250000_500.mtx 50 1E-08 400 1
# 
# python3 run_cg.py FULL_FEM_250000_500.mtx 50 1E-08 0 2
# python3 run_cg.py FULL_FEM_250000_500.mtx 50 1E-08 50 2
# python3 run_cg.py FULL_FEM_250000_500.mtx 50 1E-08 100 2
# python3 run_cg.py FULL_FEM_250000_500.mtx 50 1E-08 200 2
# python3 run_cg.py FULL_FEM_250000_500.mtx 50 1E-08 300 2
# python3 run_cg.py FULL_FEM_250000_500.mtx 50 1E-08 400 2
# 
# python3 run_cg.py FULL_FEM_250000_500.mtx 50 1E-08 0 3
# python3 run_cg.py FULL_FEM_250000_500.mtx 50 1E-08 50 3
# python3 run_cg.py FULL_FEM_250000_500.mtx 50 1E-08 100 3
# python3 run_cg.py FULL_FEM_250000_500.mtx 50 1E-08 200 3



python3 run_cg.py FULL_FEM3D_250047_63.mtx 50 1E-08 0 1
python3 run_cg.py FULL_FEM3D_250047_63.mtx 50 1E-08 50 1
python3 run_cg.py FULL_FEM3D_250047_63.mtx 50 1E-08 100 1
python3 run_cg.py FULL_FEM3D_250047_63.mtx 50 1E-08 200 1
python3 run_cg.py FULL_FEM3D_250047_63.mtx 50 1E-08 300 1
python3 run_cg.py FULL_FEM3D_250047_63.mtx 50 1E-08 400 1

python3 run_cg.py FULL_FEM3D_250047_63.mtx 50 1E-08 0 2
python3 run_cg.py FULL_FEM3D_250047_63.mtx 50 1E-08 50 2
python3 run_cg.py FULL_FEM3D_250047_63.mtx 50 1E-08 100 2
python3 run_cg.py FULL_FEM3D_250047_63.mtx 50 1E-08 200 2
python3 run_cg.py FULL_FEM3D_250047_63.mtx 50 1E-08 300 2
python3 run_cg.py FULL_FEM3D_250047_63.mtx 50 1E-08 400 2
# 
# python3 run_cg.py FULL_FEM3D_250047_63.mtx 50 1E-08 0 3
# python3 run_cg.py FULL_FEM3D_250047_63.mtx 50 1E-08 50 3
# python3 run_cg.py FULL_FEM3D_250047_63.mtx 50 1E-08 100 3
# python3 run_cg.py FULL_FEM3D_250047_63.mtx 50 1E-08 200 3



python3 run_cg.py FULL_FEM3D_216000_60.mtx 50 1E-08 0 1
python3 run_cg.py FULL_FEM3D_216000_60.mtx 50 1E-08 50 1
python3 run_cg.py FULL_FEM3D_216000_60.mtx 50 1E-08 100 1
python3 run_cg.py FULL_FEM3D_216000_60.mtx 50 1E-08 200 1
python3 run_cg.py FULL_FEM3D_216000_60.mtx 50 1E-08 300 1
python3 run_cg.py FULL_FEM3D_216000_60.mtx 50 1E-08 400 1

python3 run_cg.py FULL_FEM3D_216000_60.mtx 50 1E-08 0 2
python3 run_cg.py FULL_FEM3D_216000_60.mtx 50 1E-08 50 2
python3 run_cg.py FULL_FEM3D_216000_60.mtx 50 1E-08 100 2
python3 run_cg.py FULL_FEM3D_216000_60.mtx 50 1E-08 200 2
python3 run_cg.py FULL_FEM3D_216000_60.mtx 50 1E-08 300 2
python3 run_cg.py FULL_FEM3D_216000_60.mtx 50 1E-08 400 2
# 
# python3 run_cg.py FULL_FEM3D_216000_60.mtx 50 1E-08 0 3
# python3 run_cg.py FULL_FEM3D_216000_60.mtx 50 1E-08 50 3
# python3 run_cg.py FULL_FEM3D_216000_60.mtx 50 1E-08 100 3
# python3 run_cg.py FULL_FEM3D_216000_60.mtx 50 1E-08 200 3











