#!/bin/bash
#SBATCH -J SolveU
#SBATCH -p cpu
#SBATCH -t 00:30:00
#SBATCH -o outfile-%J
#SBATCH --cpus-per-task=4
#SBATCH -N 1 
#SBATCH -q debug
#SBATCH --ntasks-per-node 1

export OMPI_MCA_btl_base_warn_component_unused=0   ## 避免运行MEME时出现警告信息

module load anaconda/3
module load gcc/9.3.0
module load intel/2019.5.281
module load openmpi/3.0.0-new
module load ngsolve/parallel

python3 SolveU.py