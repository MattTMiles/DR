#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --job-name=J1103_mode
#SBATCH --mem=5gb
#SBATCH --tmp=5gb

#cd /fred/oz002/users/mmiles/SinglePulse/two-mode-timing-slurm-uniform-pool
touch ${1}".J1103_mode"

srun python /home/mmiles/soft/DR/J1103_two_mode.py $1

rm -f ${1}".J1103_mode"

echo done