#!/bin/bash
#SBATCH --job-name=playGod_trial
#SBATCH --output=/fred/oz002/users/mmiles/MSP_DR/clock_correction_work/playGod_noise
#SBATCH --ntasks=4
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=10GB
#SBATCH --mail-type=FAIL --mail-user=matthewmiles@swin.edu.au
#SBATCH --tmp=10GB

#cd /fred/oz002/users/mmiles/MSP_DR/clock_correction_work/playGod_noise

#psrdir=/fred/oz002/users/mmiles/MSP_DR/clock_correction_work/no_redDM_partim

#psr=$1
#parfile=MTMSP-${psr}-.par
#timfile=${psr}.tim

#### REALLY TEMPORARY VERSION DO NOT USE FOR ALL
srun --ntasks=4 --mem-per-cpu=10GB --time=12:00:00 python ~/soft/DR/clock_recover.py ${psr}
#srun python ~/soft/DR/clock_recover.py J1909-3744 /fred/oz002/users/mmiles/MSP_DR/clock_correction_work/playGod_noise/test_1909/1909_sim.tim /fred/oz002/users/mmiles/MSP_DR/clock_correction_work/playGod_noise/test_1909/test-clock-J1909-3744-.par 1909_sim_residuals.dat

