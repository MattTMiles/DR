#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --job-name=monopole_best_10
#SBATCH --mem=1gb
#SBATCH --tmp=1gb



touch "1909_monopole"


#python ~/soft/DR/clock_recover.py $1

#srun --cpus-per-task=4 --time=48:00:00 python /fred/oz002/users/mmiles/MSP_DR/clock_correction_work/enterprise_version/nanograv_clones/run_enerprise_1909_mono.py /fred/oz002/users/mmiles/MSP_DR/clock_correction_work/enterprise_version/full_search/data /fred/oz002/users/mmiles/MSP_DR/clock_correction_work/enterprise_version/full_search/full_search_noise.json /fred/oz002/users/mmiles/MSP_DR/clock_correction_work/enterprise_version/nanograv_clones/results/1909_mono/

srun --cpus-per-task=4 --time=48:00:00 python ~/soft/DR/run_enterprise.py /fred/oz002/users/mmiles/MSP_DR/clock_correction_work/enterprise_version/best_10_clock/data /fred/oz002/users/mmiles/MSP_DR/clock_correction_work/enterprise_version/best_10_clock/best_10_noise.json /fred/oz002/users/mmiles/MSP_DR/clock_correction_work/enterprise_version/nanograv_clones/results/best_10/CRN_Monopole_noDMRN_chain_$1

#srun --cpus-per-task=1 --time=03:00:00 python ~/soft/DR/clock_recover.py -data /fred/oz002/users/mmiles/MSP_DR/clock_correction_work/enterprise_version/nanograv_clones/clocks/best_10_data -results /fred/oz002/users/mmiles/MSP_DR/clock_correction_work/enterprise_version/nanograv_clones/clocks/best_10 -solo $1

rm -f "1909_monopole"


echo done
