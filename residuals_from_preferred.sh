#!/bin/bash

# This script quickly makes the new residual and pdf for a pulsar in the preferred noise model directory

pulsar=$1 ;

pref_model_dir=/fred/oz002/users/mmiles/MSP_DR/noise_modelling/new_models/preferred_model_ephs ;
pref_par=${pref_model_dir}/MTMSP-${pulsar}-.par ;

residual_dir=/fred/oz002/users/mmiles/MSP_DR/noise_subtracted_residuals ;

cd ${residual_dir} ;

ave_with_noise="${residual_dir}/partim/ave" ;
subband_with_noise="${residual_dir}/partim/16ch" ;

cd ${ave_with_noise} ;
cp ${pref_par} ${ave_with_noise} ;

echo AVERAGERES -f 0.1 >> MTMSP-${pulsar}-.par ;

echo first_tempo2

/apps/users/pulsar/skylake/software/tempo2/0759584/bin/tempo2 -f MTMSP-${pulsar}-.par ${pulsar}_16ch.tim; mv avg.dat ${pulsar}_avg_residuals.dat ;


cd ${subband_with_noise} ;
cp ${pref_par} ${subband_with_noise} ;

echo second_tempo

/apps/users/pulsar/skylake/software/tempo2/0759584/bin/tempo2 -f MTMSP-${pulsar}-.par -output general2 -s "{sat} {freqSSB} {posttn} {err}\n" -outfile ${pulsar}_16ch_residuals.dat ${pulsar}_16ch.tim ;

ave_no_noise="${residual_dir}/F_residuals/partim" ;
subband_no_noise="${residual_dir}/f16_residuals/partim" ;

cd ${ave_no_noise} ;
cp ${pref_par} ${ave_no_noise} ;

echo AVERAGERES -f 0.1 >> MTMSP-${pulsar}-.par ; 
echo TNsubtractDM 1 >> MTMSP-${pulsar}-.par ;
echo TNsubtractRed 1 >> MTMSP-${pulsar}-.par ;

echo third_tempo

/apps/users/pulsar/skylake/software/tempo2/0759584/bin/tempo2 -f MTMSP-${pulsar}-.par ${pulsar}_16ch.tim; mv avg.dat ../${pulsar}_avg_residuals.dat ;

cd ${subband_no_noise} ;
cp ${pref_par} ${subband_no_noise} ;

echo TNsubtractDM 1 >> MTMSP-${pulsar}-.par ;
echo TNsubtractRed 1 >> MTMSP-${pulsar}-.par ;

echo fourth_tempo

/apps/users/pulsar/skylake/software/tempo2/0759584/bin/tempo2 -f MTMSP-${pulsar}-.par -output general2 -s "{sat} {freqSSB} {posttn} {err}\n" -outfile ../${pulsar}_16ch_residuals.dat ${pulsar}_16ch.tim ;

#python code for creating the residual files
python ~/soft/DR/res_plotter.py ${ave_with_noise} ave ;
python ~/soft/DR/res_plotter.py ${subband_with_noise} ;

python ~/soft/DR/res_plotter.py ${ave_no_noise} ave ;
python ~/soft/DR/res_plotter.py ${subband_no_noise} ;