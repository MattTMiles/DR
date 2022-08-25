#!/bin/sh

cd /fred/oz002/users/mmiles/MSP_DR/github_ephs

for par in $(ls J*par); do

    tempo2 -f ${par} ../notebooks/notebook_tims/${par%.*} -fit DM1 -newpar
    dm1=$(grep "^DM1" new.par | awk '{print $1}')
    dm1_val = $

