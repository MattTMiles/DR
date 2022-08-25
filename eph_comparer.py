import numpy as np
import os
import glob


git_dir = '/fred/oz002/users/mmiles/MSP_DR/github_ephs'
psrcat_dir = '/fred/oz002/users/mmiles/MSP_DR/psrcat_ephs'
pulsar_list = '/fred/oz002/users/mmiles/MSP_DR/pulsar_list.txt'


with open(pulsar_list,'r') as pulsarfile:
    for pulsar in pulsarfile:
        print('Pulsar: {}'.format(pulsar))
        
        github_par = np.loadtxt(git_dir+'/*'+pulsar+'*',usecols=0)
        psrcat_par = np.loadtxt(psrcat_dir+'/*'+pulsar+'*',usecols=0)
        


