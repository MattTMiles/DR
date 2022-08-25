import numpy as np
import psrchive
import pandas as pd
import matplotlib.pyplot as plt 
import scipy.integrate
import glob
import sys
import os
import pickle

pulsar = sys.argv[1]

os.chdir('/fred/oz002/users/mmiles/MSP_DR/notebooks')

Weffs={}

for subband in sorted(glob.glob('*'+pulsar+'*0000.16chan')):
    print(subband)

    snr = os.popen('psrstat -c snr=pdmp -c snr -Q '+subband)
    snr = snr.strip()
    snr = snr.split(' ')[1]
    snr = float(snr)
    arch = psrchive.Archive_load(subband)
    arch.remove_baseline()
    period = arch[0].get_folding_period()

    binlength = period/1024 #in seconds

    archdata = arch.get_data()
    archdata = archdata[0,0,0,:]

    #Normalise to unity otherwise this won't work
    archscaled = archdata/archdata.max()

    diffs = []
    for i in range(1, len(archscaled)-1):
        diff = archscaled[i+1]-archscaled[i]
        diffs.append(diff)
    
    diffs = np.array(diffs)
    diff2 = diffs**2
    summed = np.sum(diff2)
    
    Weff = (binlength)/(summed)

    Weffs[subband]=Weff

dict_dir = '/fred/oz002/users/mmiles/MSP_DR/Weff_dicts'

with open(dict_dir+'/'+pulsar+'_weff_dict.pkl', 'wb') as f:
    pickle.dump(Weffs, f)




