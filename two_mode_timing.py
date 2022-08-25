#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 21:48:07 2021
​
@author: dreardon
"""

import numpy as np
import matplotlib.pyplot as plt
import bilby
import psrchive
import os


def fft_rotate(data, bins):
    """Return data rotated by 'bins' places to the left. The
       rotation is done in the Fourier domain using the Shift Theorem.
​
       Inputs:
            data: A 1-D numpy array to rotate.
            bins: The (possibly fractional) number of phase bins to rotate by.
​
       Outputs:
            rotated: The rotated data.
    """
    freqs = np.arange(data.size/2+1, dtype=np.float)
    phasor = np.exp(complex(0.0, 2.0*np.pi) * freqs * bins / float(data.size))
    return np.fft.irfft(phasor*np.fft.rfft(data))


"""
Block of code to generate simulated data and templates for two modes
"""

# Pulse phase array

nbins = 1024
x = np.linspace(0, 1, nbins)


strong_width = 0.02
weak_width = 0.05
noise = 0.0012


#datas = psrchive.Archive_load('/fred/oz002/users/mmiles/MSP_DR/archs/grand.p')
datas = psrchive.Archive_load('/fred/oz002/users/mmiles/MSP_DR/J1103_study/archs/grand.dly')

datas.remove_baseline()
datas.dedisperse()
datas = datas.get_data()
datas = datas[:,0,0,:]
data = datas[0]


for i, data in enumerate(datas):
# Bilby likelihood
    
    os.system('sbatch ~/soft/DR/two_mode_slurm.sh '+str(i))


