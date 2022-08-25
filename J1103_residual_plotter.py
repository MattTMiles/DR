# This is exclusively for comparing profiles and creating residuals of J1103 in the J1103 directory 

import numpy as np
import psrchive
import matplotlib.pyplot as plt
import bilby
import os
import re

data = psrchive.Archive_load('/fred/oz002/users/mmiles/MSP_DR/J1103_study/archs/grand.dly')
data.remove_baseline()
data.dedisperse()

data = data.get_data()
data = data[:,0,0,:]

nbins = 1024
x = np.linspace(0, 1, nbins)


early_template = np.load('/fred/oz002/users/mmiles/MSP_DR/J1103_study/archs/smoothed_vearly_scaled_standard.npy')
late_template = np.load('/fred/oz002/users/mmiles/MSP_DR/J1103_study/archs/smoothed_vearly_diff_standard.npy')

rfft_data = np.fft.rfft(data)
rfft_early_template = np.fft.rfft(early_template)
rfft_late_template = np.fft.rfft(late_template)

def phasor_scale_fft(rfft_data, bins):
    """
    Add a phase gradient to the rotated FFT input
    """
    freqs = np.arange(rfft_data.size, dtype=np.float)
    phasor = np.exp(complex(0.0, 2.0*np.pi) * freqs * bins / float(2*(rfft_data.size - 1)))
    return phasor*rfft_data

def avg_profile_model_fdm_fast_second_model(x, early_amp, late_amp, phase):
    """
    Model for the average profile given two modes
    NOTE: The templates must already be defined in this script
    """
    early_mode = phasor_scale_fft(early_amp*rfft_early_template, phase)
    late_mode = phasor_scale_fft(late_amp*rfft_late_template, phase)

    return np.concatenate((np.real(late_mode + early_mode), np.imag(late_mode + early_mode))), early_mode, late_mode


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

pulsar_dir = '/fred/oz002/users/mmiles/MSP_DR/J1103_study/bilby_run_19112021'


rfft_early_template = np.fft.rfft(early_template)
rfft_late_template = np.fft.rfft(late_template)
data_fdm = np.concatenate((np.real(rfft_data), np.imag(rfft_data)))

dir_files = os.listdir(pulsar_dir)  
dir_files = [ x for x in dir_files if x.startswith('J1103') ]
dir_files = [ x for x in dir_files if not x.endswith('png') ]
dir_files.sort(key=lambda f: int(re.sub('\D', '', f)))  

for i, subintdir in enumerate(dir_files):  
    if subintdir.startswith('J1103'):
        if not subintdir.endswith('png'):
            os.chdir(subintdir)
            
            results =bilby.result.read_in_result(filename='dynesty_result.json')
            late_amp = results.get_one_dimensional_median_and_error_bar('late_amp').median 
            early_amp = results.get_one_dimensional_median_and_error_bar('early_amp').median
            phase = results.get_one_dimensional_median_and_error_bar('phase').median

            active_data = data[i]
            late = late_amp*late_template
            late = fft_rotate(late, phase)
            early = early_amp*early_template
            early = fft_rotate(early, phase)
            average = late+early
            residual = active_data-average
            plt.plot(residual)
            os.chdir(pulsar_dir)

'''        fig,ax =plt.subplots(2,1)
        ax[0].plot(active_data,label='Observation')
        ax[0].plot(average,label='Average mode')
        ax[0].plot(early,label='Early mode')
        ax[0].plot(late,label='Late_residual component')
        ax[0].legend()

        ax[1].plot(residual,label='Residual (data-average)')
        ax[1].legend()
        
        os.chdir(pulsar_dir)
        fig.savefig(subintdir+'_res.png')'''
        
        

plt.show()
        







