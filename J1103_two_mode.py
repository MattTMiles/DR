import numpy as np
import matplotlib.pyplot as plt
import bilby
import psrchive
import os
import sys

i = sys.argv[1]

nbins = 1024
x = np.linspace(0, 1, nbins)

noise = 0.0012


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


#psradd versions
early_template = np.load('/fred/oz002/users/mmiles/MSP_DR/J1103_study/archs/smoothed_vearly_scaled_standard.npy')
late_template = np.load('/fred/oz002/users/mmiles/MSP_DR/J1103_study/archs/smoothed_vearly_diff_standard.npy')

#Rotate towards strong
#early_template = fft_rotate(early_template,early_template.argmax()-late_template.argmax())


datas = psrchive.Archive_load('/fred/oz002/users/mmiles/MSP_DR/J1103_study/archs/grand.dly')
datas.remove_baseline()
datas.dedisperse()
datas = datas.get_data()
datas = datas[:,0,0,:]

def avg_profile_model(x, early_amp, early_phase, late_amp, late_phase):
 
    early_mode = fft_rotate(early_amp * late_amp * early_template / np.max(early_template), late_phase + early_phase)
    late_mode = fft_rotate(late_amp * late_template / np.max(late_template), late_phase)

    return late_mode + early_mode



data = datas[int(i)]
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

def avg_profile_model_fdm_fast(x, early_amp, early_phase, late_amp, late_phase):
    """
    Model for the average profile given two modes
    NOTE: The templates must already be defined in this script
    """
    early_mode = phasor_scale_fft(early_amp*rfft_early_template, early_phase)
    late_mode = phasor_scale_fft(late_amp*rfft_late_template, late_phase)
    return np.concatenate((np.real(late_mode + early_mode), np.imag(late_mode + early_mode)))


def avg_profile_model_fdm_fast_second_model(x, early_amp, late_amp, phase):
    """
    Model for the average profile given two modes
    NOTE: The templates must already be defined in this script
    """
    early_mode = phasor_scale_fft(early_amp*rfft_early_template, phase)
    late_mode = phasor_scale_fft(late_amp*rfft_late_template, phase)

    return np.concatenate((np.real(late_mode + early_mode), np.imag(late_mode + early_mode)))

data_fdm = np.concatenate((np.real(rfft_data), np.imag(rfft_data)))
noise_fdm = noise*np.sqrt(len(data_fdm)/2)

likelihood = bilby.likelihood.GaussianLikelihood(x, data_fdm,
                                                avg_profile_model_fdm_fast_second_model)

priors = dict()
#priors['weak_amp'] = bilby.core.prior.Gaussian(mu=0.116654865, sigma=0.031803366, name='weak_amp')
priors['early_amp'] = bilby.core.prior.Uniform(0, 100, 'early_amp')
#priors['early_phase'] = bilby.core.prior.Uniform(-nbins/2, nbins/2,
#                                                'early_phase')
#priors['weak_phase'] = bilby.core.prior.Gaussian(mu=6.434, sigma=2.205031088120973, name='weak_phase')
                                             
priors['late_amp'] = bilby.core.prior.Uniform(0, 100, 'late_amp')
priors['phase'] = bilby.core.prior.Uniform(-nbins/2, nbins/2,
                                                'phase')
priors['sigma'] = bilby.core.prior.Uniform(0, 100, 'sigma')



results = bilby.core.sampler.run_sampler(
    likelihood, priors=priors, sampler='dynesty', label='dynesty',
    nlive=1000, verbose=True, resume=False, npool=4,
    outdir='J1103_subint_{}'.format(i))

