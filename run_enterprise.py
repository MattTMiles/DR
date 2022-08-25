## This script is to mimic the notebook version so that this can be efficiently slurmed

### TO RUN:
# python ~/soft/DR/run_enterprise.py <data_directory> <noisefile> <outdir>

from __future__ import division

import os, glob, json, pickle
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as sl
import sys

import enterprise
from enterprise.pulsar import Pulsar
import enterprise.signals.parameter as parameter
from enterprise.signals import utils
from enterprise.signals import signal_base
from enterprise.signals import selections
from enterprise.signals.selections import Selection
from enterprise.signals import white_signals
from enterprise.signals import gp_signals
from enterprise.signals import deterministic_signals
import enterprise.constants as const

import corner
import multiprocessing
from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc

import enterprise_extensions
from enterprise_extensions import models, model_utils, hypermodel, blocks


psrlist = None # define a list of pulsar name strings that can be used to filter.
# set the data directory
datadir = sys.argv[1]
if not os.path.isdir(datadir):
    datadir = '../data'
print(datadir)

# for the entire pta
parfiles = sorted(glob.glob(datadir + '/*par'))
timfiles = sorted(glob.glob(datadir + '/*tim'))

# filter
if psrlist is not None:
    parfiles = [x for x in parfiles if x.split('/')[-1].split('.')[0] in psrlist]
    timfiles = [x for x in timfiles if x.split('/')[-1].split('.')[0] in psrlist]

# Make sure you use the tempo2 parfile for J1713+0747!!
# ...filtering out the tempo parfile... 
#parfiles = [x for x in parfiles if 'J1713+0747_NANOGrav_12yv3.gls.par' not in x]

psrs = []
ephemeris = 'DE438'
for p, t in zip(parfiles, timfiles):
    psr = Pulsar(p, t, ephem=ephemeris)
    psrs.append(psr)

## Get parameter noise dictionary
noise_file = sys.argv[2]

params = {}
with open(noise_file, 'r') as fp:
    params.update(json.load(fp))

tmin = [p.toas.min() for p in psrs]
tmax = [p.toas.max() for p in psrs]
Tspan = np.max(tmax) - np.min(tmin)

selection = selections.Selection(selections.by_backend)
##############################################################

###### Setting up model prior things #######
###### This part is majoritally commented out as we're doing something very specific, but the outlay is still useful so keeping it ######

# white noise parameters
efac = parameter.Constant()
equad = parameter.Constant()
ecorr = parameter.Constant() # we'll set these later with the params dictionary

# red noise parameters
log10_A_red = parameter.Uniform(-20, -11)
gamma_red = parameter.Uniform(0, 7)

### Commenting the below out as we know these pulsars don't have any DM and we're really only concerned with the red noise recovery

# dm-variation parameters
log10_A_dm = parameter.Uniform(-20, -11)
gamma_dm = parameter.Uniform(0, 7)

# GW parameters (initialize with names here to use parameters in common across pulsars)
#log10_A_gw = parameter.Uniform(-18,-14)('log10_A_gw')
#gamma_gw = parameter.Constant(4.33)('gamma_gw')

# Clock parameters (initialize with names here to use parameters in common across pulsars)
log10_A_clock = parameter.Uniform(-20,-11)('log10_A_clock')
gamma_clock = parameter.Uniform(0,10)('gamma_clock')

def dm_noise(log10_A,gamma,Tspan,components=30,option="powerlaw"):
    """
    A term to account for stochastic variations in DM. It is based on spin
    noise model, with Fourier amplitudes depending on radio frequency nu
    as ~ 1/nu^2.
    """
    nfreqs = 30
    if option=="powerlaw":
      pl = utils.powerlaw(log10_A=log10_A, gamma=gamma, components=components)
    #elif option=="turnover":
    #  fc = parameter.Uniform(self.params.sn_fc[0],self.params.sn_fc[1])
    #  pl = powerlaw_bpl(log10_A=log10_A, gamma=gamma, fc=fc,
    #                    components=components)
    dm_basis = utils.createfourierdesignmatrix_dm(nmodes = components,
                                                  Tspan=Tspan)
    dmn = gp_signals.BasisGP(pl, dm_basis, name='dm_gp')

    return dmn

###############################################################################################################

##### Setting up models to use ######

# white noise
ef = white_signals.MeasurementNoise(efac=efac,log10_t2equad=equad, selection=selection)
#eq = white_signals.EquadNoise(log10_equad=equad, selection=selection)
ec = white_signals.EcorrKernelNoise(log10_ecorr=ecorr, selection=selection)

# red noise (powerlaw with 30 frequencies)
pl = utils.powerlaw(log10_A=log10_A_red, gamma=gamma_red)
rn = gp_signals.FourierBasisGP(spectrum=pl, components=30, Tspan=Tspan)

#dm_noise
dm = dm_noise(log10_A=log10_A_dm,gamma=gamma_dm,Tspan=Tspan,components=30,option="powerlaw")


# clock signal
#cpl = utils.powerlaw(log10_A=log10_A_clock, gamma=gamma_clock)
#clocknoise = gp_signals.FourierBasisGP(spectrum=cpl, components=30, Tspan=Tspan, name='clock')



## This is for the monopole signature
clocknoise = blocks.common_red_noise_block(psd='powerlaw', prior='log-uniform',
                           Tspan=Tspan, components=30, combine=True,
                           log10_A_val=None, gamma_val=None, delta_val=None,
                           logmin=None, logmax=None,
                           orf="monopole", orf_ifreq=0, leg_lmax=5,
                           name='clock', coefficients=False,
                           pshift=False, pseed=None)


# gwb (no spatial correlations)
#cpl = utils.powerlaw(log10_A=log10_A_gw, gamma=gamma_gw)
#gw = gp_signals.FourierBasisGP(spectrum=cpl, components=30, Tspan=Tspan, name='gw')

# for spatial correlations you can do...
# spatial correlations are covered in the hypermodel context later
# orf = utils.hd_orf()
# crn = gp_signals.FourierBasisCommonGP(cpl, orf,
#                                       components=30, Tspan=Tspan, name='gw')

# to add solar system ephemeris modeling...
bayesephem=False
if bayesephem:
    eph = deterministic_signals.PhysicalEphemerisSignal(use_epoch_toas=True)

# timing model
tm = gp_signals.TimingModel(use_svd=True)



# full model
if bayesephem:
    s = ef + ec + clocknoise + tm + eph
else:
    s = ef + ec + clocknoise + tm

# intialize PTA (this cell will take a minute or two to run)
models = []
        
for p in psrs:    
    models.append(s(p))
    
pta = signal_base.PTA(models)

# set white noise parameters with dictionary
pta.set_default_params(params)

# set initial parameters drawn from prior
x0 = np.hstack([p.sample() for p in pta.params])
ndim = len(x0) 

# set up the sampler:
# initial jump covariance matrix
cov = np.diag(np.ones(ndim) * 0.01**2)
outDir = sys.argv[3]

sampler = ptmcmc(ndim, pta.get_lnlikelihood, pta.get_lnprior, cov, 
                 outDir=outDir, resume=False)

# sampler for N steps
N = int(1e5)  # normally, we would use 5e6 samples (this will save time)
x0 = np.hstack([p.sample() for p in pta.params])
sampler.sample(x0, N, SCAMweight=30, AMweight=15, DEweight=50, )

