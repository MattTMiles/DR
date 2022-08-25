from __future__ import division

import os, glob, json, pickle
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as sl

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


psrs = []
ephemeris = 'DE438'
for p, t in zip(parfiles, timfiles):
    psr = Pulsar(p, t, ephem=ephemeris)
    psrs.append(psr)

## Get parameter noise dictionary
#noise = '/fred/oz002/users/mmiles/MSP_DR/clock_correction_work/enterprise_version/best_10_clock/best_10_noise.json'
noise = sys.argv[2]
params = {}
with open(noise, 'r') as fp:
    params.update(json.load(fp))

# find the maximum time span to set GW frequency sampling
tmin = [p.toas.min() for p in psrs]
tmax = [p.toas.max() for p in psrs]
Tspan = np.max(tmax) - np.min(tmin)

# define selection by observing backend
selection = selections.Selection(selections.by_backend)

# white noise parameters
efac = parameter.Constant() 
equad = parameter.Constant() 
ecorr = parameter.Constant() # we'll set these later with the params dictionary

# red noise parameters
log10_A_red = parameter.Uniform(-20, -11)
gamma_red = parameter.Uniform(0, 7)

# dm-variation parameters
log10_A_dm = parameter.Uniform(-20, -11)
gamma_dm = parameter.Uniform(0, 7)

# GW parameters (initialize with names here to use parameters in common across pulsars)
#log10_A_gw = parameter.Uniform(-18,-14)('log10_A_gw')
#gamma_gw = parameter.Constant(4.33)('gamma_gw')

# Clock parameters (initialize with names here to use parameters in common across pulsars)
log10_A_clock = parameter.Uniform(-20,-11)('log10_A_clock')
gamma_clock = parameter.Uniform(0,10)('gamma_clock')

#log10_A_clock = parameter.Uniform(-20,-11)
#gamma_clock = parameter.Uniform(0,10)

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

# white noise
ef = white_signals.MeasurementNoise(efac=efac,log10_t2equad=equad, selection=selection)
#eq = white_signals.EquadNoise(log10_equad=equad, selection=selection)
ec = white_signals.EcorrKernelNoise(log10_ecorr=ecorr, selection=selection)

# red noise (powerlaw with 30 frequencies)
pl = utils.powerlaw(log10_A=log10_A_red, gamma=gamma_red)
rn = gp_signals.FourierBasisGP(spectrum=pl, components=30, Tspan=Tspan)

# gwb (no spatial correlations)
cpl = utils.powerlaw(log10_A=log10_A_clock, gamma=gamma_clock)
clock = gp_signals.FourierBasisGP(spectrum=cpl, components=30, Tspan=Tspan, name='clock')

dm = dm_noise(log10_A=log10_A_dm,gamma=gamma_dm,Tspan=Tspan,components=30,option="powerlaw")

clock_block = blocks.common_red_noise_block(psd='powerlaw', prior='log-uniform',
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
    s = ef + ec + rn + tm + eph + gw
else:
    #s = ef + ec + rn + tm + gw
    s = ef + ec + tm + dm + rn

if bayesephem:
    s2 = ef + ec + rn + tm + eph + gw
else:
    #s = ef + ec + rn + tm + gw
    s2 = ef + ec + tm + dm + rn + clock

# intialize PTA (this cell will take a minute or two to run)
models = []
        
for p in psrs:    
    models.append(s(p))
    
pta_1 = signal_base.PTA(models)

models2 = []

for p in psrs:
    models2.append(s2(p))

pta_2 = signal_base.PTA(models2)
    
nmodels = 2
mod_index = np.arange(nmodels)

# Make dictionary of PTAs.
pta_all = dict.fromkeys(mod_index)

#pta is not including the clock
#pta2 is including the clock
pta_all[0] = pta_1
pta_all[0].set_default_params(params)
pta_all[1] = pta_2
pta_all[1].set_default_params(params)

super_model = hypermodel.HyperModel(pta_all, log_weights=[0, 0])

outDir = sys.argv[3]

sampler = super_model.setup_sampler(resume=False, outdir=outDir, sample_nmodel=True)

# This will take about an hour to sample
# To sample it properly will take 10+ hours
# sampler for N steps
N = int(1e6)  # 5e6 is a good number for a real analysis 
x0 = super_model.initial_sample()

# sample
sampler.sample(x0, N, SCAMweight=30, AMweight=15, DEweight=50, )

chain = np.loadtxt(outDir + '/chain_1.txt')
burn = int(0.25*chain.shape[0])
pars = np.loadtxt(outDir + '/pars.txt', dtype=np.unicode_)

pp = model_utils.PostProcessing(chain, pars)

ind_model = list(pars).index('nmodel')

print(pars)

chain_burn = chain[burn:,:]

ind_model = list(pars).index('nmodel')
ind_gwamp = list(pars).index('log10_A_clock')

bf, unc = model_utils.odds_ratio(chain_burn[:, ind_model], models=[0,1])
print("Bayes factor and uncertainty are:".format(bf, unc))

os.system("echo BF: "+bf+" >> "+outDir+"/evidence.txt")
os.system("echo BF_unc: "+unc+" >> "+outDir+"/evidence.txt")

log10bf = np.log10(bf)
print(log10bf)

os.system("echo log10BF: "+log10bf+" >> "+outDir+"/evidence.txt")