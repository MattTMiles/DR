# This code grabs the waveforms, and creates an average waveform out of them

import numpy as np 
import psrchive
import matplotlib.pyplot as plt
import bilby
import scipy.integrate as integrate
import sys
import os 
import glob
import argparse
import pandas as pd
import matplotlib._color_data as mcd
import matplotlib.patches as mpatch

parser = argparse.ArgumentParser(description="Subband comparison")
parser.add_argument("-data", dest="data", help="data directory to use", required = True)
parser.add_argument("-results", dest="results", help="result directory to use", required = True)
args = parser.parse_args()

data_dir = str(args.data)
results_dir = str(args.results)

if not os.path.exists(results_dir):
    try:
        os.makedirs(results_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

os.chdir(data_dir)

npt = 500
mjdpt = np.zeros(npt)

#Use 1909 for the global MJDs. We know it already works so should be fine

dir_other = "/fred/oz002/users/mmiles/MSP_DR/clock_correction_work/enterprise_version/nanograv_clones/clocks/best_10_data"

timfile_1909 = glob.glob(dir_other+"/*J1909-3744*.tim")[0]
global_mjds = np.loadtxt(timfile_1909,skiprows=1,usecols=2)

mjdmin = np.min(global_mjds)
mjdmax = np.max(global_mjds)

for i in range(npt):
    mjdpt[i] = mjdmin + ((mjdmax-mjdmin)/npt)*(i+0.5)

#Find the unique pulsars
list_init = glob.glob("J*")
list_dos = [i.split("_",1)[0] for i in list_init]
list_tres = [i.split(".",1)[0] for i in list_dos]
setlist = set(list_tres)
pulsar_list = list(setlist)


final_waveform = []

for i in range(len(np.load(pulsar_list[0]+"_clk_sim_waveform.npy"))):
    
    total_pulsar_contribution = []
    total_error_waveforms = []
    
    for pulsar in pulsar_list:
        if pulsar not in ["J1024-0719","J0125-2327","J1730-2304","J1017-7156"]:
            pulsar_waveform = np.load(pulsar+"_clk_sim_waveform.npy")
            pulsar_waveform = pulsar_waveform - np.mean(pulsar_waveform)
            error_waveform = np.load(pulsar+"_clk_sim_error_array.npy")

            pulsar_contribution = (pulsar_waveform[i]/error_waveform[i])
            total_error_waveforms.append(1/error_waveform[i])
            total_pulsar_contribution.append(pulsar_contribution)

    summed_pulsar_contribution = np.nansum(np.array(total_pulsar_contribution))
    summed_error_waveforms = np.nansum(np.array(total_error_waveforms))

    final_waveform_index = summed_pulsar_contribution/summed_error_waveforms
    final_waveform.append(final_waveform_index)

fwave = np.array(final_waveform)

mk2utc = np.loadtxt('/fred/oz002/users/mmiles/MSP_DR/clock_correction_work/enterprise_version/nanograv_clones/clocks/mk2utc.clk',skiprows=7)

plt.plot(mjdpt,fwave,label='simulated_average_waveform',color='xkcd:green')
plt.plot(mk2utc[:,0],-mk2utc[:,1] - np.mean(-mk2utc[:,1]),label='mk2utc')
plt.title("Weighted average waveform")
plt.legend()
plt.figure()


colors = mcd.XKCD_COLORS
for i, pulsar in enumerate(pulsar_list):
    if pulsar not in ["J1024-0719","J0125-2327","J1730-2304","J1017-7156"]:
        pulsar_waveform = np.load(pulsar+"_clk_sim_waveform.npy")
        error_waveform = np.load(pulsar+"_clk_sim_error_array.npy")
        plt.plot(mjdpt,pulsar_waveform,label='simulated_clock_{}'.format(pulsar),color = "tab:green")
        plt.fill_between(mjdpt, pulsar_waveform - error_waveform, pulsar_waveform + error_waveform,alpha=0.2,color="tab:green")
        plt.plot(mk2utc[:,0],-mk2utc[:,1],label='mk2utc',color='black')
        plt.legend()
        plt.figure()

plt.show()


