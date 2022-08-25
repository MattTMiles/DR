import numpy as np
import psrchive
import pandas as pd
import matplotlib.pyplot as plt 
import scipy.integrate
import glob
import sys
import os
import pickle
from scipy.optimize import curve_fit

res_folder = '/fred/oz002/users/mmiles/MSP_DR/notebooks/notebook_residuals'

weff_dir = '/fred/oz002/users/mmiles/MSP_DR/Weff_dicts/'

tim_folder = '/fred/oz002/users/mmiles/MSP_DR/notebooks/notebook_tims/'

os.chdir(weff_dir)

pulsar = sys.argv[1]

pulsar_extend = []

pulsar_file = glob.glob(pulsar+'*pkl')[0]
pulsar_dict = pd.read_pickle(pulsar_file)

eff_widths = np.array(list(pulsar_dict.values()))

cfreqs = [916,968,1016,1064,1114,1158,1210,1256,1307,1356,1404,1453,1501,1550,1600,1648]

res = glob.glob(res_folder+'/'+pulsar+'*dat')[0]

mjds = np.loadtxt(res,usecols=0)
freqs = np.loadtxt(res,usecols=1)
resids = np.loadtxt(res,usecols=2)
unc = np.loadtxt(res,usecols=3)

timfile = glob.glob(tim_folder+'*'+pulsar+'*')[0]
snrtim = np.loadtxt(timfile,skiprows=1,usecols=-1)

templist = []
snrlist = []

unc1 = []
unc2 = []
unc3 = []
unc4 = []
unc5 = []
unc6 = []
unc7 = []
unc8 = []
unc9 = []
unc10 = []
unc11 = []
unc12 = []
unc13 = []
unc14 = []
unc15 = []
unc16 = []

snr1 = []
snr2 = []
snr3 = []
snr4 = []
snr5 = []
snr6 = []
snr7 = []
snr8 = []
snr9 = []
snr10 = []
snr11 = []
snr12 = []
snr13 = []
snr14 = []
snr15 = []
snr16 = []

pulsar_extend = []

templist = [unc1,unc2,unc3,unc4,unc5,unc6,unc7,unc8,unc9,unc10,unc11,\
    unc12,unc13,unc14,unc15,unc16]

snrlist = [snr1,snr2,snr3,snr4,snr5,snr6,snr7,snr8,snr9,snr10,snr11,\
    snr12,snr13,snr14,snr15,snr16]

for i in range(len(resids)):
    
    closest_freq = min(cfreqs, key=lambda x:abs(x-freqs[i]))
    freqindex = cfreqs.index(closest_freq)

    templist[freqindex].append(unc[i])
    snrlist[freqindex].append(snrtim[i])

for j in range(len(templist)):
    pulsar_extend.append([[eff_widths[j]]*len(templist[j]),templist[j], snrlist[j],[cfreqs[j]]*len(templist[j])])


extended_array = np.hstack(pulsar_extend)

eff_width_df = pd.DataFrame(extended_array.T, columns=['eff_width','uncertainty','snr','freq'])
eff_width_df['weighted_eff_width'] = eff_width_df['eff_width']/eff_width_df['snr']


def power_law(x, a, b):
    return a*np.power(x, b)

pars, cov = curve_fit(f=power_law, xdata=eff_width_df['freq'], ydata=eff_width_df['eff_width'])

weff_fit = power_law(eff_width_df['freq'], *pars)

ax = eff_width_df.plot.scatter('weighted_eff_width','uncertainty',c='freq', colormap='viridis')

ax.set_xlabel('Weighted effective width (s)')
ax.set_ylabel(r'Uncertainty ($\mu$s)')

ax.set_title('Timing uncertainty vs weighted effective width: '+pulsar)

plt.tight_layout()
plt.savefig('{0}_uncertainty_vs_weighted_weff.pdf'.format(pulsar),format='pdf',dpi=80)
#plt.figure()


ax = eff_width_df.plot.scatter('weighted_eff_width','freq')

#ax.legend()
ax.set_xlabel('Weighted effective width (s)')
ax.set_ylabel(r'Frequency (MHz)')

ax.set_title('Frequency vs weighted effective width: '+pulsar)

plt.tight_layout()
plt.savefig('{0}_frequency_vs_weighted_weff.pdf'.format(pulsar),format='pdf',dpi=80)
#plt.show()

ax = eff_width_df.plot.scatter('freq','eff_width',label='data')

ax.plot(eff_width_df['freq'],weff_fit, linestyle = '--',color='k',label = 'power law fit; a={0:.4f}, b={1:.4f}'.format(pars[0],pars[1]))
ax.legend()
ax.set_ylabel('Effective width (s)')
ax.set_xlabel(r'Frequency (MHz)')

ax.set_title('Frequency vs effective width: '+pulsar)

plt.tight_layout()
plt.savefig('{0}_weff_vs_frequency.pdf'.format(pulsar),format='pdf',dpi=80)
#plt.show()


