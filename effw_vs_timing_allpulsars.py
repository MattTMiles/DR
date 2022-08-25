import numpy as np
import psrchive
import pandas as pd
import matplotlib.pyplot as plt 
import scipy.integrate
import glob
import sys
import os
import pickle


res_folder = '/fred/oz002/users/mmiles/MSP_DR/notebooks/notebook_residuals'

weff_dir = '/fred/oz002/users/mmiles/MSP_DR/Weff_dicts/'

pulsar_list = '/fred/oz002/users/mmiles/MSP_DR/pulsar_list.txt'

tim_folder = '/fred/oz002/users/mmiles/MSP_DR/notebooks/notebook_tims/'

os.chdir(weff_dir)

cfreqs = [916,968,1016,1064,1114,1158,1210,1256,1307,1356,1404,1453,1501,1550,1600,1648]

pulsar_extend = []
freq_total = []

with open(pulsar_list,'r') as pulsarfile:
    for pulsar in pulsarfile:
        pulsar = pulsar.strip()
        print(pulsar)
        pulsar_file = glob.glob(pulsar+'*pkl')[0]
        pulsar_dict = pd.read_pickle(pulsar_file)
        timfile = glob.glob(tim_folder+'*'+pulsar+'*')[0]

        try:
            
            eff_widths = np.array(list(pulsar_dict.values()))

            res = glob.glob(res_folder+'/'+pulsar+'*dat')[0]

            mjds = np.loadtxt(res,usecols=0)
            freqs = np.loadtxt(res,usecols=1)
            resids = np.loadtxt(res,usecols=2)
            unc = np.loadtxt(res,usecols=3)

            print(timfile)
            snrtim = np.loadtxt(timfile,skiprows=1,usecols=-1)

            #weighted_eff = eff_widths/snrtim
            
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

            templist = [unc1,unc2,unc3,unc4,unc5,unc6,unc7,unc8,unc9,unc10,unc11,\
                unc12,unc13,unc14,unc15,unc16]

            snrlist = [snr1,snr2,snr3,snr4,snr5,snr6,snr7,snr8,snr9,snr10,snr11,\
                snr12,snr13,snr14,snr15,snr16]


            for i in range(len(resids)):
            
                closest_freq = min(cfreqs, key=lambda x:abs(x-freqs[i]))
                freqindex = cfreqs.index(closest_freq)

                templist[freqindex].append(unc[i])
                snrlist[freqindex].append(snrtim[i])
                freq_total.append(closest_freq)

            for j in range(len(templist)):
                pulsar_extend.append([[eff_widths[j]]*len(templist[j]),templist[j], snrlist[j],[cfreqs[j]]*len(templist[j])])
            
            print('yes')
        except:
            print('no: {}'.format(pulsar))
            continue

extended_array = np.hstack(pulsar_extend)
        
eff_width_df = pd.DataFrame(extended_array.T, columns=['eff_width','uncertainty','snr','freq'])
eff_width_df['weighted_eff_width'] = eff_width_df['eff_width']/eff_width_df['snr']

eff_width_df.to_pickle("/fred/oz002/users/mmiles/MSP_DR/master_weff_dict.pkl")

eff_width_df.plot.scatter('weighted_eff_width','uncertainty',c='freq', colormap='viridis')

plt.yscale('log')
plt.xscale('log')
plt.show()

'''
eff_width_df.plot.scatter('weighted_eff_width','freq')
plt.xscale('log')
plt.show()
'''

'''
fig, ax = plt.subplots() 

#for i in range(len(extended_array[0])):
#    ax.scatter(extended_array[0,i],extended_array[1,i],marker='x',label=extended_array[2,i]) 

ax.scatter(extended_array[0], extended_array[1], marker='x')

ax.legend(extended_array[2])

#fig.legend(symbols, labels, fontsize = 80, bbox_transform = plt.gcf().transFigure, ncol=3)

ax.set_xlabel('Effective width (s)') 
ax.set_ylabel(r'Uncertainty ($\mu$s)') 

ax.set_title('Relationship between timing uncertainty and effective width: All pulsars') 
  
fig.tight_layout() 
fig.show()
''' 



