import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import sys

import pandas as pd

res_folder = sys.argv[1]
#res_folder = '/fred/oz002/users/mmiles/MSP_DR/notebooks/notebook_residuals'
try:
    avg = sys.argv[2]
except IndexError:
    print("Avg not selected")
    avg = None

os.chdir(res_folder)

total_data = []

for res in glob.glob('J*dat'):
    print(res)
    if avg is None:
        mjds = np.loadtxt(res,usecols=0)
        resids = np.loadtxt(res,usecols=2)
        unc = np.loadtxt(res,usecols=3)
        unc = unc*(10**-6)

    pulsar = res.split('_')[0]

    if avg is not None:
        mjds = np.loadtxt(res,usecols=0)
        resids = np.loadtxt(res,usecols=1)
        unc = np.loadtxt(res,usecols=2)
        #unc = unc*(10**-6)

        mjds,unique_ind = np.unique(mjds,return_index=True)
        resids = resids[unique_ind]
        unc = unc[unique_ind]


    weights = 1/(unc**2)  
    xbar = np.sum(resids*weights)/np.sum(weights)
    wrms = np.sqrt(np.sum(weights*((resids-xbar)**2))/np.sum(weights))
    wrms_micro = wrms*10**6
    #print(rms)

    data = [pulsar,wrms_micro]
    total_data.append(data)

    fig, ax = plt.subplots(figsize=(15,5))
    
    ax.errorbar(mjds,resids,yerr=unc, fmt='x',color='black',ecolor='black')
    ax.set_title(pulsar+r' frequency averaged; wrms = {0:.4f} $\mu$s'.format(wrms_micro))
    ax.set_xlabel('MJD')
    ax.set_ylabel(r'$\mu$s')

    fig.tight_layout()
    fig.savefig(pulsar+'_16ch_residuals.pdf', format='pdf',dpi=1000)
    plt.close()

df = pd.DataFrame(total_data,columns=["Pulsar","wrms"])
