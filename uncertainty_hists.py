import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import pandas as pd


res_folder = '/fred/oz002/users/mmiles/MSP_DR/notebooks/notebook_residuals'

os.chdir(res_folder)

total_uncertainty = []

tunc1 = []
tunc2 = []
tunc3 = []
tunc4 = []
tunc5 = []
tunc6 = []
tunc7 = []
tunc8 = []
tunc9 = []
tunc10 = []
tunc11 = []
tunc12 = []
tunc13 = []
tunc14 = []
tunc15 = []
tunc16 = []

cfreqs = [916,968,1016,1064,1114,1158,1210,1256,1307,1356,1404,1453,1501,1550,1600,1648]

for res in glob.glob('J*dat'):
    print(res)
    mjds = np.loadtxt(res,usecols=0)
    freqs = np.loadtxt(res,usecols=1)
    resids = np.loadtxt(res,usecols=2)
    unc = np.loadtxt(res,usecols=3)
    #unc = unc*(10**-6)

    total_uncertainty.append(unc)

    pulsar = res.split('_')[0]

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

    templist = [unc1,unc2,unc3,unc4,unc5,unc6,unc7,unc8,unc9,unc10,unc11,\
        unc12,unc13,unc14,unc15,unc16]

    for i in range(len(resids)):
        
        closest_freq = min(cfreqs, key=lambda x:abs(x-freqs[i]))
        freqindex = cfreqs.index(closest_freq)

        templist[freqindex].append(unc[i])

    '''  
    unc1 = unc[::16]
    unc2 = unc[1::16]
    unc3 = unc[2::16]
    unc4 = unc[3::16]
    unc5 = unc[4::16]
    unc6 = unc[5::16]
    unc7 = unc[6::16]
    unc8 = unc[7::16]
    unc9 = unc[8::16]
    unc10 = unc[9::16]
    unc11 = unc[10::16]
    unc12 = unc[11::16]
    unc13 = unc[12::16]
    unc14 = unc[13::16]
    unc15 = unc[14::16]
    unc16 = unc[15::16]
    '''

    tunc1.append(unc1)
    tunc2.append(unc2)
    tunc3.append(unc3)
    tunc4.append(unc4)
    tunc5.append(unc5)
    tunc6.append(unc6)
    tunc7.append(unc7)
    tunc8.append(unc8)
    tunc9.append(unc9)
    tunc10.append(unc10)
    tunc11.append(unc11)
    tunc12.append(unc12)
    tunc13.append(unc13)
    tunc14.append(unc14)
    tunc15.append(unc15)
    tunc16.append(unc16)


    unc_series = [unc1,unc2,unc3,unc4,unc5,unc6,unc7,unc8,unc9,unc10,unc11,unc12,unc13,unc14,unc15,unc16]


    #Total pulsar subplots
    fig,ax = plt.subplots(figsize = (15,5))

    ax.hist(unc,bins=50)
    ax.set_xlabel(r'Uncertainty ($\mu$s)')
    ax.set_ylabel('Count')
    ax.set_title('Total TOA uncertainty distribution for '+pulsar)
    fig.savefig(pulsar+'_total_uncertainties.pdf',format='pdf',dpi=1000)


    #Subbanded subplots
    fig,ax = plt.subplots(figsize = (15,5))

    for unc_s in unc_series:
        ax.hist(unc_s,bins=50)
    
    ax.set_xlabel(r'Uncertainty ($\mu$s)')
    ax.set_ylabel('Count')
    ax.set_title('Subbanded TOA uncertainty distribution for '+pulsar)

    ax.legend(['Subband_1','Subband_2','Subband_3','Subband_4','Subband_5','Subband_6','Subband_7','Subband_8','Subband_9','Subband_10','Subband_11','Subband_12','Subband_13','Subband_14','Subband_15','Subband_16'])

    fig.savefig(pulsar+'_subbanded_uncertainties.pdf',format='pdf',dpi=1000)


tunc_series = [np.concatenate(tunc1),np.concatenate(tunc2),np.concatenate(tunc3),np.concatenate(tunc4),np.concatenate(tunc5),np.concatenate(tunc6),np.concatenate(tunc7),\
    np.concatenate(tunc8),np.concatenate(tunc9),np.concatenate(tunc10),np.concatenate(tunc11),np.concatenate(tunc12),np.concatenate(tunc13),np.concatenate(tunc14),np.concatenate(tunc15),np.concatenate(tunc16)]

total_uncertainty = np.array(total_uncertainty)
total_uncertainty = np.concatenate(total_uncertainty)

plt.close()

# All pulsar uncertainties
fig,ax = plt.subplots(figsize = (15,5))

ax.hist(total_uncertainty,bins=100)
ax.set_title('Total TOA uncertainty distribution for MSP DR')
ax.set_xlabel(r"Uncertainty ($\mu$s)")
ax.set_ylabel('Count')

fig.savefig('All_pulsars_total_uncertainties.pdf',format='pdf',dpi=1000)


fig,ax = plt.subplots(figsize = (15,5))

for tunc_s in tunc_series:
    ax.hist(tunc_s,bins=50)

ax.set_xlabel(r'Uncertainty ($\mu$s)')
ax.set_ylabel('Count')
ax.set_title('Subbanded TOA uncertainty distribution for MSP DR')

ax.legend(['Subband_1','Subband_2','Subband_3','Subband_4','Subband_5','Subband_6','Subband_7','Subband_8','Subband_9','Subband_10','Subband_11','Subband_12','Subband_13','Subband_14','Subband_15','Subband_16'])

fig.savefig('All_pulsars_subbanded_uncertainties.pdf',format='pdf',dpi=1000)



fig = plt.figure(figsize=(15,2*16))

ax = fig.subplots(16,2,sharex=True)
j=0
for i, tunc_s in enumerate(tunc_series):
    ax[i,j].hist(tunc_s,bins=50)

    ax[i,j].set_xlabel(r'Uncertainty ($\mu$s)')
    ax[i,j].set_ylabel('Count')
    ax[i,j].set_title('Centre Frequency: {0}'.format(cfreqs[i]))

    if i ==7:
        j=j+1
    #ax.legend(['Subband_1','Subband_2','Subband_3','Subband_4','Subband_5','Subband_6','Subband_7','Subband_8','Subband_9','Subband_10','Subband_11','Subband_12','Subband_13','Subband_14','Subband_15','Subband_16'])
fig.suptitle('ToA uncertainties per subband')
fig.tight_layout()
fig.savefig('All_pulsars_subbanded_uncertainties_2.pdf',format='pdf',dpi=1000)



#tunc_df = pd.concat(tunc_series)







