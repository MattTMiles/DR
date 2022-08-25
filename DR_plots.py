from __future__ import division
import numpy as np
from numpy.lib.npyio import loadtxt
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cbook
from matplotlib import ticker
from matplotlib.ticker import MaxNLocator, NullLocator
from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import matplotlib.gridspec as gridspec
import bilby

from scipy import signal
from scipy.signal import convolve, deconvolve
import sys
from pylab import hist, diag
import scipy.integrate as integrate
from scipy.special import gamma, factorial
from astropy.timeseries import LombScargle

import os
import seaborn as sns
import glob

from astropy import units as u
from astropy.coordinates import SkyCoord
import psrchive

from scipy.optimize import curve_fit

import ligo.skymap.plot

from subprocess import call


## Set consistent fonts
font = FontProperties()
font.set_family('serif')
font.set_name('Times New Roman')
font.set_size(30)


### Plot of subbanded uncertainties with subbands beneath ###
def uncertainty_subbands():
    res_folder = "/fred/oz002/users/mmiles/MSP_DR/noise_subtracted_residuals/f16_residuals"
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
    
    cfreqs=[]
    cfreqs = [916,968,1016,1064,1114,1158,1210,1256,1307,1356,1404,1453,1501,1550,1600,1648]

    total_snr = []

    for res in glob.glob('J*dat'):
        print(res)

        mjds = np.loadtxt(res,usecols=0)
        freqs = np.loadtxt(res,usecols=1)/1000000
        resids = np.loadtxt(res,usecols=2)
        unc = np.loadtxt(res,usecols=3)
        #unc = unc*(10**-6)

        total_uncertainty.append(unc)

        pulsar = res.split('_')[0]

        partim = res_folder+"/partim"
        timfile = glob.glob(partim+"/"+pulsar+"*tim")
        snr = np.loadtxt(timfile[0],usecols=-1,skiprows=1)

        total_snr.append(snr)

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
            snrlist[freqindex].append(snr[i])


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

    total_uncertainty = np.concatenate(total_uncertainty)
    total_snr = np.concatenate(total_snr)

    hp_total_uncertainty = [ unc for i, unc in enumerate(total_uncertainty) if total_snr[i] > 15 ]
    hp_total_uncertainty = np.array(hp_total_uncertainty)

    hp_total_uncertainty_20 = [ unc for i, unc in enumerate(total_uncertainty) if total_snr[i] > 20 ]
    hp_total_uncertainty_20 = np.array(hp_total_uncertainty_20)   

    tunc1 = np.concatenate(tunc1)
    tunc2 = np.concatenate(tunc2)
    tunc3 = np.concatenate(tunc3)
    tunc4 = np.concatenate(tunc4)
    tunc5 = np.concatenate(tunc5)
    tunc6 = np.concatenate(tunc6)
    tunc7 = np.concatenate(tunc7)
    tunc8 = np.concatenate(tunc8)
    tunc9 = np.concatenate(tunc9)
    tunc10 = np.concatenate(tunc10)
    tunc11 = np.concatenate(tunc11)
    tunc12 = np.concatenate(tunc12)
    tunc13 = np.concatenate(tunc13)
    tunc14 = np.concatenate(tunc14)
    tunc15 = np.concatenate(tunc15)
    tunc16 = np.concatenate(tunc16)


    fig, axs = plt.subplots(nrows = 9, ncols=1, sharex=True, gridspec_kw={'height_ratios': [4, 1,1,1,1,1,1,1,1], 'hspace':0},figsize=(15,20))
    tag='viridis'

    axs[0].set_xlim(-2,12.95054032)
    axs[0].set_ylim(0,8000)

    sns.histplot(total_uncertainty, bins=100, linewidth=3,fill=False,label="Total ToA uncertainties", ax=axs[0], color=sns.color_palette()[0])
    #sns.histplot(hp_total_uncertainty, bins=100, linewidth=3,fill=False,label="Total ToA uncertainties > snr 15", ax=axs[0], color=sns.color_palette()[1])
    #sns.histplot(hp_total_uncertainty_20, bins=100, linewidth=3,fill=False,label="Total ToA uncertainties > snr 20", ax=axs[0], color=sns.color_palette()[2])
    sns.distplot(tunc16, hist=False, kde=True,kde_kws = {'shade': True, 'linewidth': 3},label=str(cfreqs[15])+" MHz", ax=axs[1],color=sns.color_palette(tag,n_colors=17)[1])
    sns.distplot(tunc15, hist=False, kde=True,kde_kws = {'shade': True, 'linewidth': 3},label=str(cfreqs[14])+" MHz" ,axlabel=False, ax=axs[1],color=sns.color_palette(tag,n_colors=17)[2]).set(xlim=0)
    sns.distplot(tunc14, hist=False, kde=True,kde_kws = {'shade': True, 'linewidth': 3},label=str(cfreqs[13])+" MHz", axlabel=False, ax=axs[2],color=sns.color_palette(tag,n_colors=17)[3]).set(xlim=0)
    sns.distplot(tunc13, hist=False, kde=True,kde_kws = {'shade': True, 'linewidth': 3},label=str(cfreqs[12])+" MHz" ,axlabel=False, ax=axs[2],color=sns.color_palette(tag,n_colors=17)[4]).set(xlim=0)
    sns.distplot(tunc12, hist=False, kde=True,kde_kws = {'shade': True, 'linewidth': 3},label=str(cfreqs[11])+" MHz", axlabel=False, ax=axs[3],color=sns.color_palette(tag,n_colors=17)[5]).set(xlim=0)
    sns.distplot(tunc11, hist=False, kde=True,kde_kws = {'shade': True, 'linewidth': 3},label=str(cfreqs[10])+" MHz", axlabel=False, ax=axs[3],color=sns.color_palette(tag,n_colors=17)[6]).set(xlim=0)
    sns.distplot(tunc10, hist=False, kde=True,kde_kws = {'shade': True, 'linewidth': 3},label=str(cfreqs[9])+" MHz", axlabel=False, ax=axs[4],color=sns.color_palette(tag,n_colors=17)[7]).set(xlim=0)
    sns.distplot(tunc9, hist=False, kde=True,kde_kws = {'shade': True, 'linewidth': 3},label=str(cfreqs[8])+" MHz", axlabel=False, ax=axs[4],color=sns.color_palette(tag,n_colors=17)[8]).set(xlim=0)
    sns.distplot(tunc8, hist=False, kde=True,kde_kws = {'shade': True, 'linewidth': 3},label=str(cfreqs[7])+" MHz", axlabel=False,ax=axs[5],color=sns.color_palette(tag,n_colors=17)[9]).set(xlim=0)
    sns.distplot(tunc7, hist=False, kde=True,kde_kws = {'shade': True, 'linewidth': 3},label=str(cfreqs[6])+" MHz", axlabel=False, ax=axs[5],color=sns.color_palette(tag,n_colors=17)[10]).set(xlim=0)
    sns.distplot(tunc6, hist=False, kde=True,kde_kws = {'shade': True, 'linewidth': 3},label=str(cfreqs[5])+" MHz", axlabel=False, ax=axs[6],color=sns.color_palette(tag,n_colors=17)[11]).set(xlim=0)
    sns.distplot(tunc5, hist=False, kde=True,kde_kws = {'shade': True, 'linewidth': 3},label=str(cfreqs[4])+" MHz", axlabel=False, ax=axs[6],color=sns.color_palette(tag,n_colors=17)[12]).set(xlim=0)
    sns.distplot(tunc4, hist=False, kde=True,kde_kws = {'shade': True, 'linewidth': 3},label=str(cfreqs[3])+" MHz", axlabel=False,ax=axs[7],color=sns.color_palette(tag,n_colors=17)[13]).set(xlim=0)
    sns.distplot(tunc3, hist=False, kde=True,kde_kws = {'shade': True, 'linewidth': 3},label=str(cfreqs[2])+" MHz", axlabel=False,ax=axs[7],color=sns.color_palette(tag,n_colors=17)[14]).set(xlim=0)
    sns.distplot(tunc2, hist=False, kde=True,kde_kws = {'shade': True, 'linewidth': 3},label=str(cfreqs[1])+" MHz", axlabel=False,ax=axs[8],color=sns.color_palette(tag,n_colors=17)[15]).set(xlim=0)
    sns.distplot(tunc1, hist=False, kde=True,kde_kws = {'shade': True, 'linewidth': 3},label=str(cfreqs[0])+" MHz", axlabel=False,ax=axs[8],color=sns.color_palette(tag,n_colors=17)[16]).set(xlim=0)

    font = 24

    for ax in axs[1:]:
        ax.set_ylabel('')
        ax.set_yticks([])
        ax.legend(fontsize=font)

        ax.axvline(linewidth=2, x=-2, color='black')
        ax.axvline(linewidth=2, x=12.95054032, color='black')

    axs[-1].axhline(linewidth=2, y=0, color='black')
    axs[0].axhline(linewidth=2, y=8000, color='black')
    axs[0].axvline(linewidth=2, x=-2, color='black')
    axs[0].axvline(linewidth=2, x=12.95054032, color='black')
    axs[0].tick_params(axis='both',labelsize=font)
    axs[-1].tick_params(axis='both',labelsize=font)


    axs[0].legend(fontsize=font)
    axs[0].set_ylabel("Counts",fontsize=font+3)
    fig.supxlabel(r"ToA Uncertainty ($\mu s$)",fontsize=font+3)
    fig.tight_layout()
    fig.show()
    fig.savefig('/fred/oz002/users/mmiles/MSP_DR/paper_plots/Uncertainty_inc_subband_density.pdf',dpi=1200)

    return total_uncertainty, hp_total_uncertainty, hp_total_uncertainty_20

def position_map(type = "galactic"):
    
    if type == "galactic":
        l_s = []
        b_s = []
        wrms_all = []

        total_data = []

        pulsar_list = "/fred/oz002/users/mmiles/MSP_DR/pulsar_list.txt"

        obs_dir = "/fred/oz002/users/mmiles/MSP_DR/reviewed_data"
        #os.chdir(eph_dir)
        list_read = open(pulsar_list,"r")

        wrms_df = pd.read_pickle("/fred/oz002/users/mmiles/MSP_DR/noise_subtracted_residuals/F_residuals/F_ave_wrms.pkl")

        for pulsar in list_read.readlines():
            pulsar = pulsar.strip("\n")
            print(pulsar)
            try:
                coords =  SkyCoord.from_name("PSR "+pulsar,parse=True)
                coords = coords.galactic
                l = coords.l.value
                b = coords.b.value

                l_s.append(l)
                b_s.append(b)

                wrms = wrms_df[wrms_df.Pulsar == pulsar].wrms.values[0]

                temp_lengths = []
                pulsar_obs_dir = obs_dir+"/"+pulsar+"/corrected_obs"
                os.chdir(pulsar_obs_dir)
                obs_list = glob.glob(pulsar_obs_dir+"/*cln")
                for obs in obs_list:
                    arch = psrchive.Archive_load(obs)
                    arint = arch[0]
                    length = arint.get_duration()
                    temp_lengths.append(length)

                lengths = np.array(temp_lengths)
                median_len = np.median(lengths)

                data = [pulsar,l,b,wrms,median_len]

                total_data.append(data)
                print("done")
            except:
                print(pulsar+" didn't find it")

    if type == "radec":
        ra_s = []
        dec_s = []
        wrms_all = []

        total_data = []

        pulsar_list = "/fred/oz002/users/mmiles/MSP_DR/pulsar_list.txt"

        obs_dir = "/fred/oz002/users/mmiles/MSP_DR/reviewed_data"
        #os.chdir(eph_dir)
        list_read = open(pulsar_list,"r")

        wrms_df = pd.read_pickle("/fred/oz002/users/mmiles/MSP_DR/noise_subtracted_residuals/F_residuals/F_ave_wrms.pkl")

        for pulsar in list_read.readlines():
            pulsar = pulsar.strip("\n")
            print(pulsar)
            try:
                coords =  SkyCoord.from_name("PSR "+pulsar,parse=True)
                coords = coords.icrs
                ra = coords.ra.hour
                dec = coords.dec.degree

                ra_s.append(ra)
                dec_s.append(dec)

                wrms = wrms_df[wrms_df.Pulsar == pulsar].wrms.values[0]

                temp_lengths = []
                pulsar_obs_dir = obs_dir+"/"+pulsar+"/corrected_obs"
                os.chdir(pulsar_obs_dir)
                obs_list = glob.glob(pulsar_obs_dir+"/*cln")
                for obs in obs_list:
                    arch = psrchive.Archive_load(obs)
                    arint = arch[0]
                    length = arint.get_duration()
                    temp_lengths.append(length)

                lengths = np.array(temp_lengths)
                median_len = np.median(lengths)

                data = [pulsar,ra,dec,wrms,median_len]

                total_data.append(data)
                print("done")
            except:
                print(pulsar+" didn't find it")




    master_df = pd.DataFrame(total_data,columns=["Pulsar","ra","dec","wrms","med_len"])
    add1 = {"Pulsar":"J1652-4838","ra":"16.52","dec":"-48.38","wrms":wrms_df[wrms_df.Pulsar == "J1652-4838"].wrms.values[0],"med_len":1135.5033205233638}
    add2 = {"Pulsar":"J2150-0326","ra":"21.50","dec":"-3.26","wrms":wrms_df[wrms_df.Pulsar == "J2150-0326"].wrms.values[0],"med_len":255.57744927102783}

    master_df = master_df.append(add1,ignore_index=True)
    master_df = master_df.append(add2,ignore_index=True)

    gal = SkyCoord(master_df[:76].ra,master_df[:76].dec,frame="galactic",unit=u.deg)
    galadd1 = SkyCoord(master_df[76:].ra,master_df[76:].dec,unit=(u.hourangle,u.deg))
    #gal = SkyCoord(master_df.ra,master_df.dec,unit=(u.hour,u.deg))
    #ra = gal.icrs.ra.hourangle
    radeg = gal.icrs.ra.deg
    radeg_extra = galadd1.icrs.ra.deg
    radeg = np.concatenate((radeg, radeg_extra))
    dec = gal.icrs.dec.deg
    decextra = galadd1.icrs.dec.deg
    dec = np.concatenate((dec,decextra))
    #data = np.array([ra,dec[:76]])
    #temp_df = pd.DataFrame(data.T,columns=["RA hr","DEC deg"])
    #master_df["ra"] = temp_df["RA hr"]
    #master_df["dec"] = temp_df["DEC deg"]

    master_df.ra = master_df.ra.astype(float)
    master_df.dec = master_df.dec.astype(float)

    deg_coords = SkyCoord(master_df.ra*15-180,master_df.dec,unit=u.deg)
    l = deg_coords.galactic.l.deg
    b = deg_coords.galactic.b.deg

    lfix, bfix = np.deg2rad(l).copy(), np.deg2rad(b).copy()
    mask = np.ma.masked_greater_equal(lfix, np.pi).mask
    lfix[mask] = lfix[mask] - 2*np.pi

    rafix, decfix = np.deg2rad(radeg).copy(), np.deg2rad(dec).copy()
    mask = np.ma.masked_greater_equal(rafix, np.pi).mask
    rafix[mask] = rafix[mask] - 2*np.pi

    max_med_len = np.max(master_df.med_len.values)

    font = 20
    cm =sns.diverging_palette(250, 30, l=65, center="dark", as_cmap=True)
    #cm = sns.color_palette("colorblind", as_cmap=True)
    fig = plt.figure(figsize=(15,10))
    
    #fig,ax = plt.subplots()
    #cm = sns.color_palette("pastel1", as_cmap=True)
    ax = fig.add_subplot(111,projection="mollweide")
    #ax = plt.axes(projection="astro hours mollweide")
    plt.grid(True)
    sc=ax.scatter(-1*rafix,decfix, c=master_df.wrms.values, norm=matplotlib.colors.LogNorm(), cmap = cm, s= 2*(master_df.med_len.values/max_med_len)*100)
    axcolor = plt.colorbar(sc,fraction=0.046, pad=0.04)
    ax.set_xlabel("Right Ascension",fontsize=font)
    ax.set_ylabel("Declination",fontsize=font)
    axcolor.set_label("Band-averaged weighted RMS uncertainty ($\mu$s)",fontsize=font)
    axcolor.ax.tick_params(labelsize=font)
    ax.set_xticks(ticks=np.radians([-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150]), labels=['10h', '8h', '6h', '4h', '2h', '0h','22h', '20h', '18h', '16h', '14h'])
    ax.tick_params(axis="both", labelsize=font)
    fig.tight_layout()
    plt.savefig('/fred/oz002/users/mmiles/MSP_DR/paper_plots/skymap_wrms_RA_DEC_FIX.pdf',dpi=1200)
    fig.show()



    return master_df

def noise_removal_example(pulsar):
    # This shows 4 panels of residuals
    # 1. 16 subbanded full residuals
    # 2. 16 subbanded without DM
    # 3. 16 subbanded without DM and Red
    # 4. ave residuals
    
    pulsar = str(pulsar)

    noise_residual_file = "/fred/oz002/users/mmiles/MSP_DR/noise_subtracted_residuals/partim/16ch/"+pulsar+"_16ch_residuals.dat"
    without_DM_file = "/fred/oz002/users/mmiles/MSP_DR/noise_subtracted_residuals/1909_noise_plots/J1909_red_present_residuals.dat"
    fully_whitened_file = "/fred/oz002/users/mmiles/MSP_DR/noise_subtracted_residuals/f16_residuals/"+pulsar+"_16ch_residuals.dat"
    ave_file = "/fred/oz002/users/mmiles/MSP_DR/noise_subtracted_residuals/F_residuals/"+pulsar+"_avg_residuals.dat"

    dm_list = "/fred/oz002/users/mmiles/MSP_DR/noise_modelling/new_models/DM_list.txt"
    red_dm_list = "/fred/oz002/users/mmiles/MSP_DR/noise_modelling/new_models/red_dm_list.txt"
    red_list = "/fred/oz002/users/mmiles/MSP_DR/noise_modelling/new_models/Red_list.txt"

    nres_mjds = np.loadtxt(noise_residual_file,usecols=0)
    nres_resids = np.loadtxt(noise_residual_file,usecols=2)
    nres_unc = np.loadtxt(noise_residual_file,usecols=3)*1e-6

    nres_weights = 1/(nres_unc**2)  
    nres_xbar = np.sum(nres_resids*nres_weights)/np.sum(nres_weights)
    nres_wrms = np.sqrt(np.sum(nres_weights*((nres_resids-nres_xbar)**2))/np.sum(nres_weights))
    nres_wrms_micro = nres_wrms*10**6

    nodm_mjds = np.loadtxt(without_DM_file,usecols=0)
    nodm_resids = np.loadtxt(without_DM_file,usecols=2)
    nodm_unc = np.loadtxt(without_DM_file,usecols=3)*1e-6

    nodm_weights = 1/(nodm_unc**2)  
    nodm_xbar = np.sum(nodm_resids*nodm_weights)/np.sum(nodm_weights)
    nodm_wrms = np.sqrt(np.sum(nodm_weights*((nodm_resids-nodm_xbar)**2))/np.sum(nodm_weights))
    nodm_wrms_micro = nodm_wrms*10**6

    white_mjds = np.loadtxt(fully_whitened_file,usecols=0)
    white_resids = np.loadtxt(fully_whitened_file,usecols=2)
    white_unc = np.loadtxt(fully_whitened_file,usecols=3)*1e-6

    white_weights = 1/(white_unc**2)  
    white_xbar = np.sum(white_resids*white_weights)/np.sum(white_weights)
    white_wrms = np.sqrt(np.sum(white_weights*((white_resids-white_xbar)**2))/np.sum(white_weights))
    white_wrms_micro = white_wrms*10**6

    ave_mjds = np.loadtxt(ave_file,usecols=0)
    ave_resids = np.loadtxt(ave_file,usecols=1)
    ave_unc = np.loadtxt(ave_file,usecols=2)

    ave_mjds,unique_ind = np.unique(ave_mjds,return_index=True)
    ave_resids = ave_resids[unique_ind]
    ave_unc = ave_unc[unique_ind]

    ave_weights = 1/(ave_unc**2)  
    ave_xbar = np.sum(ave_resids*ave_weights)/np.sum(ave_weights)
    ave_wrms = np.sqrt(np.sum(ave_weights*((ave_resids-ave_xbar)**2))/np.sum(ave_weights))
    ave_wrms_micro = ave_wrms*10**6
    
    font = 20
    if pulsar in [line.rstrip() for line in open(red_dm_list)]:
        fig, axs = plt.subplots(nrows = 4, ncols=1, sharex=True, gridspec_kw={'height_ratios': [1,1,1,1], 'hspace':0.1},figsize=(25,15))
        
        font = 20
        
        axs[0].errorbar(nres_mjds,(nres_resids-np.mean(nres_resids))*1e6,yerr=nres_unc*1e6,fmt="o",color=sns.color_palette()[2],ecolor=sns.color_palette()[0],alpha=0.25,label = "$\sigma_\mathrm{}$: {:.3f} $\mu$s".format(r"{rms}",nres_wrms_micro))
        axs[1].errorbar(nodm_mjds,(nodm_resids-np.mean(nodm_resids))*1e6,yerr=nodm_unc*1e6,fmt="o",color=sns.color_palette()[2],ecolor=sns.color_palette()[0],alpha=0.25,label = "$\sigma_\mathrm{}$: {:.3f} $\mu$s".format(r"{rms}",nodm_wrms_micro))
        axs[2].errorbar(white_mjds,(white_resids-np.mean(white_resids))*1e6,yerr=white_unc*1e6,fmt="o",color=sns.color_palette()[2],ecolor=sns.color_palette()[0],alpha=0.25,label = "$\sigma_\mathrm{}$: {:.3f} $\mu$s".format(r"{rms}",white_wrms_micro))
        axs[3].errorbar(ave_mjds,(ave_resids-np.mean(ave_resids))*1e6,yerr=ave_unc*1e6,fmt="o",color=sns.color_palette()[2],ecolor=sns.color_palette()[0],alpha=0.5,label = "$\sigma_\mathrm{}$: {:.3f} $\mu$s".format(r"{rms}",ave_wrms_micro))

        for ax in axs:
            ax.legend(fontsize=font, loc="upper right")
            #ax.ticklabel_format(axis='y', style='sci', scilimits=(-6,-6))
            ax.yaxis.set_major_locator(plt.MaxNLocator(6))
            ax.tick_params(axis='both',labelsize=font)
            ax.set_xlim(ave_mjds[0],ave_mjds[-1])

            ax.axvline(linewidth=2, x=-ave_mjds[0], color='black')
            ax.axvline(linewidth=2, x=ave_mjds[-1], color='black')
            
        #axs[0].set_ylim(-8,8)
        #axs[1].set_ylim(-8,8)
        #axs[2].set_ylim(-8,8)
        #axs[3].set_ylim(-0.8,0.8)

        #axs[0].axhline(-8e-6,8e-6)
        #axs[1].axhline(-8e-6,8e-6)
        #axs[2].axhline(-8e-6,8e-6)
        #axs[3].axhline(-8e-7,8e-7)

        axs[0].set_ylabel("Full \n residuals ($\mu$s)",fontsize=font)
        axs[1].set_ylabel("DM subtracted \n residuals ($\mu$s)",fontsize=font)
        axs[2].set_ylabel("Whitened \n residuals ($\mu$s)",fontsize=font)
        axs[3].set_ylabel("Whitened ave. \n residuals ($\mu$s)",fontsize=font)

        axs[3].set_xlabel("MJD",labelpad=4,fontsize=font)
        #fig.supylabel("Timing residuals (s)")
        
        fig.tight_layout()
        fig.savefig('/fred/oz002/users/mmiles/MSP_DR/paper_plots/pulsar_resids/{}_noise_reduce_steps.pdf'.format(pulsar),dpi=1200)
        #fig.show()

    elif pulsar in [line.rstrip() for line in open(dm_list)]:

        fig, axs = plt.subplots(nrows = 3, ncols=1, sharex=True, gridspec_kw={'height_ratios': [1,1,1], 'hspace':0.1},figsize=(25,15))
        
        #font = 14
        
        axs[0].errorbar(nres_mjds,(nres_resids-np.mean(nres_resids))*1e6,yerr=nres_unc*1e6,fmt="o",color=sns.color_palette()[2],ecolor=sns.color_palette()[0],alpha=0.25,label = "$\sigma_\mathrm{}$: {:.3f} $\mu$s".format(r"{rms}",nres_wrms_micro))
        #axs[1].errorbar(nodm_mjds,(nodm_resids-np.mean(nodm_resids))*1e6,yerr=nodm_unc*1e6,fmt="o",color=sns.color_palette()[2],ecolor=sns.color_palette()[0],alpha=0.25,label = "$\sigma_\mathrm{rms}$: {:.3f} $\mu$s".format(nodm_wrms_micro))
        axs[1].errorbar(white_mjds,(white_resids-np.mean(white_resids))*1e6,yerr=white_unc*1e6,fmt="o",color=sns.color_palette()[2],ecolor=sns.color_palette()[0],alpha=0.25,label = "$\sigma_\mathrm{}$: {:.3f} $\mu$s".format(r"{rms}",white_wrms_micro))
        axs[2].errorbar(ave_mjds,(ave_resids-np.mean(ave_resids))*1e6,yerr=ave_unc*1e6,fmt="o",color=sns.color_palette()[2],ecolor=sns.color_palette()[0],alpha=0.5,label = "$\sigma_\mathrm{}$: {:.3f} $\mu$s".format(r"{rms}",ave_wrms_micro))

        for ax in axs:
            ax.legend(fontsize=font, loc="upper right")
            #ax.ticklabel_format(axis='y', style='sci', scilimits=(-6,-6))
            ax.yaxis.set_major_locator(plt.MaxNLocator(6))
            ax.tick_params(axis='both',labelsize=font)
            ax.set_xlim(ave_mjds[0],ave_mjds[-1])

            ax.axvline(linewidth=2, x=-ave_mjds[0], color='black')
            ax.axvline(linewidth=2, x=ave_mjds[-1], color='black')
            
        #axs[0].set_ylim(-8,8)
        #axs[1].set_ylim(-8,8)
        #axs[1].set_ylim(-8,8)
        #axs[2].set_ylim(-0.8,0.8)

        #axs[0].axhline(-8e-6,8e-6)
        #axs[1].axhline(-8e-6,8e-6)
        #axs[2].axhline(-8e-6,8e-6)
        #axs[3].axhline(-8e-7,8e-7)

        axs[0].set_ylabel("Full \n residuals ($\mu$s)",fontsize=font)
        #axs[1].set_ylabel("DM noise subtracted residuals ($\mu$s)",fontsize=font)
        axs[1].set_ylabel("Whitened \n residuals ($\mu$s)",fontsize=font)
        axs[2].set_ylabel("Whitenned ave. \n residuals ($\mu$s)",fontsize=font)

        axs[2].set_xlabel("MJD",labelpad=4,fontsize=font)
        #fig.supylabel("Timing residuals (s)")
        
        fig.tight_layout()
        fig.savefig('/fred/oz002/users/mmiles/MSP_DR/paper_plots/pulsar_resids/{}_noise_reduce_steps.pdf'.format(pulsar),dpi=1200)
        #fig.show()

    elif pulsar in [line.rstrip() for line in open(red_list)]:

        fig, axs = plt.subplots(nrows = 3, ncols=1, sharex=True, gridspec_kw={'height_ratios': [1,1,1], 'hspace':0.1},figsize=(25,15))
        
        #font = 14
        
        axs[0].errorbar(nres_mjds,(nres_resids-np.mean(nres_resids))*1e6,yerr=nres_unc*1e6,fmt="o",color=sns.color_palette()[2],ecolor=sns.color_palette()[0],alpha=0.25,label = "$\sigma_\mathrm{}$: {:.3f} $\mu$s".format(r"{rms}",nres_wrms_micro))
        #axs[1].errorbar(nodm_mjds,(nodm_resids-np.mean(nodm_resids))*1e6,yerr=nodm_unc*1e6,fmt="o",color=sns.color_palette()[2],ecolor=sns.color_palette()[0],alpha=0.25,label = "$\sigma_\mathrm{rms}$: {:.3f} $\mu$s".format(nodm_wrms_micro))
        axs[1].errorbar(white_mjds,(white_resids-np.mean(white_resids))*1e6,yerr=white_unc*1e6,fmt="o",color=sns.color_palette()[2],ecolor=sns.color_palette()[0],alpha=0.25,label = "$\sigma_\mathrm{}$: {:.3f} $\mu$s".format(r"{rms}",white_wrms_micro))
        axs[2].errorbar(ave_mjds,(ave_resids-np.mean(ave_resids))*1e6,yerr=ave_unc*1e6,fmt="o",color=sns.color_palette()[2],ecolor=sns.color_palette()[0],alpha=0.5,label = "$\sigma_\mathrm{}$: {:.3f} $\mu$s".format(r"{rms}",ave_wrms_micro))

        for ax in axs:
            ax.legend(fontsize=font, loc="upper right")
            #ax.ticklabel_format(axis='y', style='sci', scilimits=(-6,-6))
            ax.yaxis.set_major_locator(plt.MaxNLocator(6))
            ax.tick_params(axis='both',labelsize=font)
            ax.set_xlim(ave_mjds[0],ave_mjds[-1])

            ax.axvline(linewidth=2, x=-ave_mjds[0], color='black')
            ax.axvline(linewidth=2, x=ave_mjds[-1], color='black')
            
        #axs[0].set_ylim(-8,8)
        #axs[1].set_ylim(-8,8)
        #axs[1].set_ylim(-8,8)
        #axs[2].set_ylim(-0.8,0.8)

        #axs[0].axhline(-8e-6,8e-6)
        #axs[1].axhline(-8e-6,8e-6)
        #axs[2].axhline(-8e-6,8e-6)
        #axs[3].axhline(-8e-7,8e-7)

        axs[0].set_ylabel("Full \n residuals ($\mu$s)",fontsize=font)
        #axs[1].set_ylabel("DM noise subtracted residuals ($\mu$s)",fontsize=font)
        axs[1].set_ylabel("Whitened \n residuals ($\mu$s)",fontsize=font)
        axs[2].set_ylabel("Whitened ave. \n residuals ($\mu$s)",fontsize=font)

        axs[2].set_xlabel("MJD",labelpad=4,fontsize=font)
        #fig.supylabel("Timing residuals (s)")
        
        fig.tight_layout()
        fig.savefig('/fred/oz002/users/mmiles/MSP_DR/paper_plots/pulsar_resids/{}_noise_reduce_steps.pdf'.format(pulsar),dpi=1200)
        #fig.show()
    
    else:
        fig, axs = plt.subplots(nrows = 2, ncols=1, sharex=True, gridspec_kw={'height_ratios': [1,1], 'hspace':0.1},figsize=(25,15))
        
        #font = 14
        
        axs[0].errorbar(nres_mjds,(nres_resids-np.mean(nres_resids))*1e6,yerr=nres_unc*1e6,fmt="o",color=sns.color_palette()[2],ecolor=sns.color_palette()[0],alpha=0.25,label = "$\sigma_\mathrm{}$: {:.3f} $\mu$s".format(r"{rms}",nres_wrms_micro))
        #axs[1].errorbar(nodm_mjds,(nodm_resids-np.mean(nodm_resids))*1e6,yerr=nodm_unc*1e6,fmt="o",color=sns.color_palette()[2],ecolor=sns.color_palette()[0],alpha=0.25,label = "$\sigma_\mathrm{rms}$: {:.3f} $\mu$s".format(nodm_wrms_micro))
        #axs[1].errorbar(white_mjds,(white_resids-np.mean(white_resids))*1e6,yerr=white_unc*1e6,fmt="o",color=sns.color_palette()[2],ecolor=sns.color_palette()[0],alpha=0.25,label = "$\sigma_\mathrm{}$: {:.3f} $\mu$s".format(r"{rms}",white_wrms_micro))
        axs[1].errorbar(ave_mjds,(ave_resids-np.mean(ave_resids))*1e6,yerr=ave_unc*1e6,fmt="o",color=sns.color_palette()[2],ecolor=sns.color_palette()[0],alpha=0.5,label = "$\sigma_\mathrm{}$: {:.3f} $\mu$s".format(r"{rms}",ave_wrms_micro))

        for ax in axs:
            ax.legend(fontsize=font, loc="upper right")
            #ax.ticklabel_format(axis='y', style='sci', scilimits=(-6,-6))
            ax.yaxis.set_major_locator(plt.MaxNLocator(6))
            ax.tick_params(axis='both',labelsize=font)
            ax.set_xlim(ave_mjds[0],ave_mjds[-1])

            ax.axvline(linewidth=2, x=-ave_mjds[0], color='black')
            ax.axvline(linewidth=2, x=ave_mjds[-1], color='black')
            
        #axs[0].set_ylim(-8,8)
        #axs[1].set_ylim(-8,8)
        #axs[1].set_ylim(-8,8)
        #axs[2].set_ylim(-0.8,0.8)

        #axs[0].axhline(-8e-6,8e-6)
        #axs[1].axhline(-8e-6,8e-6)
        #axs[2].axhline(-8e-6,8e-6)
        #axs[3].axhline(-8e-7,8e-7)

        axs[0].set_ylabel("Full \n residuals ($\mu$s)",fontsize=font)
        #axs[1].set_ylabel("DM noise subtracted residuals ($\mu$s)",fontsize=font)
        #axs[1].set_ylabel("Whitened residuals ($\mu$s)",fontsize=font)
        axs[1].set_ylabel("Ave. \n residuals ($\mu$s)",fontsize=font)

        axs[1].set_xlabel("MJD",labelpad=4,fontsize=font)
        #fig.supylabel("Timing residuals (s)")
        
        fig.tight_layout()
        fig.savefig('/fred/oz002/users/mmiles/MSP_DR/paper_plots/pulsar_resids/{}_noise_reduce_steps.pdf'.format(pulsar),dpi=1200)
        #fig.show()


def effective_width_unc(select=1):

    if select == 1:
        font = 24

        res_folder = '/fred/oz002/users/mmiles/MSP_DR/notebooks/notebook_residuals'

        weff_dir = '/fred/oz002/users/mmiles/MSP_DR/Weff_dicts/'

        tim_folder = '/fred/oz002/users/mmiles/MSP_DR/notebooks/notebook_tims/'

        os.chdir(weff_dir)

        master_dict = pd.read_pickle("/fred/oz002/users/mmiles/MSP_DR/master_weff_dict.pkl")
        
        fig = plt.figure(figsize=(15,10))
        gs = fig.add_gridspec(2, 9)
        ax1 = fig.add_subplot(gs[1, :9])
        axs = [fig.add_subplot(gs[0, :4]), fig.add_subplot(gs[0, 4:8])]

        for ax in fig.get_axes():
            ax.tick_params(axis='both',labelsize=font-2)
        
        #
        #plt.subplot(2,1,2)
        master_dict.plot.scatter('weighted_eff_width', 'uncertainty', c='freq', colormap='viridis',ax=ax1)
        cax = fig.get_axes()[-1].set_ylabel("Frequency (MHz)",fontsize=font)
        fig.get_axes()[-1].tick_params(axis='both',labelsize=font-2)

        ax1.set_xlabel("$\mathrm{W}_{\mathrm{eff}}$ weighted by S/N (s)",fontsize=font)
        ax1.set_ylabel("ToA uncertainty ($\mu$s)",fontsize=font)

        ax1.set_yscale("log")
        ax1.set_xscale("log")

        for p, pulsar in enumerate(["J1909-3744","J1327-0755"]):

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

            def parabola(x,a,b,c):

                return a*(x**2) + b*x + c

            pars, cov = curve_fit(f=power_law, xdata=eff_width_df['freq'], ydata=eff_width_df['eff_width'])
            pars2, cov2 = curve_fit(f=parabola, xdata=eff_width_df['freq'], ydata=eff_width_df['eff_width'])


            weff_fit = power_law(eff_width_df['freq'], *pars)
            
            weff_fit2 = parabola(eff_width_df['freq'], *pars2)


            if p == 0:
                axs[p].scatter(eff_width_df['freq'],eff_width_df['eff_width']*1000000,color = sns.color_palette("colorblind")[2], label="Effective width")
                axs[p].plot(eff_width_df['freq'],weff_fit*1000000, linestyle = '--',color='k',label = 'Power-law')
            else:
                axs[p].scatter(eff_width_df['freq'],eff_width_df['eff_width']*1000000,color = sns.color_palette("colorblind")[2], label="Effective width")
                axs[p].plot(eff_width_df['freq'],weff_fit2*1000000, linestyle = '--',color='k',label = 'Parabola')

            #axs[p].legend(fontsize=font-12)
            axs[p].set_ylabel('$\mathrm{W}_{\mathrm{eff}}$ ($\mu$s)',fontsize=font)
            axs[p].set_xlabel(r'Frequency (MHz)',fontsize=font)

            axs[p].set_title(pulsar,fontsize=font)
            
            fig.tight_layout()
            fig.savefig('/fred/oz002/users/mmiles/MSP_DR/paper_plots/eff_width_unc.pdf',dpi=50)
            fig.show()

    if select == 2:
        font = 24

        res_folder = '/fred/oz002/users/mmiles/MSP_DR/notebooks/notebook_residuals'

        weff_dir = '/fred/oz002/users/mmiles/MSP_DR/Weff_dicts/'

        tim_folder = '/fred/oz002/users/mmiles/MSP_DR/notebooks/notebook_tims/'

        os.chdir(weff_dir)

        master_dict = pd.read_pickle("/fred/oz002/users/mmiles/MSP_DR/master_weff_dict.pkl")
        
        fig,axs = plt.subplots(nrows=1,ncols=2,figsize=(15,10))

        for ax in axs:
            ax.tick_params(axis='both',labelsize=font-2)
        

        for p, pulsar in enumerate(["J1909-3744","J1327-0755"]):

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

            def parabola(x,a,b,c):

                return a*(x**2) + b*x + c

            pars, cov = curve_fit(f=power_law, xdata=eff_width_df['freq'], ydata=eff_width_df['eff_width'])
            pars2, cov2 = curve_fit(f=parabola, xdata=eff_width_df['freq'], ydata=eff_width_df['eff_width'])


            weff_fit = power_law(eff_width_df['freq'], *pars)
            
            weff_fit2 = parabola(eff_width_df['freq'], *pars2)


            if p == 0:
                axs[p].scatter(eff_width_df['freq'],eff_width_df['eff_width']*1000000,color = sns.color_palette("colorblind")[2], label="Effective width")
                axs[p].plot(eff_width_df['freq'],weff_fit*1000000, linestyle = '--',color='k',label = 'Power-law')
            else:
                axs[p].scatter(eff_width_df['freq'],eff_width_df['eff_width']*1000000,color = sns.color_palette("colorblind")[2], label="Effective width")
                axs[p].plot(eff_width_df['freq'],weff_fit2*1000000, linestyle = '--',color='k',label = 'Parabola')

            #axs[p].legend(fontsize=font-12)
            axs[p].set_ylabel('$\mathrm{W}_{\mathrm{eff}}$ ($\mu$s)',fontsize=font)
            axs[p].set_xlabel(r'Frequency (MHz)',fontsize=font)

            axs[p].set_title(pulsar,fontsize=font)
            
            fig.tight_layout()
            fig.savefig('/fred/oz002/users/mmiles/MSP_DR/paper_plots/eff_width_unc.pdf',dpi=50)
            fig.show()

def wrms_hists():

    pulsar_list = "/fred/oz002/users/mmiles/MSP_DR/pulsar_list.txt"

    list_read = open(pulsar_list,"r")

    total_data = []

    font = 24

    for pulsar in list_read:

        pulsar = pulsar.strip("\n")

        fully_whitened_file = "/fred/oz002/users/mmiles/MSP_DR/noise_subtracted_residuals/f16_residuals/"+pulsar+"_16ch_residuals.dat"
        ave_file = "/fred/oz002/users/mmiles/MSP_DR/noise_subtracted_residuals/F_residuals/"+pulsar+"_avg_residuals.dat"
        
        
        white_mjds = np.loadtxt(fully_whitened_file,usecols=0)
        white_resids = np.loadtxt(fully_whitened_file,usecols=2)
        white_unc = np.loadtxt(fully_whitened_file,usecols=3)*1e-6

        white_weights = 1/(white_unc**2)  
        white_xbar = np.sum(white_resids*white_weights)/np.sum(white_weights)
        white_wrms = np.sqrt(np.sum(white_weights*((white_resids-white_xbar)**2))/np.sum(white_weights))
        white_wrms_micro = white_wrms*10**6

        ave_mjds = np.loadtxt(ave_file,usecols=0)
        ave_resids = np.loadtxt(ave_file,usecols=1)
        ave_unc = np.loadtxt(ave_file,usecols=2)

        ave_mjds,unique_ind = np.unique(ave_mjds,return_index=True)
        ave_resids = ave_resids[unique_ind]
        ave_unc = ave_unc[unique_ind]

        ave_weights = 1/(ave_unc**2)  
        ave_xbar = np.sum(ave_resids*ave_weights)/np.sum(ave_weights)
        ave_wrms = np.sqrt(np.sum(ave_weights*((ave_resids-ave_xbar)**2))/np.sum(ave_weights))
        ave_wrms_micro = ave_wrms*10**6

        data = [pulsar, white_wrms_micro, ave_wrms_micro]

        total_data.append(data)

    df = pd.DataFrame(total_data,columns=["Pulsar","subband","ave"])

    sns.histplot(df.subband, bins=30,label="Sub-banded wRMS",color=sns.color_palette("colorblind")[2])
    sns.histplot(df.ave, bins=30,label="Band-averaged wRMS",color=sns.color_palette("colorblind")[0])

    plt.axvline(1,linestyle="--",color="black",label="1 $\mu$s")

    plt.legend(fontsize=font)
    plt.xlabel("RMS ($\mu$s)",fontsize=font)
    plt.ylabel("Pulsar count",fontsize=font)
    plt.show()


def clock_recover_img():

    npt = 500
    mjdpt = np.zeros(npt)
    
    mk2utc = np.loadtxt('/fred/oz002/rshannon/tempo2/clock/mk2utc.clk',skiprows=7)
    mk2utc_mjds = mk2utc[:,0]
    global_mjds = mk2utc_mjds[:-1]

    mjdmin = np.min(global_mjds)
    mjdmax = np.max(global_mjds)
    for i in range(npt):
        mjdpt[i] = mjdmin + ((mjdmax-mjdmin)/npt)*(i+0.5)

    master_clk_to_save = np.load("/fred/oz002/users/mmiles/MSP_DR/clock_correction_work/enterprise_version/nanograv_clones/clocks/best_10_data_updated_pars/results_all_extended/clk_sim_waveform.npy")
    master_error_to_save = np.load("/fred/oz002/users/mmiles/MSP_DR/clock_correction_work/enterprise_version/nanograv_clones/clocks/clk_sim_error.npy")

    mk2utc = np.loadtxt('/fred/oz002/users/mmiles/MSP_DR/clock_correction_work/enterprise_version/nanograv_clones/clocks/mk2utc.clk',skiprows=8)
    mk2utc_error = 4.925*1e-9

    arg_ind = np.argwhere(mjdpt>mk2utc[-1,0])[0]
    subset_clk = master_clk_to_save[:arg_ind[0]]
    mean_to_take = np.mean(subset_clk)

    font = 20
    fig, axs = plt.subplots(2, 1, sharex=True,gridspec_kw={'height_ratios': [3, 1]})
    fig.subplots_adjust(hspace=0)
    axs[0].plot(mjdpt,(master_clk_to_save-mean_to_take)*1e6, label='Recovered Clock Signal',color='xkcd:green')
    axs[0].fill_between(mjdpt, ((master_clk_to_save - master_error_to_save) - mean_to_take)*1e6, ((master_clk_to_save + master_error_to_save) - mean_to_take)*1e6,alpha=0.2,color="xkcd:green")
    axs[0].plot(mk2utc[:,0],((-mk2utc[:,1]) - np.mean(-mk2utc[:,1]))*1e6,label='mk2utc', color='tab:blue')
    axs[0].fill_between(mk2utc[:,0], ((-mk2utc[:,1] - mk2utc_error) - np.mean(-mk2utc[:,1]))*1e6, ((-mk2utc[:,1] + mk2utc_error) - np.mean(-mk2utc[:,1]))*1e6, alpha=0.2, color='tab:blue')
    #axs[1].scatter(master_mjds_active, master_res_active*1e6, marker='x', color='black', alpha=0.2)
    #plt.fill_between(mjdpt, (master_clk_to_save - placehold_error) - np.mean(master_clk_to_save), (master_clk_to_save + placehold_error) - np.mean(master_clk_to_save),alpha=0.2,color="xkcd:green")
    #plt.fill_between(mjdpt, (master_clk_to_save - master_error_to_save) - np.mean(master_clk_to_save), (master_clk_to_save + master_error_to_save) - np.mean(master_clk_to_save),alpha=0.2,color="xkcd:green")
    #plt.plot(mk2utc[:,0],(-mk2utc[:,1]) - np.mean(-mk2utc[:,1]),label='mk2utc')
    axs[0].set_title("Recovered clock signal",fontsize=font)
    #axs[1].ticklabel_format(axis="y",style="sci",scilimits=(-6,-6))
    axs[1].yaxis.offsetText.set_fontsize(16)
    axs[0].yaxis.offsetText.set_fontsize(16)
    axs[1].tick_params(axis='both', which='major', labelsize=font)
    #axs[0].ticklabel_format(axis="y",style="sci",scilimits=(-6,-6))
    axs[0].tick_params(axis='both', which='major', labelsize=font)
    axs[1].set_xlabel("MJD",fontsize=font)
    axs[1].set_ylabel("Residuals ($\mu$s)",fontsize=font)
    axs[0].set_ylabel("Clock signals ($\mu$s)",fontsize=font)
    fig.legend(fontsize=font)
    fig.tight_layout()
    #ig.savefig('/fred/oz002/users/mmiles/MSP_DR/paper_plots/clock_recovery.pdf',dpi=150)
    fig.show()


def DM_compare_plots():
    
    pulsar_list = "/fred/oz002/users/mmiles/MSP_DR/pulsar_list.txt"

    pref_eph_dir = "/fred/oz002/users/mmiles/MSP_DR/noise_modelling/new_models/preferred_model_ephs"
    
    model_list = "/fred/oz002/users/mmiles/MSP_DR/noise_modelling/new_models/chosen_models.list"
    tnest_dirs = "/fred/oz002/users/mmiles/MSP_DR/noise_modelling/new_models/tnest"


    total_data = []

    for pulsar in open(pulsar_list,"r").readlines():
        pulsar = pulsar.strip("\n")

        pulsar_eph = glob.glob(pref_eph_dir+"/*"+pulsar+"*")[0]
        
        openpar = open(pulsar_eph,"r")
        parcontent = openpar.read()
        openpar.seek(0)
        if "TNDMAmp" in parcontent:
            parlines = openpar.readlines()
            for line in parlines:
                if line.startswith("DM "):
                    dm = float(line.split()[1])
                #if line.startswith("TNDMAmp"):
                #    Adm = float(line.split()[-1].strip("\n"))

            openmodel = open(model_list,"r")
            openmodel.seek(0)
            for line in openmodel.readlines():
                if pulsar in line:
                    model = line.split()[-1]
            
            chosen_model = model
            models_pulsar = os.path.join(tnest_dirs,pulsar)
            cmodel_dir = os.path.join(models_pulsar,chosen_model)

            posts = np.loadtxt(cmodel_dir+"/MTMSP-"+pulsar+"-post_equal_weights.dat")


            if chosen_model == "results_ecorr":
                efacposts = posts[:,0]
                equadposts = posts[:,1]
                ecorrposts = posts[:,2]
            
            if chosen_model == "results_red":
                efacposts = posts[:,0]
                equadposts = posts[:,1]
                redampposts = posts[:,2]
                redslopeposts = posts[:,3]
            
            if chosen_model == "results_dm":
                efacposts = posts[:,0]
                equadposts = posts[:,1]
                dmampposts = posts[:,2]
                dmslopeposts = posts[:,3]

            if chosen_model == "results_adv1":
                efacposts = posts[:,0]
                equadposts = posts[:,1]
                ecorrposts = posts[:,2]
                redampposts = posts[:,3]
                redslopeposts = posts[:,4]

            if chosen_model == "results_adv2":
                efacposts = posts[:,0]
                equadposts = posts[:,1]
                ecorrposts = posts[:,2]
                dmampposts = posts[:,3]
                dmslopeposts = posts[:,4]

            if chosen_model == "results_adv3":
                efacposts = posts[:,0]
                equadposts = posts[:,1]
                ecorrposts = posts[:,2]
                redampposts = posts[:,3]
                redslopeposts = posts[:,4]
                dmampposts = posts[:,5]
                dmslopeposts = posts[:,6]

            if chosen_model == "results_adv4":
                efacposts = posts[:,0]
                equadposts = posts[:,1]
                redampposts = posts[:,2]
                redslopeposts = posts[:,3]
                dmampposts = posts[:,4]
                dmslopeposts = posts[:,5]

            Adm = np.median(dmampposts)
            Aunc = np.std(dmampposts)

            slope = np.median(dmslopeposts)
            slope_unc = np.std(dmslopeposts)

                
            data = [pulsar,dm,Adm,Aunc,slope,slope_unc]
            total_data.append(data)
    
    df = pd.DataFrame(total_data,columns=["Pulsar","DM","A_DM","A_unc","slope","slope_unc"])
    
    fig,ax = plt.subplots()
    sns.regplot(x=df.DM,y=df.A_DM)
    ax.errorbar(df.DM,df.A_DM,yerr=df.A_unc,fmt="None",capsize=5,zorder=1,color="C0")
    plt.show()

    return df
            

def example_plots(pulsars):
    for pulsar in [pulsars]:
        
        notebook_dir = '/fred/oz002/users/mmiles/MSP_DR/notebooks'

        #pulsar = 'J1909-3744'

        pulsar_dir = '/fred/oz002/users/mmiles/MSP_DR/subband_comps/highsnr_psradd_check/'+pulsar

        og_datafile = pulsar_dir+'/grand.dly'

        std = '/fred/oz002/users/mmiles/MSP_DR/github_templates/'+pulsar+'.std'   

        os.chdir(pulsar_dir)

        os.system("python ~/soft/timing/subband_check.py -portrait {0}/2D.{1}.notebook_version.ar -archive {2}/grand.dly -ephemeris {2}/{1}.par -subbands 16 -save -trim".format(notebook_dir,pulsar,pulsar_dir))
        #~/soft/timing/subband_check.py -portrait {notebook_dir}/2D.{pulsar}.notebook_version.ar -archive {pulsar_dir}/grand.dly -ephemeris {pulsar_dir}/{pulsar}.par -subbands 16 -save
        #plt.savefig('/fred/oz002/users/mmiles/MSP_DR/paper_plots/{0}_trimmed_notebook_comparison.pdf'.format(pulsar),format='pdf',dpi=1000)












'''
def appendix_plots():
    #This script will create plots for the appendix
    #One pdf per pulsar
    #Includes grand averaged profile, subbanded residuals, DM noise taken out, all noise taken out
'''



    

    




    






            


