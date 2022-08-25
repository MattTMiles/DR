import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import pandas as pd

partim_folder = "/fred/oz002/users/mmiles/MSP_DR/noise_subtracted_residuals/f16_residuals/partim"

pulsar_list = open("/fred/oz002/users/mmiles/MSP_DR/pulsar_list.txt","r")

pkl_dir = "/fred/oz002/users/mmiles/MSP_DR/Weff_dicts/"

latex_file_position = "/fred/oz002/users/mmiles/MSP_DR/latex_table.txt"



for pulsar in pulsar_list:
    print(pulsar)
    pulsar = pulsar.strip("\n")
    os.chdir(partim_folder)
    timfile = glob.glob("*"+pulsar+"*tim")[0]
    unc = np.loadtxt(timfile,usecols=3,skiprows=1)
    med_unc = np.median(unc)
    mean_unc = np.mean(unc)
    std_unc = np.std(unc)

    mjds = np.loadtxt(timfile,usecols=2,skiprows=1)
    total_days = mjds[-1] - mjds[0]
    years = total_days/365.25
    
    print("median: {}; mean: {}; std: {}".format(med_unc,mean_unc,std_unc))
    print("years: {}".format(years))

    period = "-"
    dm = "-"
    pb = "-"
    W_ave = "-"

    parfile = glob.glob("*"+pulsar+"*par")[0]
    openpar = open(parfile,"r")
    for line in openpar.readlines():
        if "F0" in line:
            spin_freq = float(line.split()[1])
            period = (1/spin_freq)*1000
            print("period (ms): {}".format(period))
        if "P0" in line:
            period = float(line.split()[1])
            print("period (ms): {}".format(period))
        if "DM " in line and "TNsubtractDM" not in line:
            dm = float(line.split()[1])
            print("DM: {}".format(dm))
        if "PB " in line:
            pb = float(line.split()[1])
            print("P_b: {}".format(pb))

    pkl_file = glob.glob(pkl_dir+"/*"+pulsar+"*dict.pkl")[0]
    pkl = pd.read_pickle(pkl_file)

    W_ave = np.mean(np.array(list(pkl.values())))*1000000
    
    print("W_eff ($\mu$s): {}".format(W_ave))

    
    #Echo the table to a text doc for a quick latex approximation
    
    period = str(round(period,2))
    dm = str(round(dm,2))
    try:
        pb = str(round(pb,2))
    except:
        pb = pb
        
    W_ave = str(round(W_ave,2))
    med_unc = str(round(med_unc,2))
    mean_unc = str(round(mean_unc,2))
    std_unc = str(round(std_unc,2))
    years = str(round(years,2))

    os.system(r"echo '{0} & {1} & {2} & {3} & {4} & {5} & {6} & {7} & {8} \\' >> ".format(pulsar,period,dm,pb,W_ave,med_unc,mean_unc,std_unc,years)+latex_file_position)
