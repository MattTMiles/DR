import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import pandas as pd
import math

partim_folder = "/fred/oz002/users/mmiles/MSP_DR/noise_subtracted_residuals/f16_residuals/partim"

pulsar_list = open("/fred/oz002/users/mmiles/MSP_DR/pulsar_list.txt","r")

model_list = "/fred/oz002/users/mmiles/MSP_DR/noise_modelling/new_models/chosen_models.list"

tnest_dirs = "/fred/oz002/users/mmiles/MSP_DR/noise_modelling/new_models/tnest"

latex_file_position = "/fred/oz002/users/mmiles/MSP_DR/latex_results.txt"

f16_residuals = "/fred/oz002/users/mmiles/MSP_DR/noise_subtracted_residuals/f16_residuals"
F_residuals = "/fred/oz002/users/mmiles/MSP_DR/noise_subtracted_residuals/F_residuals"

full_f16 = "/fred/oz002/users/mmiles/MSP_DR/noise_subtracted_residuals/partim/16ch"
full_F = "/fred/oz002/users/mmiles/MSP_DR/noise_subtracted_residuals/partim/ave"

os.chdir("/fred/oz002/users/mmiles/MSP_DR/noise_subtracted_residuals")

#pulsar_list = ["J1909-3744"]

for pulsar in pulsar_list:
    print(pulsar)
    pulsar = pulsar.strip("\n")
    os.chdir(partim_folder)
    timfile = glob.glob("*"+pulsar+"*tim")[0]
    unc = np.loadtxt(timfile,usecols=3,skiprows=1)

    mjds = np.loadtxt(timfile,usecols=2,skiprows=1)
    num_toas = len(mjds)

    fd1 = "-"
    fd2 = "-"
    dm2 = "-"
    redamp = "-"
    redamp_unc=""
    dmamp = "-"
    dmamp_unc=""
    dmgam = "-"
    dmgam_unc=""
    redgam = "-"
    redgam_unc=""

    band_ave_white = "-"
    subband_white = "-"

    subband_full = "-"
    band_ave_full = "-"

    red_chi_white = "-"

    parfile = glob.glob("*"+pulsar+"*par")[0]
    openpar = open(parfile,"r")
    for line in openpar.readlines():
        if "FD1" in line:
            fd1 = "$\checkmark$"
        if "FD2" in line:
            fd2 = "$\checkmark$"
        if "DM2" in line:
            dm2 = "$\checkmark$"
        if "TNDMAmp" in line:
            dmamp = float(line.split()[-1])
        if "TNDMGam" in line:
            dmgam = float(line.split()[-1])
        if "TNRedAmp" in line:
            redamp = float(line.split()[-1])
        if "TNRedGam" in line:
            redgam = float(line.split()[-1])


    openpar.seek(0)
    # whitened residuals
    if any("TNDMAmp" in line for line in openpar) or any("TNRedAmp" in line for line in openpar):
        os.chdir(F_residuals)
        for res in glob.glob(pulsar+"*dat"):
            unc = np.loadtxt(res,usecols=2)
            mjds = np.loadtxt(res,usecols=0)
            mjds_unique,unique_ind = np.unique(mjds,return_index=True)

            resids = np.loadtxt(res,usecols=1)
            resids_unique = resids[unique_ind]
            unc_unique = unc[unique_ind] 
            weights = 1/(unc_unique**2) 
            xbar = np.sum(resids_unique*weights)/np.sum(weights)
            wrms = np.sqrt(np.sum(weights*((resids_unique-xbar)**2))/np.sum(weights))
            band_ave_white = wrms*10**6

            openpar.seek(0)
            fitted=0
            for line in openpar.readlines():
                split=line.split() 
                if len(split) > 2:
                    try:
                        if float(split[2])==1:
                            fitted=fitted+1
                    except ValueError:
                        continue

                        
            chi_square_white = np.sum(((resids_unique-np.mean(resids_unique))/(unc_unique))**2)
            nu_white = len(resids_unique)-fitted
            red_chi_white = chi_square_white/nu_white

        os.chdir(f16_residuals)
        for res in glob.glob(pulsar+"*dat"):
            unc = np.loadtxt(res,usecols=3)
            unc = unc*(10**-6)
            resids = np.loadtxt(res,usecols=2)
            weights = 1/(unc**2)  
            xbar = np.sum(resids*weights)/np.sum(weights)
            wrms = np.sqrt(np.sum(weights*((resids-xbar)**2))/np.sum(weights))
            subband_white = wrms*10**6

    #Full residuals
    os.chdir(full_F)
    for res in glob.glob(pulsar+"*dat"):
        unc = np.loadtxt(res,usecols=2)
        mjds = np.loadtxt(res,usecols=0)
        mjds_unique,unique_ind = np.unique(mjds,return_index=True)

        resids = np.loadtxt(res,usecols=1)
        resids_unique = resids[unique_ind]
        unc_unique = unc[unique_ind] 
        weights = 1/(unc_unique**2) 
        xbar = np.sum(resids_unique*weights)/np.sum(weights)
        wrms = np.sqrt(np.sum(weights*((resids_unique-xbar)**2))/np.sum(weights))
        band_ave_full = wrms*10**6

        openpar.seek(0)
        fitted=0
        for line in openpar.readlines():
            split=line.split() 
            if len(split) > 2:
                try:
                    if float(split[2])==1:
                        fitted=fitted+1
                except ValueError:
                    continue

                    
        chi_square_full = np.sum(((resids_unique-np.mean(resids_unique))/(unc_unique))**2)
        nu_full = len(resids_unique)-fitted
        red_chi_full = chi_square_full/nu_full
        #os.system("echo {} {} >> /fred/oz002/users/mmiles/MSP_DR/noise_subtracted_residuals/full_red_chi.txt".format(pulsar,red_chi))
    
    os.chdir(full_f16)
    for res in glob.glob(pulsar+"*dat"):
        unc = np.loadtxt(res,usecols=3)
        unc = unc*(10**-6)
        resids = np.loadtxt(res,usecols=2)
        weights = 1/(unc**2)  
        xbar = np.sum(resids*weights)/np.sum(weights)
        wrms = np.sqrt(np.sum(weights*((resids-xbar)**2))/np.sum(weights))
        subband_full = wrms*10**6

    
    #Echo the table to a text doc for a quick latex approximation
        
    num_toas = str(round(num_toas))
    fd1 = str(fd1)
    fd2 = str(fd2)
    dm2 = str(dm2)

    openmodel = open(model_list,"r")
    openmodel.seek(0)
    for line in openmodel.readlines():
        if pulsar in line:
            model = line.split()[-1]
    

    model_bin = ["results_bm","results_ecorr","results_red","results_dm","results_adv1","results_adv2","results_adv3","results_adv4"]
    model_num = model_bin.index(model) +1

    chosen_model = model
    models_pulsar = os.path.join(tnest_dirs,pulsar)
    cmodel_dir = os.path.join(models_pulsar,chosen_model)

    posts = np.loadtxt(cmodel_dir+"/MTMSP-"+pulsar+"-post_equal_weights.dat")

    if chosen_model == "results_bm":
        efacposts = posts[:,0]
        equadposts = posts[:,1]
    
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

    def round_on_error(value,error):
        significant_digits = 10**math.floor(math.log(error, 10))
        return value // significant_digits * significant_digits

    def return_first_decimal(error):
        significant_digits = 10**math.floor(math.log(error, 10))
        return round(error/significant_digits)
    
    if redamp != "-":
        redamp = np.median(redampposts)
        redamp_unc = np.std(redampposts)
        redamp = round_on_error(redamp,redamp_unc)
        redamp_unc = return_first_decimal(redamp_unc)
        redamp_unc = "("+str(redamp_unc)+")"
        #redamp = str(round(redamp,3))
        redamp = str(redamp)
    if redgam != "-":
        redgam = np.median(redslopeposts)
        redgam_unc = np.std(redslopeposts)
        redgam = round_on_error(redgam,redgam_unc)
        redgam_unc = return_first_decimal(redgam_unc)
        redgam_unc = "("+str(redgam_unc)+")"
        #redgam = str(round(redgam,3))
        redgam = str(redgam)
    if dmamp != "-":
        dmamp = np.median(dmampposts)
        dmamp_unc = np.std(dmampposts)
        dmamp = round_on_error(dmamp, dmamp_unc)
        dmamp_unc = return_first_decimal(dmamp_unc)
        dmamp_unc = "("+str(dmamp_unc)+")"
        #dmamp = str(round(dmamp,3))
        dmamp = str(dmamp)
    if dmgam != "-":
        dmgam=np.median(dmslopeposts)
        dmgam_unc = np.std(dmslopeposts)
        dmgam = round_on_error(dmgam, dmgam_unc)
        dmgam_unc = return_first_decimal(dmgam_unc)
        dmgam_unc = "("+str(dmgam_unc)+")"
        #dmgam = str(round(dmgam,3))
        dmgam = str(dmgam)
    

    subband_full = str(round(subband_full,3))
    if subband_white != "-":
        subband_white = str(round(subband_white,3))
    band_ave_full = str(round(band_ave_full,3))
    if band_ave_white != "-":
        band_ave_white = str(round(band_ave_white,3))
    
    red_chi_full = str(round(red_chi_full,3))
    if red_chi_white != "-":
        red_chi_white = str(round(red_chi_white,3))
    
    print("full reduced chi: {}".format(red_chi_full))
    print("whitened reduced chi: {} \n".format(red_chi_white))

    os.system(r"echo '{0} & {1} & {2} & {3}{4} & {5}{6} & {7}{8} & {9}{10} & {11} & {12} & {13} & {14} & {15} & {16} \\' >> ".format(pulsar,num_toas,model_num,dmamp,dmamp_unc,dmgam,dmgam_unc,redamp,redamp_unc,redgam,redgam_unc,subband_full,subband_white,band_ave_full,red_chi_full,band_ave_white,red_chi_white)+latex_file_position)
