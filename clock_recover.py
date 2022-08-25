# The purpose of this script is to estimate the clock corrections from pulsar residuals where the clock corrections are turned off
# This will return the difference between the estimated clock corrections, the difference between the clock correction file and the estimations, and maybe more

# Comparing to BIPM2019
# Using DE438
# Final method will need to be using noise subtracted residuals, for now (12/04/2022) we will make this for normal residuals until the noise subtracted are ready

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

# Run with python ~/soft/DR/clock_recovery.py <data_dir> <results_dir>

#Assuming the parent directories of the tim, par, and res are:

parser = argparse.ArgumentParser(description="Subband comparison")
parser.add_argument("-data", dest="data", help="data directory to use", required = True)
parser.add_argument("-results", dest="results", help="result directory to use", required = True)
parser.add_argument("-solo", dest="solo", help="(Optional) Adding this flag will make this work on a singular pulsar", required=False)
parser.add_argument("-residuals", dest="residuals", help="(Optional) Adding this flag will make this work on a particular residual file (Not working yet)", required=False)
parser.add_argument("-extended",dest="extended", action="store_true", help="(Optional) Adding this flag result in a more advanced nested matrix method", required=False)
args = parser.parse_args()

data_dir = str(args.data)
results_dir = str(args.results)
solo = str(args.solo)
residuals = str(args.residuals)
ext = args.extended

def pwr_spec_clock(f,clk_amp,clk_gamma):
    #These are filled in based on noise results

    A = clk_amp 
    gamma = clk_gamma 
    f_c = 1 #years
    return ((10**A)**2/(12*(np.pi**2)))*((f/f_c)**(-gamma))

def pwr_spec_clock2(f,clk_gamma):
    #These are filled in based on noise results

    #A = clk_amp 
    gamma = clk_gamma 
    f_c = 1 #years
    return ((f/f_c)**(-gamma))

def pwr_spec_spin(f,spin_amp,spin_gamma):
    #These are filled in based on noise results

    A = spin_amp 
    gamma = spin_gamma 
    f_c = 1 #years
    return ((10**A)**2/(12*(np.pi**2)))*((f/f_c)**(-gamma))

def pwr_spec_dm(f,dm_amp,dm_gamma):
    #These are filled in based on noise results

    A = dm_amp 
    gamma = dm_gamma 
    f_c = 1 #years
    return ((10**A)**2)*((f/f_c)**(-gamma))

def pwr_spec_dm2(f,dm_gamma):
    #These are filled in based on noise results

    #A = dm_amp 
    gamma = dm_gamma 
    f_c = 1 #years
    return ((f/f_c)**(-gamma))

# Specify the paramaters of the simulated signal

npt = 500
mjdpt = np.zeros(npt)

#Use 1909 for the global MJDs. We know it already works so should be fine

#timfile_1909 = glob.glob("/fred/oz002/users/mmiles/MSP_DR/clock_correction_work/enterprise_version/nanograv_clones/clocks/best_10_data/*J1909-3744*.tim")[0]
resfile_1909 = "/fred/oz002/users/mmiles/MSP_DR/clock_correction_work/enterprise_version/nanograv_clones/clocks/ryans_residuals/J1909-3744_residuals.dat"
mk2utc = np.loadtxt('/fred/oz002/rshannon/tempo2/clock/mk2utc.clk',skiprows=7)
mk2utc_mjds = mk2utc[:,0]
global_mjds = mk2utc_mjds[:-1]


mjdmin = np.min(global_mjds)
mjdmax = np.max(global_mjds)

for i in range(npt):
    mjdpt[i] = mjdmin + ((mjdmax-mjdmin)/npt)*(i+0.5)

#define observing span
Tdays = mjdmax - mjdmin

nday=int(Tdays+1)
# Convert days to years
T = Tdays/365.25

#Start making directories if they don't exist 

if not os.path.exists(results_dir):
    try:
        os.makedirs(results_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

total_df = []

def master_function(pulsar,residuals=residuals):
    pulsar_dir = results_dir+"/"+pulsar

    if not os.path.exists(pulsar_dir):
        try:
            os.makedirs(pulsar_dir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    pulsar_cov = pulsar_dir+"/covariance/"

    if not os.path.exists(pulsar_cov):
        try:
            os.makedirs(pulsar_cov)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise


    os.system("cd "+pulsar_cov)


    #timfile = parent_tim+"*"+pulsar+"*.tim"
    #timfile = glob.glob(data_dir+"/*"+pulsar+"*.tim")[0]
    #timfile = sys.argv[2]
    

    

    parfile = glob.glob(data_dir+"/*"+pulsar+"*.par")[0]
    #parfile = sys.argv[3]
    if residuals=="None":
        residuals = glob.glob(data_dir+"/*"+pulsar+"_residuals.dat")[0]


    #residuals = sys.argv[4]
    t_res = np.loadtxt(residuals,usecols=2)
    freqs = np.loadtxt(residuals,usecols=1)
    mjds = np.loadtxt(residuals,usecols=0)
    uncs = np.loadtxt(residuals,usecols=3)

    num_toas = len(mjds)

    mjd_1 = mjds[0]
    mjd_end = mjds[-1]

    #VALUES FROM ENTERPRISE
    clk_amp = -13.189181408304862
    clk_gamma = 2.4067997669710457

    mjd_dummy=np.zeros(nday)

    for i in range(len(mjd_dummy)):
        mjd_dummy[i] = i
    
    #Initialize the 1D covariance matrices
    
    cov_clock = np.zeros(nday)
    cov_spin = np.zeros(nday)
    cov_dm = np.zeros(nday)
    cov_ecorr = np.zeros(nday)
    cov_white = np.zeros(nday)


    dm_amp = 0
    dm_gamma = 0
    red_amp = 0
    red_gamma = 0
    ecorr = 0
    equad = 0
    efac = 0

    openpar = open(parfile,"r")
    openpar.seek(0)
    
    for line in openpar.readlines():
        if "TNDMAmp" in line:
            dm_amp = float(line.split(" ")[-1])
        if "TNDMGam" in line: 
            dm_gamma = float(line.split(" ")[-1])
        if "TNRedAmp" in line:
            #red_amp = float(line.split(" ")[-1])
            red_amp = 0
        if "TNRedGam" in line:
            #red_gamma = float(line.split(" ")[-1])
            red_gamma = 0
        if "TNECORR" in line:
            ecorr = float(line.split(" ")[-1])
            ecorr = (1e-6*ecorr)/(365.25*86400)
        if "TNSECORR" in line:
            ecorr = float(line.split(" ")[-1])
            ecorr = ecorr/np.sqrt(256/3600)
            ecorr = (1e-6*ecorr)/(365.25*86400)
        if "TNGLobalEQ" in line or "TNGlobalEQ" in line:
            equad = float(line.split(" ")[-1])
        if "TNGlobalEF" in line:
            efac = float(line.split(" ")[-1])


    for i in range(len(mjd_dummy)):
        
        K = 2.410*10**-16

        #time step in years
        time_step = (mjd_dummy[i])/365.25
        cov_dm[i] = ((1/(K**2))*((10**dm_amp)**2)*integrate.quad(lambda f: \
            (pwr_spec_dm2(f,dm_gamma)*np.cos(2*np.pi*f*time_step))\
                ,1/T,365.25/28)[0])

        cov_clock[i] = ((10**clk_amp)**2/(12*(np.pi**2)))*integrate.quad(lambda f: \
                (pwr_spec_clock2(f,clk_gamma)*np.cos(2*np.pi*f*time_step))\
                    ,1/T,365.25/28)[0]

        cov_spin[i] = ((10**red_amp)**2/(12*(np.pi**2)))*integrate.quad(lambda f: \
            (pwr_spec_clock2(f,red_gamma)*np.cos(2*np.pi*f*time_step))\
                ,1/T,365.25/28)[0]

    
    #Initialize the covariance matrices
    

    C_clock = np.zeros((num_toas,num_toas))
    C_spin = np.zeros((num_toas,num_toas))
    C_dm = np.zeros((num_toas,num_toas))
    C_ecorr = np.zeros((num_toas,num_toas))
    C_white = np.zeros((num_toas,num_toas))
        
    for i in range(num_toas):
        #print("{}".format(i))
        #equad = -8
        C_white[i,i] =(((efac)*((uncs[i]*10**-6)/(365.25*86400)))**2) + (((10**(equad))/(365.25*86400))**2)
        
        for j in range(num_toas):
            
            interp_i = abs(int(mjds[j]-mjds[i]))
            
            C_spin[i,j] = 0
            
            #C_spin[i,j] = integrate.quad(lambda f: \
            #    (pwr_spec_spin(f,red_amp,red_gamma)*np.cos(2*np.pi*f*time_step))\
            #        ,1/T,num_toas/T)[0]

            if dm_amp != 0:
                C_dm[i,j] = cov_dm[interp_i]/((freqs[i]**2)*(freqs[j]**2))
            else:
                C_dm[i,j] = 0

            
            if abs(mjds[i]-mjds[j]) < 0.001:
                C_ecorr[i,j] = (ecorr)**2
                #C_ecorr[i,j] = 0
            else:
                C_ecorr[i,j] = 0


            C_clock[i,j] = cov_clock[interp_i]

    np.save(pulsar_cov+"/"+pulsar+"_C_clock",C_clock)
    np.save(pulsar_cov+"/"+pulsar+"_C_spin",C_spin)
    np.save(pulsar_cov+"/"+pulsar+"_C_dm",C_dm)
    np.save(pulsar_cov+"/"+pulsar+"_C_ecorr",C_ecorr)
    np.save(pulsar_cov+"/"+pulsar+"_C_white",C_white)

    print(pulsar)
    print("C_clock: {}; C_spin: {}; C_dm: {}; C_ecorr: {}; C_white: {}".format(C_clock[0][0],C_spin[0][0],C_dm[0][0],C_ecorr[0][0],C_white[0][0]))
    #Cov_total = C_clock + C_spin + C_dm + C_ecorr + C_white 
    Cov_total = C_clock + C_spin + C_ecorr + C_white 
    Cov_total_no_clock = C_spin + C_dm + C_ecorr + C_white 

    np.save(pulsar_cov+"/"+pulsar+"_Cov_total",Cov_total)

    Cov_inv = np.linalg.inv(Cov_total)

    ngrand = num_toas + npt
    Cov_simulated = np.zeros((ngrand,ngrand))

    for i in range(num_toas):
        for j in range(num_toas):
            Cov_simulated[i][j] = Cov_total_no_clock[i][j]
    
    Cov_clock_simulated = np.zeros((ngrand,ngrand))

    for i in range(ngrand):
        for j in range(ngrand):

            if i < num_toas:
                mjd1 = mjds[i]
            else:
                mjd1 = mjdpt[i-num_toas]
            
            if j < num_toas:
                mjd2 = mjds[j]
            else:
                mjd2 = mjdpt[j-num_toas]
            
            imjd = abs(int(mjd2-mjd1))

            Cov_clock_simulated[i][j] = cov_clock[imjd]

    Cov_simulated = Cov_simulated + Cov_clock_simulated
    
    Cov_inv_simulated = np.linalg.inv(Cov_simulated)

    t_res_sim = np.zeros(ngrand)

    for i in range(num_toas):
        t_res_sim[i] = t_res[i]


    error_temp = np.matmul(Cov_inv_simulated,Cov_clock_simulated)
    sim_clock_matrix = Cov_clock_simulated - np.matmul(Cov_clock_simulated,error_temp)
    sim_clock_error = np.diagonal(sim_clock_matrix)
    sim_clock_seconds = np.sqrt(sim_clock_error)*(365.25*86400)

    #Okay so this should be the clock waveform
    #if 
    temp_waveform = np.matmul(Cov_inv,t_res)
    clk_waveform = np.matmul(C_clock,temp_waveform)

    np.save(pulsar_cov+"/"+pulsar+"_clk_waveform",clk_waveform)

        
    temp_sim = np.matmul(Cov_inv_simulated,t_res_sim)
    clk_sim = np.matmul(Cov_clock_simulated,temp_sim)

    clk_sim_tosave = clk_sim[num_toas:]
    print("length is: {}".format(len(clk_sim_tosave)))
    
    np.save(pulsar_cov+"/"+pulsar+"_clk_sim_waveform",clk_sim_tosave)

    clk_error_tosave = sim_clock_seconds[num_toas:]
    np.save(pulsar_cov+"/"+pulsar+"_clk_sim_error_array",clk_error_tosave)

    
    mk2utc = np.loadtxt('/fred/oz002/users/mmiles/MSP_DR/clock_correction_work/enterprise_version/nanograv_clones/clocks/mk2utc.clk',skiprows=7)

    plt.plot(mjds,clk_waveform,label='reconstructed_clock_{}'.format(pulsar))
    plt.plot(mjdpt,clk_sim[num_toas:],label='simulated_clock_{}'.format(pulsar),color='xkcd:green')
    plt.fill_between(mjdpt, clk_sim[num_toas:] - sim_clock_seconds[num_toas:], clk_sim[num_toas:] + sim_clock_seconds[num_toas:],alpha=0.2,color="xkcd:green")
    plt.plot(mk2utc[:,0],-mk2utc[:,1],label='mk2utc')
    plt.title("How {} sees the clock".format(pulsar))
    plt.legend()
    plt.savefig(pulsar_cov+"/"+pulsar+"_clock_vs_MKT")
    plt.close()

    return Cov_total, cov_clock, t_res, mjds

if ext==False:
    if args.solo is None:
        os.chdir(data_dir)
        for par in glob.glob("J*par"):
            pulsar = par.split(".")[0]
            print(pulsar)

            master_function(pulsar,residuals)

    if args.solo is not None:
        os.chdir(data_dir)
        pulsar = solo
            #pulsar = par.split(".")[0]
        print(pulsar)
        
        master_function(pulsar,residuals)

if ext==True:
    if args.solo is not None:
        print("Extended is only for where there are multiple pulsars. It won't do anything if there's only one.")

    if args.solo is None:
        os.chdir(data_dir)
        num_pulsars = len(glob.glob("J*par"))

        npt = 500
        
        master_cov = np.zeros((num_pulsars,num_pulsars),dtype=object)
        master_cov_clock = np.zeros((num_pulsars,num_pulsars),dtype=object)
        master_res = np.zeros(num_pulsars,dtype=object)
        master_mjds = np.zeros(num_pulsars,dtype=object)

        #npt = 500 
        total_length = []
        for i, par in enumerate(glob.glob("J*par")):
            pulsar = par.split(".")[0]
            print(pulsar)

            Cov_total, cov_clock, t_res, mjds = master_function(pulsar,residuals)
            #plt.show()
            
            #Cov_total_in = np.zeros((2837, 2837),dtype=np.float64)
            #Cov_total_in[:Cov_total.shape[0],:Cov_total.shape[1]] = Cov_total
            master_cov[i][i] = Cov_total
            master_cov_clock = cov_clock
            master_res[i] = t_res
            length = len(t_res)
            total_length.append(length)
            master_mjds[i] = mjds
        
        total_length = np.sum(total_length)
        total_length_npt = total_length+npt
        
        #master_cov_inv = np.linalg.inv(master_cov)

        master_cov_inv_sim = np.zeros((total_length_npt,total_length_npt),dtype=np.float64)
        master_cov_sim = np.zeros((total_length_npt,total_length_npt),dtype=np.float64)
        master_cov_clock_sim = np.zeros((total_length_npt,total_length_npt),dtype=np.float64)
        eval_master_cov = np.zeros((total_length,total_length),dtype=np.float64)
        #raise a
        ntoa_master = 0
        
        for i in range(len(glob.glob("J*par"))):
            ntoa_int = len(master_res[i])
            
            print("Pulsar {}".format(i+1))
            for j in range(ntoa_int):
                for k in range(ntoa_int):
                    #To account for the 0's not being of matrix form do this try/except loop
                    #try:
                    #    master_cov_inv_sim[j+ntoa_master][k+ntoa_master] = master_cov_inv[i][i][j][k]
                    #except TypeError:
                    #    master_cov_inv_sim[j+ntoa_master][k+ntoa_master] = 0
                    try:
                        eval_master_cov[j+ntoa_master][k+ntoa_master] = master_cov[i][i][j][k]
                    except TypeError:
                        eval_master_cov[j+ntoa_master][k+ntoa_master] = 0
                    
            ntoa_master = ntoa_master + ntoa_int

        Cov_clock = np.zeros((total_length_npt,total_length_npt))

        master_mjds_active = np.hstack(master_mjds)
        master_res_active = np.hstack(master_res)

        for i in range(ntoa_master):
            for j in range(ntoa_master):

                interp_i = abs(int(master_mjds_active[j]-master_mjds_active[i]))

                Cov_clock[i,j] = master_cov_clock[interp_i]


        print("Creating master simulated clock")
        for i in range(total_length_npt):
            for j in range(total_length_npt):

                if i < len(master_res_active):
                    mjd1 = master_mjds_active[i]
                else:
                    mjd1 = mjdpt[i-len(master_res_active)]
                
                if j < len(master_res_active):
                    mjd2 = master_mjds_active[j]
                else:
                    mjd2 = mjdpt[j-len(master_res_active)]
                
                imjd = abs(int(mjd2-mjd1))
                master_cov_clock_sim[i][j] = master_cov_clock[imjd]
        
        #master_cov_with_clock = eval_master_cov + Cov_clock

        master_cov_inv = np.linalg.inv(eval_master_cov)

        for i in range(total_length):
            for j in range(total_length):
                master_cov_inv_sim[i][j] = master_cov_inv[i][j]

        #master_cov_total_sim = master_cov_sim + master_cov_clock_sim

        #master_cov_inv_sim = np.linalg.inv(master_cov_total_sim)
        
        master_t_res_sim = np.zeros(total_length_npt,dtype=np.float64)


        for i in range(ntoa_master):
            master_t_res_sim[i] = master_res_active[i]

        print("Creating master waveform")
        master_temp_waveform_sim = np.matmul(master_cov_inv_sim,master_t_res_sim)
        master_clk_sim = np.matmul(master_cov_clock_sim,master_temp_waveform_sim)

        master_clk_to_save = master_clk_sim[ntoa_master:]
        np.save(results_dir+"/MEGAmatrix_best10_clk_sim_waveform",master_clk_to_save)
        
        
        print("Creating error array")
        master_error_temp = np.matmul(master_cov_inv_sim,master_cov_clock_sim)
        master_sim_clock_matrix = master_cov_clock_sim - np.matmul(master_cov_clock_sim,master_error_temp)
        master_sim_clock_error = np.diagonal(master_sim_clock_matrix)
        master_sim_clock_seconds = np.sqrt(np.abs(master_sim_clock_error))*(365.25*86400)

        master_error_to_save = master_sim_clock_seconds[ntoa_master:]
        #np.save("/fred/oz002/users/mmiles/MSP_DR/clock_correction_work/enterprise_version/nanograv_clones/clocks/MEGAmatrix_best10_clk_sim_error",master_error_to_save)
        
        #print("Using 1909 error placeholder")
        #placehold_error = np.load('/fred/oz002/users/mmiles/MSP_DR/clock_correction_work/enterprise_version/nanograv_clones/clocks/best_10_updated_pars_results/J1909-3744/covariance/J1909-3744_clk_sim_error_array.npy')

        mk2utc = np.loadtxt('/fred/oz002/users/mmiles/MSP_DR/clock_correction_work/enterprise_version/nanograv_clones/clocks/mk2utc.clk',skiprows=7)
        plt.plot(mjdpt,master_clk_to_save-np.mean(master_clk_to_save),label='simulated_clock_MegaMatrix',color='xkcd:green')
        #plt.fill_between(mjdpt, (master_clk_to_save - placehold_error) - np.mean(master_clk_to_save), (master_clk_to_save + placehold_error) - np.mean(master_clk_to_save),alpha=0.2,color="xkcd:green")
        plt.fill_between(mjdpt, (master_clk_to_save - master_error_to_save) - np.mean(master_clk_to_save), (master_clk_to_save + master_error_to_save) - np.mean(master_clk_to_save),alpha=0.2,color="xkcd:green")
        plt.plot(mk2utc[:,0],(-mk2utc[:,1]) - np.mean(-mk2utc[:,1]),label='mk2utc')
        plt.title("How MegaMatrix sees the clock")
        plt.legend()
        plt.show()


        
        

