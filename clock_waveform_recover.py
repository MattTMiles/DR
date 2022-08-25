# This is only to recover the waveform, it assumes that the covariance matrices and the residuals are already prepared

import numpy as np 
import psrchive
import matplotlib.pyplot as plt
import bilby
import scipy.integrate as integrate
import sys
import os 
import glob


data_dir = sys.argv[1]
results_dir = sys.argv[2]

if not os.path.exists(results_dir):
    try:
        os.makedirs(results_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

os.chdir(data_dir)

for parfile in glob.glob("*par"):
    pulsar = parfile.split(".")[0]

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


    residuals = glob.glob(data_dir+"/*"+pulsar+"_residuals.dat")[0]
    t_res = np.loadtxt(residuals,usecols=2)

    C_clock = np.load(pulsar_cov+"/"+pulsar+"_C_clock.npy")
    Cov_total = np.load(pulsar_cov+"/"+pulsar+"_Cov_total.npy")
    Cov_inv = np.linalg.inv(Cov_total)

    
    clk_waveform = np.dot(C_clock,np.dot(Cov_inv,t_res))

    np.save(pulsar_cov+"/"+pulsar+"_clk_waveform",clk_waveform)

    