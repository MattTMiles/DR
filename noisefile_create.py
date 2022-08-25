#Takes the values of the white noise from the par files and creates .json noisefiles

import os
import glob
import json
import sys
import numpy as np

parfile_list = sys.argv[1]
parfiles = os.popen("cat "+parfile_list).read().split()

outfile_name = sys.argv[2]

red = sys.argv[3]

total_dict = {}

for parfile in parfiles:
    openpar = open(parfile,"r")
    for line in openpar.readlines():
        if line.startswith("PSR"):
            pulsar = line.split(" ")[-1].strip()
            print(pulsar)
        elif line.startswith("#TN"):
            if "GlobalEF" in line:
                print(line)
                total_dict[pulsar+"_KAT_MKBF_efac"] = float(line.split(" ")[-1])
            if "GLobalEQ" in line:
                print(line)
                total_dict[pulsar+"_KAT_MKBF_log10_t2equad"] = float(line.split(" ")[-1])
            if "ECORR" in line:
                total_dict[pulsar+"_KAT_MKBF_log10_ecorr"] = np.log10(1e-6*float(line.split(" ")[-1]))
            if red == "yes":
                if "DMAmp" in line:
                    total_dict[pulsar+"_dm_gp_log10_A"] = float(line.split(" ")[-1])
                if "DMGam" in line:
                    total_dict[pulsar+"_dm_gp_gamma"] = float(line.split(" ")[-1])
                if "RedAmp" in line:
                    total_dict[pulsar+"_KAT_MKBF_red_noise_log10_A"] = float(line.split(" ")[-1])
                if "RedGam" in line:
                    total_dict[pulsar+"_KAT_MKBF_red_noise_gamma"] = float(line.split(" ")[-1])
    
    openpar.seek(0)
    if total_dict.get(pulsar+"_KAT_MKBF_efac") is None:
        total_dict[pulsar+"_KAT_MKBF_efac"] = 0
    if total_dict.get(pulsar+"_KAT_MKBF_log10_t2equad") is None:
        total_dict[pulsar+"_KAT_MKBF_log10_t2equad"] = -9
    if total_dict.get(pulsar+"_KAT_MKBF_log10_ecorr") is None:
        total_dict[pulsar+"_KAT_MKBF_log10_ecorr"] = -9


    




            

with open(outfile_name,"a+") as outfile:
    json.dump(total_dict,outfile,indent=4)
