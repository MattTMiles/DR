import numpy as np
import os
import subprocess as sproc 

noise_dir = "/fred/oz005/users/aparthas/MeerTime_MSP_NoiseAnalysis/partim"
results_dir = "/fred/oz005/users/aparthas/MeerTime_MSP_NoiseAnalysis/partim/tnest"
model_bin = "/fred/oz002/users/mmiles/MSP_DR/bestmodels"

bestmodels_list = "/fred/oz002/users/mmiles/MSP_DR/bestmodels.list"

bestmodels = open(bestmodels_list,"r")

for line in bestmodels:
    line = line.split()
    pulsar = line[0]
    model = line[1]
    
    pulsar_dir = os.path.join(results_dir,pulsar)
    model_dir = os.path.join(pulsar_dir,model)
    os.system("cp "+model_dir+"/*par "+model_bin)

    #if len(line) >2:
    #    model2 = line[2]
    #    model2_dir = os.path.join(pulsar_dir,model2)
    #    os.system("cp "+model2_dir+"/*par "+model_bin)

    
    