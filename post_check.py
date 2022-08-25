import numpy as np
import os
import glob
import math
import corner
import matplotlib.pyplot as plt

pulsar_list = open("/fred/oz002/users/mmiles/MSP_DR/pulsar_list.txt","r")

model_list = "/fred/oz002/users/mmiles/MSP_DR/noise_modelling/new_models/chosen_models.list"

tnest_dirs = "/fred/oz002/users/mmiles/MSP_DR/noise_modelling/new_models/tnest"

parent_dir = "/fred/oz002/users/mmiles/MSP_DR/noise_modelling/new_models"

post_dir = "/fred/oz002/users/mmiles/MSP_DR/noise_modelling/new_models/post_pdfs"

os.chdir(parent_dir)

#pulsar_list = ["J1909-3744"]

for pulsar in pulsar_list:
    print(pulsar)
    pulsar = pulsar.strip("\n")

    openmodel = open(model_list,"r")
    openmodel.seek(0)
    for line in openmodel.readlines():
        if pulsar in line:
            model = line.split()[-1]
    pulsar_dir = os.path.join(tnest_dirs,pulsar)
    
    os.chdir(pulsar_dir)
    model_dir = os.path.join(pulsar_dir,model)
    os.chdir(model_dir)

    paramfile = "MTMSP-"+pulsar+"-.paramnames"
    posts = np.loadtxt("MTMSP-"+pulsar+"-post_equal_weights.dat")

    with open(paramfile,"r") as params:
        lines = params.readlines()

    lines = [ x.split()[-1].strip("\n") for x in lines]
    
    figure = corner.corner(posts[:,:-1], labels=lines,
                       quantiles=[0.16, 0.5, 0.84],
                       show_titles=True, title_kwargs={"fontsize": 12})
    
    figure.suptitle(pulsar,x=0.8)

    plt.savefig(post_dir+"/"+pulsar+"_posts.pdf")
    plt.close()
