import numpy as np 
import subprocess as sproc 
import os
import glob

os.chdir("/fred/oz002/users/mmiles/MSP_DR/github_ephs")

for par in sorted(glob.glob("J*par")):
    
    try:
        os.chdir("rm new.par")
    except:
        FileNotFoundError()

    pulsar = par.split(".")[0]
    os.system("tempo2 -f "+par+" ../notebooks/notebook_tims/*"+pulsar+"* -fit DM1 -newpar")
    dm1_line = os.popen("cat new.par | grep DM1").read().split()
    dm1_val, dm1_unc = dm1_line[1], dm1_line[-1]
    
    if abs(float(dm1_val)/float(dm1_unc)) >= 3:
        os.system("mv new.par "+par)

        os.system("tempo2 -f "+par+" ../notebooks/notebook_tims/*"+pulsar+"* -fit DM2 -newpar")
        dm2_line = os.popen("cat new.par | grep DM2").read().split()
        dm2_val, dm2_unc = dm2_line[1], dm2_line[-1]

        if abs(float(dm2_val)/float(dm2_unc)) >= 3:
            os.system("mv new.par "+par)

    try:
        os.chdir("rm new.par")
    except:
        FileNotFoundError()

    os.system("tempo2 -f "+par+" ../notebooks/notebook_tims/*"+pulsar+"* -fit FD1 -newpar")
    fd1_line = os.popen("cat new.par | grep FD1").read().split()
    fd1_val, fd1_unc = fd1_line[1], fd1_line[-1]
    
    if abs(float(fd1_val)/float(fd1_unc)) >= 3:
        os.system("mv new.par "+par)

        os.system("tempo2 -f "+par+" ../notebooks/notebook_tims/*"+pulsar+"* -fit FD2 -newpar")
        fd2_line = os.popen("cat new.par | grep FD2").read().split()
        fd2_val, fd2_unc = fd2_line[1], fd2_line[-1]

        if abs(float(fd2_val)/float(fd2_unc)) >= 3:
            os.system("mv new.par "+par)

    os.system("cp "+par+" /fred/oz002/users/mmiles/MSP_DR/all_lights_ephs/")
    




