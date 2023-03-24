from tabulate import tabulate
import numpy as np 
import json
import glob

partim = "/fred/oz002/users/mmiles/MPTA_GW/partim/"

tables = []
for par in sorted(glob.glob(partim+"/*par")):
    parfile = list(open(par).readlines())

    freq = "*"
    freq_unc = "*"
    freq_dot = "*"
    freq_dot_unc = "*"
    dm = "*"
    dm_unc = "*"
    dm_dot = "*"
    dm_dot_unc = "*"
    
    for parline in parfile:
        if parline.startswith("PSRJ"):
            jname = parline.split()[1]
            #print(jname)
        if parline.startswith("F0"):
            freq = float(parline.split()[1])
            freq_unc = float(parline.split()[3])
        if parline.startswith("F1"):
            freq_dot = float(parline.split()[1])
            freq_dot_unc = float(parline.split()[3])
        if parline.startswith("DM "):
            dm = float(parline.split()[1])
            dm_unc = float(parline.split()[3])
        if parline.startswith("DM1"):
            dm_dot = float(parline.split()[1])
            dm_dot_unc = float(parline.split()[3])

    #print(tabulate([jname, freq+" "+freq_unc, freq_dot+" "+freq_dot_unc, dm+" "+dm_unc, dm_dot+" "+dm_dot_unc], \
    #["PSRJ", "F0\n(Hz)","F1\n(s^-2)","DM\n(cm^-3 pc)","DM1\n(cm^-3 pc/yr)"], tablefmt="grid"))

    table = [jname, freq, freq_unc, freq_dot, freq_dot_unc, dm, dm_unc, dm_dot, dm_dot_unc]
    tables.append(table)

with open("MPTA_psrcat.txt", "w") as f:
    print(tabulate(tables, headers=["PSRJ", "F0 (Hz)"," ","F1 (s^-2)"," ","DM (cm^-3 pc)"," ","DM1 (cm^-3 pc/yr)"," "], tablefmt="outline"), file=f)



            
