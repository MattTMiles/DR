import numpy as np
import os
import glob

path = "/fred/oz002/users/mmiles/MSP_DR/noise_modelling/new_models/"
tnest_dir = "/fred/oz002/users/mmiles/MSP_DR/noise_modelling/new_models/tnest"

psrs = np.loadtxt(os.path.join(path,"psr.list"),dtype=np.str)

pref_models = "/fred/oz002/users/mmiles/MSP_DR/noise_modelling/new_models/preferred_model_ephs"

chosen_list = "/fred/oz002/users/mmiles/MSP_DR/noise_modelling/new_models/chosen_models.list"

result_dirs = ["results_adv1","results_adv2","results_adv3","results_adv4","results_bm","results_dm","results_ecorr","results_red",]

#psrs = ['J0437-4715']

#print("Not copying ephemerides at the moment")
os.chdir(path)
for psr in psrs:
    newpath = os.path.join(tnest_dir,psr)
    os.chdir(newpath)
    results = {}
    for result_dir in result_dirs:
        resultpath = os.path.join(newpath,result_dir)
        os.chdir(resultpath)

        try:
            stats_file = glob.glob("*stats.dat")[0]
            openstats = open(stats_file).readlines()

            evidence = openstats[1].split(":")[1]
            evidence = float(evidence.split("+/-")[0])
            results[result_dir] = evidence
        except:
            print("No evidence, {}".format(results_dir))
            results[result_dir] = 0

    sort_res = sorted(((v,k) for k,v in results.items()))
    sort_res.reverse()
    
    print(psr)
    for number, data in enumerate(sort_res):
        print(number,data)
    
    i = int(input("input number: "))
    chosen = sort_res[i][1]
    
    os.chdir(os.path.join(newpath,chosen))
    par = glob.glob(newpath+"/"+chosen+"/*par")[0]
    os.system("cp "+par+" "+pref_models)

    os.system(" echo {} {} >> {}".format(psr,chosen,chosen_list))



    

        

