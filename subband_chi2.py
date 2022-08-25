# Per frequency band code to check how for residual RFI in the offpulse regions of the pulsar
# Run in a directory with J*ar files


import matplotlib.pyplot as plt 
from matplotlib.pyplot import text
import numpy as np
import pandas as pd
import sys
import os
import subprocess as sproc 
from scipy.stats import chisquare
from scipy.optimize import minimize, curve_fit
import psrchive
import filecmp
import glob
import argparse

parser = argparse.ArgumentParser(description="Subband reduced chi2")
parser.add_argument("-portrait", dest="portrait", help="Portrait to use", required = True)
parser.add_argument("-ephemeris", dest="ephemeris", help="Ephemeris file to use", required = True)
parser.add_argument("-subbands", dest="subbands", help="Specify the number of subbands to compare", required = True)

args = parser.parse_args()

portrait = str(args.portrait)
eph = str(args.ephemeris)

subbands = int(args.subbands)

pulsar = eph.split('J')[-1]
pulsar = 'J'+pulsar
try:
    pulsar = pulsar.split('.')[0]
except:
    pass
#Load in the portrait
ppfile = psrchive.Archive_load(portrait)

#Define standard for shifting profiles later
std = ppfile
std.fscrunch()
std.pscrunch()
std.dedisperse()
std.remove_baseline()

#Dedisperse portrait file
p = sproc.Popen('pam -mD {0}'.format(portrait),shell=True)
p.wait()

#Scrunch to required number of subbands
p = sproc.Popen('pam -e {0}chan --setnchn {0} {1}'.format(subbands, portrait),shell=True)
p.wait()

port_trim = os.path.splitext(portrait)[0]
port_sub = os.path.splitext(portrait)[0]+'.'+str(subbands)+'chan'


if not os.path.isfile('grand.paz'):

    for files in glob.glob("J*ar"):
        if not os.path.isfile(files.split('.ar')[0]+'.paz'):
            p = sproc.Popen('paz -r '+files+' -e paz', shell=True)
            p.wait()
    

    #Initialise all the J*paz so that they only update when they see a new file

    p = sproc.Popen('echo "`ls J*paz`" > to_add',shell=True)
    p.wait()

    if os.path.isfile('added'):
        if not filecmp.cmp('to_add','added'):
            p = sproc.Popen('pam -mE {0} J*paz'.format(eph),shell=True)
            p.wait()

            p = sproc.Popen('psradd -P J*paz -o grand.paz',shell=True)
            p.wait()

            p = sproc.Popen('pam -mp grand.paz',shell=True)
            p.wait()

            p = sproc.Popen('pam -mD grand.paz',shell=True)
            p.wait()

            p = sproc.Popen('pam -e {0}chan --setnchn {0} grand.paz'.format(subbands),shell=True)
            p.wait()

            p = sproc.Popen('echo "`ls J*paz`" > added',shell=True)
            p.wait()

            p = sproc.Popen('psrsplit grand.{0}chan -c 1'.format(subbands),shell=True)
            p.wait()

    else:
        p = sproc.Popen('pam -mE {0} J*paz'.format(eph),shell=True)
        p.wait()

        p = sproc.Popen('psradd -P J*paz -o grand.paz',shell=True)
        p.wait()

        p = sproc.Popen('pam -mp grand.paz',shell=True)
        p.wait()

        p = sproc.Popen('pam -mD grand.paz',shell=True)
        p.wait()

        p = sproc.Popen('pam -e {0}chan --setnchn {0} grand.paz'.format(subbands),shell=True)
        p.wait()

        p = sproc.Popen('echo "`ls J*paz`" > added',shell=True)
        p.wait()

        p = sproc.Popen('psrsplit grand.{0}chan -c 1'.format(subbands),shell=True)
        p.wait()

sprof = std.get_Profile(0,0,0)
#grandfile should be full time resolution archive with x amount of subbands depending on how many 

grandfile = psrchive.Archive_load('grand.{}chan'.format(subbands)) 
grandfile.remove_baseline()
grandfile.dedisperse()

ppfile.remove_baseline()
ppfile.dedisperse()

if not os.path.isfile('{0}.0000_0000.{1}chan'.format(port_trim,subbands)):
    p = sproc.Popen('psrsplit {0} -c 1'.format(port_sub),shell=True)
    p.wait()

#Define the profile shift off of the standard 
ssf = psrchive.ProfileShiftFit()
ssf.set_standard(sprof)

granddata = grandfile.get_data()
ppdata = ppfile.get_data()

#initialise the reduced chi buckets
total_redchi = []
for i in range(subbands):
    locals()['offpulse_redchi2_subband{0}'.format(i)] = []
    locals()['total_redchi2_subband{0}'.format(i)] = []
    locals()['RFI_ratio_subband{0}'.format(i)] = []


vars_all = []
redchis_all = []
offpulse_redchis_all = []
freqs = []

#Subint should be of the shape (1,subints,1024)
for j, subint in enumerate(grandfile):
    for i in range(0,subbands):
        
        print('Subint: {0}; Subband: {1}'.format(j,i))

        locals()['grandfile{0}'.format(i)] = psrchive.Archive_load('grand.000{0}_0000.{1}chan'.format(i,subbands))
        activefile = locals()['grandfile{0}'.format(i)]
        activefile.remove_baseline()

        freq = activefile.get_centre_frequency()
        freqs.append(freq)
        
        locals()['gf{0}_prof'.format(i)] = activefile[j].get_Profile(0,0)
        active_prof = locals()['gf{0}_prof'.format(i)]

        locals()['ppfile{0}'.format(i)] = psrchive.Archive_load('{0}.000{1}_0000.{2}chan'.format(port_trim,i,subbands))
        ppfile = locals()['ppfile{0}'.format(i)]
        ppfile.remove_baseline()
        
        locals()['ppf{0}_prof'.format(i)] = ppfile[0].get_Profile(0,0)
        ppf_prof = locals()['ppf{0}_prof'.format(i)]

        ssf.apply_scale_and_shift(ppf_prof)
        psf = psrchive.ProfileShiftFit()

        psf.set_standard(ppf_prof)
        psf.set_Profile(active_prof)

        scale = psf.get_scale()[0]

        psf.apply_scale_and_shift(active_prof)


        locals()['ppdata{0}'.format(i)] = ppf_prof.get_amps()
        ppdata = locals()['ppdata{0}'.format(i)]
        locals()['granddata{0}'.format(i)] = active_prof.get_amps()
        granddata = locals()['granddata{0}'.format(i)]

        '''def func(temp_profile, amp, phs):
            return amp*fft_rotate(temp_profile, phs)

        amplitude = 1
        phase =  0

        locals()['ppdata_fitted{0}'.format(i)] = func(ppdata,amplitude,phase)
        ppdata_fitted = locals()['ppdata_fitted{0}'.format(i)]'''

        locals()['offpulse_{0}'.format(i)] = np.where(ppdata<0.01*ppdata.max())[0]
        offpulse_int = locals()['offpulse_{0}'.format(i)]

        locals()['offgrand_{0}'.format(i)] = granddata[offpulse_int]
        offgrand = locals()['offgrand_{0}'.format(i)]

        locals()['offport_{0}'.format(i)] = ppdata[offpulse_int]
        offport = locals()['offport_{0}'.format(i)]

        offres = offport-offgrand
        
        offres_sub = offres - np.mean(offres)

        #locals()['res{0}'.format(i)] = ppdata_fitted-granddata
        #res = locals()['res{0}'.format(i)]

        locals()['res{0}'.format(i)] = ppdata-granddata
        res = locals()['res{0}'.format(i)]

        locals()['var{0}'.format(i)] = np.var(offres)
        var = locals()['var{0}'.format(i)]
        #print("Off Pulse Variance: {}".format(var))

        varsub = np.var(offres_sub)

        ratio = var/varsub
        #print('RFI: {}'.format(ratio))
        locals()['RFI_ratio_subband{0}'.format(i)].append(ratio)

        locals()['chi{0}'.format(i)] = sum(((res)**2)/var)
        chi = locals()['chi{0}'.format(i)]

        chi_sub = sum(((offres)**2)/varsub)
        locals()['redchisub{0}'.format(i)] = chi_sub/(len(offpulse_int)-2)
        redchisub = locals()['redchisub{0}'.format(i)]
        locals()['offpulse_redchi2_subband{0}'.format(i)].append(locals()['redchisub{0}'.format(i)])

        
        #print('redchisub: {}'.format(redchisub))
        dof = 1022
        locals()['redchi{0}'.format(i)] = chi/dof
        redchi = locals()['redchi{0}'.format(i)]
        locals()['total_redchi2_subband{0}'.format(i)].append(locals()['redchi{0}'.format(i)])

        print('Total redchi: {0}; Offpulse redchi: {1}; Potential RFI: {2}'.format(redchi,redchisub,ratio))

        vars_all.append(var)
        redchis_all.append(redchi)
        offpulse_redchis_all.append(redchisub)

for i in range(subbands):

    np.save('offpulse_redchi_subband{0}'.format(i), locals()['offpulse_redchi2_subband{0}'.format(i)])
    np.save('total_redchi_subband{0}'.format(i), locals()['total_redchi2_subband{0}'.format(i)])
    np.save('RFI_ratio_subband{0}'.format(i), locals()['RFI_ratio_subband{0}'.format(i)])

np.save('all_vars',vars_all)
np.save('all_redchis',redchis_all)
np.save('all_redchis_offpulse',offpulse_redchis_all)


fig, ax = plt.subplots()
for i in range(subbands):
    ax.hist(locals()['offpulse_redchi2_subband{0}'.format(i)], bins=100, label = "subint {0}: {1}".format(i+1,freqs[i]))
ax.set_title(r'{0}: Off-Pulse $\chi_r^2$'.format(pulsar))
ax.set_xlabel(r'Offpulse $\chi_r^2$')
ax.set_ylabel('Count')
ax.legend()
fig.savefig('{0}_off_pulse_chi2.pdf'.format(pulsar),format='pdf',dpi=60)

fig, ax = plt.subplots()
for i in range(subbands):
    ax.hist(locals()['total_redchi2_subband{0}'.format(i)], bins=100, label = "subint {0}: {1}".format(i+1,freqs[i]))
ax.set_title(r'{0}: Total $\chi_r^2$'.format(pulsar))
ax.set_xlabel(r'Total $\chi_r^2$')
ax.set_ylabel('Count')
ax.legend()
fig.savefig('{0}_total_chi2.pdf'.format(pulsar),format='pdf',dpi=60)