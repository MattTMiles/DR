import numpy as np
import pandas as pd
import psrchive
import matplotlib.pyplot as plt

datafile = psrchive.Archive_load('/fred/oz002/users/mmiles/MSP_DR/J1103_study/archs/grand.dly')
datafile.remove_baseline()
datafile.dedisperse()

early = psrchive.Archive_load('/fred/oz002/users/mmiles/MSP_DR/J1103_study/archs/early.dlyT.sm')
late = psrchive.Archive_load('/fred/oz002/users/mmiles/MSP_DR/J1103_study/archs/late.dlyT.sm')

esf = psrchive.ProfileShiftFit()
lsf = psrchive.ProfileShiftFit()

early.remove_baseline()
early.dedisperse()
late.remove_baseline()
late.dedisperse()
earlyprof = early.get_Profile(0,0,0)
lateprof = late.get_Profile(0,0,0)

esf.set_standard(earlyprof)
lsf.set_standard(lateprof)

early_resids=[]
late_resids=[]
consistent = []
for i, subint in enumerate(datafile):
    subint.remove_baseline()
    subint.dedisperse()
    subintdata = subint.get_Profile(0,0)
    latesubint = subintdata
    earlysubint = subintdata
    esf.set_Profile(earlysubint)
    lsf.set_Profile(latesubint)
    latedata = latesubint.get_amps()
    earlydata = earlysubint.get_amps()
    late_prof_data = lateprof.get_amps()
    late_prof_data = late_prof_data/late_prof_data.max()
    early_prof_data = earlyprof.get_amps()
    early_prof_data = early_prof_data/early_prof_data.max()
    latedata = latedata/latedata.max()
    earlydata = earlydata/earlydata.max()
    late_res = late_prof_data-latedata
    early_res = early_prof_data-earlydata
    early_resids.append(early_res)
    late_resids.append(late_res)

    for j in [25,26,27,28,30,31,33]:
        if i==j:
            consistent.append(latedata)


late_resids = np.array(late_resids)
early_resids = np.array(early_resids)

total_early = np.sum(early_resids,axis=0)
total_late = np.sum(late_resids,axis=0)

consistent = np.array(consistent)

consistent_standard = np.sum(consistent,axis=0)
consistent_standard = consistent_standard/consistent_standard.max()

consistent_diff_early = early_prof_data-consistent_standard
consistent_diff_late = late_prof_data-consistent_standard

data =datafile.get_data()
data = data[:,0,0,:]

target = data[53]
target2 = data[29]
fig,ax = plt.subplots()
scaled_target = target/target.max()
scaled_target *= 0.85

scaled_target2 = target2/target2.max()
scaled_target2 *= 0.85

v_early = psrchive.Archive_load('v_early.dly')
v_early.remove_baseline()
v_early.dedisperse()

vearly_data = v_early.get_data()
vearly_data = vearly_data[:,0,0,:]

vearly_sum = np.sum(vearly_data,axis=0)
vearly_scaled = vearly_sum/vearly_sum.max()

vearly_scaled *= 0.7

vearly_diff = consistent_standard - vearly_scaled



fig,ax = plt.subplots()

ax.plot(consistent_standard,label='consistent')
ax.plot(vearly_scaled,label='vearly')
ax.plot(vearly_diff,label='res_diff')
ax.legend()
plt.show()

'''
ax.plot(scaled_target2)
ax.plot(consistent_standard)
plt.show()
'''



