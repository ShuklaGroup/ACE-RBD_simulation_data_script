import matplotlib.pyplot as plt 
import glob 
import numpy as np
import sys
from matplotlib import rc


####################################################################################parameter definition and initialization##################################################

hfont = {'fontname':'Helvetica'}
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

systems = ['ACE_holo','MUT3_holo']

color = {'ACE_holo':'cyan','MUT3_holo':'orange'}
label = {'ACE_holo':'wild','MUT3_holo':'v2.4'}

fig,axs = plt.subplots(1,1,figsize=(10,7),constrained_layout=True)

sys_rmsf = []


##################################################################################rmsf calculation and plotting for bootstrapped samples#####################################
for system in systems:
    original_resid = np.array([i for i in range(438,510)])
    resid = np.array([i+391 for i in original_resid])
    rmsf_arr = np.zeros([len(resid),20])
    count  = 0
    for filename in glob.glob('./bt_traj/' +system+'_bt_*_fluct.agr'):
        f = open(filename,'r')
        data = f.readlines()
        rmsf = []
        r  = []
        for line in data:
            words = line.split()
            if words[0][0] != "@":
                if float(words[0]) in resid:
                    rmsf.append(float(words[1]))           

        rmsf_arr[:,count] = np.array(rmsf)
        count += 1        

    sys_rmsf.append(rmsf_arr)
    rmsf_std_err = np.std(rmsf_arr,axis=1)
    rmsf_mean = np.mean(rmsf_arr,axis=1)
    plt.plot(original_resid,rmsf_mean,color=color[system],label=label[system],linewidth=3)
    plt.fill_between(original_resid, rmsf_mean+rmsf_std_err, rmsf_mean-rmsf_std_err, color=color[system], alpha=0.3)

#print(rmsf_mean)
diff_rmsf = abs(sys_rmsf[0] - sys_rmsf[1])

diff_rmsf_std_err = np.std(diff_rmsf,axis=1)
diff_rmsf_mean = np.mean(diff_rmsf,axis=1)

plt.plot(original_resid,diff_rmsf_mean,color='black',label='difference', linewidth=3)
plt.fill_between(original_resid, diff_rmsf_mean+diff_rmsf_std_err, diff_rmsf_mean-diff_rmsf_std_err, color='black', alpha=0.3)


plt.axvspan(477, 485, facecolor='grey', alpha=0.3)
plt.axvspan(497, 505, facecolor='grey', alpha=0.3)
    
axs.set_xlim(438,507)
axs.set_xticks(range(int(440),int(510)+1,20))
axs.set_xticklabels(range(int(440),int(510)+1,20))

axs.set_ylim(0,5.0)
axs.set_yticks(range(int(0),int(5)+1,1))
axs.set_yticklabels(range(int(0),int(5)+1,1))

plt.xticks(fontsize=22)
plt.yticks(fontsize=22)

axs.legend(fontsize=18)

plt.ylabel('RMSF (\AA)',**hfont,fontsize=30,fontweight='bold')
plt.xlabel('RBD Residue', **hfont,fontsize=30,fontweight='bold')

plt.savefig('rmsf.png',dpi=500)
