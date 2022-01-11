import matplotlib.pyplot as plt 
import glob 
import numpy as np
import sys
from matplotlib import rc

systems = ['MUT3_holo']

#################################################################Distribution calculation for bootstraped sample##########################################
for system in systems:
    count = 0
    bins = [i/10 for i in range(0,201,1)]
    prob_arr = np.zeros([len(bins)-1,2,20])
    for i in range(20):
        filename = './bt_traj/' +system+'_bt_' + str(i) + '.npy'
        dist = np.load(filename)
        prob = np.zeros([len(bins)-1,2])
        newbins = np.zeros([len(bins)-1,2])
        cou = 0
        for feature in range(1,4,2):
            nSD, binsSD= np.histogram(dist[:,feature]*10, bins=bins, density=True)
            averageSD = [(binsSD[j]+binsSD[j+1])/2 for j in range(len(binsSD)-1)]
            prob[:,cou] = nSD
            newbins[:,cou] = averageSD
            cou = cou + 1
        prob_arr[:,:,count] = prob 
        count +=1 
    prob_std_err = np.std(prob_arr,axis=2)
    prob_mean = np.mean(prob_arr,axis=2)
    print(prob_mean.shape,prob_std_err.shape)
    

##################################################################Figure Specifications##################################################################

fig,axs = plt.subplots(1,1,figsize=(10,7),constrained_layout=True)

hfont = {'fontname':'Helvetica'}
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

color = ['r','g']
label = ['TYR27-TYR473','TYR330-PRO499']

print(prob_mean)
for i in range(2):
    plt.plot(newbins[:,i],prob_mean[:,i],color=color[i],label=label[i], linewidth=3)
    plt.plot(newbins[:,i],prob_mean[:,i],color=color[i], linewidth=3)
    plt.fill_between(newbins[:,i], prob_mean[:,i]+prob_std_err[:,i], prob_mean[:,i]-prob_std_err[:,i], color=color[i], alpha=0.3)
    
axs.set_xlim(0,12)
axs.set_xticks(range(int(0),int(12)+1,2))
axs.set_xticklabels(range(int(0),int(12)+1,2))

axs.set_ylim(0,1)
axs.set_yticks([0,0.2,0.4,0.6,0.8,1])
axs.set_yticklabels([0,0.2,0.4,0.6,0.8,1])

plt.xticks(fontsize=22)
plt.yticks(fontsize=22)

plt.xlabel('Distance (\AA)',**hfont,fontsize=30,fontweight='bold')
plt.ylabel('Probability Density', **hfont,fontsize=30,fontweight='bold')

plt.savefig('dist_distribution.png',dpi=500)
