import matplotlib.pyplot as plt 
import glob 
import numpy as np
import sys
from matplotlib import rc


systems = ['MUT3_holo']

#################################################################Hydrogen bond calculation for bootstraped sample##########################################

for system in systems:
    count = 0
    count_hydrogen = np.zeros([4,20])
    for i in range(20):
        filename = './bt_traj/' +system+'_bt_'+ str(i) + '.npy' 
        dist = np.load(filename)
        bin_hydrogen = np.where(dist<0.4,1,0)
        count_hydrogen[:,count] = np.sum(bin_hydrogen,axis=0)
        count +=1
        
    
prob_YY = (count_hydrogen[0,:]+count_hydrogen[1,:])/80000
prob_YY_mean = np.mean(prob_YY)
prob_YY_std = np.std(prob_YY)


prob_YP = count_hydrogen[3,:]/40000
prob_YP_mean = np.mean(prob_YP)
prob_YP_std = np.std(prob_YP)

print(prob_YY_std,prob_YP_std)


##################################################################Figure Specifications#######################################################################


fig,axs = plt.subplots(1,1,figsize=(10,7),constrained_layout=True)

hfont = {'fontname':'Helvetica'}
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

color = ['r','g']


axs.bar([1.5],[prob_YY_mean],yerr=[prob_YY_std],color = 'r',label='TYR27-TYR473',align='center',alpha=0.5,ecolor='black',capsize=10)
#axs.scatter(np.random.uniform(low=1.45, high=1.55, size=(20,)),prob_YY,s=24,c='k')


axs.bar([3.5],[prob_YP_mean],yerr=[prob_YP_std],color = 'g',label='TYR330-PRO499',align='center',alpha=0.5,ecolor='black',capsize=10)
#axs.scatter(np.random.uniform(low=3.45, high=3.55, size=(20,)),prob_YP,s=24,c='k')


axs.set_xticklabels(['Loop 1', 'Loop 2'])
axs.set_xticks([1.5, 3.5])


plt.ylabel('Hydrogen Bond Probability',**hfont,fontsize=30)
    
axs.set_xlim(0,5)
axs.set_ylim(0,1)

axs.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])
axs.set_yticklabels([0.0,0.2,0.4,0.6,0.8,1.0])

plt.xticks(fontsize=22)
plt.yticks(fontsize=22)

axs.legend(fontsize=18)

plt.savefig('hbond_prob.png',dpi=500)
