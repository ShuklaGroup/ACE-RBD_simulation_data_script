import glob
import numpy as np
import pickle
import pyemma
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from matplotlib import rc

#######################################Parameters and hyperparameters definition##############################################
bins = 100        #number of bins used to divide the landscape

R = 0.001987        #Boltzman's constant (kcal/mol/K)
T = 300             #Temperature (K)

#######################################MSM loading##############################################
V24_msm = pickle.load(open("./../V2.4_holo_MSM_obj.pkl","rb")) #loading of MSM object
V24_weights = np.concatenate(V24_msm.trajectory_weights())     #weights(probability density of each frames)

wild_msm = pickle.load(open("./../wild_holo_MSM_obj.pkl","rb")) #loading of MSM object
wild_weights = np.concatenate(wild_msm.trajectory_weights())     #weights(probability density of each frames)

#######################################Figure Specification##############################################
hfont = {'fontname':'Helvetica','fontweight':'bold'}
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})   #Figure font definition 

fig_wid = 10        #Width of the genarated figure 
fig_hig = 7         #length of the genarated figure


fig, axs = plt.subplots(1,1,figsize=(10,7),constrained_layout=True)
label = ['WT complex','v2.4 complex']
color = ['cyan','orange']


#######################################Data loading##############################################
wild_feature = pickle.load(open('./../wild_holo_loop_dist.pkl','rb'))  
wild_txx = np.concatenate(wild_feature)[:,2]

V24_feature = pickle.load(open('./../V2.4_holo_loop_dist.pkl','rb'))
V24_txx = np.concatenate(V24_feature)[:,2]


#######################################Histogram plot for wild type##############################################
mean = np.dot(wild_txx,wild_weights)
nSD, binsSD= np.histogram(wild_txx, bins, density=True, weights=wild_weights)
averageSD = [(binsSD[i]+binsSD[i+1])/2 for i in range(len(binsSD)-1)]
highest =nSD[np.where(binsSD>mean)[0][0]]

plt.plot([mean for _ in range(100)],np.linspace(0,highest,100), '--',color=color[0],linewidth=3)
plt.plot(averageSD,nSD,linewidth=4,color=color[0],label=label[0])


#######################################Histogram plot for V2.4 mutant##############################################
mean = np.dot(V24_txx,V24_weights)
nSD, binsSD= np.histogram(V24_txx, bins, density=True, weights=V24_weights)
averageSD = [(binsSD[i]+binsSD[i+1])/2 for i in range(len(binsSD)-1)]
highest =nSD[np.where(binsSD>mean)[0][0]]

plt.plot([mean for _ in range(100)],np.linspace(0,highest,100), '--',color=color[1],linewidth=3)
plt.plot(averageSD,nSD,linewidth=4,color=color[1],label=label[1])

#######################################Figure Modification##############################################
axs.legend(prop={'size': 18})                                                            #legend Size

axs.set_xlim([0,3])                                                                      #min and max limit of x-axis
axs.set_ylim([0,15])                                                                     #min and max limit of y-axis

axs.set_xticks(np.linspace(0, 3, 3))                                                     #x-axis ticks
axs.set_xticklabels(np.linspace(0, 3, 3))                                                #x-axis tick labels

axs.set_yticks(np.linspace(0, 15, 5))                                                    #y-axis ticks
axs.set_yticklabels(np.linspace(0, 15, 5))                                               #y-axis tick labels

plt.xticks(fontsize=20)                                                                  #x-axis tick size
plt.yticks(fontsize=20)                                                                  #y-axis tick size

axs.set_xlabel('RBD loop 1 distance' + ' (nm)', **hfont,fontsize=24)                     #x-axis label
axs.set_ylabel('Probability Density', **hfont,fontsize=24)                               #y-axis label

plt.tight_layout()
plt.savefig('Loop_2_dist.png',dpi=300)
plt.close()







 




