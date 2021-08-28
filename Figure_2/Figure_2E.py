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
msm = pickle.load(open("./../V2.4_holo_MSM_obj.pkl","rb")) #loading of MSM object 
weights=np.concatenate(msm.trajectory_weights())           #weights(probability density of each frames)

#######################################Figure Specification##############################################
hfont = {'fontname':'Helvetica','fontweight':'bold'}
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})   #Figure font definition 

fig_wid = 10        #Width of the genarated figure 
fig_hig = 7         #length of the genarated figure
cmap = mpl.cm.jet   #color bar used in the figure

Max_energy =5       #maximum energy projected in color bar 
    
fig, axs = plt.subplots(1,1,figsize=(10,7),constrained_layout=True)
label = ['TYR27-TYR473','THR79-GLY485','TYR330-PRO499']
line = ['-','--','-.']

#######################################Data loading##############################################
dist_features = pickle.load(open('../V2.4_holo_polar_bond.pkl','rb'))       #loading the distance features 
txx_dist = np.concatenate(dist_features)                                    #concatenating the distance features  

#######################################Probability Density Calculation##############################################
for i in range(3):
    x_data = txx_dist[:,i]                                                               
    nSD, binsSD= np.histogram(x_data, bins, density=True, weights=weights)
    averageSD = [(binsSD[j]+binsSD[j+1])/2 for j in range(len(binsSD)-1)]
    plt.plot(averageSD,nSD,line[i],color='Black', label = label[i], linewidth=4)


#######################################Figure Modification##############################################
axs.legend(prop={'size': 18})                                                  #legend Size 

axs.set_xlim([0,1])                                                            #min and max limit of x-axis
axs.set_ylim([0,10])                                                           #min and max limit of y-axis

axs.set_xticks(np.linspace(0, 1, 5))                                           #x-axis ticks
axs.set_xticklabels(np.linspace(0, 1, 5))                                      #x-axis tick labels

axs.set_yticks(np.linspace(0, 10, 5))                                          #y-axis ticks 
axs.set_yticklabels(np.linspace(0, 10, 5))                                     #y-axis tick labels  

plt.xticks(fontsize=20)                                                        #x-axis tick size                 
plt.yticks(fontsize=20)                                                        #y-axis tick size                

axs.set_xlabel('distance' + ' (nm)', **hfont,fontsize=24)    #x-axis label
axs.set_ylabel('Probability Density', **hfont,fontsize=24)   #y-axis label

plt.tight_layout()
plt.savefig('polar_bond_dist.png',dpi=300)                                     



