from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
from TIBD_PopD_2017 import model

data = np.genfromtxt('ExpData.csv', dtype=None, delimiter=',', names=True, encoding='UTF-8')

print(data)
# print(data.dtype.names)

Drug = '1D11'
Cond = 'TUMOR'

tumor = np.array([[d[0],d[4]] for d in data 
                        if d[1]=='Tumor' and d[2]==Drug and d[3]==Cond])
osteoblast = np.array([[d[0],d[4]] for d in data 
                        if d[1]=='Osteoblasts' and d[2]==Drug and d[3]==Cond])
osteoclast = np.array([[d[0],d[4]] for d in data 
                        if d[1]=='Osteoclasts' and d[2]==Drug and d[3]==Cond])
bone = np.array([[d[0],d[4]] for d in data 
                        if d[1]=='Bone' and d[2]==Drug and d[3]==Cond])
# ob_control = np.array([[d[0],d[3]] for d in data if d[1]=='Osteoblasts' and d[2]=='CONTROL'])
# oc_control = np.array([[d[0],d[3]] for d in data if d[1]=='Osteoclasts' and d[2]=='CONTROL'])

# tumor_cell_density = np.array([[6.7, 2.85], [14, 7.60], [21, 89.62]])
# bv_tv = np.array([[0, 0.121], [7, 0.122], [14, 0.091], [21, 0.033]])

# colors = ['blue', 'red', 'orange', 'purple', 'green', 'cyan']

param_values = np.load("most_likely_par_leonard_prior.npy")
    
par_names = [m.name for m in model.parameters]   
idx_drug0 = par_names.index('Drug_0')
idx_kdrugB = par_names.index('kdrug_B')
idx_kdrugT = par_names.index('kdrug_T')

# Drug_0 = 1
param_values[idx_drug0] = 1
# kdrug_B = 0
param_values[idx_kdrugB] = 0
# kdrug_T = ADJUSTABLE (gammaT = 0.060220613103437244 from fit)
param_values[idx_kdrugT] = 0.1

for i,r21_val in enumerate([0.87]): #[-0.5]): #np.arange(-1,0.2,0.2)):
#     print(r21_val)
    
#     param_values = {'C_0' : 10, #15, #10,
#                     'B_0' : 10, #316, #10,
#                     'T_0' : 1,
#                     'Z_0' : 100,
#                     'LT' : 100,
#                     'alpha1' : 6,
#                     'alpha2' : 4,
#                     'beta1' : 0.2,
#                     'beta2' : 0.002, #0.02, #0.002,
#                     'g11' : 0.5, #1.1, 0.5
#                     'g12' : 1,
#                     'g21' : -0.5, 
#                     'g22' : 0,
#                     'r11' : 0.071, #0.005, #0.071,
#                     'r21' : r21_val, #0.87, 
#                     'r12' : 0,
#                     'r22' : 0.2,
#                     'k1' : 0.24,
#                     'k2' : 0.0017,
#                     'gammaT' : 0.63} 

#     tspan = np.linspace(7,25,(25-7)*10+1)
    tspan = np.linspace(0,32,321)
    sim = ScipyOdeSimulator(model, tspan, verbose=True) 
    x = sim.run(param_values=param_values)

    plt.figure('tumor')
    plt.plot(tumor[:,0], tumor[:,1], '.k', ms=10, label='Experiment')
    plt.plot(tspan, x.observables['Obs_T'], '-b', lw=2, label='Model')
    plt.xlabel('time (days)')
    plt.ylabel('tumor cell density (%)')
    plt.legend(loc=0)
    
    plt.figure('bone')
    plt.plot(bone[:,0], bone[:,1], '.k', ms=10, label='Experiment')
    plt.plot(tspan, x.observables['Obs_Z'], '-r', lw=2, label='Model')
    plt.xlabel('time (days)')
    plt.ylabel('bone density (%)')
    plt.legend(loc=0) 
'''
    plt.figure('Oc')
    plt.plot(osteoclast[:,0], osteoclast[:,1], '.k', ms=10, label='Experiment')
    plt.plot(tspan, x.observables['Obs_C'], '-r', lw=2, label='Model')
    plt.xlabel('time (days)')
    plt.ylabel('osteoclasts (count)')
    plt.legend(loc=0)
    
    plt.figure('Ob')
    plt.plot(osteoblast[:,0], osteoblast[:,1], '.k', ms=10, label='Experiment')
    plt.plot(tspan, x.observables['Obs_B'], '-b', lw=2, label='Model')
    plt.xlabel('time (days)')
    plt.ylabel('osteoblasts (count)')
    plt.legend(loc=0)
'''
plt.show()
