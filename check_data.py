# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 10:31:01 2020

@author: crowd
"""
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 14:35:53 2020

@author: crowd
"""


import numpy as np
import os.path
import math

out_dir='/scratch/lac468/bb84_pns_results/'
#out_dir=''
eve_present=True
#keylengths=np.arange(100,300,20)
keylengths=np.array([280])
transmittance=np.arange(1,11)*0.1
drift=np.arange(0,1)
mean_photon_number=np.arange(1,11)*0.1

keylength_false=1000 #only used one keylength since it doesn't affect probabilities above a minimum length

#analytically derived distribution for expected photon number in pulse, based on poisson dist. and channel loss
def prob(n_mean,rbp,n):
    return rbp**n * np.exp(-n_mean*rbp)*n_mean**n/math.factorial(n)
ns=np.array([0,1,2])

for k in range(len(keylengths)):
    
    for m in range(len(mean_photon_number)):
        for d in range(len(drift)):
            for t in range(len(transmittance)):
                results_file=out_dir+'results_'+str(eve_present)+'_'+str(keylengths[k])+'_'+'{:.2f}'.format(transmittance[t])+'_'+'{:.2f}'.format(drift[d])+'_'+'{:.1f}'.format(mean_photon_number[m])+'.txt'
                if os.path.exists(results_file)==True:
                    mat=np.loadtxt(results_file)
                    print(mean_photon_number[m])
                    print(transmittance[t])
                    print("Eve "+str(len(mat[0,:])))
 
                    
                #No eve simulation:
                file=out_dir+'results_False_'+str(keylength_false)+'_'+'{:.2f}'.format(transmittance[t])+'_'+'{:.2f}'.format(drift[d])+'_'+'{:.1f}'.format(mean_photon_number[m])+'.txt'
                if os.path.exists(file)==True:
                    mat=np.loadtxt(results_file)
                    print("No Eve "+str(len(mat[0,:])))
                    