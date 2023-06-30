# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 17:05:05 2020

@author: crowd
"""

import numpy as np
import matplotlib.pyplot as plt
import os.path

out_dir='/scratch/lac468/bb84_binom_results/'
eve_present=np.array([0,1])
keylengths=10*np.arange(1,16)
loss=np.arange(0,1)
drift=3.14/20*np.arange(0,11)

alpha1=0.25
alpha2=0.5
alpha3=0.75

for a in range(len(eve_present)):
    for b in range(len(keylengths)):
        for c in range(len(loss)):
            
            rejections=np.zeros((3,len(drift)))
            
            for d in range(len(drift)):
                
                extension='_'+str(bool(eve_present[a]))+'_'+str(keylengths[b])+'_'+'{:.2f}'.format(drift[d])
                results_file=out_dir+'binom'+extension+'.txt'
                
                if os.path.exists(results_file)==True:
                    print(extension)
                    
                    p_vals=np.loadtxt(results_file,ndmin=1)
                    """
                    bins=keylengths[b]
                    fig=plt.hist(p_vals,bins,align='mid')
                    plt.title('p-values '+extension)
                    plt.xlabel('p-value')
                    plt.ylabel('Frequency')
                    plt.savefig(out_dir+'dist_graphs/pval'+extension+'.png')
                    plt.close()
                    """

                    a1_count=0
                    a2_count=0
                    a3_count=0
                    
                    for i in range(len(p_vals)):
                        if p_vals[i]<alpha1:
                            a1_count=a1_count+1
                        if p_vals[i]<alpha2:
                            a2_count=a2_count+1
                        if p_vals[i]<alpha3:
                            a3_count=a3_count+1    
                
                    print("alpha="+str(alpha1)+", we rejected the secure key hypothesis "+str(a1_count)+"/"+str(len(p_vals))+" times.")
                    print("alpha="+str(alpha2)+", we rejected the secure key hypothesis "+str(a2_count)+"/"+str(len(p_vals))+" times.")
                    print("alpha="+str(alpha3)+", we rejected the secure key hypothesis "+str(a3_count)+"/"+str(len(p_vals))+" times.")
                    
                    rejections[0,d]=a1_count/len(p_vals)
                    rejections[1,d]=a2_count/len(p_vals)
                    rejections[2,d]=a3_count/len(p_vals)
                    
            fig=plt.figure()
            plt.plot(drift,rejections[0,:])
            plt.plot(drift,rejections[1,:])
            plt.plot(drift,rejections[2,:])
            plt.legend((alpha1,alpha2,alpha3),title='alpha values')
            plt.title('Proportion of Null Hypothesis Rejections '+str(bool(eve_present[a]))+'_'+str(keylengths[b]))
            plt.xlabel('Drift (radians)')
            plt.ylabel('% Rejected')
            fig.savefig(out_dir+'rejection_graphs/rejection_drift'+'_'+str(bool(eve_present[a]))+'_'+str(keylengths[b]))
            plt.close()
                    
                    
                    
                    