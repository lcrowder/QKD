# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 12:42:21 2020

@author: crowd
"""
import numpy as np
import matplotlib.pyplot as plt
import os.path
import math

drift=np.arange(0,10)*3.14/40
transmittance=np.array([1])
keylengths=np.arange(1,2)*1000
eve_present=np.array([1])
#out_dir='./scratch/bb84_results_inr/'
out_dir='./scratch/bb84_results_inr_test/'

graph_dir='./scratch/bb84_inr_graphs_test/'

for e in range(len(eve_present)):
    for k in range(len(keylengths)):
        for t in range(len(transmittance)):
            
            QBER=np.zeros(len(drift))
            QBER_std=np.zeros(len(drift))
            EER=np.zeros(len(drift))
            eve_percent=np.zeros(len(drift))
            eve_percent_std=np.zeros(len(drift))
            avg_sifted_kl=np.zeros(len(drift))
            
            EI=np.zeros(len(drift))
            
            for d in range(len(drift)):
                
                if drift[d] <= np.pi/6:
                    EI[d]= 100 * 4*math.sin(drift[d])**2 * 2/3
                else:
                    EI[d]= 100 *(1-1/(4*math.cos(drift[d])**2))
                
                EER[d]=(math.sin(drift[d]))**2
            
                results_file=out_dir+'results_'+str(bool(eve_present[e]))+'_'+str(keylengths[k])+'_'+'{:.2f}'.format(transmittance[t])+'_'+'{:.2f}'.format(drift[d])+'.txt'
                if os.path.exists(results_file)==True:
                    
                    data=np.loadtxt(results_file)
                    bitmatch=data[:,1]
                    
                    avg_sifted_kl[d]=np.mean(bitmatch)
                    
                    QBER[d]=1-avg_sifted_kl[d]/keylengths[k]
                    QBER_std[d]=np.std(bitmatch/keylengths[k])
                    
                    eve_percent[d]=sum(data[:,3])/len(data[:,3])
                    eve_percent_std[d]=np.std(data[:,3])
        
        fig=plt.figure()
        plt.plot(drift,EER,'b')
        plt.errorbar(drift,QBER,yerr=QBER_std,color='r')
        plt.legend(('Expected','Simulated'))
        plt.title('QBER for eve='+str(bool(eve_present[e]))+', key length='+str(int(keylengths[k])))
        plt.xlabel('Drift (radians)')
        plt.ylabel('Qubit Error Rate')
        fig.savefig(graph_dir+'QBER_'+str(bool(eve_present[e]))+'_'+str(keylengths[k])+'.png')
        plt.close()
        
        if eve_present[e]==1:
            fig=plt.figure()
            plt.plot(drift,EI)
            plt.errorbar(drift,eve_percent,yerr=eve_percent_std)
            plt.legend(('Theory','Actual'))
            plt.title("Eve's Knowledge of Key with initial length "+str(int(keylengths[k])))
            plt.xlabel('Drift (radians)')
            plt.ylabel('% of key eve knows')
            fig.savefig(graph_dir+'eve_percent_'+str(keylengths[k])+'.png')
            plt.close()