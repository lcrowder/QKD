# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 15:41:02 2021

@author: crowd
"""

import numpy as np
import matplotlib.pyplot as plt
import os.path

drift=np.arange(0,1)*3.14/40
transmittance=np.arange(1,11)*0.1
keylengths=np.arange(1,2)*1000
eve_present=np.array([1])
out_dir='./scratch/bb84_results_c/'
graph_dir='./scratch/bb84_c_graphs/'

chance=transmittance
d=0

def renyi_theory(p): #Theoretical renyi info of slutsky-brandt attack
    return np.log2( 1+ (4*p*(1-2*p))/((1-p)**2 ))

for e in range(len(eve_present)):
    for k in range(len(keylengths)):
        
        QBER=np.zeros(len(chance))
        QBER_std=np.zeros(len(chance))
        EER=np.zeros(len(chance))
        eve_percent=np.zeros(len(chance))
        eve_percent_std=np.zeros(len(chance))
        avg_sifted_kl=np.zeros(len(chance))
        
        EI=np.zeros(len(chance)) #expected eve info
        
        shannon_info=np.zeros(len(chance))
        shannon_info_std=np.zeros(len(chance))
        renyi_info=np.zeros(len(chance))
        renyi_info_std=np.zeros(len(chance))
        
        renyi_SB=np.zeros(len(chance))
            
        for t in range(len(chance)):
        
            results_file=out_dir+'results_'+str(bool(eve_present[e]))+'_'+str(keylengths[k])+'_'+'{:.2f}'.format(chance[t])+'_'+'{:.2f}'.format(drift[d])+'.txt'
            if os.path.exists(results_file)==True:
                
                print(results_file)
                
                
                data=np.loadtxt(results_file,max_rows=500)
                bitmatch=data[:,1]
                
                avg_sifted_kl[t]=np.mean(bitmatch)
                
                QBER[t]=1-avg_sifted_kl[t]/keylengths[k]
                QBER_std[t]=np.std(bitmatch/keylengths[k])
                
                eve_percent[t]=sum(data[:,3])/len(data[:,3])
                eve_percent_std[t]=np.std(data[:,3])
                
                shannon_info[t]=sum(data[:,4])/len(data[:,4])
                shannon_info_std[t]=np.std(data[:,4])
                
                renyi_info[t]=sum(data[:,5])/len(data[:,5])
                renyi_info_std[t]=np.std(data[:,5])
                
                renyi_SB[t]=renyi_theory(QBER[t])
                
                
        """
        fig=plt.figure()
        plt.plot(drift,EER,'b')
        plt.errorbar(drift,QBER,yerr=QBER_std,color='r')
        plt.legend(('Expected','Simulated'))
        plt.title('QBER for eve='+str(bool(eve_present[e]))+', key length='+str(int(keylengths[k])))
        plt.xlabel('Drift (radians)')
        plt.ylabel('Qubit Error Rate')
        #fig.savefig('./scratch/bb84_inr_graphs/QBER_'+str(bool(eve_present[e]))+'_'+str(keylengths[k])+'.png')
        fig.savefig(graph_dir+'QBER_'+str(bool(eve_present[e]))+'_'+str(keylengths[k])+'.png')
        plt.close()
        """
     
        if eve_present[e]==1:
       
            """
            fig=plt.figure()
            plt.plot(drift,EI)
            plt.errorbar(drift,eve_percent,yerr=eve_percent_std)
            plt.legend(('Theory','Actual'))
            plt.title("Eve's Knowledge of Key with initial length "+str(int(keylengths[k])))
            plt.xlabel('Drift (radians)')
            plt.ylabel('% of key eve knows')
            #fig.savefig('./scratch/bb84_inr_graphs/eve_percent_'+str(keylengths[k])+'.png')
            fig.savefig(graph_dir+'eve_percent_'+str(keylengths[k])+'.png')
            plt.close()
            """
            
            fig=plt.figure()
            plt.plot(QBER,renyi_info)
            plt.plot(QBER,renyi_SB)
            plt.legend(('Direct Measurement','SB attack'))
            plt.title('Eve/Bob mutual information')
            plt.xlabel('error rate')
            plt.ylabel('Renyi Information')
            fig.savefig(graph_dir+'information_'+str(keylengths[k])+'.png')
            plt.close()
            