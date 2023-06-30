# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 19:21:04 2020

@author: crowd
"""
"""
Analyze data from combined pns and intercept-and-resend attacks
"""
import numpy as np
import matplotlib.pyplot as plt
import os.path
import math
from scipy.interpolate import CubicSpline

out_dir='/scratch/lac468/bb84_combined_new_results/'
out_dir_inr='./scratch/bb84_results_inr_test/'
eve_chances_file='/scratch/lac468/eve_chances.txt'
#out_dir=''
eve_present=True
#keylengths=np.arange(100,300,20)
keylengths=np.array([1000])
transmittance=np.arange(1,11)*0.1
drift=np.arange(0,10)*3.14/40
mean_photon_number=np.arange(1,11)*0.1

#analytically derived distribution for expected photon number in pulse, based on poisson dist. and channel loss
def prob(n_mean,rbp,n):
    return rbp**n * np.exp(-n_mean*rbp)*n_mean**n/math.factorial(n)

for k in range(len(keylengths)):
    
    Pb0=np.zeros((len(mean_photon_number),len(drift),len(transmittance))) #array for P(n_b=0), num Bob actually receives
    Pb0_std=np.zeros((len(mean_photon_number),len(drift),len(transmittance))) #standard dev
    eve_percent=np.zeros((len(mean_photon_number),len(drift),len(transmittance))) #What percent of the key does eve know?
    eve_percent_std=np.zeros((len(mean_photon_number),len(drift),len(transmittance))) #standard dev of eve_percent
    QBER=np.zeros((len(mean_photon_number),len(drift),len(transmittance)))
    QBER_std=np.zeros((len(mean_photon_number),len(drift),len(transmittance)))
    
    #analytical derived expectation for what Bob would receive from lossy channel w/ no Eve:
    P0=np.zeros((len(mean_photon_number),len(transmittance)))
    
    #no eve simulation results:
    P0_false=np.zeros((len(mean_photon_number),len(drift),len(transmittance)))
    P0_false_std=np.zeros((len(mean_photon_number),len(drift),len(transmittance)))
    QBER_false=np.zeros((len(mean_photon_number),len(drift),len(transmittance)))
    QBER_false_std=np.zeros((len(mean_photon_number),len(drift),len(transmittance)))
    
    EER=np.zeros(len(drift))
    EI=np.zeros((len(mean_photon_number),len(drift),len(transmittance)))  #eves information, in theory
    
    "inr data:"
    eve_percent_inr=np.zeros(len(drift))
    eve_percent_inr_std=np.zeros(len(drift))
    QBER_inr=np.zeros(len(drift))
    QBER_inr_std=np.zeros(len(drift))
    
    print('INR '+str(eve_percent_inr[0]))
    
    for m in range(len(mean_photon_number)):
        for d in range(len(drift)):
            
            EER[d]=math.sin(drift[d])**2
            
            results_file_inr=out_dir_inr+'results_'+str(eve_present)+'_'+str(keylengths[k])+'_'+'{:.2f}'.format(transmittance[-1])+'_'+'{:.2f}'.format(drift[d])+'.txt'
            if os.path.exists(results_file_inr)==True:
                basismatch_temp=np.loadtxt(results_file_inr,usecols=0)
                bitmatch_temp=np.loadtxt(results_file_inr,usecols=1)
                eve_percent_temp=np.loadtxt(results_file_inr,usecols=3)
                
                QBER_inr[d]=np.sum(1-bitmatch_temp/basismatch_temp)/len(bitmatch_temp)
                QBER_inr_std[d]=np.std(1-bitmatch_temp/basismatch_temp)
                eve_percent_inr[d]=np.sum(eve_percent_temp)/len(eve_percent_temp)
                eve_percent_inr_std[d]=np.std(eve_percent_temp)
                
            
            for t in range(len(transmittance)):
                results_file=out_dir+'results_'+str(eve_present)+'_'+str(keylengths[k])+'_'+'{:.2f}'.format(transmittance[t])+'_'+'{:.2f}'.format(drift[d])+'_'+'{:.1f}'.format(mean_photon_number[m])+'.txt'
                if os.path.exists(results_file)==True:
                    
                    basismatch_temp=np.loadtxt(results_file,usecols=0)
                    bitmatch_temp=np.loadtxt(results_file,usecols=1)
                    eve_percent_temp=np.loadtxt(results_file,usecols=2)
                    Pb0_temp=np.loadtxt(results_file,usecols=3)
                    
                    QBER[m,d,t]=np.sum(1-bitmatch_temp/basismatch_temp)/len(bitmatch_temp)
                    QBER_std[m,d,t]=np.std(1-bitmatch_temp/basismatch_temp)
                    eve_percent[m,d,t]=np.sum(eve_percent_temp)/len(eve_percent_temp)
                    eve_percent_std[m,d,t]=np.std(eve_percent_temp)
                    Pb0[m,d,t]=np.sum(Pb0_temp)/len(Pb0_temp)
                    Pb0_std[m,d,t]=np.std(Pb0_temp)
                    
                    P0[m,t]=prob(mean_photon_number[m],transmittance[t],0) #theory
                    
                
                    """
                    fig=plt.figure(1)
                    plt.plot(np.zeros(len(Pb0)),Pb0,'*')
                    plt.plot(ns[0],P0,'r+')
                    fig=plt.figure(2)
                    plt.plot(np.ones(len(Pb1)),Pb1,'*')
                    plt.plot(ns[1],P1,'r+')
                    """
                    """
                    print(keylengths[k])
                    print('transmittance '+str(transmittance[t]))
                    print('mpn '+str(mean_photon_number[m]))
                    print('drift '+str(drift[d]))
                    print('')
                    """
                #No eve simulation:
                file=out_dir+'results_False_'+str(keylengths[k])+'_'+'{:.2f}'.format(transmittance[t])+'_'+'{:.2f}'.format(drift[d])+'_'+'{:.1f}'.format(mean_photon_number[m])+'.txt'
                if os.path.exists(file)==True:
                    
                    basismatch_temp=np.loadtxt(file,usecols=0)
                    bitmatch_temp=np.loadtxt(file,usecols=1)
                    Pb0_temp=np.loadtxt(file,usecols=3)
                    
                    QBER_false[m,d,t]=np.sum(1-bitmatch_temp/basismatch_temp)/len(bitmatch_temp)
                    QBER_false_std[m,d,t]=np.std(1-bitmatch_temp/basismatch_temp)
                    P0_false[m,d,t]=np.sum(Pb0_temp)/len(Pb0_temp)
                    P0_false_std[m,d,t]=np.std(Pb0_temp)
                    
                
                #Use pre-generated data to decide the proportion of n=1,2,3 photon pulses Eve will delete.
                #interpolation with eve_chances.txt works for 0.025 <= <n> <= 1.025, 0.1 <= transmittance <= 1.0
                eve_chances_data=np.loadtxt(eve_chances_file)
                t_points=eve_chances_data[0:20,1] #transmittance data points solved for directly with eves_chances algorithm
                n_spline=round(2*mean_photon_number[m],1)/2 #round <n> to closest known data point 
                
                #Load data for specific <n>, interpolate C[i] as function of transmittance with cubic spline.  
                C1_points=eve_chances_data[:,2][eve_chances_data[:,0]==n_spline]
                C2_points=eve_chances_data[:,3][eve_chances_data[:,0]==n_spline]
                C3_points=eve_chances_data[:,4][eve_chances_data[:,0]==n_spline]
                cs1=CubicSpline(t_points,C1_points)
                C1=CubicSpline.__call__(cs1,transmittance[t])
                cs2=CubicSpline(t_points,C2_points)
                C2=CubicSpline.__call__(cs2,transmittance[t])
                cs3=CubicSpline(t_points,C3_points)
                C3=CubicSpline.__call__(cs3,transmittance[t])
                
                #If Eve is not present, she does not kill any qubits
                if eve_present==False:
                    C1=0
                    C2=0
                    C3=0
                
                #convevnient poisson dist for particular mpn.
                def pois(x):
                    return prob(mean_photon_number[m],1,x)
                   
                prop_single_photons= (pois(1)*(1-C1)) / (1-(pois(0) + pois(1)*C1 + pois(2)*C2 + pois(3)*C3 ))
                print(prop_single_photons)
                if prop_single_photons==0.0:
                    temp_chance=1.0
                else:
                    temp_chance=(4*math.sin(drift[d])**2)/prop_single_photons
                
                if temp_chance <= 1:
                    EI[m,d,t]=1-prop_single_photons + 2/3*4*math.sin(drift[d])**2 
                else: 
                    EI[m,d,t]=1 - prop_single_photons/(1+2*math.cos(drift[d])**2/(1-prop_single_photons/4))
                
                EI[m,d,t]=100*EI[m,d,t]
                
               
        """
        "Eve's knowledge vs transmittance, drift, mpn"
        T,D=np.meshgrid(transmittance,drift)
        fig=plt.figure()
        plt.contourf(T,D,eve_percent[m,:,:])
        plt.colorbar(label='% of key Eve learned')
        plt.title("Eve's knowledge of key (<n>="+'{:.1f}'.format(mean_photon_number[m])+')')
        plt.xlabel('Transmittance')
        plt.ylabel('Drift')
        fig.savefig('/scratch/lac468/bb84_combined_new_graphs/eve_percent_trans_drift_'+str(keylengths[k])+'_'+'{:.1f}'.format(mean_photon_number[m])+'.png')
        plt.close()
        
        
        "QBER vs transmittance, drift, mpn"
        T,D=np.meshgrid(transmittance,drift)
        fig=plt.figure()
        plt.contourf(T,D,QBER[m,:,:])
        plt.colorbar(label='QBER')
        plt.title("QBER (<n>="+'{:.1f}'.format(mean_photon_number[m])+')')
        plt.xlabel('Transmittance')
        plt.ylabel('Drift')
        fig.savefig('/scratch/lac468/bb84_combined_new_graphs/QBER_trans_drift_'+str(keylengths[k])+'_'+'{:.1f}'.format(mean_photon_number[m])+'.png')
        plt.close()
        """
    
    #intermediate index values
    mpn_q1=len(mean_photon_number)//4
    drift_q1=len(drift)//4
    trans_q1=len(transmittance)//4
    mpn_mid=len(mean_photon_number)//2
    drift_mid=len(drift)//2
    trans_mid=len(transmittance)//2
    mpn_q3=mpn_q1*3
    drift_q3=drift_q1*3
    trans_q3=trans_q1*3
    """
    "P(n=0) vs transmittance graph"
    fig=plt.figure()
    plt.plot(transmittance,P0[mpn_mid,:],'b')    
    plt.errorbar(transmittance,Pb0[mpn_mid,drift_mid,:],yerr=Pb0_std[mpn_mid,drift_mid,:],color='r')
    plt.errorbar(transmittance,P0_false[mpn_mid,drift_mid,:],yerr=P0_false_std[mpn_mid,drift_mid,:],color='g')
    plt.title('Combined InR & PNS: Expected vs simulated n=0 pulses (<n>='+'{:.1f}'.format(mean_photon_number[mpn_mid])+")")
    plt.xlabel('Channel Transmittance')
    plt.ylabel('Proportion of pulses Bob receives with n=0')
    plt.legend(('Expected outcome: No Eve','Simulated PNS attack','Simulated with no eve'))
    fig.savefig('/scratch/lac468/bb84_combined_new_graphs/P0_trans_'+str(keylengths[k])+'_'+'{:.2f}'.format(drift[drift_mid])+'_'+'{:.1f}'.format(mean_photon_number[mpn_mid])+'.png')
    plt.close()
    
    "P(n=0) vs drift"
    fig=plt.figure()
    plt.errorbar(drift,Pb0[mpn_mid,:,trans_mid],yerr=Pb0_std[mpn_mid,:,trans_mid],color='r')
    plt.errorbar(drift,P0_false[mpn_mid,:,trans_mid],yerr=P0_false_std[mpn_mid,:,trans_mid],color='g')
    plt.errorbar(drift,Pb0[mpn_q1,:,trans_mid],yerr=Pb0_std[mpn_q1,:,trans_mid])
    plt.errorbar(drift,P0_false[mpn_q1,:,trans_mid],yerr=P0_false_std[mpn_q1,:,trans_mid])
    plt.errorbar(drift,Pb0[mpn_q3,:,trans_mid],yerr=Pb0_std[mpn_q3,:,trans_mid])
    plt.errorbar(drift,P0_false[mpn_q3,:,trans_mid],yerr=P0_false_std[mpn_q3,:,trans_mid])
    plt.legend(('True 1','False 1','True 2','False 2','True 3','False 3'))
    plt.title('Combined InR & PNS: Expected vs simulated n=0 pulses (<n>='+'{:.1f}'.format(mean_photon_number[mpn_q1])+","+'{:.1f}'.format(mean_photon_number[mpn_mid])+","+'{:.1f}'.format(mean_photon_number[mpn_q3])+")")
    plt.xlabel('Channel Drift')
    plt.ylabel('Proportion of pulses Bob receives with n=0')
    #plt.legend(('Expected outcome: No Eve','Simulated PNS attack','Simulated with no eve'))
    fig.savefig('/scratch/lac468/bb84_combined_new_graphs/P0_drift_'+str(keylengths[k])+'_'+'{:.2f}'.format(drift[drift_mid])+'_'+'{:.1f}'.format(mean_photon_number[mpn_mid])+'.png')
    plt.close()
    
    "QBER vs drift"
    fig=plt.figure()
    plt.errorbar(drift,QBER[mpn_q1,:,trans_mid],yerr=QBER_std[mpn_q1,:,trans_mid])
    plt.errorbar(drift,QBER_false[mpn_q1,:,trans_mid],yerr=QBER_false_std[mpn_q1,:,trans_mid])
    plt.errorbar(drift,QBER[mpn_mid,:,trans_mid],yerr=QBER_std[mpn_mid,:,trans_mid])
    plt.errorbar(drift,QBER_false[mpn_mid,:,trans_mid],yerr=QBER_false_std[mpn_mid,:,trans_mid])
    plt.errorbar(drift,QBER[mpn_q3,:,trans_mid],yerr=QBER_std[mpn_q3,:,trans_mid])
    plt.errorbar(drift,QBER_false[mpn_q3,:,trans_mid],yerr=QBER_false_std[mpn_q3,:,trans_mid])
    legend=('True,<n>='+'{:.1f}'.format(mean_photon_number[mpn_q1]),
            'False,<n>='+'{:.1f}'.format(mean_photon_number[mpn_q1]),
            'True,<n>='+'{:.1f}'.format(mean_photon_number[mpn_mid]),
            'False,<n>='+'{:.1f}'.format(mean_photon_number[mpn_mid]),
            'True,<n>='+'{:.1f}'.format(mean_photon_number[mpn_q3]),
            'False,<n>='+'{:.1f}'.format(mean_photon_number[mpn_q3]),)
    plt.legend(legend)
    plt.title('Combined InR & PNS: QBER')
    plt.xlabel('Channel Drift')
    plt.ylabel('QBER')
    fig.savefig('/scratch/lac468/bb84_combined_new_graphs/QBER_drift_'+str(keylengths[k])+'_'+'{:.2f}'.format(drift[drift_mid])+'_'+'{:.1f}'.format(mean_photon_number[mpn_mid])+'.png')
    plt.close()
    
    
    "QBER vs mpn"
    fig=plt.figure()
    plt.errorbar(mean_photon_number,QBER[:,drift_q1,trans_mid],yerr=QBER_std[:,drift_q1,trans_mid])
    plt.errorbar(mean_photon_number,QBER_false[:,drift_q1,trans_mid],yerr=QBER_false_std[:,drift_q1,trans_mid])
    plt.errorbar(mean_photon_number,QBER[:,drift_mid,trans_mid],yerr=QBER_std[:,drift_mid,trans_mid])
    plt.errorbar(mean_photon_number,QBER_false[:,drift_mid,trans_mid],yerr=QBER_false_std[:,drift_mid,trans_mid])
    plt.errorbar(mean_photon_number,QBER[:,drift_q3,trans_mid],yerr=QBER_std[:,drift_q3,trans_mid])
    plt.errorbar(mean_photon_number,QBER_false[:,drift_q3,trans_mid],yerr=QBER_false_std[:,drift_q3,trans_mid])
    legend=('True,drift='+'{:.2f}'.format(drift[drift_q1]),
            'False,drift='+'{:.2f}'.format(drift[drift_q1]),
            'True,drift='+'{:.2f}'.format(drift[drift_mid]),
            'False,drift='+'{:.2f}'.format(drift[drift_mid]),
            'True,drift='+'{:.2f}'.format(drift[drift_q3]),
            'False,drift='+'{:.2f}'.format(drift[drift_q3]),)
    plt.legend(legend,loc='lower right')
    plt.title('Combined InR & PNS: QBER')
    plt.xlabel('Mean Photon Number')
    plt.ylabel('QBER')
    fig.savefig('/scratch/lac468/bb84_combined_new_graphs/QBER_mpn_'+str(keylengths[k])+'_'+'{:.2f}'.format(transmittance[trans_mid])+'.png')
    plt.close()
    
    
    "QBER vs transmittance"
    fig=plt.figure()
    plt.errorbar(transmittance,QBER[mpn_mid,drift_mid,:],yerr=QBER_std[mpn_mid,drift_mid,:])
    plt.errorbar(transmittance,QBER_false[mpn_mid,drift_mid,:],yerr=QBER_false_std[mpn_mid,drift_mid,:])
    plt.errorbar(transmittance,QBER[mpn_q1,drift_mid,:],yerr=QBER_std[mpn_q1,drift_mid,:])
    plt.errorbar(transmittance,QBER_false[mpn_q1,drift_mid,:],yerr=QBER_false_std[mpn_q1,drift_mid,:])
    plt.errorbar(transmittance,QBER[mpn_q3,drift_mid,:],yerr=QBER_std[mpn_q3,drift_mid,:])
    plt.errorbar(transmittance,QBER_false[mpn_q3,drift_mid,:],yerr=QBER_false_std[mpn_q3,drift_mid,:])
    legend=('True,<n>='+'{:.1f}'.format(mean_photon_number[mpn_q1]),
            'False,<n>='+'{:.1f}'.format(mean_photon_number[mpn_q1]),
            'True,<n>='+'{:.1f}'.format(mean_photon_number[mpn_mid]),
            'False,<n>='+'{:.1f}'.format(mean_photon_number[mpn_mid]),
            'True,<n>='+'{:.1f}'.format(mean_photon_number[mpn_q3]),
            'False,<n>='+'{:.1f}'.format(mean_photon_number[mpn_q3]),)
    plt.legend(legend)
    plt.title('Combined InR & PNS: QBER')
    plt.xlabel('Channel Transmittance')
    plt.ylabel('QBER')
    fig.savefig('/scratch/lac468/bb84_combined_new_graphs/QBER_trans_'+str(keylengths[k])+'_'+'{:.2f}'.format(drift[drift_mid])+'_'+'{:.1f}'.format(mean_photon_number[mpn_mid])+'.png')
    plt.close()
    """
    
    "Eve's information vs drift"
    fig=plt.figure()
    plt.errorbar(drift,eve_percent[mpn_q1,:,trans_mid],yerr=eve_percent_std[mpn_q1,:,trans_mid])
    plt.errorbar(drift,eve_percent[mpn_mid,:,trans_mid],yerr=eve_percent_std[mpn_mid,:,trans_mid])
    plt.errorbar(drift,eve_percent[mpn_q3,:,trans_mid],yerr=eve_percent_std[mpn_q3,:,trans_mid])
    legend=('{:.1f}'.format(mean_photon_number[mpn_q1]),'{:.1f}'.format(mean_photon_number[mpn_mid]),'{:.1f}'.format(mean_photon_number[mpn_q3]))
    plt.legend(legend,title='<n>')
    plt.xlabel('Drift')
    plt.ylabel('% of key Eve learned')
    plt.title("Eve's knowledge of key (T="+'{:.1f}'.format(transmittance[trans_mid])+')')
    fig.savefig('/scratch/lac468/bb84_combined_new_graphs/eve_percent_drift_'+str(keylengths[k])+'.png')
    plt.close()
    
    
    fig=plt.figure()
    plt.plot(drift,EI[mpn_mid,:,trans_mid])
    #plt.errorbar(drift,eve_percent[mpn_q1,:,trans_mid],yerr=eve_percent_std[mpn_q1,:,trans_mid])
    plt.errorbar(drift,eve_percent[mpn_mid,:,trans_mid],yerr=eve_percent_std[mpn_mid,:,trans_mid])
    #plt.errorbar(drift,eve_percent[mpn_q3,:,trans_mid],yerr=eve_percent_std[mpn_q3,:,trans_mid])
    #legend=('{:.1f}'.format(mean_photon_number[mpn_q1]),'{:.1f}'.format(mean_photon_number[mpn_mid]),'{:.1f}'.format(mean_photon_number[mpn_q3]))
    legend=('Theory ('+'{:.1f}'.format(mean_photon_number[mpn_mid])+')','Actual ('+'{:.1f}'.format(mean_photon_number[mpn_mid])+')')
    plt.legend(legend,title='<n>')
    plt.xlabel('Drift')
    plt.ylabel('% of key Eve learned')
    plt.title("Eve's knowledge of key (T="+'{:.1f}'.format(transmittance[trans_mid])+')')
    fig.savefig('/scratch/lac468/bb84_combined_new_graphs/eve_percent_drift_temp_'+str(keylengths[k])+'.png')
    plt.close()
    
    "Eve's informataion vs mpn"
    fig=plt.figure()
    plt.errorbar(mean_photon_number,eve_percent[:,drift_mid,trans_mid],yerr=eve_percent_std[:,drift_mid,trans_mid])
    plt.xlabel('<n>')
    plt.ylabel('% of key Eve learned')
    plt.title("Eve's knowledge of key (T="+'{:.1f}'.format(transmittance[trans_mid])+',D='+'{:.1f}'.format(drift[drift_mid])+')')
    fig.savefig('/scratch/lac468/bb84_combined_new_graphs/eve_percent_mpn_'+str(keylengths[k])+'.png')
    plt.close()
    
    """
    "QBER InR and combined vs drift"
    fig=plt.figure()
    #expected
    plt.plot(drift,EER,'r')
    #inr
    plt.errorbar(drift,QBER_inr,yerr=QBER_inr_std,color='g')
    #combined
    plt.errorbar(drift,QBER[mpn_mid,:,-1],yerr=QBER_std[mpn_mid,:,-1],color='b')
    legend=('Expected', 'Intercept & Resend', 'I&R and PNS' )
    plt.legend(legend)
    plt.title('QBER of I&R, combined I&R/PNS attacks (transmittance='+'{:.2f}'.format(transmittance[-1])+', <n>=' +'{:.1f}'.format(mean_photon_number[mpn_mid])+')')
    plt.xlabel('Channel Drift')
    plt.ylabel('QBER')
    fig.savefig('/scratch/lac468/bb84_combined_new_graphs/QBER_inr_drift_'+str(keylengths[k])+'_'+'{:.2f}'.format(transmittance[-1])+'_'+'{:.1f}'.format(mean_photon_number[mpn_mid])+'.png')
    plt.close()  
    """
    