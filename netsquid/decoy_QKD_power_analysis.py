# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 15:11:48 2021

@author: crowd

This code performs power analyses on variations of QKD with decoys.
Assumptions of this code:
    -Non-decoy states are sent with equal probability (As in Vince/Brit's two-freq decoy BB84 scheme)
    -Eve and Bob choose bases randomly with 1/2 probability
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

def power_2sided(n,p_0,p_a,alpha):
    #n=number of objects
    #p_0, p_a = proportions for null and alternative hypotheses
    #Solving for beta (prob of a type 2 error)
    
    if p_0==0:
        #find probability that, given Eve is eavesdropping:
            #-no qubits are flipped (if checking QBER), or 
            #-no qubits are decoys (if checking decoy rate)
        #This equals the probability Eve going undetected = power = 1-beta
        beta=stats.binom.pmf(0,n,p_a) #1-P(Eve goes undetected)
    else:
        k1=stats.binom.ppf(alpha/2,n,p_0)
        #print('stats k1 ',k1)
        k2=stats.binom.isf(alpha/2,n,p_0)
        #print('stats k2 ',k2)
        #Now to find power
        beta=stats.binom.cdf(k2,n,p_a)-stats.binom.cdf(k1-1,n,p_a)
    return(1-beta)

def power_1sided(n,p_0,p_a,alpha): #for when the hypothesis test is p_a > p_0
    if p_0==0:
        #find probability that, given Eve is eavesdropping:
            #-no qubits are flipped (if checking QBER), or 
            #-no qubits are decoys (if checking decoy rate)
        #This equals the probability Eve going undetected = power = 1-beta
        power=1-stats.binom.pmf(0,n,p_a) #1-P(Eve goes undetected)
    else:
        k=stats.binom.isf(alpha,n,p_0)
        power=stats.binom.sf(k,n,p_a)
    return power


log2_num_check_bits=np.arange(1,12)
error_rates=np.arange(0.0,0.25,0.025) #error rate due to depolarization
alpha=0.1 #probability of falsely concluding that eve is eavesdropping
min_power=0.99 #minimum probability of falsely concluding that eve is NOT present

##############################################################################################

#arrays for probability of detecting eve vs key length 
detection_probs_trad=np.zeros(len(log2_num_check_bits))
detection_probs_raw_decoy=np.zeros(len(log2_num_check_bits))
detection_probs_sifted_decoy=np.zeros(len(log2_num_check_bits))

#arrays for key length vs error rate
ncb_trad1=np.zeros(len(error_rates)) #check bits for traditional bb84 (1 sided)
ncb_trad2=np.zeros(len(error_rates)) #check bits for traditional bb84 (2 sided)
ncb_raw_decoy1=np.zeros(len(error_rates)) #check bits for decoy bb84 raw key comparison (1 sided)
ncb_raw_decoy2=np.zeros(len(error_rates)) #check bits for decoy bb84 raw key comparison (2 sided)
ncb_sifted_decoy1=np.zeros(len(error_rates)) #check bits for decoy bb84 sifted key comparison (1 sided)

#QBERS and DRs when no errors:
QBER_0_noerr=0 #error rate when there is no Eve and no error (null hypothesis)
QBER_a=0.25 #error rate when there is Eve and no error (alternate hypothesis)
raw_DR_0_noerr=1/6 #raw key decoy rate with no errors and no eve (null hypothesis)
raw_DR_a=5/24 #raw key decoy rate with no errors and eve (alternate hypothesis)
sifted_DR_0_noerr=0 #sifted key decoy rate with no errors and no eve (null hypothesis)
sifted_DR_a=1/12 #sifted key decoy rate with no errors and eve (alternate hypothesis)

#calculate power per keylength for different protocols
for i in range(len(log2_num_check_bits)):
    num_check_bits=int(2**log2_num_check_bits[i])
    detection_probs_trad[i]=power_1sided(num_check_bits,QBER_0_noerr,QBER_a,alpha)
    detection_probs_raw_decoy[i]=power_1sided(num_check_bits, raw_DR_0_noerr, raw_DR_a, alpha)
    detection_probs_sifted_decoy[i]=power_1sided(num_check_bits, sifted_DR_0_noerr, sifted_DR_a, alpha)
       
#Traditional: calculate number of bits to check in order to have power >= 0.99
for x in range(len(error_rates)):
    #For the code below, I want to iterate over keylengths, and test if it is sufficient to have a power above some value.
    QBER_0=error_rates[x] #error rate when there is no Eve (null hypothesis)
    #QBER_a=1/2*EER + 1/4 #error rate when there is Eve and error (alternate hypothesis)
    
    #1 sided power analysis
    for keylengths1 in range(100,5500,100): #Checking in multiples of 100
        power=power_1sided(keylengths1,QBER_0,QBER_a,alpha)
        if power >= min_power:
            k_temp1=keylengths1
            for keylengths1_refined in range(k_temp1-99,k_temp1+1): #Checking one-by-one in the range where k is.
                power=power_1sided(keylengths1_refined,QBER_0,QBER_a,alpha)
                if power >= min_power:
                    ncb_trad1[x]=keylengths1_refined
                    break  
            break
    
    #2 sided power analysis
    for keylengths2 in range(100,5500,100): #Checking in multiples of 100
        power=power_2sided(keylengths2,QBER_0,QBER_a,alpha)
        if power >= min_power:
            k_temp=keylengths2
            for keylengths2_refined in range(k_temp-99,k_temp+1): #Checking one-by-one in the range where k is.
                power=power_2sided(keylengths2_refined,QBER_0,QBER_a,alpha)
                if power >= min_power:
                    ncb_trad2[x]=keylengths2_refined
                    break  
            break


###############################################################################
#Decoy states

#Decoy, raw key comparison: calculate number of bits to check in order to have power >= 0.99
for x in range(len(error_rates)):
    DR_0=1/6+1/6*error_rates[x] #expected decoy rate when there is no Eve (null hypothesis)
    
    for keylengths2 in range(1000,200000,1000): #Checking in multiples of 100
        power=power_2sided(keylengths2,DR_0,raw_DR_a,alpha)
        if power >= min_power:
            k_temp=keylengths2
            for keylengths2_refined in range(k_temp-990,k_temp+10,10): #Checking one-by-one in the range where k is.
                power=power_2sided(keylengths2_refined,DR_0,raw_DR_a,alpha)
                if power >= min_power:
                    ncb_raw_decoy2[x]=keylengths2_refined
                    break  
            break
    for keylengths1 in range(1000,200000,1000): #Checking in multiples of 100
        power=power_1sided(keylengths1,DR_0,raw_DR_a,alpha)
        if power >= min_power:
            k_temp1=keylengths1
            for keylengths1_refined in range(k_temp1-990,k_temp1+10,10): #Checking one-by-one in the range where k is.
                power=power_1sided(keylengths1_refined,DR_0,raw_DR_a,alpha)
                if power >= min_power:
                    ncb_raw_decoy1[x]=keylengths1_refined
                    break  
            break

#Decoy, sifted key comparison: calculate number of bits to check in order to have power >= 0.99
for x in range(len(error_rates)):
    #For the code below, I want to iterate over keylengths, and test if it is sufficient to have a power above some value.
    DR_0=1/3*error_rates[x] #expected decoy rate when there is no Eve (null hypothesis)
    
    for keylengths1 in range(100,15000,100): #Checking in multiples of 100
        power=power_1sided(keylengths1,DR_0,sifted_DR_a,alpha)
        if power >= min_power:
            k_temp1=keylengths1
            for keylengths1_refined in range(k_temp1-99,k_temp1+1): #Checking one-by-one in the range where k is.
                power=power_1sided(keylengths1_refined,DR_0,sifted_DR_a,alpha)
                if power >= min_power:
                    ncb_sifted_decoy1[x]=keylengths1_refined
                    break  
            break

#plot probability of detecting eve vs key length 
fig=plt.figure()
plt.plot(2**log2_num_check_bits,detection_probs_trad)
plt.plot(2**log2_num_check_bits,detection_probs_sifted_decoy)
plt.plot(2**log2_num_check_bits,detection_probs_raw_decoy)
plt.legend(('Traditional BB84','Decoy BB84 (sifted key)','Decoy BB84 (raw key)'))
plt.xlabel("# of bits to check")
plt.ylabel('probability of detecting eve')
plt.xscale("log",base=2)
plt.title('Decoy BB84: Required length of key to check (no system errors)')

#plot traditional bb84 key length vs error rate, 1 and 2 sided
fig=plt.figure()
plt.plot(error_rates,ncb_trad1)
plt.plot(error_rates,ncb_trad2)
plt.legend(('1 sided','2 sided'))
plt.yscale("log",base=2)
plt.xlabel('Error rate due to depolarization')
plt.ylabel("# of bits to check")
plt.title('Trad BB84: Required length of sifted key subset to compare for 99% confidence')

fig=plt.figure()
plt.plot(error_rates,ncb_raw_decoy1)
plt.plot(error_rates,ncb_raw_decoy2)
plt.legend(('1 sided','2 sided'))
plt.yscale("log",base=2)
plt.xlabel('Error rate due to depolarization')
plt.ylabel("# of bits to check")
plt.title('Decoy BB84: Required length of raw key to check for 99% confidence')

fig=plt.figure()
plt.plot(error_rates,ncb_trad1)
plt.plot(error_rates,ncb_sifted_decoy1)
plt.plot(error_rates,ncb_raw_decoy1)
plt.legend(('Traditional BB84','Decoy BB84 (sifted key)','Decoy BB84 (raw key)'))
plt.yscale("log",base=2)
plt.xlabel('Error rate due to depolarization')
plt.ylabel("# of bits to compare/check")
plt.title('Traditional vs Decoy BB84: # of bits needed for 99% confidence')


######################################################################################
"""
Calculate the raw key length required for each protocol to produce a given final 
sifted error_free key length with 99% confidence, plotted against system error rates.
"""
#n = sample size for hypothesis test (# of check-bits)
#k = raw key length
#f = final sifted error_free key length
#p_E = probability of error

rkl_trad1=np.zeros(len(error_rates))
rkl_raw_decoy1=np.zeros(len(error_rates))
rkl_sifted_decoy1=np.zeros(len(error_rates))
final_sifted_keylength=2**8
for k in range(len(error_rates)):
    #Raw decoy check: Bob checks all qubits regardless of basis or errors, so sample size = length of raw key
    #eqn: n+f=1/2*(1-p_E)*k
    rkl_raw_decoy1[k]=max(ncb_raw_decoy1[k],2*final_sifted_keylength/(1-error_rates[k]))
    
    #sifted decoy check: half of the raw key is removed from sample size due to sifting. But errors are included into sample size
    #eqn: n=1/2*k  ,  f=1/2*(1-p_E)*k = (1-p_E)*n
    rkl_sifted_decoy1[k]=2*max(ncb_sifted_decoy1[k],final_sifted_keylength/(1-error_rates[k]))
    
    #trad bb84: sample size is separate from final sifted key since bits are publicly compared
    #note that a subset of the sifted key is removed and compared, and THEN errors in final key are tossed out
    #eqn: f = (1/2*k - n)*(1-p_E)  =>  k = 2*(n + f/(1-p_E))
    rkl_trad1[k]=2*(ncb_trad1[k]+final_sifted_keylength/(1-error_rates[k]))
    

fig=plt.figure()
plt.plot(error_rates,rkl_trad1)
plt.plot(error_rates,rkl_sifted_decoy1)
plt.plot(error_rates,rkl_raw_decoy1)
plt.legend(('Traditional BB84','Decoy BB84 (Bob checks sifted key)','Decoy BB84 (Bob checks raw key)'))
plt.yscale("log",base=2)
plt.xlabel('Error rate due to depolarization')
plt.ylabel("raw key length")
plt.title('Raw key length needed for 99% confidence in '+str(final_sifted_keylength)+'-bit, errorless sifted key')

########################

raw_keylengths=2**np.arange(6,13)
dim=(len(raw_keylengths),len(error_rates))
fsl_raw_decoy1=np.zeros(dim)
fsl_sifted_decoy1=np.zeros(dim)
fsl_trad1=np.zeros(dim)

for k in range(len(raw_keylengths)):
    for p in range(len(error_rates)):
        if ncb_raw_decoy1[p] > raw_keylengths[k]:
            fsl_raw_decoy1[k,p]=0 #if sample size required is larger than available bits, no key can be made with 99% confidence
        else:
            fsl_raw_decoy1[k,p]=raw_keylengths[k]*(1-error_rates[p])
        
        if ncb_sifted_decoy1[p] > 1/2*raw_keylengths[k]:
            fsl_sifted_decoy1[k,p]=0
        else:
            fsl_sifted_decoy1[k,p]=1/2*raw_keylengths[k]*(1-error_rates[p])
        
        if ncb_trad1[p] > 1/2*raw_keylengths[k]:
            fsl_trad1[k,p]=0
        else:
            fsl_trad1[k,p]=(1/2*raw_keylengths[k]-ncb_trad1[p]) *(1-error_rates[p]) 

e_index=4
fig=plt.figure()
plt.plot(raw_keylengths,fsl_trad1[:,e_index])
plt.plot(raw_keylengths,fsl_sifted_decoy1[:,e_index])
#plt.plot(raw_keylengths,fsl_raw_decoy1[:,e_index])
#plt.legend(('Traditional BB84','Decoy BB84 (Bob checks sifted key)','Decoy BB84 (Bob checks raw key)'))
plt.legend(('Traditional BB84','Decoy BB84 (Bob checks sifted key)'))
plt.yscale("log",base=2)
plt.xscale("log",base=2)
plt.xlabel('raw key length')
plt.ylabel('longest possible secret key length')
plt.title('Best Secret Key Length with 99% security (error rate='+str(error_rates[e_index])+')')

rk_index=-1
fig=plt.figure()
plt.plot(error_rates,fsl_trad1[rk_index,:])
plt.plot(error_rates,fsl_sifted_decoy1[rk_index,:])
#plt.plot(error_rates,fsl_raw_decoy1[rk_index,:])
#plt.legend(('Traditional BB84','Decoy BB84 (Bob checks sifted key)','Decoy BB84 (Bob checks raw key)'))
plt.legend(('Traditional BB84','Decoy BB84 (Bob checks sifted key)'))
plt.yscale("log",base=2)
plt.xlabel('error rate')
plt.ylabel('longest possible secret key length')
plt.title('Best Secret Key Length with 99% security (raw key length='+str(raw_keylengths[rk_index])+')')

