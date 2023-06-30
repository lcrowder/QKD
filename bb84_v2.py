# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 09:34:30 2020

@author: crowd
"""
"""
This code simulates the basic BB84 protocol, with intercept and resend attack. Eve has an errorless
channel and she can manipulate the qubits how she pleases (cause fake drift/loss)
"""


import numpy as np
import random as random
import secrets as secrets
import math
import scipy.stats as stat
import sys as sys
"""
myself = sys.argv[0]  #this is reading in the first argument in the driver file, which is name of the python file that is being executed
out_dir = sys.argv[1]  #this is the 2nd argument which is output directory 
in_dir = sys.argv[2]   #this is the 3rd argument which is the input directory 
trials=int(sys.argv[3])
eve_present=bool(int(sys.argv[4]))
keylengths=int(sys.argv[5])  #this is the 4th argument which is the keylength, which is an integer so need to convert to integer 
transmittance=float(sys.argv[6])
drift=float(sys.argv[7])

print("myself ",myself)  #just testing that the readin went correct 
print("out_dir ",out_dir)
print("in_dir ",in_dir)
print("keylengths ",keylengths)
"""

#local PC parameters
trials=1
eve_present=True
keylengths=1000
transmittance=0.1
drift=0.0


#chance_eve_measures=min(4*math.sin(drift)**2,1.0)
chance_eve_measures=transmittance #Being lazy here and using transmittance b/c its already built in to monsoon bash scripts
transmittance=1.0
eve_drift=drift

#if chance_eve_measures==1.0:
    #eve_drift=math.asin(math.sqrt(2*math.sin(drift)**2-1/2))

print('drift='+str(drift))
print('C='+str(chance_eve_measures))
print('eve_drift='+str(eve_drift))

comp_1=np.array([0.,1.])
comp_0=np.array([1.,0.])
hadamard_0= 1./math.sqrt(2.)*np.array([1.,1.])
hadamard_1= 1./math.sqrt(2.)*np.array([1.,-1.])
H=(1/math.sqrt(2))*np.array([[1., 1. ],[1. ,-1. ]])
zero=np.array([0.,0.])
tol=1e-08
k=keylengths

def measure(basis,qub):
    if basis == 0 and (abs(qub[0])>tol or abs(qub[1])>tol):
        measured=random.choices( [comp_0,comp_1] , [qub[0]**2,qub[1]**2] )
        qub=measured[0]
        "Using random.choices creates a list with one element which is an array. Because of this I had to "
        "make qubit just the first element of the list. There is probably a MUCH better way to do this."
        
        if abs(qub[0]-comp_0[0])<=tol and abs(qub[1]-comp_0[1])<=tol:
            bits=0
        elif abs(qub[0]-comp_1[0])<=tol and abs(qub[1]-comp_1[1])<=tol:
            bits=1   
    elif basis==1 and (abs(qub[0])>tol or abs(qub[1])>tol):
        q_hadamard=np.matmul(H, qub)
        measured=random.choices( [hadamard_0, hadamard_1] , [q_hadamard[0]**2 , q_hadamard[1]**2] )
        qub=measured[0]
        if abs(qub[0]-hadamard_0[0])<=tol and abs(qub[1]-hadamard_0[1])<=tol:
            bits=0
        elif abs(qub[0]-hadamard_1[0])<=tol and abs(qub[1]-hadamard_1[1])<=tol:
            bits=1
    elif abs(qub[0])<=tol and abs(qub[1])<=tol:
        bits=-1
    return [bits,qub]
def plog(p1,p2):
    if p1==0.0 or p2==0.0:
        return 0.0
    else:
        return p1*np.log2(p2)

"""
results_file=out_dir+'results_'+str(eve_present)+'_'+str(keylengths)+'_'+'{:.2f}'.format(chance_eve_measures)+'_'+'{:.2f}'.format(drift)+'.txt'
with open(results_file, 'w') as f:  
    f.write("#sifted key length, num matching bits, p-value, % sifted key eve knows \n")
f.close()
"""
for j in range(trials):
    
    basismatch=0
    bitmatch=0
    
    alice_bits=np.zeros(k)
    alice_basis=np.zeros(k)
    bob_bits=np.zeros(k)
    bob_basis=np.zeros(k)
    eve_basis=np.zeros(k)
    eve_bits=np.zeros(k)-1
    
    for i in range(k):
        
        #print('qubit #'+str(i+1))
        
        alice_bits[i]=secrets.choice([0,1])
        alice_basis[i]=secrets.choice([0,1])
        if alice_bits[i]==0 and alice_basis[i]==0:
            qubit=comp_0
        elif alice_bits[i]==1 and alice_basis[i]==0:
            qubit=comp_1
        elif alice_bits[i]==0 and alice_basis[i]==1:
            qubit= hadamard_0
        elif alice_bits[i]==1 and alice_basis[i]==1:
            qubit= hadamard_1
        
        #print('qubit is '+str(qubit))
        
        "Eve eavesdrops"
        temp_transmittance=transmittance
        rand_drift=drift
        if eve_present==True:
            #print('eve is here')
            does_eve_measure=random.choices([0,1],[1-chance_eve_measures,chance_eve_measures])[0]
            if does_eve_measure==1:
                #print('eve is gonna measure')
                eve_basis[i]=secrets.choice([0,1])
                [eve_bits[i],qubit]=measure(eve_basis[i],qubit)
                #print('eve measured and qubit = '+str(qubit))
            elif does_eve_measure==0:
                #print('eve is NOT gonna measure')
                eve_bits[i]==-1
            temp_transmittance=1.0
            rand_drift=eve_drift
            
        "ERROR: loss and drift"
        qubit=random.choices( [qubit,zero] , [temp_transmittance,1-temp_transmittance] )
        qubit=qubit[0]
        
        #rand_drift=random.gauss(0,drift)
        driftmat=np.array([[math.cos(rand_drift),-1*math.sin(rand_drift)],[math.sin(rand_drift),math.cos(rand_drift)]])
        qubit=np.matmul(driftmat,qubit)
        
        "Bob receives information"
        bob_basis[i]=alice_basis[i]
        #bob_basis[i]=secrets.choice([0,1])
        [bob_bits[i],qubit]=measure(bob_basis[i],qubit)
        #print('bob measures qubit and it = '+str(qubit))
  
    "Sift key"
    lost_bits=bob_bits==-1
    bob_bits=bob_bits[lost_bits==False]
    alice_bits=alice_bits[lost_bits==False]
    alice_basis=alice_basis[lost_bits==False]
    bob_basis=bob_basis[lost_bits==False]
    eve_bits=eve_bits[lost_bits==False]
    eve_basis=eve_basis[lost_bits==False]

    basismatch=sum((alice_basis==bob_basis))
    bitmatch=sum((alice_basis==bob_basis) & (alice_bits==bob_bits))
    eves_info=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (eve_basis==alice_basis) & (eve_bits==alice_bits))
    
    ESR=(math.cos(drift))**2 #expected success rate
    p_val_2=stat.binom_test(bitmatch,basismatch,ESR,alternative='two-sided')
    eve_key_percent=100*eves_info/bitmatch
    
    """
    print('sifted key length = '+str(bitmatch))
    print('expected error = '+'{:.3f}'.format(1-ESR))
    print('actual QBER '+'{:.3f}'.format(1-bitmatch/basismatch))
    
    print('eve knows '+'{:.2f}'.format(eve_key_percent)+'% of the key')
    """

    Pb0=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (bob_bits==0))/bitmatch
    Pb1=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (bob_bits==1))/bitmatch
    Pe0=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (eve_bits==0))/bitmatch
    Pe1=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (eve_bits==1))/bitmatch
    Pbe00=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (eve_bits==bob_bits) & (eve_bits==0))/bitmatch
    Pbe01=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (bob_bits==0) & (eve_bits==1))/bitmatch
    Pbe10=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (bob_bits==1) & (eve_bits==0))/bitmatch
    Pbe11=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (eve_bits==bob_bits) & (eve_bits==1))/bitmatch
    
    renyi_info=-np.log2(Pb0**2+Pb1**2)+ Pe0*np.log2((Pbe00/Pe0)**2+(Pbe10/Pe0)**2) + Pe1*np.log2((Pbe01/Pe1)**2 + (Pbe11/Pe1)**2)
    shannon_info=-plog(Pb0,Pb0) - plog(Pb1,Pb1) + (plog(Pbe00,Pbe00/Pe0) + plog(Pbe10,Pbe10/Pe0)) + (plog(Pbe01,Pbe01/Pe1) + plog(Pbe11,Pbe11/Pe1))
    
    """
    with open(results_file, 'a') as f:  
        f.write("%d %d %.3f %.3f %.3f %.3f\n" % (basismatch,bitmatch,p_val_2,eve_key_percent,shannon_info,renyi_info))
    f.close()
    """
