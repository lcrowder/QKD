# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 13:13:34 2020

@author: crowd
"""
"""
This code performs basic BB84 protocol, excluding the loss of qubits and the basis reconciliation
between Alice and Bob, since these two factors only reduce the keylength. Thus the keylength
variable here is the portion of the sifted keylength that Alice and Bob compare. Using the expected
success rate (ESR), the code performs a two-sided binomial test and outputs the p-value.
"""

import scipy.stats as stat
import numpy as np
import random as random
import secrets as secrets
import math
import sys as sys

myself = sys.argv[0]  #this is reading in the first argument in the driver file, which is name of the python file that is being executed
out_dir = sys.argv[1]  #this is the 2nd argument which is output directory 
in_dir = sys.argv[2]   #this is the 3rd argument which is the input directory 
trials=int(sys.argv[3])
eve_present=bool(int(sys.argv[4]))
keylengths=int(sys.argv[5])  #this is the 4th argument which is the keylength, which is an integer so need to convert to integer 
loss=float(sys.argv[6])
drift=float(sys.argv[7])

print("myself ",myself)  #just testing that the readin went correct 
print("out_dir ",out_dir)
print("in_dir ",in_dir)
print("keylengths ",keylengths)

comp_1=np.array([0.,1.])
comp_0=np.array([1.,0.])
hadamard_0= 1./math.sqrt(2.)*np.array([1.,1.])
hadamard_1= 1./math.sqrt(2.)*np.array([1.,-1.])
H=(1/math.sqrt(2))*np.array([[1., 1. ],[1. ,-1. ]])
zero=np.array([0.,0.])
tol=1e-08
received_bit_prop=1/(10.0**(loss/10.0))
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


results_file=out_dir+'binom_'+str(eve_present)+'_'+str(keylengths)+'_'+'{:.2f}'.format(drift)+'.txt'
with open(results_file, 'w') as f:  
    f.write("#num matching bits, p-value for 2-tailed test\n")
f.close()

for j in range(trials):
    
    print('trial '+str(j))
    
    basismatch=0
    bitmatch=0
    
    alice_bits=np.zeros(k)
    alice_basis=np.zeros(k)
    bob_bits=np.zeros(k)
    bob_basis=np.zeros(k)
    eve_basis=np.zeros(k)
    eve_bits=np.zeros(k)
    
    for i in range(k):
        
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
        
        "Eve eavesdrops"
        if eve_present==True:
            eve_basis[i]=secrets.choice([0,1])
            [eve_bits[i],qubit]=measure(eve_basis[i],qubit)
        
        "ERROR: loss and drift"
        #qubit=random.choices( [qubit,zero] , [received_bit_prop,1-received_bit_prop] )
        #qubit=qubit[0]
        
        #rand_drift=random.gauss(0,drift)
        rand_drift=drift
        driftmat=np.array([[math.cos(rand_drift),-1*math.sin(rand_drift)],[math.sin(rand_drift),math.cos(rand_drift)]])
        qubit=np.matmul(driftmat,qubit)
        
    
        "Bob receives information"
        #bob_basis[i]=secrets.choice([0,1])
        bob_basis[i]=alice_basis[i] #Since we're working with an already sifted key.
        [bob_bits[i],qubit]=measure(bob_basis[i],qubit)
    """
    "Sift key"
    lost_bits=bob_bits==-1
    bob_bits=bob_bits[lost_bits==False]
    alice_bits=alice_bits[lost_bits==False]
    alice_basis=alice_basis[lost_bits==False]
    bob_basis=bob_basis[lost_bits==False]
    """
    for n in range(len(bob_bits)):
        if (alice_basis[n]==bob_basis[n]) and (alice_bits[n] != bob_bits[n]):
            basismatch=basismatch+1
        elif (alice_basis[n]==bob_basis[n]) and (alice_bits[n] == bob_bits[n]):
            basismatch=basismatch+1       
            bitmatch=bitmatch+1
    

    #In this case we are only looking at an already sifted key. 
    #We want an exact, known sifted keylength to perform a binomial test.
    #Therefore, we don't do any loss simulation and bob automatically matches alice's basis.
    #All that we do is let Eve manipulate the qubit going from alice to bob.    
    """
    print(alice_basis)
    print(bob_basis)
    print(eve_basis)
    print(alice_bits)
    print(bob_bits)    
    print(bitmatch)
    print(basismatch)
    """
    
    ESR=(math.cos(rand_drift))**2 #expected success rate
    p_val_2=stat.binom_test(bitmatch,basismatch,ESR,alternative='two-sided')
    
    
    with open(results_file, 'a') as f:  
        f.write("%d %.3f\n" % (bitmatch, p_val_2))
    f.close()
    