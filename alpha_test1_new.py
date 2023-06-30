# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 17:54:47 2020

@author: crowd
"""

import numpy as np
import random as random
import secrets as secrets
import math
"import matplotlib.pyplot as plt"
import scipy.stats as stat
import sys as sys

myself = sys.argv[0]  #this is reading in the first argument in the driver file, which is name of the python file that is being executed
out_dir = sys.argv[1]  #this is the 2nd argument which is output directory 
in_dir = sys.argv[2]   #this is the 3rd argument which is the input directory 

keylengths=int(sys.argv[3])  #this is the 4th argument which is the keylength, which is an integer so need to convert to integer 
alpha=sys.argv[4].split(',')
alpha=np.array(alpha)
alpha=alpha.astype(np.float)
trials=int(sys.argv[5])
eve_present=bool(int(sys.argv[6]))
loss=float(sys.argv[7])
loss_dev=float(sys.argv[8])
drift=float(sys.argv[9])

print("myself ",myself)  #just testing that the readin went correct 
print("out_dir ",out_dir)
print("in_dir ",in_dir)
print("keylengths ",keylengths)
print("alpha ",alpha)

with open(out_dir+'bb84_hyp_results.txt', 'w') as f:
    
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
    
    comp_1=np.array([0.,1.])
    comp_0=np.array([1.,0.])
    hadamard_0= 1./math.sqrt(2.)*np.array([1.,1.])
    hadamard_1= 1./math.sqrt(2.)*np.array([1.,-1.])
    H=(1/math.sqrt(2))*np.array([[1., 1. ],[1. ,-1. ]])
    zero=np.array([0.,0.])
    tol=1e-08
    received_bit_prop=1/(10.0**(loss/10.0))
    
    QBERtotals=np.zeros((keylengths,trials,len(alpha)))
    QBERstd=np.zeros((keylengths,len(alpha)))
    QBERavg=np.zeros((keylengths,len(alpha)))
    safe=np.zeros((keylengths,len(alpha)))
    
    for a in range(len(alpha)):
        print("alpha is ",alpha[a])
        
        for k in range(keylengths):
            real_k=k+5
            "start keylengths at 5. No reason to do extremely short key lengths."
            print("k is ",real_k)
            for j in range(trials):
                
                basismatch=0
                err=0
                
                alice_bits=np.zeros(real_k)
                alice_basis=np.zeros(real_k)
                bob_bits=np.zeros(real_k)
                bob_basis=np.zeros(real_k)
                eve_basis=np.zeros(real_k)
                eve_bits=np.zeros(real_k)
                
                for i in range(real_k):
                    
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
                    
                    "ERROR: loss and drift"
                    randlossprop=random.gauss(received_bit_prop,loss_dev*received_bit_prop)
                    qubit=random.choices( [qubit,zero] , [randlossprop,1-randlossprop] )
                    qubit=qubit[0]
                    
                    rand_drift=random.gauss(0,drift)
                    driftmat=np.array([[math.cos(rand_drift),-1*math.sin(rand_drift)],[math.sin(rand_drift),math.cos(rand_drift)]])
                    qubit=np.matmul(driftmat,qubit)
                    
                    "Eve eavesdrops"
                    if eve_present==True:
                        eve_basis[i]=secrets.choice([0,1])
                        [eve_bits[i],qubit]=measure(eve_basis[i],qubit)
                
                    "Bob receives information"
                    bob_basis[i]=secrets.choice([0,1])
                    [bob_bits[i],qubit]=measure(bob_basis[i],qubit)
              
                "Sift key"
                lost_bits=bob_bits==-1
                bob_bits=bob_bits[lost_bits==False]
                alice_bits=alice_bits[lost_bits==False]
                alice_basis=alice_basis[lost_bits==False]
                bob_basis=bob_basis[lost_bits==False]
                
                for n in range(len(bob_bits)):
                    if (alice_basis[n]==bob_basis[n]) and (alice_bits[n] != bob_bits[n]):
                        basismatch=basismatch+1
                        err=err+1
                    elif (alice_basis[n]==bob_basis[n]) and (alice_bits[n] == bob_bits[n]):
                        basismatch=basismatch+1   
                
                if basismatch !=0:
                    QBERtotals[k,j,a]=err/basismatch
            
                "Hypothesis Test"
                stderror=math.sqrt(QBERtotals[k,j,a]*(1-QBERtotals[k,j,a])/float(real_k))
                if stderror==0.0:
                    stderror=tol
                z=(QBERtotals[k,j,a]-(math.sin(drift))**2)/stderror
                pval=1-stat.norm.cdf(z)
                if pval > alpha[a]:
                    safe[k,a]=safe[k,a]+1
            
            QBERstd[k,a]=np.std(QBERtotals[k,:,a])
            QBERavg[k,a]=np.mean(QBERtotals[k,:,a])
            
            print(real_k,QBERavg[k,a],QBERstd[k,a])
            
        saferate=safe/trials
        safe_avg=np.mean(saferate[:,a])
        safe_std=np.std(saferate[:,a])
        print(alpha[a],safe_avg,safe_std)
    
f.close()
