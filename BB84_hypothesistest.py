# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 09:04:29 2019

@author: crowd
"""
import numpy as np
import random as random
import secrets as secrets
import math
import matplotlib.pyplot as plt
import scipy.stats as stat

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


"VARIABLES TO MANIPULATE"
"-----------------------------"
keylengths=np.array([20,30,40,50,60,70,80,90,100])
eve_present=False
loss=0
loss_dev=0.05
drift=0.2  
trials=1000
alpha=np.array([0.01,0.05,0.1])
"-----------------------------"

comp_1=np.array([0.,1.])
comp_0=np.array([1.,0.])
hadamard_0= 1./math.sqrt(2.)*np.array([1.,1.])
hadamard_1= 1./math.sqrt(2.)*np.array([1.,-1.])
H=(1/math.sqrt(2))*np.array([[1., 1. ],[1. ,-1. ]])
zero=np.array([0.,0.])
tol=1e-08
received_bit_prop=1/(10.0**(loss/10.0))


QBERtotals=np.zeros((len(keylengths),trials,len(alpha)))
QBERstd=np.zeros((len(keylengths),len(alpha)))
QBERavg=np.zeros((len(keylengths),len(alpha)))
safe=np.zeros((len(keylengths),len(alpha)))

for a in range(len(alpha)):
    
    for k in range(len(keylengths)):
    
        for j in range(trials):
            
            basismatch=0
            err=0
            
            alice_bits=np.zeros(keylengths[k])
            alice_basis=np.zeros(keylengths[k])
            bob_bits=np.zeros(keylengths[k])
            bob_basis=np.zeros(keylengths[k])
            eve_basis=np.zeros(keylengths[k])
            eve_bits=np.zeros(keylengths[k])
            
            for i in range(keylengths[k]):
                
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
            stderror=math.sqrt(QBERtotals[k,j,a]*(1-QBERtotals[k,j,a])/keylengths[k])
            if stderror==0.0:
                stderror=tol
            z=(QBERtotals[k,j,a]-(math.sin(drift))**2)/stderror
            pval=1-stat.norm.cdf(z)
            if pval > alpha[a]:
                safe[k,a]=safe[k,a]+1
        
        QBERstd[k,a]=np.std(QBERtotals[k,:,a])
        QBERavg[k,a]=np.mean(QBERtotals[k,:,a])
    
saferate=safe/trials

plt.figure(1)
plt.plot(keylengths,saferate[:,0],'o')
plt.plot(keylengths,saferate[:,1],'*')
plt.plot(keylengths,saferate[:,2],'s')
plt.xlabel("Key Length")
plt.ylabel("Rate of determined Safety")
plt.legend(['alpha=0.01','0.05','0.1'])

plt.figure(2)
plt.errorbar(keylengths,QBERavg[:,0],QBERstd[:,0],fmt='-o')
plt.errorbar(keylengths,QBERavg[:,1],QBERstd[:,1],fmt='-o')
plt.errorbar(keylengths,QBERavg[:,2],QBERstd[:,2],fmt='-o')
plt.xlabel("Key Length")
plt.ylabel("QBER")
plt.legend(['alpha=0.01','0.05','0.1'])