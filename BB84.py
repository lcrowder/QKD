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
keylengths=np.arange(10,501,35)
eve_present=False
loss=np.array([0,5])
loss_dev=0.05
drift=0.1  
trials=100
QBERmax=0.02
"-----------------------------"

comp_1=np.array([0.,1.])
comp_0=np.array([1.,0.])
hadamard_0= 1./math.sqrt(2.)*np.array([1.,1.])
hadamard_1= 1./math.sqrt(2.)*np.array([1.,-1.])
H=(1/math.sqrt(2))*np.array([[1., 1. ],[1. ,-1. ]])
zero=np.array([0.,0.])

detection=np.zeros((len(loss),len(keylengths)))
detection_stddev=np.zeros((len(loss),len(keylengths)))
QBER=np.zeros((len(loss),len(keylengths)))
QBER_stddev=np.zeros((len(loss),len(keylengths)))
tol=1e-08

for x in range(len(loss)):
    received_bit_prop=1/(10.0**(loss[x]/10.0))
    
    for k in range(len(keylengths)):
    
        alice_bits=np.zeros(keylengths[k])
        alice_basis=np.zeros(keylengths[k])
        bob_bits=np.zeros(keylengths[k])
        bob_basis=np.zeros(keylengths[k])
        eve_basis=np.zeros(keylengths[k])
        eve_bits=np.zeros(keylengths[k])
     
        gotchas=np.zeros(trials)
        QBERtotals=np.zeros(trials)
        
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
              
            "Remove lost bits"
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
                QBERtotals[j]=err/basismatch
            if basismatch !=0 and QBERtotals[j] > QBERmax:
                gotchas[j]=1
        
        
        QBER_stddev[x,k]=np.std(QBERtotals)
        QBER[x,k]= np.mean(QBERtotals)
        detection[x,k]= np.mean(gotchas)
        detection_stddev[x,k]=np.std(gotchas)


plt.figure(1)
plt.errorbar(keylengths,detection[0,:],detection_stddev[0,:],fmt='-o')
plt.errorbar(keylengths,detection[1,:],detection_stddev[1,:],fmt='-o')
"plt.errorbar(keylengths,detection[2,:],detection_stddev[2,:],fmt='-o')"
"plt.legend(['1db','5dB','10dB'])  "
plt.xlabel('Initial Key Length')
plt.ylabel('Probability of Creating Secure Key (QBER< )')
plt.title('BB84: Secure Key Generation at Given Key Lengths')

plt.figure(2)
plt.errorbar(keylengths,QBER[0,:],QBER_stddev[0,:],fmt='-o')
plt.errorbar(keylengths,QBER[1,:],QBER_stddev[1,:],fmt='-o')
"plt.errorbar(keylengths,QBER[2,:],QBER_stddev[2,:],fmt='-o')"
"plt.legend(['1db','5dB','10dB'])  "

plt.xlabel('Initial Key Length')
plt.ylabel('QBER')
plt.title('BB84: QBER vs Initial Key Length')


"Ultralow loss fiber: 0.164 dB/km"