# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 08:40:51 2021

@author: crowd
"""
"""
Simulation of basic BB84 using the qutip package. qutip package in temp_env anaconda environment.
"""


import numpy as np
import matplotlib.pyplot as plt
import random as random
import secrets as secrets
import qutip

I=qutip.qeye(2)
H=qutip.Qobj(1/np.sqrt(2)*np.array([[1., 1.],[1., -1.]]))

zero=qutip.basis(2,0)
one=qutip.basis(2,1)
plus=zero.transform(H)
minus=one.transform(H)

def measure(basis,qubit):
    projection_0=basis.dag()*zero.proj()*basis #0 in measurement basis
    projection_1=basis.dag()*one.proj()*basis #1 in measurement basis
    qubit=random.choices([zero.transform(basis),one.transform(basis)],[qutip.expect(projection_0,qubit),qutip.expect(projection_1,qubit)])
    bit=abs(round(np.real(one.overlap((basis*qubit)[0])))) #project the qubit onto the 1 state in basis. If Eve measured 1 => <1|1>=1, if Eve measured 0 => <0|1> = 0.
    return [qubit,bit]

def rot_oper(theta): #rotate state about Y-axis on the bloch sphere by angle theta.
    return qutip.Qobj(np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]]))

def plog(p1,p2):
    if p1==0.0 or p2==0.0:
        return 0.0
    else:
        return p1*np.log2(p2)


drift=np.arange(0,10)*np.pi/40
#drift=np.array([np.pi/16])
k=1000 #keylength

QBER=np.zeros(len(drift))
EER=np.zeros(len(drift)) #expected error rate
EEP=np.zeros(len(drift)) #expected eve percent
eve_percent=np.zeros(len(drift))
renyi_info=np.zeros(len(drift))
shannon_info=np.zeros(len(drift))

for d in range(len(drift)):
    
    EEP[d]=np.cos(drift[d])**2 / (1/2+np.cos(drift[d])**2)
    EER[d]=1/2*np.sin(drift[d])**2 + 1/4
    
    drift_operator=rot_oper(drift[d])
    print(drift_operator)
    
    alice_bits=np.zeros(k)
    alice_basis=np.zeros(k)
    bob_bits=np.zeros(k)
    bob_basis=np.zeros(k)
    eve_basis=np.zeros(k)
    eve_bits=np.zeros(k)-1
    
    for i in range(k):
    
        #Alice generates random qubit    
        alice_bits[i]=secrets.choice([0,1])
        alice_basis[i]=secrets.choice([0,1])
        
        if alice_basis[i]==0:
            current_basis=I
        elif alice_basis[i]==1:
            current_basis=H
        qubit=qutip.basis(2,int(alice_bits[i])).transform(current_basis)
        #print(qubit)
        
        #drift
        qubit=qubit.transform(drift_operator)
        #print(qubit)
        
        #Eve eavesdrops
        eve_basis[i]=secrets.choice([0,1]) #randomly select basis to use for measurement
        if eve_basis[i]==0:
            current_eve_basis=I
        elif eve_basis[i]==1:
            current_eve_basis=H
        
        #Eve measures qubit. project qubit onto orthogonal states in Eve basis, randomly choose w/ corresponding probability.
        [qubit,eve_bits[i]]=measure(current_eve_basis,qubit)
        
        #Bob receives qubit
        bob_basis[i]=secrets.choice([0,1])
        if bob_basis[i]==0:
            current_bob_basis=I
        elif bob_basis[i]==1:
            current_bob_basis=H
        
        #Bob measures qubit. project qubit onto orthogonal states in Bob basis, randomly choose w/ corresponding probability.
        [qubit,bob_bits[i]]=measure(current_bob_basis,qubit)
        
    basismatch=sum((alice_basis==bob_basis))
    bitmatch=sum((alice_basis==bob_basis) & (alice_bits==bob_bits))
    
    eve_percent[d]=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (eve_basis==alice_basis) & (eve_bits==alice_bits))/bitmatch
    QBER[d]=1-bitmatch/basismatch
    
    Pb0=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (bob_bits==0))/bitmatch
    Pb1=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (bob_bits==1))/bitmatch
    Pe0=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (eve_bits==0))/bitmatch
    Pe1=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (eve_bits==1))/bitmatch
    Pbe00=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (eve_bits==bob_bits) & (eve_bits==0))/bitmatch
    Pbe01=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (bob_bits==0) & (eve_bits==1))/bitmatch
    Pbe10=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (bob_bits==1) & (eve_bits==0))/bitmatch
    Pbe11=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (eve_bits==bob_bits) & (eve_bits==1))/bitmatch
    
    renyi_info[d]=-np.log2(Pb0**2+Pb1**2)+ Pe0*np.log2((Pbe00/Pe0)**2+(Pbe10/Pe0)**2) + Pe1*np.log2((Pbe01/Pe1)**2 + (Pbe11/Pe1)**2)
    shannon_info[d]=-plog(Pb0,Pb0) - plog(Pb1,Pb1) + (plog(Pbe00,Pbe00/Pe0) + plog(Pbe10,Pbe10/Pe0)) + (plog(Pbe01,Pbe01/Pe1) + plog(Pbe11,Pbe11/Pe1))
    
    print(basismatch)
    print(QBER[d])
    print(eve_percent[d])

fig=plt.figure(1)
plt.title('QBER')
plt.ylabel('QBER')
plt.xlabel('drift')
plt.plot(drift,QBER,'.')
plt.plot(drift,EER)
plt.legend(('Simulated','Theory'))

fig=plt.figure(2)
plt.title('Eve percent of key')
plt.ylabel('percent of key eve knows')
plt.xlabel('drift')
plt.plot(drift,eve_percent,'.')
plt.plot(drift,EEP)
plt.legend(('simulated','theory'))

fig=plt.figure(3)
plt.title('information')
plt.ylabel('entropy')
plt.xlabel('drift')
plt.plot(drift,shannon_info,'.')
plt.plot(drift,renyi_info,'*')
plt.legend(('Shannon','Renyi'))


