# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 17:25:27 2021

@author: crowd
"""

"""
Entangling probe attack (Slutsky-Brandt attack) with qutip package. qutip package in temp_env anaconda environment.
"""

import numpy as np
import random as random
import secrets as secrets
import qutip
import sys as sys

myself = sys.argv[0]  #this is reading in the first argument in the driver file, which is name of the python file that is being executed
out_dir = sys.argv[1]  #this is the 2nd argument which is output directory 
in_dir = sys.argv[2]   #this is the 3rd argument which is the input directory 
trials=int(sys.argv[3])
eve_present=bool(int(sys.argv[4]))
keylengths=int(sys.argv[5])  #this is the 4th argument which is the keylength, which is an integer so need to convert to integer 
transmittance=float(sys.argv[6])
drift=float(sys.argv[7])
prob_error=float(sys.argv[8]) #error probability Eve is willing to create.

print("myself ",myself)  #just testing that the readin went correct 
print("out_dir ",out_dir)
print("in_dir ",in_dir)
print("keylengths ",keylengths)
k=keylengths

I=qutip.qeye(2)
H=qutip.Qobj(1/np.sqrt(2)*np.array([[1., 1.],[1., -1.]]))
zero=qutip.basis(2,0)
one=qutip.basis(2,1)
plus=zero.transform(H)
minus=one.transform(H)
CNOT=qutip.tensor(zero.proj(),I)+qutip.tensor(one.proj(),qutip.sigmax())

def measure(basis,qubit):
    projection_0=basis.dag()*zero.proj()*basis #0 in measurement basis
    projection_1=basis.dag()*one.proj()*basis #1 in measurement basis
    qubit=random.choices([zero.transform(basis),one.transform(basis)],[qutip.expect(projection_0,qubit),qutip.expect(projection_1,qubit)])
    bit=abs(round(np.real(one.overlap((basis*qubit)[0])))) #project the qubit onto the 1 state in basis. If Eve measured 1 => <1|1>=1, if Eve measured 0 => <0|1> = 0.
    return [qubit,bit]

def measure_entangled(state): #measures in computational basis (must transform if for different basis)
    state_00=qutip.tensor(zero,zero)
    state_01=qutip.tensor(zero,one)
    state_10=qutip.tensor(one,zero)
    state_11=qutip.tensor(one,one)
    
    state=random.choices([state_00,state_01,state_10,state_11],
                         [qutip.expect(qutip.tensor(zero.proj(),zero.proj()),state),
                          qutip.expect(qutip.tensor(zero.proj(),one.proj()),state),
                          qutip.expect(qutip.tensor(one.proj(),zero.proj()),state),
                          qutip.expect(qutip.tensor(one.proj(),one.proj()),state)])[0]
    bit_0=abs(round(np.real(state_10.overlap((state)))))+abs(round(np.real(state_11.overlap((state)))))
    bit_1=abs(round(np.real(state_01.overlap((state)))))+abs(round(np.real(state_11.overlap((state)))))
    return [state,bit_0,bit_1]

def rot_oper(theta): #rotate state about Y-axis on the bloch sphere by angle theta.
    return qutip.Qobj(np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]]))

results_file=out_dir+'results_'+str(eve_present)+'_'+str(keylengths)+'_'+'{:.2f}'.format(transmittance)+'_'+'{:.2f}'.format(drift)+'_'+'{:.2f}'.format(prob_error)+'.txt'
with open(results_file, 'w') as f:  
    f.write("#basismatch, bitmatch, proportion of key eve knows, renyi information\n")
f.close()

for j in range(trials):
    
    print("trial "+str(j+1))
    
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
       
        #Eve eavesdrops
        control_basis=rot_oper(np.pi/8)
        target_basis=I
        C=np.sqrt(1-2*prob_error)
        S=np.sqrt(2*prob_error)
        target_qubit=((C+S)*zero + (C-S)*one).unit()  #AKA probe qubit
        
        #Eve applies transformations tot entangle her ancillary qubit
        qubit=qubit.transform(control_basis.dag()) #rotate control qubit before CNOT gate to create entanglement
        state=qutip.tensor(qubit,target_qubit) #tensor product of two qubits
        state=state.transform(CNOT) #apply CNOT gate to qubits
        state=state.transform(qutip.tensor(control_basis,I)) #transform contol qubit back to original alice basis
        
        #Errors: qubit drifts in Alice-Bob channel.
        state=state.transform(qutip.tensor(rot_oper(drift),I))
        
        #Bob receives qubit
        #bob_basis[i]=secrets.choice([0,1])
        bob_basis[i]=alice_basis[i]
        if bob_basis[i]==0:
            current_bob_basis=I
        elif bob_basis[i]==1:
            current_bob_basis=H
        
        state=state.transform(qutip.tensor(current_bob_basis.dag(),I)) #transform into Bob's measurement basis.
        
        [state,bob_bits[i],eve_bits[i]]=measure_entangled(state) #Eve always projects onto her computational basis
        state=state.transform(qutip.tensor(current_bob_basis,I))
        
        eve_basis[i]=alice_basis[i] #eve learns basis choice after alice and bob reconcile on public channel
    
    
    basismatch=sum((alice_basis==bob_basis))
    bitmatch=sum((alice_basis==bob_basis) & (alice_bits==bob_bits))
    eves_info=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (eve_basis==alice_basis) & (eve_bits==alice_bits))/bitmatch
    QBER=1-bitmatch/basismatch
    
    print(basismatch)
    print(QBER)
    print(eves_info)
    
    #calculate probabilities of eve and bob bit combinations, used to find eve's information gain.
    Pb0=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (bob_bits==0))/bitmatch
    Pb1=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (bob_bits==1))/bitmatch
    Pe0=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (eve_bits==0))/bitmatch
    Pe1=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (eve_bits==1))/bitmatch
    Pbe00=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (eve_bits==bob_bits) & (eve_bits==0))/bitmatch
    Pbe01=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (bob_bits==0) & (eve_bits==1))/bitmatch
    Pbe10=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (bob_bits==1) & (eve_bits==0))/bitmatch
    Pbe11=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (eve_bits==bob_bits) & (eve_bits==1))/bitmatch
    
    info=-np.log2(Pb0**2+Pb1**2)+ Pe0*np.log2((Pbe00/Pe0)**2+(Pbe10/Pe0)**2) + Pe1*np.log2((Pbe01/Pe1)**2 + (Pbe11/Pe1)**2)
    
    with open(results_file, 'a') as f:  
        f.write("%d %d %.4f %.4f \n" % (basismatch, bitmatch, eves_info, info))
    f.close()