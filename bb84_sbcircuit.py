# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 09:54:33 2021

@author: crowd
"""

import numpy as np
import random as random
import secrets as secrets
import itertools
import matplotlib.pyplot as plt
import qutip

I=qutip.qeye(2)
H=qutip.Qobj(1/np.sqrt(2)*np.array([[1., 1.],[1., -1.]]))

zero=qutip.basis(2,0)
one=qutip.basis(2,1)
plus=zero.transform(H)
minus=one.transform(H)

CNOT=qutip.tensor(zero.proj(),I)+qutip.tensor(one.proj(),qutip.sigmax())
CNOT_12=qutip.tensor(CNOT,I)
CNOT_13=qutip.tensor(zero.proj(),I,I)+qutip.tensor(one.proj(),I,qutip.sigmax())

def rot_oper(theta): #rotate state about Y-axis on the bloch sphere by angle theta.
    return qutip.Qobj(np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]]))

def measure_general(state,num): #num=number of qubits
    eigenstates=[]
    probabilities=[]
    bit_permutations=list(itertools.product((0,1),repeat=num))
    index=[]
    for x in range(len(bit_permutations)):
        index.append(x)
        bits=bit_permutations[x]
        eigenstate=qutip.basis(2,bits[0])
        for i in range(len(bits)-1):
            eigenstate=qutip.tensor(eigenstate,qutip.basis(2,bits[i+1]))
        eigenstates.append(eigenstate)
        probabilities.append(qutip.expect(eigenstate.proj(),state))  
        #print(eigenstate)
        #print(bits)
    measured_index=random.choices(index,probabilities)[0]
    state=eigenstates[measured_index]
    out_bits=bit_permutations[measured_index]
    return [state,out_bits]

"""
drift=np.pi/6
qubit=qutip.basis(2,0)
prob_error=1/4
C=np.sqrt(1-2*prob_error)
S=np.sqrt(2*prob_error)
target_qubit=((C+S)*zero + (C-S)*one).unit()  #AKA probe qubit

state=qutip.tensor(qubit,target_qubit)

qc=qutip.qip.circuit.QubitCircuit(3)
#qc.add_gate('RY',targets=0 ,arg_value=-np.pi/4)
#qc.add_gate('CNOT',controls=0,targets=4)
#qc.add_gate('RY',targets=0,arg_value=np.pi/4)
#Basic error correction
qc.add_gate('CNOT',controls=0,targets=1)
qc.add_gate('CNOT',controls=0,targets=2)
qc.add_gate('RY',targets=0,arg_value=2*drift) #drift
qc.add_gate('RY',targets=1,arg_value=2*drift) #drift
qc.add_gate('RY',targets=2,arg_value=2*drift) #drift
sb_circuit=qutip.qip.operations.gate_sequence_product(qc.propagators())
#print(state.transform(sb_circuit))


q0=one
q1=zero
q2=zero

q=qutip.tensor(q0,q1,q2)
print(q)
q=q.transform(sb_circuit)
print(q)

e_bits=[]
indices=[]
bit_permutations=list(itertools.product((0,1),repeat=3))
for n in range(1000):
    q=qutip.tensor(q0,q1,q2)
    q=q.transform(sb_circuit)
    [q,temp_bits]=measure_general(q,3)
    e_bits.append(temp_bits)
    for m in range(len(bit_permutations)):
        if bit_permutations[m]==temp_bits:
            indices.append(m)

print('000: '+str(sum(np.array(indices)==0)))
print('001: '+str(sum(np.array(indices)==1)))
print('010: '+str(sum(np.array(indices)==2)))
print('011: '+str(sum(np.array(indices)==3)))
print('100: '+str(sum(np.array(indices)==4)))
print('101: '+str(sum(np.array(indices)==5)))
print('110: '+str(sum(np.array(indices)==6)))
print('111: '+str(sum(np.array(indices)==7)))
"""
##########################################################################################3


def renyi_theory(p): #Theoretical renyi info of slutsky-brandt attack
    return np.log2( 1+ (4*p*(1-2*p))/((1-p)**2 ))

def renyi_theory_ec_drift(p,d): #Theoretical renyi info of slutsky-brandt attack with ec and drift
    return 1+np.log2(1/2+2*p*(1-2*p)*np.cos(2*d)**2)

def eve_percent_theory(p): #I dont remember deriving this, and I have no record of it.
    return ((1-p)/2 + np.sqrt(p-2*p**2) )/(1-p)

def qber_sbec(p,d):
    return 1-p-(1-2*p)*np.cos(d)**2

p=np.arange(0,0.55,0.05)
#p=np.array([0])

drift=np.pi/8*0

QBER=np.zeros(len(p))
QBER_sbec=np.zeros(len(p))
eves_info=np.zeros(len(p))
renyi_info=np.zeros(len(p))
eves_info_theory=np.zeros(len(p))
test_info=np.zeros(len(p))
ec_drift_info=np.zeros(len(p))

for j in range(len(p)):
    
    k=1000 #keylength
    alice_bits=np.zeros(k)
    alice_basis=np.zeros(k)
    bob_bits=np.zeros(k)
    bob_basis=np.zeros(k)
    eve_basis=np.zeros(k)
    eve_bits=np.zeros(k)-1
    
    for i in range(k):
    
        "Alice generates random qubit"    
        alice_bits[i]=secrets.choice([0,1])
        alice_basis[i]=secrets.choice([0,1])
        #alice_basis[i]=0
        if alice_basis[i]==0:
            current_basis=I
        elif alice_basis[i]==1:
            current_basis=H
        qubit=qutip.basis(2,int(alice_bits[i])).transform(current_basis)
        
        "Eve eavesdrops"
        prob_error=p[j] #error probability Eve is willing to create.
        
        control_basis=rot_oper(np.pi/8)
        target_basis=I
        C=np.sqrt(1-2*prob_error)
        S=np.sqrt(2*prob_error)
        target_qubit=((C+S)*zero + (C-S)*one).unit()  #AKA probe qubit, for SB attack
        
        "ancillary error correction bits"
        q1=zero
        q2=zero 
        state=qutip.tensor(qubit,q1,q2,target_qubit) #state of entire system of qubits

        "Quantum circuit for SB attack"
        qc=qutip.qip.circuit.QubitCircuit(4)
        qc.add_gate('RY',targets=0 ,arg_value=-np.pi/4)
        qc.add_gate('CNOT',controls=0,targets=3)
        qc.add_gate('RY',targets=0,arg_value=np.pi/4)
        "Basic error correction"
        qc.add_gate('CNOT',controls=0,targets=1)
        qc.add_gate('CNOT',controls=0,targets=2)
        qc.add_gate('RY',targets=0,arg_value=2*drift) #drift
        qc.add_gate('RY',targets=1,arg_value=2*drift) #drift
        qc.add_gate('RY',targets=2,arg_value=2*drift) #drift
        sb_ec_circuit=qutip.qip.operations.gate_sequence_product(qc.propagators())
        
        state=state.transform(sb_ec_circuit) #pass state through quantum circuit
        
        "Bob measures"
        #bob_basis[i]=secrets.choice([0,1])
        bob_basis[i]=alice_basis[i]
        if bob_basis[i]==0:
            current_bob_basis=I
        elif bob_basis[i]==1:
            current_bob_basis=H
        state=state.transform(qutip.tensor(current_bob_basis.dag(),I,I,I)) #transform into Bob's measurement basis.
        
        #print(i)
        #print(state)
        bit_perm_4=list(itertools.product((0,1),repeat=4))
        #for a in range(len(bit_perm_4)):
            #print('P('+str(bit_perm_4[a])+')=', '{:.4f}'.format(abs(state[a,0])**2))
        [state,bits]=measure_general(state,4) #Eve always projects onto her computational basis
        state=state.transform(qutip.tensor(current_bob_basis,I,I,I))
        bob_bits[i]=bits[0]
        

        "Eve's measurement"
        eve_basis[i]=alice_basis[i] #eve learns basis choice after alice and bob reconcile on public channel
        state=state.transform(qutip.tensor(I,current_bob_basis.dag(),current_bob_basis.dag(),I))
        [state,bits]=measure_general(state,4)
        check_bits=bits[1:len(bits)] #q1,q2, and probe
        if sum(np.array(check_bits)) <= 1: #all 0s or two 0s and one 1
            eve_bits[i]=0
        else: #at least two 1s
            eve_bits[i]=1        
        
        #eve_bits[i]=bits[len(bits)-1] 
        
        """
        print('alice bit ',alice_bits[i])
        print('bob bit ',bob_bits[i])
        print('check bits ',check_bits)
        print('eve bit ',eve_bits[i])
        """
        
    basismatch=sum((alice_basis==bob_basis))
    bitmatch=sum((alice_basis==bob_basis) & (alice_bits==bob_bits))
    eves_info[j]=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (eve_basis==alice_basis) & (eve_bits==alice_bits))/bitmatch
    QBER[j]=1-bitmatch/basismatch
    
    renyi_info[j]=renyi_theory(prob_error)
    eves_info_theory[j]=eve_percent_theory(prob_error)
    
    #print(basismatch)
    #print(QBER)
    #print(eves_info)
    
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
    test_info[j]=info
    ec_drift_info[j]=renyi_theory_ec_drift(prob_error,drift)
    QBER_sbec[j]=qber_sbec(prob_error,drift)

fig1=plt.figure()
plt.title('QBER')
plt.plot(p,QBER,'*')
#plt.plot(p,QBER_sbec)

fig2=plt.figure()
plt.title('percent of key')
plt.plot(p,eves_info_theory)
plt.plot(p,eves_info,'o')
plt.legend(('eve percent theory','eve percent'))

fig3=plt.figure()
plt.title('renyi info')
plt.plot(p,renyi_info)
#plt.plot(p,ec_drift_info)
plt.plot(p,test_info,'o')
plt.legend(('plain SB no drift','simulated'))

renyi_qber=np.zeros(len(QBER))
for f in range(len(QBER)):
    renyi_qber[f]=renyi_theory(QBER[f])

fig4=plt.figure()
plt.title('info vs qber')
plt.plot(QBER,test_info,'*')
plt.plot(QBER,renyi_qber)

