# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 13:30:44 2021

@author: crowd
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 08:40:51 2021

@author: crowd
"""
"""
Entangling probe attack (Slutsky-Brandt attack) with qutip package. qutip package in temp_env anaconda environment.
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

def renyi_theory(p): #Theoretical renyi info of slutsky-brandt attack
    return np.log2( 1+ (4*p*(1-2*p))/((1-p)**2 ))

def eve_percent_theory(p): #I dont remember deriving this, and I have no record of it.
    return ((1-p)/2 + np.sqrt(p-2*p**2) )/(1-p)

p=np.arange(0,0.55,0.05)
#p=np.array([0.25])

drift_ab=0.2 #drift in alice-bob quantum channel
drift_e=0.0 #drift in eve's quantum channel

QBER=np.zeros(len(p))
eves_info=np.zeros(len(p))
renyi_info=np.zeros(len(p))
eves_info_theory=np.zeros(len(p))
test_info=np.zeros(len(p))



for j in range(len(p)):
    
    k=1000 #keylength
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
        """
        print("alice basis "+str(alice_basis[i]))
        print('alice bit '+str(alice_bits[i]))
        """
        if alice_basis[i]==0:
            current_basis=I
        elif alice_basis[i]==1:
            current_basis=H
        qubit=qutip.basis(2,int(alice_bits[i])).transform(current_basis)
        """
        print()
        print("qubit")
        print(qubit)
        """
        #Eve eavesdrops
        prob_error=p[j] #error probability Eve is willing to create.
        control_basis=rot_oper(np.pi/8)
        target_basis=I
        C=np.sqrt(1-2*prob_error)
        S=np.sqrt(2*prob_error)
        target_qubit=((C+S)*zero + (C-S)*one).unit()  #AKA probe qubit
        """
        print()
        print("target_qubit")
        print(target_qubit)
        """
        qubit=qubit.transform(control_basis.dag()) #rotate control qubit before CNOT gate to create entanglement
        """
        print()
        print('rotated qubit')
        print(qubit)
        """
        state=qutip.tensor(qubit,target_qubit) #tensor product of two qubits
        """
        print()
        print('2 qubit state')
        print(state)
        """
        state=state.transform(CNOT) #apply CNOT gate to qubits
        """
        print()
        print('state after cnot')
        print(state)
        """
        state=state.transform(qutip.tensor(control_basis,I)) #transform contol qubit back to original alice basis
        """
        print()
        print('state transformed back to alice basis')
        print(state)
        """
        #Bob receives qubit
        
        #bob_basis[i]=secrets.choice([0,1])
        bob_basis[i]=alice_basis[i]
        if bob_basis[i]==0:
            current_bob_basis=I
        elif bob_basis[i]==1:
            current_bob_basis=H
        """
        print()
        print('bob basis '+str(bob_basis[i]))
        """
        state=state.transform(qutip.tensor(current_bob_basis.dag(),I)) #transform into Bob's measurement basis.
        """
        print()
        print('state transformed to bob basis')
        print(state)
        """
        [state,bob_bits[i],eve_bits[i]]=measure_entangled(state) #Eve always projects onto her computational basis
        state=state.transform(qutip.tensor(current_bob_basis,I))
        
        eve_basis[i]=alice_basis[i] #eve learns basis choice after alice and bob reconcile on public channel
        """
        print()
        print('bob measured state')
        print(state)
        print(' bob bit '+str(bob_bits[i]))
    
        print()
        print('eve measured state')
        print(state)
        print('eve bit '+str(eve_bits[i]))
        """
    
    
    basismatch=sum((alice_basis==bob_basis))
    bitmatch=sum((alice_basis==bob_basis) & (alice_bits==bob_bits))
    eves_info[j]=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (eve_basis==alice_basis) & (eve_bits==alice_bits))/bitmatch
    QBER[j]=1-bitmatch/basismatch
    
    renyi_info[j]=renyi_theory(prob_error)
    eves_info_theory[j]=eve_percent_theory(prob_error)
    
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
    test_info[j]=info

fig1=plt.figure()
plt.title('QBER')
plt.plot(p,QBER,'*')

fig2=plt.figure()
plt.title('percent of key')
plt.plot(p,eves_info_theory)
plt.plot(p,eves_info,'o')
plt.legend(('eve percent theory','eve percent'))

fig3=plt.figure()
plt.title('renyi info')
plt.plot(p,renyi_info)
plt.plot(p,test_info,'o')
plt.legend(('theory','simulated'))