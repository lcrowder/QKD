"""
SB and PNS
"""

import numpy as np
import random as random
import secrets as secrets
import qutip
from scipy.interpolate import CubicSpline
import math
import sys as sys
import itertools

I=qutip.qeye(2)
H=qutip.Qobj(1/np.sqrt(2)*np.array([[1., 1.],[1., -1.]]))
zero=qutip.basis(2,0)
one=qutip.basis(2,1)
plus=zero.transform(H)
minus=one.transform(H)
CNOT=qutip.tensor(zero.proj(),I)+qutip.tensor(one.proj(),qutip.sigmax())
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

def rot_oper(theta): #rotate state about Y-axis on the bloch sphere by angle theta.
    return qutip.Qobj(np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]]))

def renyi_theory(p): #Theoretical renyi info of slutsky-brandt attack
    return np.log2( 1+ (4*p*(1-2*p))/((1-p)**2 ))

def eve_percent_theory(p): #I dont remember deriving this, and I have no record of it.
    return ((1-p)/2 + np.sqrt(p-2*p**2) )/(1-p)

#analytically derived distribution for expected photon number in pulse, based on poisson dist. and channel loss
def prob(n_mean,rbp,n):
    return rbp**n * np.exp(-n_mean*rbp)*n_mean**n/math.factorial(n)


#########################################################################################
#PARAMETERS

#Monsoon parameters
myself = sys.argv[0]  #this is reading in the first argument in the driver file, which is name of the python file that is being executed
out_dir = sys.argv[1]  #this is the 2nd argument which is output directory 
in_dir = sys.argv[2]   #this is the 3rd argument which is the input directory 
trials=int(sys.argv[3])
eve_present=bool(int(sys.argv[4]))
keylengths=int(sys.argv[5])  #this is the 4th argument which is the keylength, which is an integer so need to convert to integer 
transmittance=float(sys.argv[6])
drift=float(sys.argv[7])
mean_photon_number=float(sys.argv[8])

k=keylengths

print("myself ",myself)  #just testing that the readin went correct 
print("out_dir ",out_dir)
print("in_dir ",in_dir)
print("keylengths ",keylengths)

eve_chances_file='/scratch/lac468/eve_chances.txt'
"""
#local params
eve_chances_file='eve_chances.txt'
transmittance=0.6
mean_photon_number=0.3
drift=0.3
k=1000 #keylength
"""

#########################################################################################
#EVE ATTACK PARAMETERS

#convevnient poisson dist for particular mpn.
def pois(x):
    return prob(mean_photon_number,1,x)

#Use pre-generated data to decide the proportion of n=1,2,3 photon pulses Eve will delete.
#interpolation with eve_chances.txt works for 0.025 <= <n> <= 1.025, 0.1 <= transmittance <= 1.0
eve_chances_data=np.loadtxt(eve_chances_file)
t_points=eve_chances_data[0:20,1] #transmittance data points solved for directly with eves_chances algorithm
n_spline=round(2*mean_photon_number,1)/2 #round <n> to closest known data point 

#Load data for specific <n>, interpolate C[i] as function of transmittance with cubic spline.  
C1_points=eve_chances_data[:,2][eve_chances_data[:,0]==n_spline]
C2_points=eve_chances_data[:,3][eve_chances_data[:,0]==n_spline]
C3_points=eve_chances_data[:,4][eve_chances_data[:,0]==n_spline]
cs1=CubicSpline(t_points,C1_points)
C1=CubicSpline.__call__(cs1,transmittance)
cs2=CubicSpline(t_points,C2_points)
C2=CubicSpline.__call__(cs2,transmittance)
cs3=CubicSpline(t_points,C3_points)
C3=CubicSpline.__call__(cs3,transmittance)

prop_single_photons= (pois(1)*(1-C1)) / (1-(pois(0) + pois(1)*C1 + pois(2)*C2 + pois(3)*C3 ))
print(prop_single_photons)

if np.sin(drift)**2 < prop_single_photons/3:
    prob_error=np.sin(drift)**2/prop_single_photons
    eve_drift=0
else:
    prob_error=1/3
    eve_drift=np.arcsin(np.sqrt((np.sin(drift)**2-prop_single_photons/3)/(1-prop_single_photons)))


#########################################################################################
##############################################################
#DO BB84


results_file=out_dir+'results_'+str(eve_present)+'_'+str(keylengths)+'_'+'{:.2f}'.format(transmittance)+'_'+'{:.2f}'.format(drift)+'_'+'{:.1f}'.format(mean_photon_number)+'.txt'
with open(results_file, 'w') as f:  
    f.write("#sifted key length, num matching bits, P(n_b=0), % of the key eve knows, renyi info\n")
f.close()

for j in range(trials):
    
    print('trial '+str(j+1))
    
    alice_bits=np.zeros(k)
    alice_basis=np.zeros(k)
    bob_bits=np.zeros(k)
    bob_basis=np.zeros(k)
    eve_basis=np.zeros(k)
    eve_bits=np.zeros(k)-1
    
    #counting pulses and photon numbers
    alice_photon_number=np.zeros(k) #for attenuated laser pulse
    current_photon_number=np.zeros(k) 
    bob_photon_number=np.zeros(k) #number bob receives, to check everything worked
    poisson_lost_pulses=0
    pulses_eve_kills=0
    all_pulses=0
    pn1_after_kills=0

    for i in range(k):
        
        send_pulse=False
        while send_pulse==False:
            alice_photon_number[i]=np.random.poisson(mean_photon_number)
            all_pulses=all_pulses+1
            #photons generated according to Poisson distribution.
            if alice_photon_number[i]==0:
                poisson_lost_pulses=poisson_lost_pulses+1
                #pulse never went through b/c of Poisson.
            if alice_photon_number[i]==1:
                does_eve_kill=random.choices([1,0],[C1,1-C1])[0]
                if does_eve_kill==0:
                    send_pulse=True
                    #Eve doesn't kill so photon goes through
                if does_eve_kill==1:
                    pulses_eve_kills=pulses_eve_kills+1
                    #One photon is sent, but Eve measures the PN and kills it.
            if alice_photon_number[i] == 2:
                does_eve_kill=random.choices([1,0],[C2,1-C2])[0]
                if does_eve_kill==0:
                    send_pulse=True
                    #Eve doesn't kill so photons go through
                if does_eve_kill==1:
                    pulses_eve_kills=pulses_eve_kills+1
                    #Two photons are sent, but Eve measures the PN and kills them.
            if alice_photon_number[i] == 3:
                does_eve_kill=random.choices([1,0],[C3,1-C3])[0]
                if does_eve_kill==0:
                    send_pulse=True
                    #Eve doesn't kill so photons go through
                if does_eve_kill==1:
                    pulses_eve_kills=pulses_eve_kills+1
                    #Three photons are sent, but Eve measures the PN and kills them.
            elif alice_photon_number[i] > 3:
                send_pulse==True
                #For mpn <=1, few pulses have n>3, and there is no need for eve to do anything with these.
    
    
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
      
        #Eve eavesdrops
        current_photon_number[i]=alice_photon_number[i]
        
        #If n=1, do SB attack.
        if current_photon_number[i]==1:
            pn1_after_kills=pn1_after_kills+1
            
            control_basis=rot_oper(np.pi/8)
            C=np.sqrt(1-2*prob_error)
            S=np.sqrt(2*prob_error)
            target_qubit=((C+S)*zero + (C-S)*one).unit()  #AKA probe qubit    
            qubit=qubit.transform(control_basis.dag()) #rotate control qubit before CNOT gate to create entanglement    
            state=qutip.tensor(qubit,target_qubit) #tensor product of two qubits    
            state=state.transform(CNOT) #apply CNOT gate to qubits    
            state=state.transform(qutip.tensor(control_basis,I)) #transform contol qubit back to original alice basis
            
            #Bob receives qubit
            bob_basis[i]=alice_basis[i]
            if bob_basis[i]==0:
                current_bob_basis=I
            elif bob_basis[i]==1:
                current_bob_basis=H
            
            state=state.transform(qutip.tensor(current_bob_basis.dag(),I)) #transform into Bob's measurement basis.
            [state,bits]=measure_general(state,2) #Eve always projects onto her computational basis
            bob_bits[i]=bits[0]
            eve_bits[i]=bits[1]
            state=state.transform(qutip.tensor(current_bob_basis,I))
            
            eve_basis[i]=alice_basis[i] #eve learns basis choice after alice and bob reconcile on public channel
      
        
        #If there is more than 1 photon, eve steals one without altering the qubit (PNS).
        elif current_photon_number[i] > 1:
            current_photon_number[i]=current_photon_number[i]-1  #Eve splits one photon off
            #Eve stores qubit in quantum memory, and waits until she knows the correct basis from Alice and Bob reconcilitation
            eve_basis[i]=alice_basis[i]
            eve_bits[i]=alice_bits[i]
            qubit=qubit.transform(rot_oper(eve_drift))
            
            #Bob receives qubit
            bob_basis[i]=alice_basis[i]
            if bob_basis[i]==0:
                current_bob_basis=I
            elif bob_basis[i]==1:
                current_bob_basis=H
            
            qubit=qubit.transform(current_bob_basis.dag())
            [qubit,bits]=measure_general(qubit,1)
            bob_bits[i]=bits[0]
            qubit=qubit.transform(current_bob_basis)
            
        bob_photon_number[i]=current_photon_number[i]
        if bob_photon_number[i]==0:
            bob_bits[i]=-1
        
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
    
    bob_num_0=sum(bob_photon_number==0)
    
    "Sift key"
    lost_bits=bob_bits==-1
    bob_bits=bob_bits[lost_bits==False]
    alice_bits=alice_bits[lost_bits==False]
    alice_basis=alice_basis[lost_bits==False]
    bob_basis=bob_basis[lost_bits==False]
    eve_basis=eve_basis[lost_bits==False]
    eve_bits=eve_bits[lost_bits == False]
    
    
    bob_num_1=sum(bob_photon_number==1)
    bob_num_2=sum(bob_photon_number==2)
    basismatch=sum((alice_basis==bob_basis))
    bitmatch=sum((alice_basis==bob_basis) & (alice_bits==bob_bits))
    
    Pbn0=(poisson_lost_pulses+pulses_eve_kills+bob_num_0)/all_pulses #prob od bob photon numbers
    Pbn1=bob_num_1/all_pulses
    Pbn2=bob_num_2/all_pulses
    """
    print('Simulation QBER: '+str(1-bitmatch/basismatch))
    print('Expected QBER from drift: '+str(np.sin(drift)**2))
    print()
    print('What Bob receives:')
    print('proportion of pulses w/ n=0: '+str(Pbn0))
    print('proportion of pulses w/ n=1: '+str(Pbn1))
    print('proportion of pulses w/ n=2: '+str(Pbn2))   
    print('sum='+str(Pbn0+Pbn1+Pbn2))
    print()
    print('If no Eve, Bob expects:') 
    print('proportion of pulses w/ n=0: '+str(prob(mean_photon_number,transmittance,0)))
    print('proportion of pulses w/ n=1: '+str(prob(mean_photon_number,transmittance,1)))
    print('proportion of pulses w/ n=2: '+str(prob(mean_photon_number,transmittance,2)))
    print('sum='+str(prob(mean_photon_number,transmittance,0)+prob(mean_photon_number,transmittance,1)+prob(mean_photon_number,transmittance,2)))
    """
    
    eve_percent=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (eve_basis==alice_basis) & (eve_bits==alice_bits))/bitmatch
    print()
    print('Prop. of key eve learned: '+str(eve_percent))
    
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
    print()
    print('Renyi info: '+str(info))
    
    with open(results_file, 'a') as f:  
        f.write("%d %d %.4f %.4f %.4f\n" % (basismatch,bitmatch,Pbn0,eve_percent,info))
    f.close()