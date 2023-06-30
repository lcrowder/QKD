# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 13:31:51 2020

@author: crowd
"""
"""
This code simulates BB84 with an attenuated laser pulse, and Eve performs the photon number
splitting attack. A low mean photon number causes significant loss in photons. So for a keylength k,
Alice sends information until k pulses all with photon number > 0 make it through.

PNS Attack 1 (current): Eve has a lossless channel and uses a statistical formula to calculate how
many single photon pulses to delete to imitate channel loss. She does NOT perform intercept-and-resend
attack on undeleted single photon pulses. On all other pulses (n>1), Eve does QND measurement to split
off one photon without altering the qubit, thus knowing that bit of the key.
*Note this method is only good for a low mean photon number such as 0.1

The way this is coded results in a sifted key exactly as long as set in the variable 'keylengths',
as long as eve_present=True. In other words, all the losses from Poisson distribution and Eve deleting
pulses takes place in the 'while send_pulse' loop. 


"""

import numpy as np
import random as random
import secrets as secrets
import math
import scipy.stats as stat
import sys as sys

##########################################################################
###################### I N P U T S ########################################

"""
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

print("myself ",myself)  #just testing that the readin went correct 
print("out_dir ",out_dir)
print("in_dir ",in_dir)
print("keylengths ",keylengths)
"""

#local PC parameters
trials=1
eve_present=True
keylengths=1000
transmittance=0.9
drift=0.0
mean_photon_number=1.0


######################################################################
########## D E F I N E   Q U B I T S   &   M E A S U R E M E N T ########

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

##########################################################################
############ E V E ' S   P A R A M E T E R S ################################################

#Use pre-generated data to decide the proportion of n=1,2,3 photon pulses Eve will delete.

eve_chances_file='eve_chances.txt'
#eve_chances_file='/scratch/lac468/eve_chances.txt'
eve_chances_data=np.loadtxt(eve_chances_file)
mask=(eve_chances_data[:,0]==mean_photon_number) & (eve_chances_data[:,1]==transmittance)
C1=eve_chances_data[:,2][mask][0]
C2=eve_chances_data[:,3][mask][0]
C3=eve_chances_data[:,4][mask][0]

#If Eve is not present, she does not kill any qubits
if eve_present==False:
    C1=0
    C2=0
    C3=0

##########################################################################
############ D O   B B 8 4 ################################################

"""
results_file=out_dir+'results_'+str(eve_present)+'_'+str(keylengths)+'_'+'{:.2f}'.format(transmittance)+'_'+'{:.2f}'.format(drift)+'_'+'{:.1f}'.format(mean_photon_number)+'.txt'
with open(results_file, 'w') as f:  
    f.write("#sifted key length, % of the key eve knows, P(n_b=0), P(n_b=1), P(n_b=2)\n")
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
    eve_bits=np.zeros(k)
    
    #counting pulses and photon numbers
    alice_photon_number=np.zeros(k) #for attenuated laser pulse
    current_photon_number=np.zeros(k) 
    bob_photon_number=np.zeros(k) #number bob receives, to check everything worked
    poisson_lost_pulses=0
    pulses_eve_kills=0
    all_pulses=0
    bob_num_0=0
    bob_num_1=0
    bob_num_2=0
    eves_info=np.zeros(k)
    
    alice_1=0
    
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
                alice_1=alice_1+1
                does_eve_kill=random.choices([1,0],[C1,1-C1])[0]
                if does_eve_kill==0:
                    send_pulse=True
                if does_eve_kill==1:
                    pulses_eve_kills=pulses_eve_kills+1
                    #One photon is sent, but Eve measures the PN and kills it.
            if alice_photon_number[i] == 2:
                does_eve_kill=random.choices([1,0],[C2,1-C2])[0]
                if does_eve_kill==0:
                    send_pulse=True
                if does_eve_kill==1:
                    pulses_eve_kills=pulses_eve_kills+1
                    #Two photons are sent, but Eve measures the PN and kills them.
            if alice_photon_number[i] == 3:
                does_eve_kill=random.choices([1,0],[C3,1-C3])[0]
                if does_eve_kill==0:
                    send_pulse=True
                if does_eve_kill==1:
                    pulses_eve_kills=pulses_eve_kills+1
                    #Three photons are sent, but Eve measures the PN and kills them.
            elif alice_photon_number[i] > 3:
                send_pulse==True
                
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
        
        "Eve eavesdrops: PNS"
        
        temp_transmittance=transmittance
        current_photon_number[i]=alice_photon_number[i]
        if eve_present==True:
            """
            #intercept-and-resend attack
            if alice_photon_number[i]==1:
                eve_basis[i]=secrets.choice([0,1])
                [eve_bits[i],qubit]=measure(eve_basis[i],qubit)
            """
            #If there is more than 1 photon, eve steals one without altering the qubit.
            if alice_photon_number[i] > 1:
                current_photon_number[i]=current_photon_number[i]-1
                eve_basis[i]=alice_basis[i]
                #Eve splits one photon off
                [eve_bits[i],qubit]=measure(eve_basis[i],qubit)
                #eves_info[i]=1
            #Assume Eve has the ability to use a lossless channel
            temp_transmittance=1.0
            drift=0.0
        
        "ERROR: loss and drift"
        #find how many photons are lost, and delete the qubit if all photons are lost
        qubit_survival=random.choices( [1,0] , [temp_transmittance,1-temp_transmittance],k=int(current_photon_number[i]))
        bob_photon_number[i]=np.sum(qubit_survival)
        if bob_photon_number[i]==0:
            qubit=zero
            bob_num_0=bob_num_0+1
            
        #rand_drift=random.gauss(0,drift)
        rand_drift=drift
        driftmat=np.array([[math.cos(rand_drift),-1*math.sin(rand_drift)],[math.sin(rand_drift),math.cos(rand_drift)]])
        qubit=np.matmul(driftmat,qubit)
        
        "Bob receives information"
        #bob_basis[i]=secrets.choice([0,1])
        bob_basis[i]=alice_basis[i]
        [bob_bits[i],qubit]=measure(bob_basis[i],qubit)
  
        if int(bob_photon_number[i])==1:
            bob_num_1=bob_num_1+1
        if int(bob_photon_number[i])==2:
            bob_num_2=bob_num_2+1
            
            
    "Sift key"
    lost_bits=bob_bits==-1
    bob_bits=bob_bits[lost_bits==False]
    alice_bits=alice_bits[lost_bits==False]
    alice_basis=alice_basis[lost_bits==False]
    bob_basis=bob_basis[lost_bits==False]
    eve_basis=eve_basis[lost_bits==False]
    eves_info=eves_info[lost_bits==False]
    eve_bits=eve_bits[lost_bits == False]
    
    for n in range(len(bob_bits)):
        if (alice_basis[n]==bob_basis[n]) and (alice_bits[n] != bob_bits[n]):
            basismatch=basismatch+1
        elif (alice_basis[n]==bob_basis[n]) and (alice_bits[n] == bob_bits[n]):
            basismatch=basismatch+1       
            bitmatch=bitmatch+1
            if eve_present==True and eve_basis[n]==alice_basis[n] and eve_bits[n]==alice_bits[n]:
                eves_info[n]=1
           
    ESR=(math.cos(rand_drift))**2 #expected success rate
    p_val_2=stat.binom_test(bitmatch,basismatch,ESR,alternative='two-sided')
    eve_key_percent=100*sum(eves_info)/bitmatch
    
    Pb0=(poisson_lost_pulses+pulses_eve_kills+bob_num_0)/all_pulses
    Pb1=bob_num_1/all_pulses
    Pb2=bob_num_2/all_pulses
    
    print(all_pulses)
    print('What Bob receives:')
    print('proportion of pulses w/ n=0: '+str(Pb0))
    print('proportion of pulses w/ n=1: '+str(Pb1))
    print('proportion of pulses w/ n=2: '+str(Pb2))
    if eve_present==True:
        print('eve knows '+str(eve_key_percent)+'% of the key')
    print('QBER '+'{:.3f}'.format(1-bitmatch/basismatch))
    
    """
    with open(results_file, 'a') as f:  
        f.write("%d %.3f %.6f %.6f %.6f\n" % (basismatch,eve_key_percent,Pb0,Pb1,Pb2))
    f.close()
    """
#############################################################################


#analytically derived distribution for expected photon number in pulse, based on poisson dist. and channel loss
def prob(n_mean,rbp,n):
    return rbp**n * np.exp(-n_mean*rbp)*n_mean**n/math.factorial(n)

print('If no Eve, Bob expects:') 
print('proportion of pulses w/ n=0: '+str(prob(mean_photon_number,transmittance,0)))
print('proportion of pulses w/ n=1: '+str(prob(mean_photon_number,transmittance,1)))
print('proportion of pulses w/ n=2: '+str(prob(mean_photon_number,transmittance,2)))
