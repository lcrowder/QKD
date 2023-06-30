# -*- coding: utf-8 -*-
"""
Created on Sun Oct 11 17:11:54 2020

@author: crowd
"""
"""
BB84, with PNS and Intercept and resend combined.
"""

import numpy as np
import random as random
import secrets as secrets
import math
import scipy.stats as stat
import sys as sys
from scipy.interpolate import CubicSpline

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

eve_chances_file='/scratch/lac468/eve_chances.txt'

"""
eve_chances_file='eve_chances.txt'
#local PC parameters
trials=1
eve_present=True
keylengths=1000
transmittance=0.7
drift=math.pi/6-0.02
mean_photon_number=0.1

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

#If Eve is not present, she does not kill any qubits
if eve_present==False:
    C1=0
    C2=0
    C3=0

#analytically derived distribution for expected photon number in pulse, based on poisson dist. and channel loss
def prob(n_mean,rbp,n):
    return rbp**n * np.exp(-n_mean*rbp)*n_mean**n/math.factorial(n)
#convevnient poisson dist for particular mpn.
def pois(x):
    return prob(mean_photon_number,1,x)
   
prop_single_photons= (pois(1)*(1-C1)) / (1-(pois(0) + pois(1)*C1 + pois(2)*C2 + pois(3)*C3 ))
print(prop_single_photons)
if prop_single_photons==0.0:
    temp_chance=1.0
else:
    temp_chance=(4*math.sin(drift)**2)/prop_single_photons

chance_eve_measures=min(temp_chance,1.0)
eve_drift=0.0

if chance_eve_measures==1.0:
    eve_drift=math.acos(math.cos(drift)/math.sqrt(1-prop_single_photons/4))
    #eve_drift=math.asin(math.sqrt((math.sin(drift)**2-1/4)/(1-prop_single_photons/2)))

print(eve_drift)

pn1_after_kills=0
n1_measures=0

##########################################################################
############ D O   B B 8 4 ################################################

"""
results_file=out_dir+'results_'+str(eve_present)+'_'+str(keylengths)+'_'+'{:.2f}'.format(transmittance)+'_'+'{:.2f}'.format(drift)+'_'+'{:.1f}'.format(mean_photon_number)+'.txt'
with open(results_file, 'w') as f:  
    f.write("#sifted key length, num matching bits, % of the key eve knows, P(n_b=0)\n")
f.close()
"""
for j in range(trials):
    
    alice_bits=np.zeros(k)
    alice_basis=np.zeros(k)
    bob_bits=np.zeros(k)
    bob_basis=np.zeros(k)
    eve_basis=np.zeros(k)
    eve_bits=np.zeros(k)-1 #Eve does not know the key bits by default
    
    #counting pulses and photon numbers
    alice_photon_number=np.zeros(k) #for attenuated laser pulse
    current_photon_number=np.zeros(k) 
    bob_photon_number=np.zeros(k) #number bob receives, to check everything worked
    poisson_lost_pulses=0
    pulses_eve_kills=0
    all_pulses=0
    
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
        
        "Eve eavesdrops: PNS and I&R"
        
        #set temporary transmittance and drift values, since Eve may change these
        temp_transmittance=transmittance
        temp_drift=drift
        current_photon_number[i]=alice_photon_number[i]
        
        if eve_present==True:
            
            #Assume Eve has the ability to use a lossless channel
            temp_transmittance=1.0
            temp_drift=eve_drift
            
            #intercept-and-resend attack
            if alice_photon_number[i]==1:
                pn1_after_kills=pn1_after_kills+1
                does_eve_measure=random.choices([0,1],[1-chance_eve_measures,chance_eve_measures])[0]
                
                if does_eve_measure==1:
                    #Eve measures so she randomly selects basis and does basic measurement of qubit
                    n1_measures=n1_measures+1
                    eve_basis[i]=secrets.choice([0,1])
                    [eve_bits[i],qubit]=measure(eve_basis[i],qubit)
            
            #If there is more than 1 photon, eve steals one without altering the qubit (PNS).
            if alice_photon_number[i] > 1:
                current_photon_number[i]=current_photon_number[i]-1  #Eve splits one photon off
                #Eve stores qubit in quantum memory, and waits until she knows the correct basis from Alice and Bob reconcilitation
                eve_basis[i]=alice_basis[i]
                [eve_bits[i],qubit]=measure(eve_basis[i],qubit)
                
        "ERROR: loss and drift"
        #find how many photons are lost, and delete the qubit if all photons are lost
        qubit_survival=random.choices( [1,0] , [temp_transmittance,1-temp_transmittance],k=int(current_photon_number[i]))
        bob_photon_number[i]=np.sum(qubit_survival)
        if bob_photon_number[i]==0:
            qubit=zero
        
        #apply drift to qubit
        driftmat=np.array([[math.cos(temp_drift),-1*math.sin(temp_drift)],[math.sin(temp_drift),math.cos(temp_drift)]])
        qubit=np.matmul(driftmat,qubit)
        
        "Bob receives information"
        #bob_basis[i]=secrets.choice([0,1])
        bob_basis[i]=alice_basis[i]
        [bob_bits[i],qubit]=measure(bob_basis[i],qubit)
            
    
    bob_num_0=sum(bob_photon_number==0)
    
    "Sift key"
    lost_bits=bob_bits==-1
    bob_bits=bob_bits[lost_bits==False]
    alice_bits=alice_bits[lost_bits==False]
    alice_basis=alice_basis[lost_bits==False]
    bob_basis=bob_basis[lost_bits==False]
    eve_basis=eve_basis[lost_bits==False]
    eve_bits=eve_bits[lost_bits == False]
           
    basismatch=sum((alice_basis==bob_basis))
    bitmatch=sum((alice_basis==bob_basis) & (alice_bits==bob_bits))
    bob_num_1=sum((bob_photon_number==1))
    bob_num_2=sum((bob_photon_number==2))
    eves_info=sum((alice_basis==bob_basis) & (alice_bits==bob_bits) & (eve_basis==alice_basis) & (eve_bits==alice_bits))
    
    ESR=(math.cos(drift))**2 #expected success rate
    p_val_2=stat.binom_test(bitmatch,basismatch,ESR,alternative='two-sided')
    
    eve_key_percent=100*eves_info/bitmatch
    
    Pb0=(poisson_lost_pulses+pulses_eve_kills+bob_num_0)/all_pulses
    Pb1=bob_num_1/all_pulses
    Pb2=bob_num_2/all_pulses
    """
    print('QBER '+'{:.3f}'.format(1-bitmatch/basismatch))
    print('EER '+'{:.3f}'.format(1-ESR))
    
    
    if eve_present==True:
        print('eve knows '+str(eve_key_percent)+'% of the key')
        print('prop_single_photons theory ' + str(prop_single_photons))
        print('prop_single_photons simulated ' + str(pn1_after_kills/keylengths))
        print('chance_eve_measures '+ str(chance_eve_measures))
        print('prop of single photons measured '+str(n1_measures/pn1_after_kills))
        print('prop of total qubits measured '+str(n1_measures/keylengths))
        print('prop of measured qubits that caused error '+str((basismatch-bitmatch)/n1_measures))
    
    print('What Bob receives:')
    print('proportion of pulses w/ n=0: '+str(Pb0))
    print('proportion of pulses w/ n=1: '+str(Pb1))
    print('proportion of pulses w/ n=2: '+str(Pb2))   
    """
    """
    with open(results_file, 'a') as f:  
        f.write("%d %d %.3f %.6f\n" % (basismatch,bitmatch,eve_key_percent,Pb0))
    f.close()
    """
#############################################################################
"""
print('If no Eve, Bob expects:') 
print('proportion of pulses w/ n=0: '+str(prob(mean_photon_number,transmittance,0)))
print('proportion of pulses w/ n=1: '+str(prob(mean_photon_number,transmittance,1)))
print('proportion of pulses w/ n=2: '+str(prob(mean_photon_number,transmittance,2)))
"""