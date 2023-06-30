# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 19:23:22 2020

@author: crowd
"""

"(loss, drift, current_siftedkey, QBERavg, QBERstd, current_safe[0], current_safe[1], current_safe[2])"
"Eve is present here"

#plot:
#Keylength vs sifted key
#surface plot: pick keylength, (1) plot QBER over drift and loss, (2) plot each saferate over drift and loss
#surface plot: pick loss, (1) QBER vs keylength and drift, (2) saferate vs keylength and drift

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

key_250=np.array([np.zeros(8)])


with open('C:\\Users\\crowd\\Documents\\Python Scripts\\QKD\\err_results_eve\\bb84_errortests_results_250.txt','r') as f:
    
    for line in f:
        new_line=np.fromstring(line, dtype=float, sep=' ')
        if len(new_line)==8:
            key_250=np.append(key_250, np.array([new_line]), axis=0)
    f.close()

key_250=key_250[1:10202,:]
loss_1=key_250[::101 ,0]
drift_1=key_250[0:101, 1]
D,L=np.meshgrid(drift_1,loss_1)

Q=np.zeros((len(loss_1),len(drift_1))) #QBERavg
Qstd=np.zeros((len(loss_1),len(drift_1))) #QBER std
SK=np.zeros((len(loss_1),len(drift_1))) #Sifted key
S1=np.zeros((len(loss_1),len(drift_1))) #saferate, alpha=0.01
S2=np.zeros((len(loss_1),len(drift_1))) #saferate, alpha=0.05
S3=np.zeros((len(loss_1),len(drift_1))) #saferate, alpha =0.10

for i in range(101):
    SK[i,:]=np.transpose( key_250[101*i:101*(i+1),2]    )
    Q[i,:]=np.transpose( key_250[101*i:101*(i+1),3]    )
    Qstd[i,:]=np.transpose( key_250[101*i:101*(i+1),4]    )
    S1[i,:]=np.transpose( key_250[101*i:101*(i+1),5]    )
    S2[i,:]=np.transpose( key_250[101*i:101*(i+1),6]    )
    S3[i,:]=np.transpose( key_250[101*i:101*(i+1),7]    )

fig1= plt.figure()
ax=plt.axes(projection='3d')
surf1=ax.plot_surface(D,L,Q, cmap='viridis')
ax.set_xlabel("Drift (radians)")
ax.set_ylabel("Loss (dB)")
ax.set_zlabel("QBER")
plt.title("Loss and drift Relation to QBER, at key length 250")
fig1.colorbar(surf1)
ax.view_init(30,230)
plt.savefig('QBER_eve_250.png')

fig2= plt.figure()
ax=plt.axes(projection='3d')
surf2=ax.plot_surface(D,L,Qstd, cmap='viridis')
ax.set_xlabel("Drift (radians)")
ax.set_ylabel("Loss (dB)")
ax.set_zlabel("QBER standard deviation")
fig2.colorbar(surf2)
plt.title("Loss and drift Relation to QBER std, at key length 250")
ax.view_init(30,210)

plt.savefig('QBERstd_eve_250.png')

fig3= plt.figure()
ax=plt.axes(projection='3d')
surf3=ax.plot_surface(D,L,SK, cmap='viridis')
ax.set_xlabel("Drift (radians)")
ax.set_ylabel("Loss (dB)")
ax.set_zlabel("sifted key length")
fig3.colorbar(surf3)
plt.title("Loss and drift Relation to sifted key length, at key length 250")
ax.view_init(30,60)

plt.savefig('siftedkeylength_eve_250.png')

fig4= plt.figure()
ax=plt.axes(projection='3d')
surf4=ax.plot_surface(D,L,S1, cmap='viridis')
ax.set_xlabel("Drift (radians)")
ax.set_ylabel("Loss (dB)")
ax.set_zlabel("safe rate")
plt.title("Loss and drift Relation to safe rate (alpha=0.01), at key length 250")
fig4.colorbar(surf4)
ax.view_init(60,210)

plt.savefig('saferate0.01_eve_250.png')

fig5= plt.figure()
ax=plt.axes(projection='3d')
surf5=ax.plot_surface(D,L,S2, cmap='viridis')
ax.set_xlabel("Drift (radians)")
ax.set_ylabel("Loss (dB)")
ax.set_zlabel("safe rate")
plt.title("Loss and drift Relation to safe rate (alpha=0.05), at key length 250")
fig5.colorbar(surf5)
ax.view_init(60,210)

plt.savefig('saferate0.05_eve_250.png')

fig6= plt.figure()
ax=plt.axes(projection='3d')
surf6=ax.plot_surface(D,L,S3, cmap='viridis')
ax.set_xlabel("Drift (radians)")
ax.set_ylabel("Loss (dB)")
ax.set_zlabel("safe rate")
plt.title("Loss and drift Relation to safe rate (alpha=0.10), at key length 250")
fig6.colorbar(surf6)
ax.view_init(60,210)
plt.savefig('saferate0.10_eve_250.png')

