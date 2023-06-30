# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 17:51:55 2020

@author: crowd
"""
import numpy as np
import matplotlib.pyplot as plt

data=np.array([[0, 0, 0, 0, 0]])

with open('bb84_hyp_results.txt','r') as f:
    
    for line in f:
        new_line=np.fromstring(line, dtype=float, sep=' ')
        if len(new_line)==5:
            data=np.append(data, np.array([new_line]), axis=0)
    f.close()

plt.figure(1)
plt.title("Average Sifted key length vs Initial key length")
plt.xlabel("Initial key length")
plt.ylabel("Average Sifted key length")
plt.plot(data[1:250,0], data[1:250,1],"b.")
plt.plot(data[251:500,0],data[251:500,1],"g.")
plt.plot(data[501:750,0],data[501:750,1],"r.")
plt.legend(["alpha =0.01","0.05","0.10"])
plt.savefig("siftedkeyplot.png")

plt.figure(2)
plt.xlabel("Initial key length")
plt.ylabel("QBER")
plt.title("QBER vs intitial key length (alpha=0.01)")
plt.errorbar(data[1:250,0], data[1:250,2], data[1:250,3],marker="o" )
plt.savefig("QBERplot1.png")

plt.figure(3)
plt.xlabel("Initial key length")
plt.ylabel("QBER")
plt.title("QBER vs intitial key length (alpha=0.05)")
plt.errorbar(data[251:500,0],data[251:500,2], data[251:500,3],marker="o")
plt.savefig("QBERplot2.png")

plt.figure(4)
plt.xlabel("Initial key length")
plt.ylabel("QBER")
plt.title("QBER vs intitial key length (alpha=0.10)")
plt.errorbar(data[501:750,0],data[501:750,2],data[501:750,3],marker="o")
plt.savefig("QBERplot3.png")

plt.figure(5)
plt.xlabel("Initial key length")
plt.ylabel("Safety rate")
plt.title("Statistically determined safety rate vs Initial key length")
plt.plot(data[1:250,0], data[1:250,4],"b.")
plt.plot(data[251:500,0],data[251:500,4],"g.")
plt.plot(data[501:750,0],data[501:750,4],"r.")
plt.legend(["alpha =0.01","0.05","0.10"])
plt.savefig("saferateplot.png")