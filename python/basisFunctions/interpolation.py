# -*- coding: utf-8 -*-
"""
Created on Thu May 13 14:41:19 2021

@author: benja
"""
import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt

plt.close('all')

nodes = 6;
sps = 100;

np.random.seed(1)
sig1 = np.random.rand(nodes);
sig2 = sig.resample(sig1,nodes*sps)[:-sps];

sig1x = np.linspace(0,nodes,nodes*sps, endpoint=False)[:-sps]

fig, ax = plt.subplots(1, figsize = (4, 2))
ax.plot(sig1x,sig2)
ax.set_title("Function")
labels = ['$x_0$','$x_1$','$x_2$','$x_3$','$x_4$','$x_5$']
ax.set_xticks(np.arange(0,6))
ax.set_xticklabels(labels)

fig.savefig("func.svg")
#ax.set_aspect(15)
fig, ax = plt.subplots(1, figsize = (4, 2))
ax.plot(sig1x,sig2)
ax.plot(sig1)
ax.set_title("Linear interpolation")
ax.set_xticks(np.arange(0,6))
ax.set_xticklabels(labels)

fig.savefig("func_interp1.svg")


#### pre-generate all the basis functions

I_list = np.arange(0,(nodes-1)*sps)
I_list = np.reshape(I_list,(nodes-1,sps))
#for each I
phi = np.zeros((nodes,(nodes-1)*sps));
for i,I in enumerate(I_list):
    phi[i,I] =  np.arange(sps,0,-1)/sps
    phi[i+1,I] =  np.arange(0,sps)/sps
    
fig2, ax2 = plt.subplots(1, figsize = (4, 2))
ax2.plot(sig1x,phi.T)
ax2.set_title("Basis functions")
ax2.set_xticks(np.arange(0,6))
ax2.set_xticklabels(labels)
fig2.savefig("Hat_functions.svg")

M = np.zeros((nodes,nodes))


#### Compute the M matrix #####
for i in range(nodes):
    for j in range(nodes):
        #integrade
        integral = np.sum(phi[i,:]*phi[j,:])/sps
        
        M[i,j] = integral
        
        
print(M)


#### Asemble the b matrix ###

b = np.zeros(nodes)
for i in range(nodes):
    b[i] = np.sum(sig2*phi[i,:])/sps
    
print(b.T)
xi = np.linalg.solve(M,b)

print(xi)

Phf = np.zeros((nodes-1)*sps)
for j in range(nodes):
    Phf =  Phf + phi[j,:]*xi[j]



fig3, ax3 = plt.subplots(1, figsize = (4, 2))
ax3.plot(sig1x,sig2)
ax3.plot(sig1x,Phf)
ax3.set_title("L^{2} projection")
ax3.set_xticks(np.arange(0,6))
ax3.set_xticklabels(labels)

#### tjek ####
print(f"Error: {np.sum(np.square(Phf-sig2))}")
print(f"Sum of distance: {np.sum((Phf-sig2))}")


fig3.savefig("func_interp2.svg")


##### demonstation of hat functions ####

fig4, ax4 = plt.subplots(1, figsize = (4,2))
ax4.plot(sig1x,sig2)
ax4.plot([0,1],sig1[0:2])
ax4.plot(sig1x,phi[0,:]*sig1[0])
ax4.plot(sig1x,phi[1,:]*sig1[1])
ax4.set_title("Hat basis functions to linear function")
ax4.set_xticks(np.arange(0,6))
ax4.set_xticklabels(labels)

fig4.savefig("hat_demo.svg")