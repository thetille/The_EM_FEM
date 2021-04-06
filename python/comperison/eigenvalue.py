# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 15:33:41 2021

@author: benja
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as cons
import scipy.io as io



###### analytical solution #######
m_list = np.arange(0,3)
n_list = np.arange(0,3)
l_list = np.arange(0,3)

a = 0.2
b = 0.2
d = 0.4

k = np.zeros((3,3,3))
for m in m_list:
    for n in n_list:
        for l in l_list:
            k[m,n,l] = np.sqrt( (((m*cons.pi)/a)**2) + (((n*cons.pi)/b)**2) + (((l*cons.pi)/d)**2))
            
            
f = np.zeros((4,4,4))
for m in m_list:
    for n in n_list:
        for l in l_list:
            f[m,n,l] = (cons.c*k[m,n,l])/(2*cons.pi)
            
            
fig, ax = plt.subplots()
f_list = np.reshape(f,(-1))
f_list = np.sort(f_list)/(1e9)



cstData = io.loadmat('cst_data.mat')['cst_data']

ax.plot(cstData, '+', label='CST')

f_list = f_list[42:-1]
#f_list = f_list[2:-1]
ax.plot(f_list, 'x',label='Analytical')
plt.legend()