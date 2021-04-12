# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 15:33:41 2021

@author: benja
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as cons
import scipy.io as io


wave_nummber = 15
###### analytical solution #######
m_list = np.arange(0,wave_nummber)
n_list = np.arange(0,wave_nummber)
l_list = np.arange(0,wave_nummber)

a = 0.2
b = 0.2
d = 0.4

k = np.zeros((wave_nummber,wave_nummber,wave_nummber))
for m in m_list:
    for n in n_list:
        for l in l_list:
            k[m,n,l] = np.sqrt( (((m*cons.pi)/a)**2) + (((n*cons.pi)/b)**2) + (((l*cons.pi)/d)**2))
            
            
f = np.zeros((wave_nummber,wave_nummber,wave_nummber))
for m in m_list:
    for n in n_list:
        for l in l_list:
            f[m,n,l] = (cons.c*k[m,n,l])/(2*cons.pi)
            
            
fig, ax = plt.subplots()
f_list = np.reshape(f,(-1))
f_list = np.sort(f_list)/(1e9)



#cstData = io.loadmat('cst_data.mat')['cst_data']

onlyEFeildData = np.reshape(io.loadmat('onlyEField.mat')['fr'][151:-1]/1e9,-1)
EAndHFeildData = np.reshape(io.loadmat('EAndHField.mat')['fr']/1e9,-1)

cstData = np.genfromtxt('CST_data.txt', delimiter='\n')

ax.plot(cstData, '+', label='CST')

f_list = f_list[5:-1]
f_list = np.insert(f_list, 8 , f_list[7])
#f_list = f_list[2:-1]
ax.plot(f_list, 'x',label='Analytical')

ax.plot(onlyEFeildData, 'x',label='only E Field data')
ax.plot(EAndHFeildData, 'x',label='E and H Field data')

plt.legend()

print("Analyticaly \t CST Data \t E Field Data \t cst,analyt diff \t E field,analyt diff")
for data in zip(f_list, cstData, onlyEFeildData, EAndHFeildData):
    print("{:.5f} \t\t {:.5f} \t\t {:.5f} \t\t {:.5f} \t\t {:.5f}".format(data[1],data[2],data[3],data[1]-data[2],data[1]-data[3]))