# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

pi = 3.14
rho = 1;
eps = 1;

nx = 8
ny = 8
x = np.linspace(-1,1,nx)
y = np.linspace(-1,1,nx)

x,y = np.meshgrid(x,y)

u = -y/np.sqrt(x**2 + y**2)
v = x/np.sqrt(x**2 + y**2)

plt.quiver(x,y,u,v,scale=15)
plt.axis('off')
plt.savefig("cur1.svg")
plt.show()



z = np.zeros((nx,ny))
u = 0
v = 0
w = 1/np.sqrt(x**2+y**2)

fig = plt.figure()
ax = fig.gca(projection='3d')

plt.quiver(x,y,z,u,v,w,length=0.3,arrow_length_ratio=0.2)
ax.set_zlim((-0.01,1))
ax.set_ylim((-1,1))
ax.set_xlim((-1,1))
#plt.axis('off')
ax.set_yticklabels([])
ax.set_ylabel('y')
ax.set_xticklabels([])
ax.set_xlabel('x')
ax.set_zticklabels([])
ax.set_xlabel('z')
ax.view_init(elev=40., azim=70)
plt.savefig("cur1_res.svg")
plt.show()
