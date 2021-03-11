import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

pi = 3.14
rho = 1;
eps = 1;

nx = 10
ny = 10
x = np.linspace(-1,1,nx)
y = np.linspace(-1,1,nx)

x,y = np.meshgrid(x,y)

u = x/np.sqrt(x**2 + y**2)
v = y/np.sqrt(x**2 + y**2)

plt.quiver(x,y,u,v,scale=15)
plt.axis('off')
plt.savefig("divergence.svg",bbox_inches = 0, transparent = True)
plt.show()


#divergance is'
nx = 50
ny = 50
x = np.linspace(-1,1,nx)
y = np.linspace(-1,1,nx)

x,y = np.meshgrid(x,y)

div = 1/np.sqrt(x**2 + y**2)

#norm = cm.colors.Normalize(vmax=abs(Z).max(), vmin=-abs(Z).max())
#div = 10*np.log10(div/10)
plt.pcolor(div,norm=colors.LogNorm(vmin=div.min()))#,extent=[-1, 1, -1, 1], norm=colors.LogNorm(vmin=div.min(), vmax=div.max()) )
plt.colorbar()
plt.axis('off')
plt.savefig("divergence_res.svg",bbox_inches = 0, transparent = True)
plt.show()
