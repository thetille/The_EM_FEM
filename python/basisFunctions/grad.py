import matplotlib.pyplot as plt
import numpy as np
#generate 3d points

def plot_lines(ax):
    ax.plot( (hhat[0],what[0]), (hhat[1],what[1]),'ro-', zs = (hhat[2],what[2]))
    ax.plot( (hhat[0],vhat[0]), (hhat[1],vhat[1]),'ro-', zs = (hhat[2],vhat[2]))
    ax.plot( (hhat[0],uhat[0]), (hhat[1],uhat[1]),'ro-', zs = (hhat[2],uhat[2]))
    
    ax.plot( (vhat[0],what[0]), (vhat[1],what[1]),'ro-', zs = (vhat[2],what[2]))
    ax.plot( (uhat[0],vhat[0]), (uhat[1],vhat[1]),'ro-', zs = (uhat[2],vhat[2]))
    ax.plot( (what[0],uhat[0]), (what[1],uhat[1]),'ro-', zs = (what[2],uhat[2]))

plt.close('all')

w_res = 0.1
v_res = 0.1
u_res = 0.1
w_mesh = []
v_mesh = []
u_mesh = []
phi1 = []
phi2 = []
phi3 = []
phi4 = []

N1 = np.array([])
N2 = np.array([])
N3 = []
N4 = []
N5 = []
N6 = []

uhat = np.array([1,0,0])
vhat = np.array([0,1,0])
what = np.array([0,0,1])
hhat = np.array([0,0,0])

w_list = np.arange(0,1,w_res)
for w in w_list:
    v_list = np.arange(0,1-w,v_res)
    for v in v_list:
        u_list = np.arange(0,1-v-w,u_res)
        for u in u_list:
            w_mesh.append(w)
            v_mesh.append(v)
            u_mesh.append(u)
            
            
N1 = np.zeros((len(w_mesh),3))
N2 = np.zeros((len(w_mesh),3))
N3 = np.zeros((len(w_mesh),3))
N4 = np.zeros((len(w_mesh),3))
N5 = np.zeros((len(w_mesh),3))
N6 = np.zeros((len(w_mesh),3))

for i,(w,v,u) in enumerate(zip(w_mesh,v_mesh,u_mesh)):         
    phi1.append(1-u-v-w)
    phi2.append(u)
    phi3.append(v)
    phi4.append(w)
    
    N1[i,:] = ( ((1-w-v)*uhat) + (u*vhat) + (u*what))
    N2[i,:] = (-v*uhat+u*vhat)
    N3[i,:] = (-v*uhat+(u+w-1)*vhat-v*what)
    N4[i,:] = (w*uhat+w*vhat+(1-v-u)*what)
    N5[i,:] = -w*uhat+ u*what
    N6[i,:] = -w*vhat + v*what
    
    #M1[i,:] = 2*(u*uhat+v*vhat+(w-1)*what)
    #M2[i,:] = 2*(u*uhat+(v-1)*vhat+w*what)
    #M3[i,:] = 2*(u*uhat+(v-1)*vhat+w*what)
    
    
#fig = plt.figure(figsize=plt.figaspect(0.25))
fig = plt.figure()
ax = fig.add_subplot(2, 2, 1, projection='3d')
ax.scatter(w_mesh,v_mesh,u_mesh,c = phi1)
plot_lines(ax)
ax = fig.add_subplot(2, 2, 2, projection='3d')
ax.scatter(w_mesh,v_mesh,u_mesh,c = phi2) 
plot_lines(ax)
ax = fig.add_subplot(2, 2, 3, projection='3d')
ax.scatter(w_mesh,v_mesh,u_mesh,c = phi3)
plot_lines(ax)
ax = fig.add_subplot(2, 2, 4, projection='3d')
ax.scatter(w_mesh,v_mesh,u_mesh,c = phi4)
plot_lines(ax)

fig = plt.figure()
ax = fig.add_subplot(2, 3, 1, projection='3d')
mag =( N1[:,0]**2+N1[:,1]**2+N1[:,2]**2 )
mag = mag/np.max(mag)
c = np.concatenate((mag, np.repeat(mag, 2)))
c = plt.cm.viridis(c)
ax.quiver(w_mesh,v_mesh,u_mesh,N1[:,0],N1[:,1],N1[:,2],mag,length=0.1,normalize=True, color = c)
plot_lines(ax)

ax = fig.add_subplot(2, 3, 2, projection='3d')
mag =( N2[:,0]**2+N2[:,1]**2+N2[:,2]**2 )
mag = mag/np.max(mag)
c = np.concatenate((mag, np.repeat(mag, 2)))
c = plt.cm.viridis(c)
ax.quiver(w_mesh,v_mesh,u_mesh,N2[:,0],N2[:,1],N2[:,2],mag,length=0.1,normalize=True, color = c)
plot_lines(ax)

ax = fig.add_subplot(2, 3, 3, projection='3d')
mag =( N3[:,0]**2+N3[:,1]**2+N3[:,2]**2 )
mag = mag/np.max(mag)
c = np.concatenate((mag, np.repeat(mag, 2)))
c = plt.cm.viridis(c)
ax.quiver(w_mesh,v_mesh,u_mesh,N3[:,0],N3[:,1],N3[:,2],mag,length=0.1,normalize=True,  color = c)
plot_lines(ax)

ax = fig.add_subplot(2, 3, 4, projection='3d')
mag =( N4[:,0]**2+N4[:,1]**2+N4[:,2]**2 )
mag = mag/np.max(mag)
c = np.concatenate((mag, np.repeat(mag, 2)))
c = plt.cm.viridis(c)
ax.quiver(w_mesh,v_mesh,u_mesh,N4[:,0],N4[:,1],N4[:,2],mag,length=0.1,normalize=True, color = c)
plot_lines(ax)

ax = fig.add_subplot(2, 3, 5, projection='3d')
mag =( N5[:,0]**2+N5[:,1]**2+N5[:,2]**2 )
mag = mag/np.max(mag)
c = np.concatenate((mag, np.repeat(mag, 2)))
c = plt.cm.viridis(c)
ax.quiver(w_mesh,v_mesh,u_mesh,N5[:,0],N5[:,1],N5[:,2],mag,length=0.1,normalize=True, color = c)
plot_lines(ax)

ax = fig.add_subplot(2, 3, 6, projection='3d')
mag =( N6[:,0]**2+N6[:,1]**2+N6[:,2]**2 )
mag = mag/np.max(mag)
c = np.concatenate((mag, np.repeat(mag, 2)))
c = plt.cm.viridis(c)
ax.quiver(w_mesh,v_mesh,u_mesh,N6[:,0],N6[:,1],N6[:,2],mag,length=0.1,normalize=True, color = c)
plot_lines(ax)




    




# fig = plt.figure()
# ax = fig.add_subplot(1, 2, 1, projection='3d')
# mag =( N1[:,0]**2+N1[:,0]**2+N1[:,0]**2 )
# ax.quiver(w_mesh,v_mesh,u_mesh,N1[:,0],N1[:,1],N1[:,2],mag,length=0.05,normalize=True, cmap=plt.cm.hsv)

# ax = fig.add_subplot(1, 2, 2, projection='3d')
# mag =( N1[:,0]**2+N1[:,0]**2+N1[:,0]**2 )
# mag = mag/np.max(mag)
# # Repeat for each body line and two head lines
# c = np.concatenate((mag, np.repeat(mag, 2)))
# c = plt.cm.viridis(c)
# ax.quiver(w_mesh,v_mesh,u_mesh,N1[:,0],N1[:,1],N1[:,2], colors=c ,length=0.05,normalize=True)