# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 14:12:14 2021

@author: benja
"""
import gmsh
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as a3
import matplotlib.colors as colors
import scipy.io as io

# If sys.argv is passed to gmsh.initialize(), Gmsh will parse the command line
# in the same way as the standalone Gmsh app:

gmsh.initialize(sys.argv)

name_list = ["cylinder_waveguide{}","waveguide_model2{}"]
name = name_list[1]
gmsh.open(name.format(".geo"))

def meshSizeCallback(dim, tag, x, y, z):
    return 0.08 - 0.15*z

gmsh.model.mesh.setSizeCallback(meshSizeCallback)

gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)



gmsh.model.mesh.generate(3)

groups = gmsh.model.getPhysicalGroups(dim=-1)

groups_name = []
for group in groups:
    groups_name.append(gmsh.model.getPhysicalName(dim = group[0], tag=group[1]))

enteisesForGroupe = []
for group in groups:
    enteisesForGroupe.append(gmsh.model.getEntitiesForPhysicalGroup(group[0],group[1]))

no2xyz = np.reshape(gmsh.model.mesh.getNodes(dim = -1, tag=-1)[1],(-1,3))

if np.size(no2xyz,0) != np.size(np.unique(no2xyz,axis = 0),0):
    raise Exception("Sorry, repeats") 


################ find all faces #########################
ellements = gmsh.model.mesh.getElements(dim = 3, tag = -1)

el2no = np.reshape(ellements[2],(-1,4))

fa2no_all = np.zeros((np.size(el2no,0)*4,3), dtype=np.int64)
i = 0
for element2no in el2no:
    fa2no_all[i,:] = element2no[[1,2,3]]
    i += 1
    fa2no_all[i,:] = element2no[[0,2,3]]
    i += 1
    fa2no_all[i,:] = element2no[[0,1,3]]
    i += 1
    fa2no_all[i,:] = element2no[[0,1,2]]
    i += 1

fa2no_all = np.sort(fa2no_all, axis=-1)
fa2no_all = np.unique(fa2no_all,axis = 0)



# ax = a3.Axes3D(plt.figure())
# ax.set_xlim(np.min(no2xyz[:,0]),np.max(no2xyz[:,0]))
# ax.set_ylim(np.min(no2xyz[:,1]),np.max(no2xyz[:,1]))
# ax.set_zlim(np.min(no2xyz[:,2]),np.max(no2xyz[:,2]))


# for i in range(len(fa2no_all)):
#     nodes = fa2no_all[i,:]
#     vtx = no2xyz[nodes-1];
#     tri = a3.art3d.Poly3DCollection([vtx])
#     tri.set_color(colors.rgb2hex(np.random.rand(3)))
#     tri.set_edgecolor('k')
#     ax.add_collection3d(tri)
#     plt.show()
#     plt.pause(0.01)


################ find all edges #########################
ed2no_all = np.zeros((len(el2no)*6,2), dtype=np.int64)

i = 0
for element2no in el2no:
    ed2no_all[i,:] = element2no[[0,1]]
    i += 1
    ed2no_all[i,:] = element2no[[0,2]]
    i += 1
    ed2no_all[i,:] = element2no[[0,3]]
    i += 1
    ed2no_all[i,:] = element2no[[1,2]]
    i += 1
    ed2no_all[i,:] = element2no[[1,3]]
    i += 1
    ed2no_all[i,:] = element2no[[2,3]]
    i += 1
    
ed2no_all = np.sort(ed2no_all, axis=-1)
ed2no_all = np.unique(ed2no_all,axis = 0)
    
# ax = a3.Axes3D(plt.figure())
# ax.set_xlim(np.min(no2xyz[:,0]),np.max(no2xyz[:,0]))
# ax.set_ylim(np.min(no2xyz[:,1]),np.max(no2xyz[:,1]))
# ax.set_zlim(np.min(no2xyz[:,2]),np.max(no2xyz[:,2]))

# for i in range(len(ed2no_all)):
#     nodes = ed2no_all[i,:]
#     vtx = no2xyz[nodes-1];
#     tri = a3.art3d.Poly3DCollection([vtx])
#     #tri.set_color(colors.rgb2hex(np.random.rand(3)))
#     tri.set_edgecolor(colors.rgb2hex(np.random.rand(3)))
#     ax.add_collection3d(tri)
# plt.show()

################ find pec faces #########################
groupidx = groups_name.index('port1')
fac2no_port1 = gmsh.model.mesh.getElements(dim = groups[groupidx][0], tag = enteisesForGroupe[groupidx][0])
fac2no_port1 = np.reshape(fac2no_port1[2],(-1,3))

groupidx = groups_name.index('port2')
fac2no_port2 = gmsh.model.mesh.getElements(dim = groups[1][0], tag = enteisesForGroupe[1][0])
fac2no_port2 = np.reshape(fac2no_port2[2],(-1,3))

groupidx = groups_name.index('bound')
fac2no_bound = np.empty((0,3),dtype = np.uint64)
for entei in enteisesForGroupe[groupidx][:]:
    nodes = gmsh.model.mesh.getElements(dim = groups[groupidx][0], tag = entei)[2]
    fac2no_bound = np.concatenate((fac2no_bound,np.reshape(nodes,(-1,3))))

fig = plt.figure()
ax = fig.add_subplot(211, projection='3d')
for nodes in fac2no_port1:
    xyz = no2xyz[nodes-1]
    ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2])
    
#ax = fig.add_subplot(211, projection='3d')   
for nodes in fac2no_port2:
    xyz = no2xyz[nodes-1]
    ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2])
    
ax = fig.add_subplot(212, projection='3d')   
for nodes in fac2no_bound:
    xyz = no2xyz[nodes-1]
    ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2])
    
################ find pec edges #########################
ed2no_port1 = np.concatenate((fac2no_port1[:,[0,1]],fac2no_port1[:,[1,2]],fac2no_port1[:,[0,2]]))
ed2no_port1 = np.sort(ed2no_port1, axis = -1)
ed2no_port1 = np.unique(ed2no_port1,axis = 0)

ed2no_port2 = np.concatenate((fac2no_port2[:,[0,1]],fac2no_port2[:,[1,2]],fac2no_port2[:,[0,2]]))
ed2no_port2 = np.sort(ed2no_port2, axis = -1)
ed2no_port2 = np.unique(ed2no_port2,axis = 0)

ed2no_bound = np.concatenate((fac2no_bound[:,[0,1]],fac2no_bound[:,[1,2]],fac2no_bound[:,[0,2]]))
ed2no_bound = np.sort(ed2no_bound, axis = -1)
ed2no_bound = np.unique(ed2no_bound,axis = 0)



ax = a3.Axes3D(plt.figure())
ax.set_xlim(np.min(no2xyz[:,0]),np.max(no2xyz[:,0]))
ax.set_ylim(np.min(no2xyz[:,1]),np.max(no2xyz[:,1]))
ax.set_zlim(np.min(no2xyz[:,2]),np.max(no2xyz[:,2]))

for i in range(len(ed2no_port1)):
    nodes = ed2no_port1[i,:]
    vtx = no2xyz[nodes-1];
    tri = a3.art3d.Poly3DCollection([vtx])
    #tri.set_color(colors.rgb2hex(np.random.rand(3)))
    tri.set_edgecolor("red")
    ax.add_collection3d(tri)
plt.show()


for i in range(len(ed2no_port2)):
    nodes = ed2no_port2[i,:]
    vtx = no2xyz[nodes-1];
    tri = a3.art3d.Poly3DCollection([vtx])
    #tri.set_color(colors.rgb2hex(np.random.rand(3)))
    tri.set_edgecolor("red")
    ax.add_collection3d(tri)
plt.show()

for i in range(len(ed2no_bound)):
    nodes = ed2no_bound[i,:]
    vtx = no2xyz[nodes-1];
    tri = a3.art3d.Poly3DCollection([vtx])
    #tri.set_color(colors.rgb2hex(np.random.rand(3)))
    tri.set_edgecolor("green")
    ax.add_collection3d(tri)
plt.show()

ed2no_pec = np.concatenate((ed2no_port1,ed2no_port2,ed2no_bound))
ed2no_pec = np.sort(ed2no_pec, axis = -1)
ed2no_pec = np.unique(ed2no_pec,axis = 0)

el2ma = np.ones((1,len(el2no)))

el2no = el2no.astype('float64') 
ed2no_port1 = ed2no_port1.astype('float64') 
ed2no_port2 = ed2no_port2.astype('float64') 
ed2no_bound = ed2no_bound.astype('float64') 
fa2no_all = fa2no_all.astype('float64')
ed2no_pec = ed2no_pec.astype('float64')

mdic = {"no2xyz": no2xyz.T, "el2no": el2no.T, "el2ma": el2ma, "ed2no_all": ed2no_all.T, "ed2no_port1": ed2no_port1.T, "ed2no_port2": ed2no_port2.T, "ed2no_bound": ed2no_bound.T, "fa2no_all": fa2no_all.T, "ed2no_pec": ed2no_pec.T}
io.savemat(name.format(".mat"), mdic)

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(no2xyz[:,0],no2xyz[:,1],no2xyz[:,2])

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# nodes = ed2no_all[2]
# nodes = np.unique(nodes)
# xyz = no2xyz[nodes-1]
# ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2])










# dim = 2
# port12elle = np.empty((0,dim+1),dtype = 'int64')
# for entities in enteisesForGroupe[0]:
#     nodes = (np.reshape(gmsh.model.mesh.getElements(dim = dim, tag = entities)[2],(-1,dim+1)))
#     port12elle = np.concatenate((port12elle ,nodes),axis = 0)

# # ax3 = fig.add_subplot(3,1,3, projection='3d')
# corrd = np.empty((0,3))
# for ellement in port12elle :
#     for node in ellement:
#         corrd = np.concatenate((corrd,[no2xyz[node]]),axis = 0)
        
#     print(corrd)
    
# ax.scatter(corrd[:,0],corrd[:,1],corrd[:,2],marker = 's')








 
    

#hej1 = gmsh.model.mesh.getNodesByElementType(1,tag = 5)



# hej2 = gmsh.model.mesh.getNodesByElementType(2)

# hej3 = gmsh.model.mesh.getElementsByType(1)

# hej4 = gmsh.model.mesh.getNodesForPhysicalGroup(dim=2, tag=55)



#hej6 = gmsh.model.mesh.getEdges([1, 29])

#A = np.reshape(hej[1],(len(hej[0]),3))



# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(A[:,0],A[:,0],A[:,0])

# if '-nopopup' not in sys.argv:
#     gmsh.fltk.run()