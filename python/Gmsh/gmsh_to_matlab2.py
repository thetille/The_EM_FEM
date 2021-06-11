# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 14:12:14 2021

@author: benja
"""
import gmsh
import sys
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as a3
import matplotlib.colors as colors
import scipy.io as io

# If sys.argv is passed to gmsh.initialize(), Gmsh will parse the command line
# in the same way as the standalone Gmsh app:

gmsh.initialize(sys.argv)

name_list = ["waveguide_flat_rod{}","waveguide_model2{}","waveguide_model3{}","waveguide_model3 - simple{}","waveguide_model3_flat{}","waveguide_model3_flat_long{}","waveguide_with_3_ports{}"]
name = name_list[0]
gmsh.open(name.format(".geo"))

def meshSizeCallback(dim, tag, x, y, z):
    return 0.03 #- 0.35*x

plot = False


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

temp =gmsh.model.mesh.getNodes(dim = -1, tag=-1)
no2xyz = np.zeros((temp[0].size,3))
no2xyz[temp[0]-1,:] =  np.reshape(temp[1],(-1,3))

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

################ find pec and port faces #########################
port_list = []
for i,group_name in enumerate(groups_name):
    if group_name.startswith('port'):
        port_list.append(i)

fac2no_port = np.empty(len(port_list), dtype=object)
for i,port in enumerate(port_list):
    fac2no_port_temp = (gmsh.model.mesh.getElements(dim = groups[port][0], tag = enteisesForGroupe[port][0]))
    #fac2no_port_temp =  np.array(fac2no_port_temp, dtype=object)
    fac2no_port_temp = np.reshape(fac2no_port_temp[2],(-1,3))
    fac2no_port[i] = fac2no_port_temp.T

groupidx = groups_name.index('pec')
fac2no_bound = np.empty((0,3),dtype = np.uint64)
for entei in enteisesForGroupe[groupidx][:]:
    nodes = gmsh.model.mesh.getElements(dim = groups[groupidx][0], tag = entei)[2]
    fac2no_bound = np.concatenate((fac2no_bound,np.reshape(nodes,(-1,3))))

if plot:

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')   
    for nodes in fac2no_bound:
        xyz = no2xyz[nodes-1]
        ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2])
    
################ find port and pec edges #########################
ed2no_port = np.empty(len(port_list), dtype=object)
for i,port in enumerate(port_list):
    ed2no_port_temp = np.concatenate((fac2no_port[i][[0,1],:].T,fac2no_port[i][[1,2],:].T,fac2no_port[i][[0,2],:].T))
    ed2no_port_temp = np.sort(ed2no_port_temp, axis = -1)
    ed2no_port_temp = np.unique(ed2no_port_temp,axis = 0)
    ed2no_port[i] = ed2no_port_temp.T

ed2no_pec = np.concatenate((fac2no_bound[:,[0,1]],fac2no_bound[:,[1,2]],fac2no_bound[:,[0,2]]))
ed2no_pec = np.sort(ed2no_pec, axis = -1)
ed2no_pec = np.unique(ed2no_pec,axis = 0)

################ find volume materials #########################

el2ma = np.zeros(len(el2no),dtype = np.int64)*-1;
for i,volume in enumerate((3,4)):
    elInMa = (gmsh.model.mesh.getElements(dim = groups[volume][0], tag = enteisesForGroupe[volume][0]))
    #fac2no_port_temp =  np.array(fac2no_port_temp, dtype=object)
    for k,el in enumerate(elInMa[1][0]):
        el2ma[np.where(ellements[1][0] == el)] = i+1
        
np.where(el2ma == 0)


if plot:
    
    # ax = a3.Axes3D(plt.figure())
    # ax.set_xlim(np.min(no2xyz[:,0]),np.max(no2xyz[:,0]))
    # ax.set_ylim(np.min(no2xyz[:,1]),np.max(no2xyz[:,1]))
    # ax.set_zlim(np.min(no2xyz[:,2]),np.max(no2xyz[:,2]))
    
    # for tri in fac2no_port[0].T:
    #     nodes = no2xyz[tri-1]
    #     ax.plot(np.append(nodes[:,0],nodes[1,0]),np.append(nodes[:,1],nodes[0,1]),np.append(nodes[:,2],nodes[0,2]))
    
    
    
    
    ax = a3.Axes3D(plt.figure())
    ax.set_xlim(np.min(no2xyz[:,0]),np.max(no2xyz[:,0]))
    ax.set_ylim(np.min(no2xyz[:,1]),np.max(no2xyz[:,1]))
    ax.set_zlim(np.min(no2xyz[:,2]),np.max(no2xyz[:,2]))
    
    for j,port in enumerate(port_list):
        for i in range(ed2no_port[j].shape[1]):
            nodes = ed2no_port[j][:,i]
            vtx = no2xyz[nodes-1];
            tri = a3.art3d.Poly3DCollection([vtx])
            #tri.set_color(colors.rgb2hex(np.random.rand(3)))
            tri.set_edgecolor("red")
            ax.add_collection3d(tri)

    for i in range(len(ed2no_pec)):
        nodes = ed2no_pec[i,:].T
        vtx = no2xyz[nodes-1];
        tri = a3.art3d.Poly3DCollection([vtx])
        #tri.set_color(colors.rgb2hex(np.random.rand(3)))
        tri.set_edgecolor("green")
        ax.add_collection3d(tri)
    plt.show()
    
    fig = plt.figure()
    colors = ((1,0,0),(0,1,0))
    ax = fig.add_subplot(111, projection='3d')   
    for i,nodes in enumerate(el2no):
        xyz = no2xyz[nodes-1]
        ax.scatter(np.mean(xyz[:,0]),np.mean(xyz[:,1]),np.mean(xyz[:,2]),color = colors[el2ma[i]-1])


el2no = el2no.astype('float64') 
ed2no_pec = ed2no_pec.astype('float64') 
fa2no_all = fa2no_all.astype('float64')

mdic = {"ed2no_pec": ed2no_pec.T, "ed2no_port":ed2no_port, "no2xyz": no2xyz.T, "el2no": el2no.T, "el2ma": el2ma, "ed2no_all": ed2no_all.T, "fa2no_all": fa2no_all.T,"port_fac2no_list": fac2no_port}

io.savemat(name.format("10.mat"), mdic)