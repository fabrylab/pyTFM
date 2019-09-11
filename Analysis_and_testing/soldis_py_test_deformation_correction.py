from __future__ import absolute_import, division, print_function
from solidspy import solids_GUI
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import copy
import os
from datetime import datetime as dt
import solidspy.preprocesor as pre
import solidspy.postprocesor as pos
import solidspy.assemutil as ass
import solidspy.solutil as sol
import matplotlib
from functions_for_cell_colonie import *

folder="/media/user/GINA1-BK/data_traktion_force_microscopy/fem_example2/"
plt.close("all")


t1=time.time()
## generating input files
# regular grid of size n
n=100
l=np.sqrt(n).astype(int)
c_x,c_y=np.meshgrid(range(int(np.sqrt(n))),range(int(np.sqrt(n))))


nodes=np.zeros((n,5))
nodes[:,0]=np.arange(n)
nodes[:,1]=c_x.flatten()
nodes[:,2]=c_y.flatten()
nodes=nodes.astype(int)

n_el=int((np.sqrt(n)-1)**2) # number of elements
elements=np.zeros((n_el,7))



elements[:,0]=np.arange(n_el) #id
elements[:,1]=1 #element type/geometry (squares)
elements[:,2]=0 # elastic properties reference

sqr=[(c_x[:-1,:-1],c_y[:-1,:-1]), (c_x[:-1,:-1]+1,c_y[:-1,:-1]), (c_x[:-1,:-1]+1,c_y[:-1,:-1]+1),(c_x[:-1,:-1],c_y[:-1,:-1]+1)] #needs "counter clockwise orientation
sqr=[(x.flatten(),y.flatten()) for x,y in sqr]



ids=np.zeros(c_x.shape)   # array where all coordinates are attributed their id
ids[c_x.flatten(),c_y.flatten()]=np.arange(n,dtype=int)


elements[:,3]=ids[sqr[0][0],sqr[0][1]]  ## zusammenfassen
elements[:,4]=ids[sqr[1][0],sqr[1][1]]
elements[:,5]=ids[sqr[2][0],sqr[2][1]]
elements[:,6]=ids[sqr[3][0],sqr[3][1]]
elements=elements.astype(int)

# writing to files
loads=np.zeros((n,3))
loads[:,0]=np.arange(n)

loads[:50, 2] = -0.001
loads[50:100, 2] = 0.001

np.savetxt(os.path.join(folder,"nodes.txt"),nodes,fmt="%i",delimiter="\t")
np.savetxt(os.path.join(folder,"loads.txt"),loads,fmt="%.3f",delimiter="\t")
np.savetxt(os.path.join(folder,"eles.txt"),elements,fmt="%i",delimiter="\t")
t2=time.time()
print("grid assembly ", t2-t1)

start_time = dt.now()
echo = False




nodes, mats, elements, loads = pre.readin(folder=folder)
if echo:
    pre.echomod(nodes, mats, elements, loads, folder=folder)
DME, IBC, neq = ass.DME(nodes, elements)  # boundary conditions asembly??
# DME
#IBC list of "equations in x y direction for each node, -1 indicates no equations of movement, as completly fixed --> related to boundary conditions
# neq total number of equations
# DME is like list of independent equations per element(square, so just rearanged IBC for elements,side node zeros are at the end because of fixed length( for other element shape=
#
#
# print("Number of nodes: {}".format(nodes.shape[0]))
print("Number of elements: {}".format(elements.shape[0]))
print("Number of equations: {}".format(neq))



# System assembly
t1=time.time()
KG = ass.assembler(elements, mats, nodes, neq, DME,sparse=False)
RHSG = ass.loadasem(loads, IBC, neq)


t2=time.time()
print("sys assembly ",t2-t1)

# System solution

t1=time.time()
UG1 = sol.static_sol(KG, RHSG)  # solver without constrainsts
UG2,rx = sol.static_sol_cond(KG, RHSG, np.ones((l,l))) #solver with constratinst to zero translation and zero rotation

if not (np.allclose(KG.dot(UG1) / KG.max(), RHSG / KG.max())):
    print("The system is not in equilibrium!")

t2=time.time()
print("system solution ", t2-t1)

end_time = dt.now()
#print('Duration for system solution: {}'.format(end_time - start_time))

# Post-processing
start_time = dt.now()
UC1 = pos.complete_disp(IBC, nodes, UG1)  # uc are x and y displacements
UC2 = pos.complete_disp(IBC, nodes, UG2)  # uc are x and y displacements

dx1,dy1=make_field(nodes,UC1,(l,l))
dx2,dy2=make_field(nodes,UC2,(l,l))
dx3,dy3,trans,angle=correct_rotation(dx1,dy1,np.ones((l,l)))  #correcting rotation and translation after solveing equation system

# note rotation correction doesnt seem to be that exact.....

UC3 = pos.complete_disp(IBC, nodes, make_solids_py_values_list(nodes,dx3,dy3,mask=np.ones((l,l))))

print(calculate_rotation(dx1,dy1, np.ones((l,l))))
print(calculate_rotation(dx2,dy2, np.ones((l,l))))
print(calculate_rotation(dx3,dy3, np.ones((l,l))))

# caclualting strains and so on
E_nodes1, S_nodes1 = pos.strain_nodes(nodes, elements, mats, UC1)
E_nodes2, S_nodes2 = pos.strain_nodes(nodes, elements, mats, UC2)
E_nodes3, S_nodes3 = pos.strain_nodes(nodes, elements, mats, UC3)


# plot loads
fig=plot_arrows(nodes,loads[:,1],loads[:,2],scale_ratio=0.06,dims=(l,l),title="loads",origin="upper")
# plotting deformation field
fig=plot_arrows(nodes,UC1[:,0],UC1[:,1],scale_ratio=0.06,title="deformation no correction",origin="upper")
fig=plot_arrows(nodes,UC2[:,0],UC2[:,1],scale_ratio=0.06,title="deformation constraint equations",origin="upper")
fig=plot_arrows(nodes,UC3[:,0],UC3[:,1],scale_ratio=0.06,title="deformation correction after solving",origin="upper")


# plotting strains
plot_fields(nodes,fields=[E_nodes1[:,0],E_nodes1[:,1],E_nodes1[:,2]],titles=["x_strain_1","y_strain_1","xy_strain_1"],cbar_str="strain")
plot_fields(nodes,fields=[E_nodes2[:,0],E_nodes2[:,1],E_nodes2[:,2]],titles=["x_strain_2","y_strain_2","xy_strain_2"],cbar_str="strain")
plot_fields(nodes,fields=[E_nodes3[:,0],E_nodes3[:,1],E_nodes3[:,2]],titles=["x_strain_3","y_strain_3","xy_strain_3"],cbar_str="strain")









