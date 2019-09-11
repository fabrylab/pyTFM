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
from grid_setup_solids_py import *
from tqdm import tqdm

plt.close("all")
#plt.ioff()
pixelsize = 6.25 / 40 # µm/pixel pixelsize of the original images

folder="/media/user/GINA1-BK/data_traktion_force_microscopy/fem_example2/"
# loading mask
db = clickpoints.DataFile(
    "/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/mask_cell_boundary3.cdb",
    "r")
mask = db.getMask(frame=0).data
# loading traction forces
t_x = np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/tx.npy")
t_y = np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/ty.npy")

# some pre clean up
mask = remove_small_holes(mask, 100)
mask = remove_small_objects(label(mask), 1000) > 0  # removing other small bits
# interpolation to size of traction force array
mask_int = interpolation(mask, t_x.shape)
# further preparatio of mask data
mask_area, mask_boundaries = prepare_mask(mask_int)

ps_new=pixelsize*np.mean(np.array(mask.shape)/np.array(mask_area.shape)) # pixelsize of fem grid in µm (?? is this ok??)


f_x=t_x*((ps_new*(10**-6))**2) # point force for each node from tractions
f_y=t_y*((ps_new*(10**-6))**2) ## this factor is just for numerical reasons.... ## alsow try to mae this more efficient
f_x[~mask_area]=np.nan # setting all values outside of maske area to zero
f_y[~mask_area]=np.nan
#f_x=f_x.astype("float128")
#f_y=f_x.astype("float128")

f_x_c1=f_x-np.nanmean(f_x) #normalizing traction force to sum up to zero (no displacement)
f_y_c1=f_y-np.nanmean(f_y)

f_x_c2,f_y_c2,p=correct_torque(f_x_c1,f_y_c1,mask_area)







# setup of the grid
nodes, elements, loads, mats = grid_setup(mask_area, f_x_c2, f_y_c2, 100, 0.5)

nodes, elements, loads2, mats = grid_setup(mask_area, f_x_c1, f_y_c1, 100, 0.5)  # without torque correction

#nodes, elements, loads3, mats = grid_setup(mask_area, f_x, f_y, 1, 0.5)  #




#nodes[594]=np.array([594,  45,  44,   -1,   -1])
nodes[433]=np.array([433,  31,  40,  0,  0])

k=np.random.choice(nodes[:,0])
#nodes[k,[3,4]]= np.array([-1,   -1])
get_torque2(nodes,loads)

# note: assume incomressibility : possion ratio=0.5, should be about the same as any other value according to tambe et al 2013
#plot_grid(nodes, elements, inverted_axis=True)
#plot_grid(nodes,elements,inverted_axis=False,symbol_size=4,arrows=True)
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
RHSG2 = ass.loadasem(loads, IBC, neq) # not torque corrected loads
t2=time.time()
print("sys assembly ",t2-t1)










# System solution

t1=time.time()
UG1 = sol.static_sol(KG, RHSG,nodes)  # solver without constrainsts
UG2,rx = sol.static_sol_cond(KG, RHSG, nodes,mask_area) #solver with constratinst to zero translation and zero rotation

UG3 = sol.static_sol(KG, RHSG2,nodes)  # solver without constrainsts not torque corrected
UG4,rx = sol.static_sol_cond(KG, RHSG2, nodes,mask_area) #solver with constratinst to zero translation and zero rotation, torque corrected



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

dx1,dy1=make_field(nodes,UC1,mask_area.shape)
dx2,dy2=make_field(nodes,UC2,mask_area.shape)
dx3,dy3,trans,angle=correct_rotation(dx1,dy1,mask_area)  #correcting rotation and translation after solveing equation system

# note rotation correction doesnt seem to be that exact.....
UC3 = pos.complete_disp(IBC, nodes, make_solids_py_values_list(nodes,dx3,dy3,mask=mask_area))

UC4 = pos.complete_disp(IBC, nodes, UG3) # not torque corrected, normal solver
UC5 = pos.complete_disp(IBC, nodes, UG4) # not torque corrected with constraints
dx4,dy4=make_field(nodes,UC3,mask_area.shape)
dx5,dy5=make_field(nodes,UC4,mask_area.shape)




print(get_torque1(f_x,f_y,mask_area)) # show torque correction
print(get_torque1(f_x_c2,f_y_c2,mask_area)) # answer the question: is this alot of torque or rather litttle??ß

print(calculate_rotation(dx1,dy1, np.ones(mask_area.shape)))
print(calculate_rotation(dx2,dy2, np.ones(mask_area.shape)))
print(calculate_rotation(dx3,dy3, np.ones(mask_area.shape)))
print(calculate_rotation(dx4,dy4, np.ones(mask_area.shape)))
print(calculate_rotation(dx5,dy5, np.ones(mask_area.shape)))



#E_nodes1, S_nodes1 = pos.strain_nodes(nodes, elements, mats, UC1) #torque correction normal solver



# caclualting strains and so on
E_nodes1, S_nodes1 = pos.strain_nodes(nodes, elements, mats, UC1) #torque correction normal solver
E_nodes2, S_nodes2 = pos.strain_nodes(nodes, elements, mats, UC2) #torque correction solver wtih constraints
E_nodes3, S_nodes3 = pos.strain_nodes(nodes, elements, mats, UC3) # rotation correction after normal solver , with torque correction
E_nodes4, S_nodes4 = pos.strain_nodes(nodes, elements, mats, UC4) # no torque correction normal solver
E_nodes5, S_nodes5 = pos.strain_nodes(nodes, elements, mats, UC5) # no torque correction solver with constraints


# plotting loads
fig=plot_arrows(nodes,f_x[mask_area],f_y[mask_area],scale_ratio=0.06,dims=mask_area.shape,title="loads_uncorrected",origin="upper")
fig=plot_arrows(nodes,f_x_c2[mask_area],f_y_c2[mask_area],scale_ratio=0.06,dims=mask_area.shape,title="loads_corrected",origin="upper")

# plotting deformation field
fig=plot_arrows(nodes,UC1[:,0],UC1[:,1],scale_ratio=0.06,dims=mask_area.shape,title="def torque correction normal solver",origin="upper")
fig=plot_arrows(nodes,UC2[:,0],UC2[:,1],scale_ratio=0.06,dims=mask_area.shape,title="def torque correction solver wtih constraints",origin="upper")
fig=plot_arrows(nodes,UC3[:,0],UC3[:,1],scale_ratio=0.06,dims=mask_area.shape,title="def rotation correction after normal solver , with torque correction",origin="upper")
fig=plot_arrows(nodes,UC4[:,0],UC4[:,1],scale_ratio=0.06,dims=mask_area.shape,title="def no torque correction normal solver",origin="upper")
fig=plot_arrows(nodes,UC5[:,0],UC5[:,1],scale_ratio=0.06,dims=mask_area.shape,title="def no torque correction solver with constraints",origin="upper")
# plotting strains
plot_fields(nodes,fields=[E_nodes1[:,0],E_nodes1[:,1],E_nodes1[:,2]],dims=mask_area.shape,titles=["x_strain_1","y_strain_1","xy_strain_1"],cbar_str="strain")
plot_fields(nodes,fields=[E_nodes2[:,0],E_nodes2[:,1],E_nodes2[:,2]],dims=mask_area.shape,titles=["x_strain_2","y_strain_2","xy_strain_2"],cbar_str="strain")
plot_fields(nodes,fields=[E_nodes3[:,0],E_nodes3[:,1],E_nodes3[:,2]],dims=mask_area.shape,titles=["x_strain_3","y_strain_3","xy_strain_3"],cbar_str="strain")
plot_fields(nodes,fields=[E_nodes4[:,0],E_nodes4[:,1],E_nodes4[:,2]],dims=mask_area.shape,titles=["x_strain_4","y_strain_4","xy_strain_4"],cbar_str="strain")
plot_fields(nodes,fields=[E_nodes5[:,0],E_nodes5[:,1],E_nodes5[:,2]],dims=mask_area.shape,titles=["x_strain_5","y_strain_5","xy_strain_5"],cbar_str="strain")










    ## todo:
    # get scaling correct and also youngsmodulus in sigma


