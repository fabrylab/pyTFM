

'''
script to show that stress free boundary conditions are automatically fullfilled if i ad a layer of zero traction force
elemetns around the cells

'''


from __future__ import absolute_import, division, print_function
#import matplotlib
#matplotlib.use('Agg')

import solidspy.postprocesor as pos
import solidspy.assemutil as ass
import solidspy.solutil as sol

from grid_setup_solids_py import *
from solids_py_stress_functions import *

pixelsize = 6.25 / 40 # µm/pixel pixelsize of the original images

folder="/media/user/GINA1-BK/data_traktion_force_microscopy/fem_example2/"
# loading mask
db = clickpoints.DataFile(
    "/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/mask_cell_boundary3.cdb",
    "r")
mask = db.getMask(frame=0).data
# loading traction forces
#t_x = np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/11tx.npy")
#t_y = np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/11ty.npy")
t_x = np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/11tx.npy")
t_y = np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/11ty.npy")
# interpolate traction force array:


#t_x_resize=cv.resize(t_x,dsize=(int(mask.shape[1]*0.7),int(mask.shape[0]*0.7)),interpolation=cv.INTER_LINEAR)
#t_y_resize=cv.resize(t_y,dsize=(int(mask.shape[1]*0.7),int(mask.shape[0]*0.7)),interpolation=cv.INTER_LINEAR) ## looks good
t_x_resize=t_x
t_y_resize=t_y

# some pre clean up
mask = remove_small_holes(mask, 100)
mask = remove_small_objects(label(mask), 1000) > 0  # removing other small bits
# interpolation to size of traction force array
mask_int = interpolation(mask, t_x_resize.shape)
# further preparatio of mask data
mask_area, mask_boundaries = prepare_mask(mask_int)
mask_area=binary_dil(mask_area,iterations=40)

mask_edge=np.logical_and(mask_area,~binary_erosion(mask_area))
mask_area_outer=binary_dil(mask_area,iterations=2)

mask_edge_outer=np.logical_and(mask_area_outer,~binary_erosion(mask_area_outer))
mask_edge_outer1=np.logical_and(~mask_area,binary_erosion(mask_area_outer))

ps_new=pixelsize*np.mean(np.array(mask.shape)/np.array(mask_area.shape)) # pixelsize of fem grid in µm (?? is this ok??)






#1 method extended grid
f_x=t_x_resize*((ps_new*(10**-6))**2) # point force for each node from tractions
f_y=t_y_resize*((ps_new*(10**-6))**2) ## this factor is just for numerical reasons.... ## alsow try to mae this more efficient
f_x[~mask_area]=np.nan # setting all values outside of maske area to zero
f_y[~mask_area]=np.nan


#f_x=f_x.astype("float128")
#f_y=f_x.astype("float128")
f_x_c1=f_x-np.nanmean(f_x) #normalizing traction force to sum up to zero (no displacement)
f_y_c1=f_y-np.nanmean(f_y)
f_x_c2,f_y_c2,p=correct_torque(f_x_c1,f_y_c1,mask_area)


f_x_c2[np.logical_and(~mask_area,mask_area_outer)]=0 # trying zero force at outer boundarie
f_y_c2[np.logical_and(~mask_area,mask_area_outer)]=0





nodes, elements, loads, mats = grid_setup(mask_area_outer, f_x_c2, f_y_c2, 1, 0.49)



#get_torque2(nodes,loads)

# note: assume incomressibility : possion ratio=0.5, should be about the same as any other value according to tambe et al 2013
#plot_grid(nodes, elements, inverted_axis=True)
#plot_grid(nodes,elements,inverted_axis=False,symbol_size=4,arrows=True)
DME, IBC, neq = ass.DME(nodes, elements)  # boundary conditions asembly??

print("Number of elements: {}".format(elements.shape[0]))
print("Number of equations: {}".format(neq))


# System assembly
KG = ass.assembler(elements, mats, nodes, neq, DME,sparse=True)
RHSG = ass.loadasem(loads, IBC, neq)

# System solution


UG_sol,rx = sol.static_sol_cond(KG, RHSG,mask_area_outer,verbose=False) #solver with constratinst to zero translation and zero rotation

#UG_sol = sol.static_sol(KG, RHSG,nodes) #solver with constratinst to zero translation and zero rotation
UG=UG_sol[0]


#UG1,rx = sol.static_sol_cond(KG1, RHSG,mask_area) #solver with constratinst to zero translation and zero rotation
#norm1=np.sqrt(np.sum((RHSG-np.dot(KG.toarray(),UG[0]))**2)) # same norm as returned by sparse solver


#norm2=np.sqrt(np.sum((RHSG-np.dot(KG1,UG1[0]))**2))# same norm as returned by sparse solver
if not (np.allclose(KG.dot(UG) / KG.max(), RHSG / KG.max())):
    print("The system is not in equilibrium!")

UC = pos.complete_disp(IBC, nodes, UG)  # uc are x and y displacements
E_node1, S_nodes = pos.strain_nodes(nodes, elements, mats, UC) # stresses and strains
stress_tensor=calculate_stress_tensor(S_nodes,nodes,dims=mask_area.shape) # assembling the stress tensor
sigma_max,sigma_min,tau_max, phi_n,phi_shear,sigma_avg=all_stress_measures(S_nodes, nodes, dims=mask_area.shape)





# line stresses for inner edge in method 1
graph,points=mask_to_graph(mask_edge_outer) # representation of boundaries as graph
points=order_points(graph,points)
n_b,n_array_b=normal_vector(points,dims=mask_area.shape)
stress_vector_b=calculate_stress_vector(n_array_b,stress_tensor)
stress_vector_norm_b=np.linalg.norm(stress_vector_b,axis=2)


#check_normal_vectors(mask_inner_edge,points,n_b)


# line stresses for cells method1
graph,points=mask_to_graph(mask_boundaries) # representation of boundaries as graph
n_l1,n_array_l1=normal_vector_from_graph(graph,points,dims=mask_area.shape) # all nomral vectors of this graph
stress_vector_l1=calculate_stress_vector(n_array_l1,stress_tensor)
stress_vector_norm_l1=np.linalg.norm(stress_vector_l1,axis=2)



### check for natural boundary condition:
# try to show sigam*n=T*n*n, where sigam:stress tensor, n normal vector towards outside, T traction force vector

################# work this out please......
T=np.append(np.expand_dims(f_x_c2,-1),np.expand_dims(f_y_c2,-1),axis=2)

bc1=np.zeros(n_array_b.shape)

## outward (along nomral vector) directed traction forces
bc1[:,:,0]=(n_array_b[:,:,0]*T[:,:,0]+n_array_b[:,:,1]*T[:,:,1])*n_array_b[:,:,0]
bc1[:,:,1]=(n_array_b[:,:,0]*T[:,:,0]+n_array_b[:,:,1]*T[:,:,1])*n_array_b[:,:,1]


bc2=np.zeros(n_array_b.shape)
bc2[:,:,0]=stress_tensor[:,:,0,0]*n_array_b[:,:,0]+stress_tensor[:,:,0,1]*n_array_b[:,:,1]
bc2[:,:,1]=stress_tensor[:,:,1,0]*n_array_b[:,:,0]+stress_tensor[:,:,1,1]*n_array_b[:,:,1]




show_quiver_mask(T[:,:,0],T[:,:,1],mask_edge_outer)
show_quiver_mask(n_array_b[:,:,0],n_array_b[:,:,1])
show_quiver_mask(bc1[:,:,0],bc1[:,:,1],mask_edge_outer)
show_quiver_mask(bc2[:,:,0],bc2[:,:,1],mask_edge_outer)
show_quiver_mask(T[:,:,0],T[:,:,1],mask_area_outer)
fig=plot_map(sigma_min,origin="upper",title="minimal principal stress",cbar_str="force in N/pixel" ,mask=mask_area_outer)
#plot_arrows(nodes,UC[:,0],UC[:,1],filter=[0,6],dims=mask_area.shape,origin="upper")
fig=plot_fields(nodes,fields=[S_nodes[:,0],S_nodes[:,1],S_nodes[:,2]],dims=mask_area.shape,titles=["absolute value of\nx_stress","absolute value of\ny_stress","absolute value of\nxy_stress"],cbar_str="stress in N/pixel",origin="upper",mask=mask_area)

#### result: doesnt seem to really work like this_---> do more research
#-> two layers or even more of zero traction forces dont help
# -> expanding the grid doesnt help either: ask ben
# ask ben about this...
