'''
script to compare two methods of handling the problem of force outside of the cell colony

# goals: compare to "projection of forces to edge, by: getting the projected forces
# and additionally get line of the boundaray in expanded case and calculate force trasnsimtion on this line
'''


from __future__ import absolute_import, division, print_function
#import matplotlib
#matplotlib.use('Agg')

import solidspy.postprocesor as pos
import solidspy.assemutil as ass
import solidspy.solutil as sol

from grid_setup_solids_py import *
from solids_py_stress_functions import *


def project_forces_on_border(fx,fy,mask_inner_edge,mask_border_region):
    '''
    function to project force from a border region on points of "mask_inner_edge".
    Finds the closest points in the border region for each edge point. And adds
    the sum off all their forces to this edge point
    :param fx: forces in x direction
    :param fy: forces in y direction
    :param mask_inner_edge: mask of the inner edge, needs to directly border the border region
    :param mask_border_region: mask of the border region
    :return:
    '''

    points1=np.array(np.where(mask_inner_edge)).T # edge between outer and innner mask
    points2=np.array(np.where(mask_border_region)).T # outer region only



    distances = np.linalg.norm(points1[None, :] - points2[:, None],
                               axis=2)  # first substracting full all values in matrix1 with all values n matrix2,
    # the taking the norm
    assign=np.argmin(distances,axis=1)
    assign_dict=defaultdict(list) # writing to dictionary with key: index of point in edge region, values: indices of points in border regions
    # one could actually scip this step..
    #check_closet_neigbours(points1, points2, assign,mask1=mask_expanded,mask2=mask_area)
    for p2_ind,p1_ind in enumerate(assign):
        assign_dict[p1_ind].append(p2_ind)



    # projecting forces on the edge region
    for p1_ind,p2_inds in assign_dict.items():
        fy[points1[p1_ind][0],points1[p1_ind][1]]+=np.sum(fy[points2[p2_inds][:,0],points2[p2_inds][:,1]])
        fx[points1[p1_ind][0], points1[p1_ind][1]] += np.sum(fx[points2[p2_inds][:, 0], points2[p2_inds][:, 1]])


    fig=show_quiver(fx,fy,filter=[0,0],scale_ratio=0.04)

    ax=fig.get_axes()[0]
    ax.imshow(mask_border_region*2,alpha=0.5)
    ax.imshow(mask_inner_edge,alpha=0.5)
    for p1_ind, p2_inds in assign_dict.items():
        ax.text(points1[p1_ind][1], points1[p1_ind][0], str(p1_ind),color="red",fontsize=20)
        for p in p2_inds:

            ax.text(points2[p][1],points2[p][0],str(p1_ind),color="blue",fontsize=20)
    #### note inner edge points that have no assignemnet are not shown


    return fx,fy




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
t_x = np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/11tx_200.npy")
t_y = np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/11ty_200.npy")
# interpolate traction force array:


#t_x_resize=cv.resize(t_x,dsize=(int(mask.shape[1]*0.1),int(mask.shape[0]*0.1)),interpolation=cv.INTER_LINEAR)
#t_y_resize=cv.resize(t_y,dsize=(int(mask.shape[1]*0.1),int(mask.shape[0]*0.1)),interpolation=cv.INTER_LINEAR) ## looks good
t_x_resize=t_x
t_y_resize=t_y

# some pre clean up
mask = remove_small_holes(mask, 100)
mask = remove_small_objects(label(mask), 1000) > 0  # removing other small bits
# interpolation to size of traction force array
mask_int = interpolation(mask, t_x_resize.shape)
# further preparatio of mask data


mask_area, mask_boundaries,mask_expanded,borders2 = prepare_mask_spread(mask_int,dil_factor=1.1)
mask_inner_edge=np.logical_and(mask_area,~binary_erosion(mask_area))
mask_outer_edge=np.logical_and(mask_expanded,~binary_erosion(mask_expanded))
ps_new=pixelsize*np.mean(np.array(mask.shape)/np.array(mask_area.shape)) # pixelsize of fem grid in µm (?? is this ok??)

######################
#mask_expanded=np.ones((int(mask_area.shape[0]*0.1),int(mask_area.shape[1]*0.1))).astype(bool)

#t_x_resize=cv.resize(t_x,dsize=(int(mask_expanded.shape[1]),int(mask_expanded.shape[0])),interpolation=cv.INTER_LINEAR)
#t_y_resize=cv.resize(t_y,dsize=(int(mask_expanded.shape[1]),int(mask_expanded.shape[0])),interpolation=cv.INTER_LINEAR)
###############################



#1 method extended grid
f_x_1=t_x_resize*((ps_new*(10**-6))**2) # point force for each node from tractions
f_y_1=t_y_resize*((ps_new*(10**-6))**2) ## this factor is just for numerical reasons.... ## alsow try to mae this more efficient
f_x_1[~mask_expanded]=np.nan # setting all values outside of maske area to zero
f_y_1[~mask_expanded]=np.nan
#f_x=f_x.astype("float128")
#f_y=f_x.astype("float128")
f_x_c1_1=f_x_1-np.nanmean(f_x_1) #normalizing traction force to sum up to zero (no displacement)
f_y_c1_1=f_y_1-np.nanmean(f_y_1)
f_x_c2_1,f_y_c2_1,p=correct_torque(f_x_c1_1,f_y_c1_1,mask_expanded)

#2 method force projection
f_x_2p=t_x_resize*((ps_new*(10**-6))**2) # point force for each node from tractions
f_y_2p=t_y_resize*((ps_new*(10**-6))**2) ## this factor is just for numerical reasons.... ## alsow try to mae this more efficient




f_x_2,f_y_2=project_forces_on_border(f_x_2p,f_y_2p,mask_inner_edge,np.logical_and(mask_expanded,~mask_area))
f_x_2[~mask_area]=np.nan # setting all values outside of maske area to zero
f_y_2[~mask_area]=np.nan

#f_x=f_x.astype("float128")
#f_y=f_x.astype("float128")
f_x_c1_2=f_x_2-np.nanmean(f_x_2) #normalizing traction force to sum up to zero (no displacement)
f_y_c1_2=f_y_2-np.nanmean(f_y_2)
f_x_c2_2,f_y_c2_2,p=correct_torque(f_x_c1_2,f_y_c1_2,mask_expanded)








nodes1, elements1, loads1, mats1 = grid_setup(mask_expanded, f_x_c2_1, f_y_c2_1, 1, 0.49)
nodes2, elements2, loads2, mats2 = grid_setup(mask_area, f_x_c2_2, f_y_c2_2, 1, 0.49)




#get_torque2(nodes,loads)

# note: assume incomressibility : possion ratio=0.5, should be about the same as any other value according to tambe et al 2013
#plot_grid(nodes, elements, inverted_axis=True)
#plot_grid(nodes,elements,inverted_axis=False,symbol_size=4,arrows=True)
DME1, IBC1, neq1 = ass.DME(nodes1, elements1)  # boundary conditions asembly??
DME2, IBC2, neq2 = ass.DME(nodes2, elements2)  # boundary conditions asembly??
print("Number of elements: {}".format(elements1.shape[0]))
print("Number of equations: {}".format(neq1))
print("Number of elements: {}".format(elements2.shape[0]))
print("Number of equations: {}".format(neq2))


# System assembly
KG1 = ass.assembler(elements1, mats1, nodes1, neq1, DME1,sparse=True)
RHSG1 = ass.loadasem(loads1, IBC1, neq1)

KG2 = ass.assembler(elements2, mats2, nodes2, neq2, DME2,sparse=True)
RHSG2 = ass.loadasem(loads2, IBC2, neq2)

# System solution


UG_sol1,rx1 = sol.static_sol_cond(KG1, RHSG1,mask_expanded,verbose=False) #solver with constratinst to zero translation and zero rotation
UG1=UG_sol1[0]

UG_sol2,rx2 = sol.static_sol_cond(KG2, RHSG2,mask_area,verbose=False) #solver with constratinst to zero translation and zero rotation
UG2=UG_sol2[0]
#UG1,rx = sol.static_sol_cond(KG1, RHSG,mask_area) #solver with constratinst to zero translation and zero rotation
#norm1=np.sqrt(np.sum((RHSG-np.dot(KG.toarray(),UG[0]))**2)) # same norm as returned by sparse solver


#norm2=np.sqrt(np.sum((RHSG-np.dot(KG1,UG1[0]))**2))# same norm as returned by sparse solver
if not (np.allclose(KG1.dot(UG1) / KG1.max(), RHSG1 / KG1.max())):
    print("The system is not in equilibrium!")
if not (np.allclose(KG2.dot(UG2) / KG2.max(), RHSG2 / KG2.max())):
    print("The system is not in equilibrium!")
UC1 = pos.complete_disp(IBC1, nodes1, UG1)  # uc are x and y displacements
E_nodes1, S_nodes1 = pos.strain_nodes(nodes1, elements1, mats1, UC1) # stresses and strains
stress_tensor1=calculate_stress_tensor(S_nodes1,nodes1,dims=mask_area.shape) # assembling the stress tensor

UC2 = pos.complete_disp(IBC2, nodes2, UG2)  # uc are x and y displacements
E_nodes2, S_nodes2 = pos.strain_nodes(nodes2, elements2, mats2, UC2) # stresses and strains
stress_tensor2=calculate_stress_tensor(S_nodes2,nodes2,dims=mask_area.shape) # assembling the stress tens






# line stresses for inner edge in method 1
graph,points=mask_to_graph(mask_inner_edge) # representation of boundaries as graph
points=order_points(graph,points)
n_b,n_array_b=normal_vector(points,dims=mask_area.shape)
stress_vector_b=calculate_stress_vector(n_array_b,stress_tensor1)
stress_vector_norm_b=np.linalg.norm(stress_vector_b,axis=2)


#check_normal_vectors(mask_inner_edge,points,n_b)


# line stresses for cells method1
graph,points=mask_to_graph(mask_boundaries) # representation of boundaries as graph
n_l1,n_array_l1=normal_vector_from_graph(graph,points,dims=mask_area.shape) # all nomral vectors of this graph
stress_vector_l1=calculate_stress_vector(n_array_l1,stress_tensor1)
stress_vector_norm_l1=np.linalg.norm(stress_vector_l1,axis=2)


# line stresses for cells method2
graph,points=mask_to_graph(mask_boundaries) # representation of boundaries as graph
n_l2,n_array_l2=normal_vector_from_graph(graph,points,dims=mask_area.shape) # all nomral vectors of this graph
stress_vector_l2=calculate_stress_vector(n_array_l2,stress_tensor2)
stress_vector_norm_l2=np.linalg.norm(stress_vector_l2,axis=2)

#plt.figure()
#plt.imshow(mask_inner_edge)
#plt.quiver(n_array_b[:,:,0],n_array_b[:,:,1],scale=0.2, scale_units="xy", angles="xy")





##
lines_spline_points=borders.lines_spline_points
lines_splines=borders.lines_splines
lines_points=borders.lines_points
# plot linestresses over border as continous curves:
lines_interpol,min_v,max_v=interpolation_for_stress_and_normal_vector(lines_splines,lines_points,stress_tensor,pixel_length=ps_new,interpol_factor=6)
fig=plot_continous_boundary_stresses(mask_int,lines_interpol,min_v,max_v,plot_n_arrows=True)





### other possible stress measures, just for a nice picture
sigma_max1,sigma_min1,tau_max1, phi_n1,phi_shear1,sigma_avg1=all_stress_measures(S_nodes1, nodes1, dims=mask_area.shape)
sigma_max2,sigma_min2,tau_max2, phi_n2,phi_shear2,sigma_avg2=all_stress_measures(S_nodes2, nodes2, dims=mask_area.shape)
sigma_max_abs1=np.maximum(np.abs(sigma_min1),np.abs(sigma_max1)) ### highest possible norm of the stress tensor
sigma_max_abs2=np.maximum(np.abs(sigma_min2),np.abs(sigma_max2)) ### highest possible norm of the stress tensor



'''
## plotting results:

#plot_stress_vectors(mask_boundaries,stress_vector*0.5*10**10,origin="upper")

#i=np.round(i,2)

# showing method 1 forces
fig=show_quiver(f_x_c2_1,f_y_c2_1,filter=[0,2])
mask_area_show=np.zeros(mask_expanded.shape)+np.nan
mask_area_show[mask_area]=1
fig.gca().imshow(mask_expanded.astype(int)+mask_area_show.astype(int),alpha=0.5)
fig.suptitle("traction forces method1")
ax_cbar=fig.get_axes()[1]
ax_cbar.set_ylabel("traction forces in Pa")
#i = str(np.round(i,2))
#plt.savefig("/media/user/GINA1-BK/data_traktion_force_microscopy/traktion_forces_from_TFM.png", dpi=300)


# showing method 2 forces
fig=show_quiver(f_x_c2_2,f_y_c2_2,filter=[0,2])
mask_area_show=np.zeros(mask_area.shape)+np.nan
mask_area_show[mask_area]=1
fig.gca().imshow(mask_area_show,alpha=0.5)
fig.suptitle("traction forces method2")
ax_cbar=fig.get_axes()[1]
ax_cbar.set_ylabel("traction forces in Pa")
#i = str(np.round(i,2))






'''


fig=show_quiver(t_x_resize,t_y_resize,filter=[0,2])
mask_area_show=np.zeros(mask_area.shape)+np.nan
mask_area_show[mask_area]=1
fig.gca().imshow(mask_area.astype(int)+mask_expanded.astype(int),alpha=0.5)



fig=show_quiver(f_x_c2_2,f_y_c2_2,filter=[0,2])

m=mask_area.astype(int)+mask_expanded.astype(int)*10
m=m.astype(float)
m[~mask_expanded]=np.nan

fig.gca().imshow(m,alpha=0.6)




fig=show_quiver(f_x_c2_1,f_y_c2_1,filter=[0,2])



fig.gca().imshow(m,alpha=0.6)


#fig=plot_arrows(nodes,loads[:,1],loads[:,2],scale_ratio=0.1,title="nodal laods",origin="upper",cbar_str="loads in N",dims=mask_area.shape,mask=mask_expanded,filter=[0,2])
fig=plot_arrows(nodes1,UC1[:,0],UC1[:,1],scale_ratio=0.1,title="deformation method1",origin="upper",dims=mask_area.shape,mask=mask_expanded,filter=[0,4],overlay_mask=mask_area)
fig=plot_arrows(nodes2,UC2[:,0],UC2[:,1],scale_ratio=0.1,title="deformation method2",origin="upper",dims=mask_area.shape,mask=mask_expanded,filter=[0,4],overlay_mask=mask_area)
#plt.savefig("/media/user/GINA1-BK/data_traktion_force_microscopy/fem_sigma_test/deformation%s.png"%str(i),dpi=300)



fig=show_quiver_mask(stress_vector_b[:,:,0],stress_vector_b[:,:,1],mask_inner_edge,filter=False) # stress vector on inner edge as vector
fig.gca().set_title("stress vector on inner edge method1")
fig=show_quiver_mask(f_x_c2_2,f_y_c2_2,mask_inner_edge,filter=False)   # forces on inner egdge on projection method
fig.gca().set_title("forces acting on inner edge method2")


#fig=plot_map(stress_vector_norm_l,origin="upper",title="norm of stress vector on cell boundaries",cbar_str="force in N/pixel" ,mask=mask_boundaries)
fig=plot_map(stress_vector_norm_l1,origin="upper",title="norm of stress vectors method 1",cbar_str="force in N/pixel" ,mask=mask_boundaries,v_range=(-0.5*10**-6,0.5*10**-6))
fig=plot_map(stress_vector_norm_l2,origin="upper",title="norm of stress vectors method 2",cbar_str="force in N/pixel" ,mask=mask_boundaries,v_range=(-0.5*10**-6,0.5*10**-6))


#fig=plot_map(stress_vector_norm_l,origin="upper",title="norm of stress vector on cell boundaries",cbar_str="force in N/pixel" ,mask=mask_boundaries)

sigma_max_abs1_part=np.zeros(sigma_max_abs1.shape)+np.nan
sigma_max_abs1_part[mask_area]=sigma_max_abs1[mask_area]
v_max=np.nanmax([np.nanmax(sigma_max_abs2),np.nanmax(sigma_max_abs1_part)])*ps_new/1.5
fig=plot_map(sigma_max_abs1_part*ps_new,origin="upper",title="maxima principal stress method 1",cbar_str="force in N/pixel" ,mask=mask_expanded,v_range=(0,v_max))
fig=plot_map(sigma_max_abs2*ps_new,origin="upper",title="maxima principal stress method 2",cbar_str="force in N/pixel" ,mask=mask_area,v_range=(0,v_max))


'''
fig=plot_fields(nodes1,fields=[np.abs(S_nodes1[:,0]),np.abs(S_nodes1[:,1]),np.abs(S_nodes1[:,2])],dims=mask_area.shape,titles=["absolute value of\nx_stress","absolute value of\ny_stress","absolute value of\nxy_stress"],cbar_str="stress in N/pixel",origin="upper",mask=mask_area)

#plt.savefig("/media/user/GINA1-BK/data_traktion_force_microscopy/fem_sigma_test/stress_vector%s.png"%str(i),dpi=300)
#plot_map(tau_max,cbar_str="angle",origin="upper",title="orientation of the maximum principle stress")
fig=plot_fields(nodes,fields=[S_nodes[:,0],S_nodes[:,1],S_nodes[:,2]],dims=mask_area.shape,titles=["x_stress","y_stress","xy_stress"],cbar_str="stress in N/pixel",origin="upper",mask=mask_area)
#plt.savefig("/media/user/GINA1-BK/data_traktion_force_microscopy/fem_sigma_test/stress%s.png"%str(i),dpi=300)
fig=plot_fields(nodes,fields=[np.abs(S_nodes[:,0]),np.abs(S_nodes[:,1]),np.abs(S_nodes[:,2])],dims=mask_area.shape,titles=["absolute value of\nx_stress","absolute value of\ny_stress","absolute value of\nxy_stress"],cbar_str="stress in N/pixel",origin="upper",mask=mask_area)
#plt.savefig("/media/user/GINA1-BK/data_traktion_force_microscopy/fem_sigma_test/stress_abs%s.png"%str(i),dpi=300)
#np.save("/media/user/GINA1-BK/data_traktion_force_microscopy/fem_sigma_test/def_data%s.npy"%str(i),UC)



        ## todo:
        # get scaling correct and also youngsmodulus in sigma
        # throughly evaluate the spares lsq solver, this can probably be improved
'''

## conclusions:
# force projection is not cut, a lot of extreme values
# maybe with some mean filter etc, but thats not nice either...
# expanded method is better, maybe also check out constrain traction force microcopy... ask ben about htis
# ausserdem christoph fragen: masss für übereinstimmung