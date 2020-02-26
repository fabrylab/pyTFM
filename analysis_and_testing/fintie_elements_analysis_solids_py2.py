from __future__ import absolute_import, division, print_function
#import matplotlib
#matplotlib.use('Agg')
import clickpoints
import numpy as np
import solidspy.postprocesor as pos
import solidspy.assemutil as ass
from andreas_TFM_package.functions_for_cell_colonie import *
from andreas_TFM_package.grid_setup_solids_py import *
from andreas_TFM_package.solids_py_stress_functions import *

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
mask_area = prepare_mask_FEM(mask_int)
#mask_area=np.ones(mask_area.shape).astype(bool)
#mask_area=binary_dil(mask_area,iterations=40)
#plt.figure();plt.imshow(mask_area)
ps_new=pixelsize*np.mean(np.array(mask.shape)/np.array(mask_area.shape)) # pixelsize of fem grid in µm (?? is this ok??)


f_x=t_x_resize*((ps_new*(10**-6))**2) # point force for each node from tractions
f_y=t_y_resize*((ps_new*(10**-6))**2) ## this factor is just for numerical reasons.... ## alsow try to mae this more efficient
f_x[~mask_area]=np.nan # setting all values outside of maske area to zero
f_y[~mask_area]=np.nan
#f_x=f_x.astype("float128")
#f_y=f_x.astype("float128")
f_x_c1=f_x-np.nanmean(f_x) #normalizing traction force to sum up to zero (no displacement)
f_y_c1=f_y-np.nanmean(f_y)
f_x_c2,f_y_c2,p=correct_torque(f_x_c1,f_y_c1,mask_area)

#get_torque1(f_y,f_x,mask_area)


# setup of the grid
nodes, elements, loads, mats = grid_setup(mask_area, f_x_c1, f_y_c1, 1, 0.49)
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
UG_sol,rx =custom_solver(KG, RHSG,mask_area,verbose=True) #solver with constratinst to zero translation and zero rotation


#UG=UG_sol[0]
#UG = sol.static_sol(KG, RHSG,nodes) #no constraints

#UG1,rx = sol.static_sol_cond(KG1, RHSG,mask_area) #solver with constratinst to zero translation and zero rotation
#norm1=np.sqrt(np.sum((RHSG-np.dot(KG.toarray(),UG[0]))**2)) # same norm as returned by sparse solver


#norm2=np.sqrt(np.sum((RHSG-np.dot(KG1,UG1[0]))**2))# same norm as returned by sparse solver
if not (np.allclose(KG.dot(UG_sol) / KG.max(), RHSG / KG.max())):
    print("The system is not in equilibrium!")
UC = pos.complete_disp(IBC, nodes, UG_sol)  # uc are x and y displacements

E_nodes, S_nodes = pos.strain_nodes(nodes, elements, mats, UC) # stresses and strains
stress_tensor=calculate_stress_tensor(S_nodes,nodes,dims=mask_area.shape) # assembling the stress tensor
#pos.fields_plot(elements, nodes, UC, E_nodes=E_nodes, S_nodes=S_nodes)
#####

# average shear and normal stress on the colony area
avg_shear,avg_normal_stress=calculate_mean_stress_measure(mask_area,stress_tensor, ps_new)


### other possible stress measures, just for a nice picture
sigma_max,sigma_min,tau_max, phi_n,phi_shear,sigma_avg=all_stress_measures(S_nodes, nodes, dims=mask_area.shape)
sigma_max_abs=np.maximum(np.abs(sigma_min),np.abs(sigma_max)) ### highest possible norm of the stress tensor



#n_array=c_l.return_n_array()  # calcualte all normal vectors n cell boundaries using the splines

# line stresses for cells method2
graph,points=mask_to_graph(mask_boundaries) # representation of boundaries as graph
n_l,n_array=normal_vector_from_graph(graph,points,dims=mask_area.shape) # all nomral vectors of this graph
stress_vector=calculate_stress_vector(n_array,stress_tensor)
stress_vector_norm=np.linalg.norm(stress_vector,axis=2)  # note: arbitarty choice for orientation at border edges





##
lines_spline_points=borders.lines_spline_points
lines_splines=borders.lines_splines
lines_points=borders.lines_points
# plot linestresses over border as continous curves:
lines_interpol,min_v,max_v=interpolation_for_stress_and_normal_vector(lines_splines,lines_points,stress_tensor,pixel_length=ps_new,interpol_factor=6)
lines_interpol=add_normal_or_shear_component(lines_interpol)

evaluate_all_stress_measures(lines_interpol,borders,norm_levels=["points","lines","cells"],types=["t_vecs","tn","ts"],show_histogramm=True)
fig=plot_continous_boundary_stresses((200,200),borders.edge_lines,lines_interpol,min_v,max_v,
                                     mask_boundaries=borders.mask_boundaries,plot_t_vecs=True,scale_ratio=0.05,arrow_filter=4)

#check_normal_vectors_graph(mask_boundaries,n,points) # plotting normal vectors
#check_normal_vectors_array(mask_boundaries,n_array)





## plotting results:

#plot_stress_vectors(mask_boundaries,stress_vector*0.5*10**10,origin="upper")

#i=np.round(i,2)

### note problem plotting large imasges, please improve
fig=show_quiver(t_x_resize,t_y_resize,filter=[0,3],scale_ratio=0.04)
mask_area_show=np.zeros(mask_area.shape)+np.nan
mask_area_show[mask_area.astype(bool)]=1
fig.gca().imshow(mask_area,alpha=0.5)
fig.suptitle("traction forces on the area of the cell colony")
ax_cbar=fig.get_axes()[1]
ax_cbar.set_ylabel("traction forces in Pa")
#i = str(np.round(i,2))

#plt.savefig("/media/user/GINA1-BK/data_traktion_force_microscopy/traktion_forces_from_TFM.png",
 #           dpi=300)

'''
fig=show_quiver(tx_h,ty_h,filter=[0,2])
mask_area_show=np.zeros(mask_area.shape)+np.nan
mask_area_show[mask_area]=1
fig.gca().imshow(mask_area,alpha=0.5)
'''

fig=plot_arrows(nodes,loads[:,1],loads[:,2],scale_ratio=0.1,title="nodal laods",origin="upper",cbar_str="loads in N",dims=mask_area.shape,mask=mask_area,filter=[0,2])
fig=plot_arrows(nodes,UC[:,0],UC[:,1],scale_ratio=0.1,title="deformation",origin="upper",dims=mask_area.shape,mask=mask_area,filter=[0,4])
#plt.savefig("/media/user/GINA1-BK/data_traktion_force_microscopy/fem_sigma_test/deformation%s.png"%str(i),dpi=300)

fig=plot_map(np.maximum(np.abs(sigma_min),np.abs(sigma_max))*ps_new,cbar_str="stress in N/µm",origin="upper",title="absolute value of\nmaximal principal stress componenet",mask=mask_area,mask_overlay=mask_int)
plt.savefig("/media/user/GINA1-BK/data_traktion_force_microscopy/max_stress.png",dpi=300)
fig=plot_map(stress_vector_norm,origin="upper",title="norm of stress vector on cell boundaries",cbar_str="force in N/pixel" ,mask=mask_boundaries)
#plt.savefig("/media/user/GINA1-BK/data_traktion_force_microscopy/fem_sigma_test/stress_vector%s.png"%str(i),dpi=300)
#plot_map(tau_max,cbar_str="angle",origin="upper",title="orientation of the maximum principle stress")
fig=plot_fields(nodes,fields=[S_nodes[:,0],S_nodes[:,1],S_nodes[:,2]],dims=mask_area.shape,titles=["x_stress","y_stress","xy_stress"],cbar_str="stress in N/pixel",origin="upper",mask=mask_area)#,mask_overlay=mask_int)
#plt.savefig("/media/user/GINA1-BK/data_traktion_force_microscopy/fem_sigma_test/stress%s.png"%str(i),dpi=300)
fig=plot_fields(nodes,fields=[np.abs(S_nodes[:,0]),np.abs(S_nodes[:,1]),np.abs(S_nodes[:,2])],dims=mask_area.shape,titles=["absolute value of\nx_stress","absolute value of\ny_stress","absolute value of\nxy_stress"],cbar_str="stress in N/pixel",origin="upper",mask=mask_area)
#plt.savefig("/media/user/GINA1-BK/data_traktion_force_microscopy/fem_sigma_test/stress_abs%s.png"%str(i),dpi=300)
#np.save("/media/user/GINA1-BK/data_traktion_force_microscopy/fem_sigma_test/def_data%s.npy"%str(i),UC)



        ## todo:
        # get scaling correct and also youngsmodulus in sigma
        # throughly evaluate the spares lsq solver, this can probably be improved

