'''
script to evaluate the influenc of poisson ratios on the resulting stresses,
according to
Monolayer Stress Microscopy: Limitations, Artifacts, and
Accuracy of Recovered Intercellular Stresses Tambe 2013
influence for sigma in [0.3,05] should be minimal

'''

from __future__ import absolute_import, division, print_function
#import matplotlib
#matplotlib.use('Agg')

import solidspy.postprocesor as pos
import solidspy.assemutil as ass

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
t_x = np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/11tx_200.npy")
t_y = np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/11ty_200.npy")
# interpolate traction force array:


#t_x_resize=cv.resize(t_x,dsize=(int(mask.shape[1]*0.3),int(mask.shape[0]*0.3)),interpolation=cv.INTER_LINEAR)
#t_y_resize=cv.resize(t_y,dsize=(int(mask.shape[1]*0.3),int(mask.shape[0]*0.3)),interpolation=cv.INTER_LINEAR) ## looks good
t_x_resize=t_x
t_y_resize=t_y

# some pre clean up
mask = remove_small_holes(mask, 100)
mask = remove_small_objects(label(mask), 1000) > 0  # removing other small bits
# interpolation to size of traction force array
mask_int = interpolation(mask, t_x_resize.shape)
# further preparatio of mask data
mask_area, mask_boundaries = prepare_mask(mask_int)


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

stress_vector_norm_list=[]
sigma_max_list=[]
for sigma in [0.2]:# iterating over sigma values
    # setup of the grid
    nodes, elements, loads, mats = grid_setup(mask_area, f_x_c2, f_y_c2, 1, sigma)
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


    UG_sol,rx = custom_solver(KG, RHSG,mask_area,verbose=True) #solver with constratinst to zero translation and zero rotation
    UG=UG_sol[0]
    #UG1,rx = sol.static_sol_cond(KG1, RHSG,mask_area) #solver with constratinst to zero translation and zero rotation
    #norm1=np.sqrt(np.sum((RHSG-np.dot(KG.toarray(),UG[0]))**2)) # same norm as returned by sparse solver


    #norm2=np.sqrt(np.sum((RHSG-np.dot(KG1,UG1[0]))**2))# same norm as returned by sparse solver
    if not (np.allclose(KG.dot(UG) / KG.max(), RHSG / KG.max())):
        print("The system is not in equilibrium!")
    UC = pos.complete_disp(IBC, nodes, UG)  # uc are x and y displacements
    E_nodes, S_nodes = pos.strain_nodes(nodes, elements, mats, UC) # stresses and strains
    stress_tensor=calculate_stress_tensor(S_nodes,nodes,dims=mask_area.shape) # assembling the stress tensor


    graph,points=mask_to_graph(mask_boundaries) # representation of boundaries as graph

    n,n_array=normal_vector_from_graph(graph,points,dims=mask_area.shape) # all nomral vectors of this graph

    #check_normal_vectors_graph(mask_boundaries,n,points) # plotting normal vectors
    #check_normal_vectors_array(mask_boundaries,n_array)

    stress_vector=calculate_stress_vector(n_array,stress_tensor)
    #check_normal_vectors_array(mask_boundaries,stress_vector)
    stress_vector_norm=np.linalg.norm(stress_vector,axis=2)


    ### other possible stress measures, just for a nice picture
    sigma_max,sigma_min,tau_max, phi_n,phi_shear,sigma_avg=all_stress_measures(S_nodes, nodes, dims=mask_area.shape)

    stress_vector_norm_list.append( stress_vector_norm)
    sigma_max_list.append(sigma_max)


## relative change between sigma=0.3 and other values (ignore sigma=0.2)
#sigma_max_rel=[(i-sigma_max_list[1])/np.mean(i-sigma_max_list[1]) for i in sigma_max_list[1:]]
sigma_max_rel=[normalizing(i) for i in sigma_max_list[1:]]
stress_norm_rel=[i-stress_vector_norm_list[1] for i in stress_vector_norm_list[1:]]



# correlation coefficinet between two data sets
# 1 if perfect correlation, 0 if no correlation
mask_boundaries=mask_boundaries.astype(bool)
c1=np.corrcoef(sigma_max_list[1][mask_area],sigma_max_list[3][mask_area])[0,1]
c2=np.corrcoef(sigma_max_list[1][mask_area],sigma_max_list[2][mask_area])[0,1]
c3=np.corrcoef(stress_vector_norm_list[1][mask_boundaries],stress_vector_norm_list[3][mask_boundaries])[0,1]
c4=np.corrcoef(stress_vector_norm_list[1][mask_boundaries],stress_vector_norm_list[2][mask_boundaries])[0,1]
#---> so this is a good measure
### nice illustration of that behavior just ofr me:
plt.figure()
plt.plot(sigma_max_list[1][mask_area],sigma_max_list[3][mask_area],"o")
plt.figure()
plt.plot(stress_vector_norm_list[1][mask_boundaries],stress_vector_norm_list[3][mask_boundaries],"o")

### just show all with same colorbar and range
v_max=np.nanmax(np.hstack(sigma_max_list[1:]))
v_min=np.nanmin(np.hstack(sigma_max_list[1:]))

fig=plot_arrows(nodes,loads[:,1],loads[:,2],origin="upper",title="loads",cbar_str="force in N/pixel" ,dims=mask_area.shape)
for i in sigma_max_list:
    fig = plot_map(i, origin="upper", title="maximal principal stress", cbar_str="force in N/pixel", mask=mask_area,v_range=(v_min,v_max))

    #fig=plot_map(i,origin="upper",title="change in maximal principal stress",cbar_str="force in N/pixel" ,mask=mask_area,v_range=0)




v_max=np.nanmax(np.hstack(stress_vector_norm_lst[1:]))
v_min=np.nanmin(np.hstack(stress_vector_norm_list[1:]))
for i in  stress_vector_norm_list:
    fig=plot_map(i,origin="upper",title="norm of stress vector",cbar_str="force in N/pixel" ,mask=mask_boundaries,v_range=(v_min,v_max))




## conclusion:
#sigma is indeed not important, gives almost the same results....















'''
other measures for similaarity that ae not so important


# rmsd between 0.3 and 0.5
IQR=np.percentile(np.hstack([sigma_max_list[1],sigma_max_list[3]]),0.75)-np.percentile(np.hstack([sigma_max_list[1],sigma_max_list[3]]),0.25)# interquantile range
rmsd=np.sqrt(np.sum((sigma_max_list[1]-sigma_max_list[3])**2)/(np.sum(mask_area)))/IQR# normalized rmsd by interquantile range

IQR=np.percentile(np.hstack([sigma_max_list[1],sigma_max_list[2]]),0.75)-np.percentile(np.hstack([sigma_max_list[1],sigma_max_list[2]]),0.25)# interquantile range
rmsd=np.sqrt(np.sum((sigma_max_list[1]-sigma_max_list[2])**2)/(np.sum(mask_area)))/IQR# normalized rmsd by interquantile range





### this would if e.g all values in array1 would be shifted by constant...
# zero if equal,, dont know if that is a good measure
max_range=np.max(np.hstack([sigma_max_list[1],sigma_max_list[3]]))-np.min(np.hstack([sigma_max_list[1],sigma_max_list[3]]))# interquantile range
er=np.sum(np.abs(sigma_max_list[1]-sigma_max_list[3]))/(np.sum(mask_area)*max_range)# normalized abs error
#should be close to zero


max_range=np.max(np.hstack([sigma_max_list[2],sigma_max_list[3]]))-np.min(np.hstack([sigma_max_list[2],sigma_max_list[3]]))# interquantile range
er=np.sum(np.abs(sigma_max_list[2]-sigma_max_list[3]))/(np.sum(mask_area)*max_range)# normalized abs error


## with a corretion for the mean (.. is that usefull?)
s1_shift=sigma_max_list[1]-np.mean(sigma_max_list[1])
s3_shift=sigma_max_list[3]-np.mean(sigma_max_list[3])
max_range=np.max(np.hstack([s1_shift,s3_shift]))-np.min(np.hstack([s1_shift,s3_shift]))
er=np.sum(np.abs(s1_shift-s3_shift))/(np.sum(mask_area)*max_range)# normalized abs error


s1_shift=sigma_max_list[1]-np.mean(sigma_max_list[1])
s2_shift=sigma_max_list[2]-np.mean(sigma_max_list[2])
max_range=np.max(np.hstack([s1_shift,s2_shift]))-np.min(np.hstack([s1_shift,s2_shift]))
er=np.sum(np.abs(s1_shift-s2_shift))/(np.sum(mask_area)*max_range)# normalized abs error
'''
