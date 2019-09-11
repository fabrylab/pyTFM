'''
answering the qquestions:
how is the behavior at the border change at when i expand the grid with a layer of zero force?
how big is this influence?
what if fix the grid at some distances form these points__> should be usefull to avoid rotation and such effects

'''
from __future__ import absolute_import, division, print_function
import solidspy.postprocesor as pos
import solidspy.assemutil as ass
import solidspy.solutil as sol
from grid_setup_solids_py import *

folder="/media/user/GINA1-BK/data_traktion_force_microscopy/fem_example2/"
plt.close("all")


shape=(300,300)
# forces
fx=np.zeros(shape)
fy=np.zeros(shape)
mask_p=np.zeros(shape)
fx[100:200,100:150]=-0.0001
fx[100:200,150:200]=0.0001
#fy[30:60]=0
mask_p[np.abs(fx)>0]=1



fig=show_quiver(fx, fy,filter=[0,2], scale_ratio=0.2)
mask_area_show=np.zeros(mask_p.shape)+np.nan
mask_area_show[mask_p.astype(bool)]=1
fig.gca().imshow(mask_area_show,alpha=0.5)

#fig=show_quiver(fx, fy,filter=[0,2], scale_ratio=0.2)

#fig.gca().imshow(mask,alpha=0.5)





#fig, axs = plt.subplots(4, 5)
#fig.set_size_inches(15, 10)
#axs = axs.flatten()
#plt.tight_layout()
v_min=0
v_max=0.0020
for j,i in enumerate([0,50]):
    if i==0:
        mask=mask_p
    else:
        mask=binary_dil(mask_p,iterations=i)

    nodes, elements, loads, mats = grid_setup(mask, fx, fy, 1, 0.49)

    DME, IBC, neq = ass.DME(nodes, elements)


    print("Number of elements: {}".format(elements.shape[0]))
    print("Number of equations: {}".format(neq))

    KG = ass.assembler(elements, mats, nodes, neq, DME,sparse=True)
    RHSG = ass.loadasem(loads, IBC, neq)

    UG_sol,rx = sol.static_sol_cond(KG, RHSG,mask,verbose=False) #solver with constratinst to zero translation and zero rotation
    UG=UG_sol[0]


    if not (np.allclose(KG.dot(UG) / KG.max(), RHSG / KG.max())):
        print("The system is not in equilibrium!")



    UC = pos.complete_disp(IBC, nodes, UG)  # uc are x and y displacements
    E_nodes, S_nodes = pos.strain_nodes(nodes, elements, mats, UC) # stresses and strains
    stress_tensor=calculate_stress_tensor(S_nodes,nodes,dims=mask.shape) # assembling the stress tensor


    ### other possible stress measures, just for a nice picture
    sigma_max,sigma_min,tau_max, phi_n,phi_shear,sigma_avg=all_stress_measures(S_nodes, nodes, dims=mask.shape)
    sigma_max_abs = np.maximum(np.abs(sigma_min), np.abs(sigma_max))  ### highest possible norm of the stress tensor

    if i==0:
        sigma_max_abs_first=copy.deepcopy(sigma_max_abs)

    mask=mask.astype(bool)
    sigma_max_abs_show=np.zeros(sigma_max_abs.shape)+np.nan
    sigma_max_abs_show[mask.astype(bool)]=sigma_max_abs[mask.astype(bool)]
    #sigma_max_show[mask_p.astype(bool)]= sigma_max_first[mask_p.astype(bool)]-sigma_max[mask_p.astype(bool)] # plotting diffrences
    plt.figure()
    im=plt.imshow(sigma_max_abs_show,vmin=v_min,vmax=v_max)
    plt.colorbar(im)

    #fig=plot_map(np.abs(sigma_max),cbar_str="stress in N/pixel",origin="upper",title="absolute value of\nmaximal principal stress componenet",mask=mask)
    #fig=plot_fields(nodes,fields=[S_nodes[:,0],S_nodes[:,1],S_nodes[:,2]],dims=mask.shape,titles=["x_stress","y_stress","xy_stress"],cbar_str="stress in N/pixel",origin="upper",mask=mask)

