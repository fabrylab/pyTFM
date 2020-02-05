# script concerned with evaluating the general accuracy, effects of geometry, pixelsize and fem grid selection.
from skimage import draw
from skimage import measure
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import copy
from contextlib import suppress
from andreas_TFM_package.utilities_TFM import round_flexible
from andreas_TFM_package.functions_for_cell_colonie import *
from andreas_TFM_package.grid_setup_solids_py import grid_setup
from andreas_TFM_package.TFM_functions_for_clickpoints import FEM_simulation
import sys
import os
sys.path.insert(0,'/home/user/Software/tracktion_force_microscopy/tracktion_force_microscopy/analysis_and_testing/')
from simulating_deformation import *



def setup_geometry(shape="rectangle", im_shape=(200,200),shape_size=100):
    if shape=="rectangle":
        start=((im_shape[0]-shape_size)/2,(im_shape[1]-shape_size)/2)
        end=((im_shape[0]+shape_size)/2,(im_shape[1]+shape_size)/2)
        if any([s<0 for s in start]) or any([end[0]>im_shape[0],end[1]>im_shape[1]]):
            raise Exception("shape is larger the image boundaries")
        x_coords,y_coords=draw.rectangle(start, end, shape=im_shape)
        matrix=np.zeros(im_shape,dtype=int)
        matrix[x_coords.astype(int),y_coords.astype(int)]=1
        return matrix



def setup_vector_field(mask,distribution="split_y",x=0,y=0):
    mask = mask.astype(bool)
    vx = np.zeros(mask.shape)
    vy = np.zeros(mask.shape)
    if distribution=="split_y":
        center=measure.regionprops(mask)[0].centroid
        ind_x,ind_y=np.where(mask) # all pixels of the mask
        ind_x_above=ind_x[ind_y>int(center[0])] # all pixels above the centroid
        ind_x_below = ind_x[ind_y < int(center[0])]  # all pixels below the centroid
        ind_y_above = ind_y[ind_y > int(center[0])]  # all pixels above the centroid
        ind_y_below = ind_y[ind_y < int(center[0])]  # all pixels below the centroid
        # loading the forces
        vx[ind_x_above,ind_y_above] = x
        vx[ind_x_below, ind_y_below] = -x
        vy[ind_x_above, ind_y_above] = y
        vy[ind_x_below, ind_y_below] = -y
    return vx,vy

def setup_stress_field(mask,distribution="uniform",sigma_n=1,sigma_shear=0):
    mask=mask.astype(bool)
    sigma_x = np.zeros(mask.shape)
    sigma_y = np.zeros(mask.shape)
    sigma_xy = np.zeros(mask.shape)

    if distribution == "uniform":
        sigma_x[mask] = sigma_n
        sigma_y[mask] = sigma_n

    stress_tensor=np.zeros((mask.shape[0],mask.shape[1],2,2))
    stress_tensor[:, :, 0,0] =sigma_x
    stress_tensor[:, :, 0,1] =sigma_xy
    stress_tensor[:, :, 1,0] =sigma_xy
    stress_tensor[:, :, 1,1] =sigma_y

    return stress_tensor

def force_and_torque_correction(f_x,f_y,mask_area):
    f_x_c1 = f_x - np.nanmean(f_x)  # normalizing traction force to sum up to zero (no displacement)
    f_y_c1 = f_y - np.nanmean(f_y)
    f_x_c2, f_y_c2, p = correct_torque(f_x_c1, f_y_c1, mask_area)
    return f_x_c2,f_y_c2,

def standard_measures(tx=None,ty=None,u=None,v=None,stress_tesnor=None,mask=None,pixelsize_tract=None,pixelsize_og=None):
    # stresses
    shear,mean_normal_stress, mean_shear,energy_points, cont_energy, contractile_force=[None for i in range(6)]
    mask=mask.astype(bool)

    # suppress is equivalent to try and expect...pass
    with suppress(TypeError): stress_tensor[:, :, 0, 1] / (pixelsize_tract * 10 ** -6)  # shear component of the stress tensor
    with suppress(TypeError): mean_normal_stress = ((stress_tensor[:, :, 0, 0] + stress_tensor[:, :, 1, 1]) / 2) / (pixelsize_tract * 10 ** -6)
    with suppress(TypeError): mean_normal_stress = np.mean(mean_normal_stress[mask])
    with suppress(TypeError): mean_shear = np.mean(shear[mask])
    # line tension
    # forces
    with suppress(TypeError):  energy_points = contractile_energy_points(u, v, tx, ty, pixelsize_og, pixelsize_tract)  # contractile energy at any point
    # needs interpolation of mask possibley
    with suppress(TypeError): cont_energy = np.sum(energy_points[mask])
    with suppress(TypeError): contractile_force, proj_x, proj_y, center = contractillity(tx, ty, pixelsize_tract, mask)
    return mean_normal_stress,mean_shear,cont_energy,contractile_force






def simple_height_correction_check(u,v,h_range): ## implemnent the more advanced one
    '''
    simple function to check the convergence of substrate-height corrected traction force calculated
    :param u:
    :param v:
    :param h_range:
    :return:
    '''

    fig=plt.figure()
    plt.title("traction force convergence")
    ys=[]
    hs=[]
    for h in tqdm(range(80,90)):
        tx, ty = ffttc_traction_finite_thickness(u, v, pixelsize1=1, pixelsize2=1, h=h, young=1, sigma=0.49,
                                                 filter="gaussian")
        sum_abs=np.sum(np.sqrt(tx**2+ty**2))
        hs.append(h)
        ys.append(sum_abs)
    plt.plot(hs,ys,"o",color="C1")


def compare_scalar_fields(f1,f2,plot_diffrence=False):
    ## correlation coefficient
    r=np.corrcoef(f1.flatten(), f2.flatten())[1,0]
    r_squared=r**2
    ##mittler absoluter fehler normalirt mit standart abweichung
    #https://de.wikipedia.org/wiki/Mittlerer_absoluter_Fehler
    #https://en.wikipedia.org/wiki/Root-mean-square_deviation#Normalized_root-mean-square_deviation
    abs_error=np.mean(np.abs(f1-f2)) # absouluter fehler
    xmin=np.nanmin(np.concatenate([f1,f2]))
    xmax=np.nanmax(np.concatenate([f1,f2]))
    xmean=np.nanmean(np.concatenate([f1,f2]))
    abs_error_norm1=abs_error/(xmax-xmin)
    abs_error_norm2=abs_error/(xmean)
    ##mean average deviation
    # custom measures tells us: "how much deviates f1 from f2 in % for every pixel
    abs_differences=np.abs(f1-f2)/f1
    abs_differences[np.isinf(abs_differences)]=np.nan
    avg_deviation=np.nanmean(abs_differences) # not a very good measures because non robust and not exactly zero values are completely ignored

    res={"r":r,"r_squared":r_squared,"abs_error":abs_error,"abs_error_norm1":abs_error_norm1,"abs_error_norm2":abs_error_norm2,"avg_deviation":avg_deviation}
    res={key:value.astype(float) for key,value in res.items()}
    if plot_diffrence:
        fig,axs=plt.subplots(2,2)
        axs[0, 0].imshow(f1)
        axs[0, 1].imshow(f2)
        axs[1, 0].imshow(abs_differences)
        axs[1, 1].table(cellText=[[(round_flexible(x,3))] for x in (res.values())],rowLabels=list(res.keys()),loc='center right',colWidths=[0.3,0.3])
        axs[1,1].axis("off")
        plt.show()


def add_title(fig,title_str):
    ax = fig.axes[0]
    ax.set_title(title_str, x=0.5, y=0.9, transform=ax.transAxes, color="white")

def general_display(plot_types=[],pixelsize=1):
    if "deformation" in plot_types or plot_types=="all":
        fig=show_quiver(u,v,scale_ratio=0.05,filter=[0,5],width=0.007,headlength=3,headwidth=4,headaxislength=2, cmap="jet")
        add_title(fig, "deformation")
    if "mask" in plot_types or plot_types=="all":
        fig=plt.figure()
        plt.imshow(mask)
        add_title(fig,"mask")
    if "forces" in plot_types or plot_types=="all":
        fig=show_quiver(fx,fy,scale_ratio=0.05,filter=[0,5],width=0.007,headlength=3,headwidth=4,headaxislength=2, cmap="jet")
        add_title(fig, "forces")
    if "shear" in plot_types or plot_types=="all":
        shear = stress_tensor[:, :, 0, 1]  # shear component of the stress tensor
        fig=show_map_clickpoints(shear,cbar_style="out")
        add_title(fig,"shear")
    if "mean_normal_stress" in plot_types or plot_types=="all":
        mean_normal_stress = ((stress_tensor[:, :, 0, 0] + stress_tensor[:, :, 1,1]) / 2)
        fig=show_map_clickpoints(mean_normal_stress,cbar_style="out")
        add_title(fig, "mean_normal_stress")

    if "full_stress_tensor" in plot_types:
        fig,axs=plt.subplots(1,3)
        im=axs[0].imshow(stress_tensor[:, :, 0, 0])
        axs[0].set_title("sig_xx", x=0.5, y=0.9, transform=axs[0].transAxes, color="white")
        #plt.colorbar(im)
        im = axs[1].imshow(stress_tensor[:, :, 1, 1])
        axs[1].set_title("sig_yy", x=0.5, y=0.9, transform=axs[1].transAxes, color="white")
        im = axs[2].imshow(stress_tensor[:, :, 1, 0])
        axs[2].set_title("sig_xy", x=0.5, y=0.9, transform=axs[2].transAxes, color="white")


###backward workflow



young=1
h=100
pixelsize=1




##forward worklow
mask=setup_geometry() # set the system geometry
# standard work flow forward

u,v= setup_vector_field(mask,x=1,y=0)
#simple_height_correction_check(u,v,150)
tx, ty = ffttc_traction_finite_thickness_wrapper(u, v, pixelsize1=pixelsize, pixelsize2=pixelsize, h=h, young=young, sigma = 0.49, filter = "gaussian")
fx,fy = tx,ty # assuming pixelsize == 1
#fx,fy=setup_vector_field(mask,x=1,y=0) # load forces

#show_quiver(fx,fy)
check_unbalanced_forces(fx,fy)
fx_corr,fy_corr=force_and_torque_correction(fx,fy,mask) # correct forces
nodes, elements, loads, mats = grid_setup(mask, -fx_corr, -fy_corr, young, sigma=0.5) # construct fem grid
# solve FEM system
UG_sol, stress_tensor=FEM_simulation(nodes, elements, loads, mats, mask, system_type="colony") # colony means the pure neumann FEM is applied
general_display(plot_types="all",pixelsize=pixelsize)
standard_measures(mask=mask)



## backwards workflow
# deformation from forces
u1,v1,z1=finite_thickness_convolution(fx, fy,pixelsize, h, young, sigma=0.5,kernel_size=None) # somwhat of an approximation
compare_scalar_fields(u,u1)
u2,v2,z1=finite_thickenss_convolution_exact(fx, fy,pixelsize, h, young, sigma=0.5,kernel_size=None) # best solution possible

# tractions from stresses
stress_tensor=setup_stress_field(mask,distribution="uniform",sigma_n=1,sigma_shear=0)

# example stress tensor
stress_tensor=np.load("/media/user/GINA1-BK/data_traction_force_microscopy/WT_vs_KO_images/KOshift/02stress_tensor.npy")
def traction_from_stress(stress_tensor):
    '''
    relation is simple "force balance"
    tx=d(sigma_xx)/dx + d(sigma_xy)/dy
    ty=d(sigma_yy)/dy + d(sigma_yx)/dx
    :param stress_tensor: 
    :return: 
    '''

    dsxx_x=np.gradient(stress_tensor[:,:,0,0],axis=1)
    dsyy_y = np.gradient(stress_tensor[:, :,1,1], axis=0)
    dsxy_x = np.gradient(stress_tensor[:, :, 0, 1], axis=1)
    dsxy_y = np.gradient(stress_tensor[:, :, 0, 1], axis=0)
    tx =  dsxx_x + dsxy_y
    ty = dsyy_y + dsxy_x

  # get a better gradient with spline interpolation ???

  #--> seriously do this....