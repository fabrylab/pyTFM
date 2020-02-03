# script concerned with evaluating the general accuracy, effects of geometry, pixelsize and fem grid selection.
from skimage import draw
from skimage import measure
import matplotlib.pyplot as plt
import numpy as np
import copy
from contextlib import suppress
from andreas_TFM_package.functions_for_cell_colonie import *
from andreas_TFM_package.grid_setup_solids_py import grid_setup
from andreas_TFM_package.TFM_functions_for_clickpoints import FEM_simulation

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



def setup_traction_force_map(mask,distribution="split_y",fx=0,fy=0):
    forces_x = np.zeros(mask.shape)
    forces_y = np.zeros(mask.shape)
    if distribution=="split_y":
        center=measure.regionprops(mask)[0].centroid
        ind_x,ind_y=np.where(mask) # all pixels of the mask
        ind_x_above=ind_x[ind_y>int(center[0])] # all pixels above the centroid
        ind_x_below = ind_x[ind_y < int(center[0])]  # all pixels below the centroid
        ind_y_above = ind_y[ind_y > int(center[0])]  # all pixels above the centroid
        ind_y_below = ind_y[ind_y < int(center[0])]  # all pixels below the centroid
        # loading the forces
        forces_x[ind_x_above,ind_y_above] = fx
        forces_x[ind_x_below, ind_y_below] = -fx
        forces_y[ind_x_above, ind_y_above] = fy
        forces_y[ind_x_below, ind_y_below] = -fy
    return forces_x,forces_y



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




def add_title(fig,title_str):
    ax = fig.axes[0]
    ax.set_title(title_str, x=0.5, y=0.9, transform=ax.transAxes, color="white")

def general_display(plot_types=[],pixelsize=1):

    if "mask" in plot_types or plot_types=="all":
        fig=plt.figure()
        plt.imshow(mask)
        add_title(fig,"mask")
    if "forces" in plot_types or plot_types=="all":
        fig=show_quiver(fx,fy,scale_ratio=0.05,filter=[0,5],width=0.01, cmap="jet")
        add_title(fig, "forces")
    if "shear" in plot_types or plot_types=="all":
        shear = stress_tensor[:, :, 0, 1]/ (pixelsize * 10 ** -6)  # shear component of the stress tensor
        fig=show_map_clickpoints(shear,cbar_style="out")
        add_title(fig,"shear")
    if "mean_normal_stress" in plot_types or plot_types=="all":
        mean_normal_stress = ((stress_tensor[:, :, 0, 0] + stress_tensor[:, :, 1,1]) / 2) / (pixelsize * 10 ** -6)
        fig=show_map_clickpoints(mean_normal_stress,cbar_style="out")
        add_title(fig, "mean_normal_stress")








###backward workflow

###forward worklow

mask=setup_geometry() # set the system geometry
# standard work flow forward
fx,fy=setup_traction_force_map(mask,fx=1,fy=0) # load forces
#show_quiver(fx,fy)
check_unbalanced_forces(fx,fy)
fx_corr,fy_corr=force_and_torque_correction(fx,fy,mask) # correct forces
nodes, elements, loads, mats = grid_setup(mask, -fx_corr, -fy_corr, 1, sigma=0.5) # construct fem grid
# solve FEM system
UG_sol, stress_tensor=FEM_simulation(nodes, elements, loads, mats, mask, system_type="colony") # colony means the pure neumann FEM is applied
general_display(plot_types="all",pixelsize=1)
standard_measures(mask=mask)