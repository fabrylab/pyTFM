# script concerned with evaluating the general accuracy, effects of geometry, pixelsize and fem grid selection.
from skimage import draw
from skimage import measure
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import copy
import clickpoints
from contextlib import suppress
from pyTFM.utilities_TFM import round_flexible,gaussian_with_nans,make_display_mask,createFolder
from pyTFM.functions_for_cell_colonie import *
from pyTFM.grid_setup_solids_py import grid_setup,interpolation
from pyTFM.TFM_functions_for_clickpoints import FEM_simulation,try_to_load_traction
from pyTFM.graph_theory_for_cell_boundaries import mask_to_graph,find_path_circular
import sys
from collections import defaultdict
from skimage.morphology import binary_erosion
from skimage.draw import circle
from scipy.ndimage import binary_dilation
from scipy.ndimage.morphology import binary_fill_holes
from itertools import chain,product
import os
sys.path.insert(0,'/home/user/Software/pyTFM/analysis_and_testing/')
from simulating_deformation import *
from playing_with_strains import *

import cv2
from scipy.ndimage.morphology  import distance_transform_edt
from itertools import product


def traction_wrapper(u, v, pixelsize, h, young,mask=None, filter="gaussian"):
    tx, ty = ffttc_traction_finite_thickness_wrapper(u, v, pixelsize1=pixelsize, pixelsize2=pixelsize, h=h,
                                                         young=young, sigma=0.49, filter=filter)
    fx, fy = tx*(pixelsize**2), ty*(pixelsize**2)
    check_unbalanced_forces(fx, fy)
    fx_corr, fy_corr = force_and_torque_correction(fx, fy, mask)  # correct forces
    return fx_corr, fy_corr

def stress_wrapper(mask, fx, fy, young, sigma=0.5):
    nodes, elements, loads, mats = grid_setup(mask, -fx, -fy, young, sigma=0.5,edge_factor=0)  # construct fem grid
    # solve FEM system
    UG_sol, stress_tensor = FEM_simulation(nodes, elements, loads, mats, mask, verbose=False,
                                             system_type="colony")  # colony means the pure neumann FEM is applied
    return UG_sol, stress_tensor


def setup_geometry(shape="rectangle", im_shape=(300,300),shape_size=100):
    if shape=="rectangle":
        start=((im_shape[0]-shape_size)/2,(im_shape[1]-shape_size)/2)
        end=((im_shape[0]+shape_size)/2,(im_shape[1]+shape_size)/2)
        if any([s<0 for s in start]) or any([end[0]>im_shape[0],end[1]>im_shape[1]]):
            raise Exception("shape is larger the image boundaries")
        x_coords,y_coords=draw.rectangle(start, end, shape=im_shape)
        matrix=np.zeros(im_shape,dtype=int)
        matrix[x_coords.astype(int),y_coords.astype(int)]=1
    if shape == "circle":
        matrix = np.zeros(im_shape, dtype=int)
        center = regionprops(matrix+1)[0].centroid
        circ=circle(center[0],center[1],shape_size)
        matrix[circ] = 1

    return matrix



def setup_vector_field(mask,distribution="split_y",x=0,y=0):
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

def setup_stress_field(mask,distribution="uniform",sigma_n=1,sigma_shear=0,sigma_gf=4,diameter_gf=4):
    mask=mask.astype(bool)
    sigma_x = np.zeros(mask.shape)
    sigma_y = np.zeros(mask.shape)
    sigma_xy = np.zeros(mask.shape)

    if distribution == "uniform":
        sigma_x[binary_erosion(mask)] = sigma_n # binary errosion because boundary effects of gradient
        sigma_y[binary_erosion(mask)] = sigma_n
    if distribution == "gaussian_flattened_circle":
        center=regionprops(mask.astype(int))[0].centroid
        shape_length=regionprops(mask.astype(int))[0].equivalent_diameter
        circ=circle(center[0],center[1],radius=shape_length/diameter_gf)
        sigma_x[circ] = 1
        sigma_y[circ] = 1
        sigma_x = gaussian_filter(sigma_x, sigma=sigma_gf)
        sigma_y = gaussian_filter(sigma_y, sigma=sigma_gf)
        sigma_x[~mask] = 0
        sigma_y[~mask] = 0
        # normalizing to get sum on mean of stress of at each pixel to 1
        sigma_x[mask] = (sigma_x[mask]/(np.sum(sigma_x[mask]))*np.sum(mask))
        sigma_y[mask] = (sigma_y[mask]/(np.sum(sigma_y[mask]))*np.sum(mask))
        #sigma_x[binary_erosion(mask)] = sigma_n  # binary errosion because boundary effects of gradient
        #sigma_y[binary_erosion(mask)] = sigma_n

    if distribution == "gaussian_flattened_rectangle":

        sigma_x[binary_erosion(mask,iterations=diameter_gf)] = 1
        sigma_y[binary_erosion(mask,iterations=diameter_gf)] = 1
        sigma_x = gaussian_filter(sigma_x, sigma=sigma_gf)
        sigma_y = gaussian_filter(sigma_y, sigma=sigma_gf)
        sigma_x[~mask] = 0
        sigma_y[~mask] = 0
        # normalizing to get sum on mean of stress of at each pixel to 1
        sigma_x[mask] = (sigma_x[mask] / (np.sum(sigma_x[mask])) * np.sum(mask))
        sigma_y[mask] = (sigma_y[mask] / (np.sum(sigma_y[mask])) * np.sum(mask))
        # sigma_x[binary_erosion(mask)] = sigma_n  # binary errosion because boundary effects of gradient
        # sigma_y[binary_erosion(mask)] = sigma_n

    if distribution == "gaussian":
        center=regionprops(mask.astype(int))[0].centroid
        sigma_x[int(center[0]),int(center[1])]=1
        sigma_y[int(center[0]), int(center[1])] = 1
        sigma_x = gaussian_filter(sigma_x, sigma=sigma_gf)
        sigma_y = gaussian_filter(sigma_y, sigma=sigma_gf)
        mask=mask.astype(bool)
        sigma_x[~mask] = 0
        sigma_y[~mask] = 0
        # normalizing to get sum on mean of stress of at each pixel to 1
        sigma_x[mask] = (sigma_x[mask]/(np.sum(sigma_x[mask]))*np.sum(mask))
        sigma_y[mask] = (sigma_y[mask]/(np.sum(sigma_y[mask]))*np.sum(mask))
        #sigma_x[binary_erosion(mask)] = sigma_n  # binary errosion because boundary effects of gradient
        #sigma_y[binary_erosion(mask)] = sigma_n

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

def standard_measures(mask,pixelsize_tract=1,pixelsize_og=1,mean_normal_list=None,stress_tensor_b=None,stress_tensor_f=None):

    mask=mask.astype(bool)
    mask_exp=binary_dilation(mask,iterations=15) #addapt to
    global cont_energy_b,cont_energy_f,contractile_force_b,contractile_force_f,\
        mean_normal_stress_b,mean_normal_stress_f,mean_shear_b,mean_shear_f,avg_normal_stress_be

    cont_energy_b, cont_energy_f, contractile_force_b, contractile_force_f,\
       mean_normal_stress_b, mean_normal_stress_f, mean_shear_b, mean_shear_f,avg_normal_stress=[None for i in range(9)]

    # suppress is equivalent to try and expect...pass
    with suppress(TypeError, NameError): shear_f=stress_tensor_f[:, :, 0, 1] / (pixelsize_tract )  # shear component of the stress tensor
    with suppress(TypeError, NameError): mean_normal_stress_f = ((stress_tensor_f[:, :, 0, 0] + stress_tensor_f[:, :, 1, 1]) / 2) / (pixelsize_tract)
    with suppress(TypeError, NameError): mean_normal_stress_f = np.mean(mean_normal_stress_f[mask])
    with suppress(TypeError, NameError): mean_shear_f = np.mean(shear_f[mask])
    # line tension
    # forces
    with suppress(TypeError, NameError):  energy_points_f = contractile_energy_points(u_b, v_b, fx_f, fy_f, pixelsize_og, pixelsize_tract)  # contractile energy at any point
    # needs interpolation of mask possibley
    with suppress(TypeError, NameError): cont_energy_f = np.sum(energy_points_f[mask_exp])
    with suppress(TypeError, NameError): contractile_force_f, proj_x, proj_y, center = contractillity(fx_f, fy_f, pixelsize_tract, mask_exp)

    with suppress(TypeError, NameError): shear_b=stress_tensor_b[:, :, 0, 1] / (pixelsize_tract)  # shear component of the stress tensor
    with suppress(TypeError, NameError): mean_normal_stress_b = ((stress_tensor_b[:, :, 0, 0] + stress_tensor_b[:, :, 1, 1]) / 2) / (pixelsize_tract)
    with suppress(TypeError, NameError): mean_normal_stress_b = np.mean(mean_normal_stress_b[mask])
    with suppress(TypeError, NameError): mean_shear_b = np.mean(shear_b[mask])
    # line tension
    # forces
    with suppress(TypeError, NameError):  energy_points_b = contractile_energy_points(u_b, v_b, fx_b, fy_b, pixelsize_og, pixelsize_tract)  # contractile energy at any point
    # needs interpolation of mask possibley
    with suppress(TypeError, NameError): cont_energy_b = np.sum(energy_points_b[mask_exp])
    with suppress(TypeError, NameError): contractile_force_b, proj_x, proj_y, center = contractillity(fx_b, fy_b, pixelsize_tract, mask_exp)

    with suppress(TypeError, NameError):
        avg_normal_stress_be=[np.mean(ms[mask]) for ms in mean_normal_list]

    return mask_exp
    #return mean_normal_stress,mean_shear,cont_energy,contractile_force






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
    r=np.corrcoef(f1.flatten(), f2.flatten())[1,0] # r gives  mostly the spatial distribution
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
    return res
def full_field_comparision():
    res_dict={}
    with suppress(NameError): res_dict["forces"]=compare_scalar_fields(np.sqrt(fx_b**2+fy_b**2),np.sqrt(fx_f**2+fy_f**2))
    with suppress(NameError): res_dict["forces_inner"] = compare_scalar_fields(np.sqrt(fx_b[mask_exp] ** 2 + fy_b[mask_exp] ** 2), np.sqrt(fx_f[mask_exp] ** 2 + fy_f[mask_exp] ** 2))
    with suppress(NameError): mean_normal_stress_b = ((stress_tensor_b[:, :, 0, 0] + stress_tensor_b[:, :, 1, 1]) / 2)
    with suppress(NameError): mean_normal_stress_f = ((stress_tensor_f[:, :, 0, 0] + stress_tensor_f[:, :, 1, 1]) / 2)
    with suppress(NameError): res_dict["mean normal stress"] = compare_scalar_fields(mean_normal_stress_b, mean_normal_stress_f)
    with suppress(NameError): res_dict["shear"] = compare_scalar_fields(stress_tensor_b[:, :, 0, 1], stress_tensor_f[:, :, 0, 1])
    return res_dict

def add_title(fig,title_str):
    ax = fig.axes[0]
    ax.set_title(title_str, x=0.5, y=0.9, transform=ax.transAxes, color="white")

def add_mask(fig,mask):
    mask_show = make_display_mask(mask)
    ax = fig.axes[0]
    ax.imshow(mask,alpha=0.4)

def filter_arrows_for_square(fx,fy,filter=6):
    mask_filtered=np.zeros(fx.shape)
    fx_filtered, fy_filtered=np.zeros(fx.shape),np.zeros(fx.shape)

    mask_uf=np.logical_or(fx!=0,fy!=0)
    out_line_graph, points = mask_to_graph(mask_uf, d=np.sqrt(2))
    circular_path = find_path_circular(out_line_graph, 0)
    circular_path=[x for i,x in enumerate(circular_path) if i%filter==0]
    circular_path.append(circular_path[0])  # to plot a fully closed loop
    mask_filtered[points[circular_path][:,0],points[circular_path][:,1]]=1
    mask_filtered=mask_filtered.astype(bool)

    fx_filtered[mask_filtered] = fx[mask_filtered]
    fy_filtered[mask_filtered] = fy[mask_filtered]
    return fx_filtered, fy_filtered

from scipy.signal import convolve2d

def custom_edge_filter(arr):
    arr_out=copy.deepcopy(arr).astype(int)
    shape1 = np.array([[0, 1, 0], [1, 1, 0], [0, 0, 0]])
    shape2 = np.array([[0, 1, 0], [0, 1, 1], [0, 0, 0]])
    shape3 = np.array([[0, 0, 0], [0, 1, 1], [0, 1, 0]])
    shape4 = np.array([[0, 0, 0], [1, 1, 0], [0, 1, 0]])
    for s in [shape1,shape2,shape3,shape4]:
        rem_mask=convolve2d(arr,s,mode="same")==3
        arr_out[rem_mask]=0
    return arr_out.astype(bool)



def display_mask(fig,mask,display_type,type=1,color="C1",d=np.sqrt(2),ax=None):

    mask=mask.astype(int)
    if display_type=="outline":
        out_line=mask-binary_erosion((mask))
        out_line=custom_edge_filter(out_line) # risky
        out_line_graph,points=mask_to_graph(out_line,d=d)
        circular_path=find_path_circular(out_line_graph,0)
        circular_path.append(circular_path[0]) # to plot a fully closed loop
        if type == 1:
            ax = fig.axes[0] if ax is None else ax
            ax.plot(points[circular_path][:,1],points[circular_path][:,0],"--",color=color,linewidth =5)
        if type == 2:
            for ax in fig.axes:
                ax.plot(points[circular_path][:, 1], points[circular_path][:, 0], "--", color=color, linewidth=5)


    if display_type=="overlay":

        mask_show = make_display_mask(mask)
        if type == 1:
            ax = fig.axes[0] if ax is None else ax
            ax.imshow(mask_show,alpha=0.4)
        if type == 2:
            for ax in fig.axes:
                ax.imshow(mask_show,alpha=0.4)

    if display_type == "windowed":
        mask_show = make_display_mask(mask)
        mask_window=copy.deepcopy(mask_show)
        mask_window[np.isnan(mask_show)]=1
        mask_window[mask_show]=np.nan
        if type == 1:
            ax = fig.axes[0] if ax is None else ax
            ax.imshow(mask_window, alpha=0.4)
        if type == 2:
            for ax in fig.axes:
                ax.imshow(mask_window, alpha=0.4)

def display_stress_tensor(stress_tensor,mask=None,title_str=""):

    fig, axs = plt.subplots(1, 3)
    plt.suptitle(title_str)
    im = axs[0].imshow(stress_tensor[:, :, 0, 0])
    axs[0].set_title("sig_xx", x=0.5, y=0.9, transform=axs[0].transAxes, color="white")
    # plt.colorbar(im)
    im = axs[1].imshow(stress_tensor[:, :, 1, 1])
    axs[1].set_title("sig_yy", x=0.5, y=0.9, transform=axs[1].transAxes, color="white")
    im = axs[2].imshow(stress_tensor[:, :, 1, 0])
    axs[2].set_title("sig_xy", x=0.5, y=0.9, transform=axs[2].transAxes, color="white")

    return fig

def general_display(plot_types=[],mask=None,pixelsize=1,display_type="outline",f_type="not circular",cmap="coolwarm",max_dict=defaultdict(lambda: None),
                    mean_normal_list=None, mask_exp_list=None,out_folder="",fx_f=None,fy_f=None,mask_exp=None,
                    border_ex_test=None,scalar_comaprisons=None,plot_gt_exp=True,stress_tensor_f=None,stress_tensor_b=None):

    '''
    plot_types=["deformation_forward","deformation_backwards","mask","forces_forward","forces_forward","shear_forward","mean_normal_stress_forward",
    "shear_backward","mean_normal_stress_backward","full_stress_tensor_forward","full_stress_tensor_backward"]

    :param plot_types:
    :param pixelsize:
    :return:
    '''
    figs={}
    createFolder(out_folder)
    types={"def_f":"deformation_forward","def_b":"deformation_backwards","mask":"mask",
          "f_f":"forces_forward","f_b":"forces_backward","st_f":"full_stress_tensor_forward",
    "st_b":"full_stress_tensor_backward","sh_f":"shear_forward","sh_b":"shear_backward",
           "norm_f":"mean_normal_stress_forward","norm_b":"mean_normal_stress_backward",
           "m":"key measures","r":"correlation","exp_test": "test for border expansion","exp_test2":"test for border expansion"}


    if types["def_f"] in plot_types or plot_types=="all":
        fig, ax = show_quiver(u_f, v_f, scale_ratio=0.05, filter=[0, 5], width=0.007, headlength=3, headwidth=4,
                                      headaxislength=2, cbar_tick_label_size=30, cmap=cmap, cbar_style="not-clickpoints",
                                      vmin=0,vmax=max_dict["def"])
        add_title(fig, types["def_f"])
        display_mask(fig, mask,display_type=display_type,d=np.sqrt(2))
        display_mask(fig, mask_exp, display_type=display_type, color="C2",d=np.sqrt(2))
        fig.savefig(os.path.join(out_folder,types["def_f"]+".png"))
    if types["def_b"] in plot_types or plot_types=="all":
        fig, ax = show_quiver(u_b,v_b,scale_ratio=0.05,filter=[0,5],width=0.007,headlength=3,headwidth=4,
                                    headaxislength=2,cbar_tick_label_size=30,
                                    cmap=cmap,cbar_style="not-clickpoints",vmin=0,vmax=max_dict["def"])
        add_title(fig, types["def_b"])
        display_mask(fig, mask,display_type=display_type)
        display_mask(fig, mask_exp, display_type=display_type,color="C2",d=np.sqrt(2))
        fig.savefig(os.path.join(out_folder,  types["def_b"] + ".png"))

    if types["mask"] in plot_types or plot_types=="all":
        fig=plt.figure()
        plt.imshow(mask)
        add_title(fig, types["mask"])
        fig.savefig(os.path.join(out_folder, types["mask"] + ".png"))

    if types["f_f"] in plot_types or plot_types=="all":
        fig, ax = show_quiver(fx_f, fy_f, scale_ratio=0.05, filter=[0, 5], width=0.007, headlength=3, headwidth=4,
                                      headaxislength=2, cbar_tick_label_size=30, cmap=cmap, cbar_style="not-clickpoints",
                                      vmin=0,vmax=max_dict["force"])
        add_title(fig, types["f_f"])
        display_mask(fig, mask,display_type=display_type)
        display_mask(fig, mask_exp, display_type=display_type, color="C2",d=np.sqrt(2))
        plt.tight_layout()
        fig.savefig(os.path.join(out_folder, types["f_f"] + ".png"))
    if types["f_b"] in plot_types or plot_types=="all":
        if f_type=="circular":
            fx_filtered, fy_filtered=filter_arrows_for_square(fx_b,fy_b,filter=4) # only works when initial forces are circular
            fig, ax = show_quiver(fx_filtered, fy_filtered, scale_ratio=0.05, filter=[0, 0], width=0.007, headlength=3, headwidth=4,
                                      headaxislength=2, cbar_tick_label_size=30, cmap=cmap, cbar_style="not-clickpoints",
                                      vmin=0,vmax=max_dict["force"])
            im=fig.axes[0].imshow(np.sqrt(fx_b**2+fy_b**2),cmap=cmap,
                                      vmin=0,vmax=max_dict["force"])
            fig.axes[1].remove()
            cb=plt.colorbar(im)
            cb.ax.tick_params(labelsize=30)
            plt.tight_layout()
        else:
            fx_filtered, fy_filtered=fx_b,fy_b
            fig, ax = show_quiver(fx_filtered, fy_filtered, scale_ratio=0.05, filter=[0, 5], width=0.007,
                                          headlength=3, headwidth=4,
                                          headaxislength=2, cbar_tick_label_size=30, cmap=cmap, cbar_style="not-clickpoints",
                                          vmin=0,vmax=max_dict["force"])
        add_title(fig, types["f_b"] )
        display_mask(fig, mask,display_type=display_type)
        display_mask(fig, mask_exp, display_type=display_type, color="C2",d=np.sqrt(2))
        fig.savefig(os.path.join(out_folder, types["f_b"] + ".png"))

    if types["sh_f"] in plot_types or plot_types=="all":
        shear = stress_tensor_f[:, :, 0, 1]  # shear component of the stress tensor
        fig,ax = show_map_clickpoints(shear, cbar_style="out", cmap=cmap, cbar_tick_label_size=30, vmin=0,vmax=max_dict["stress"])
        add_title(fig,types["sh_f"])
        display_mask(fig, mask,display_type=display_type)
        fig.savefig(os.path.join(out_folder, types["sh_f"] + ".png"))
    if types["norm_f"] in plot_types or plot_types=="all":
        mean_normal_stress = ((stress_tensor_f[:, :, 0, 0] + stress_tensor_f[:, :, 1, 1]) / 2)
        fig,ax = show_map_clickpoints(mean_normal_stress,cbar_style="out",cmap=cmap,cbar_tick_label_size=30,vmin=0,vmax=max_dict["stress"])
        add_title(fig, types["norm_f"])
        display_mask(fig, mask,display_type=display_type)
        fig.savefig(os.path.join(out_folder, types["norm_f"] + ".png"))

    if types["sh_b"] in plot_types or plot_types=="all":
        shear = stress_tensor_b[:, :, 0, 1]  # shear component of the stress tensor
        fig,ax = show_map_clickpoints(shear,cbar_style="out",cmap=cmap,cbar_tick_label_size=30,vmin=0,vmax=max_dict["stress"])
        display_mask(fig, mask,display_type=display_type)
        add_title(fig,types["sh_b"])
        fig.savefig(os.path.join(out_folder, types["sh_b"] + ".png"))
    if types["norm_b"] in plot_types or plot_types=="all":
        mean_normal_stress = ((stress_tensor_b[:, :, 0, 0] + stress_tensor_b[:, :, 1,1]) / 2)
        fig,ax = show_map_clickpoints(mean_normal_stress, cbar_style="out", cmap=cmap, cbar_tick_label_size=30, vmin=0,vmax=max_dict["stress"])
        add_title(fig, types["norm_b"])
        display_mask(fig, mask,display_type=display_type)
        fig.savefig(os.path.join(out_folder, types["norm_b"] + ".png"))

    if types["st_f"] in plot_types or plot_types=="all":
        fig = display_stress_tensor(stress_tensor_f,mask, title_str=types["st_f"])
        display_mask(fig, mask,display_type=display_type,type=2)
        fig.savefig(os.path.join(out_folder, types["st_f"] + ".png"))
    if types["st_b"] in plot_types or plot_types=="all":
        fig = display_stress_tensor(stress_tensor_b,mask, title_str=types["st_b"])
        display_mask(fig, mask,display_type=display_type,type=2)
        fig.savefig(os.path.join(out_folder, types["st_b"] + ".png"))

    if types["m"] in plot_types or plot_types=="all":
        values = [cont_energy_b,cont_energy_f,contractile_force_b,contractile_force_f,
                mean_normal_stress_b,mean_normal_stress_f,mean_shear_b,mean_shear_f]
        values_r = list(chain.from_iterable([[values[i]/values[i],values[i+1]/values[i]] for i in range(0,len(values),2)]))
        lables=["strain energy","strain energy","contractility","contractility",
                "mean normal stress","mean normal stress","mean shear stress","mean shear stress"]
        pos = list(chain.from_iterable([[p-0.2,p+0.2] for p in range(int(len(values)/2))]))
        fig = plt.figure()
        plt.bar(pos[::2],values_r[::2],width=0.4, color="C1",label="backwards")
        plt.bar(pos[1::2], values_r[1::2], width=0.4, color="C2",label="forwards")
        plt.xticks(pos,lables,rotation="70",fontsize=15)
        plt.yticks(fontsize=20)
        plt.title(types["m"])
        plt.tight_layout()
        fig.savefig(os.path.join(out_folder, types["m"] + ".png"))
    if types["r"] in plot_types or plot_types=="all":
        rs={key:value['r_squared'] for key,value in scalar_comaprisons.items() if not np.isnan(value['r_squared'])}
        fig = plt.figure(figsize=(2,5))
        plt.bar(list(rs.keys()), list(rs.values()), width=0.6, color="C5")
        plt.xticks(rotation="70",fontsize=15)
        plt.yticks(fontsize=20)
        plt.ylim((0,1))
        plt.tight_layout()
        fig.savefig(os.path.join(out_folder, types["r"] + ".png"))
    if types["exp_test"] in plot_types or plot_types=="all" and len(mean_normal_list)>0:
        fig=show_exp_test(mean_normal_list=mean_normal_list,max_dict=max_dict,mask_exp_list=mask_exp_list)
        fig.savefig(os.path.join(out_folder, types["exp_test"] + ".png"))

    if types["exp_test2"] in plot_types or plot_types=="all" and len(mean_normal_list)>0:

        fig =plt.figure()
        avg_normal_stress_be_rel=[x/mean_normal_stress_b for x in avg_normal_stress_be]
        if plot_gt_exp:
            plt.plot(border_ex_test, avg_normal_stress_be_rel,color="C3")
            plt.plot(border_ex_test, [1]*(len(border_ex_test)),color="C4")
            plt.ylim((0, mean_normal_stress_b * 1.2))
        else:
            plt.plot(border_ex_test,  [x/np.max(avg_normal_stress_be) for x in avg_normal_stress_be],
                     color="C3")
            plt.ylim((0,1.1))
        plt.xticks(fontsize=20,rotation="70")
        plt.yticks(fontsize=20)
        plt.title(types["exp_test"])
        plt.tight_layout()
        fig.savefig(os.path.join(out_folder, types["exp_test"] + "2.png"))
    if types["exp_test2"] in plot_types or plot_types=="all" and len(mean_normal_list)>0:
        sub_folder_stress = createFolder(os.path.join(out_folder,"stresses"))
        sub_folder_force = createFolder(os.path.join(out_folder, "forces"))
        for i,(ms,mask_expand) in enumerate(zip(mean_normal_list,mask_exp_list)):
            fig, ax = show_map_clickpoints(ms, cbar_style="out", cmap=cmap, cbar_tick_label_size=30, vmin=0,
                                   vmax=max_dict["stress"])
            add_title(fig, types["norm_b"])
            display_mask(fig, mask, display_type=display_type, color="#FFEF00")
            display_mask(fig, mask_expand, display_type=display_type,color="C3")
            fig.savefig(os.path.join(sub_folder_stress, types["exp_test"] + "%s.png"%str(i)))

            fig, ax = show_quiver(fx_f, fy_f, scale_ratio=0.3, filter=[0, 10], width=0.007, headlength=3,
                              headwidth=4,
                              headaxislength=2, cbar_tick_label_size=30, cmap=cmap, cbar_style="not-clickpoints",
                              vmin=0, vmax=1600)
            add_title(fig, types["f_f"])
            display_mask(fig, mask, display_type=display_type, color="#FFEF00")
            display_mask(fig, mask_expand, display_type=display_type, color="C3")
            plt.tight_layout()
            fig.savefig(os.path.join(sub_folder_force, types["exp_test"] + "%s.png" % str(i)))








def show_exp_test(mean_normal_list=None,max_dict=None,mask_exp_list=None):
    n=int(np.ceil(np.sqrt(len(mean_normal_list))))
    fig,axs=plt.subplots(n,n)
    axs=axs.flatten()
    for i in range(len(mean_normal_list)):
        axs[i].imshow(mean_normal_list[i],vmax=max_dict["stress"])
        axs[i].set_title("mean normal stress #%s"%str(i), x=0.5, y=0.9, transform=axs[i].transAxes, color="white")
        display_mask(fig, mask_exp_list[i], display_type="outline", type=1, color="C1", d=np.sqrt(2), ax=axs[i])

    return fig
def plot_gradient_normal_stress(stress_tensor):

    dsxx_x = np.gradient(stress_tensor[:, :, 0, 0], axis=1)
    dsyy_y = np.gradient(stress_tensor[:, :, 1, 1], axis=0)
    dsxx_y = np.gradient(stress_tensor[:, :, 0, 0], axis=0)
    dsyy_x = np.gradient(stress_tensor[:, :, 1, 1], axis=1)

    fig, axs = plt.subplots(2, 2)
    plt.suptitle("gradient")
    im = axs[0,0].imshow(dsxx_x)
    axs[0,0].set_title("dsxx_x", x=0.5, y=0.9, transform=axs[0,0].transAxes, color="white")
    im = axs[0,1].imshow(dsyy_y)
    axs[0,1].set_title("dsyy_y", x=0.5, y=0.9, transform=axs[0,1].transAxes, color="white")
    im = axs[1,0].imshow(dsxx_y)
    axs[1,0].set_title("dsxx_y", x=0.5, y=0.9, transform=axs[1,0].transAxes, color="white")
    im = axs[1,1].imshow(dsyy_x)
    axs[1,1].set_title("dsyy_x", x=0.5, y=0.9, transform=axs[1,1].transAxes, color="white")




def gaussian_stress_tensor(stress_tensor,sigma):
    stress_tensor_filtered=np.zeros(stress_tensor.shape)
    stress_tensor_filtered[:, :, 0, 0] = gaussian_with_nans(stress_tensor[:, :, 0, 0],sigma)
    stress_tensor_filtered[:, :, 0, 1] = gaussian_with_nans(stress_tensor[:, :, 0, 1],sigma)
    stress_tensor_filtered[:, :, 1, 0] = gaussian_with_nans(stress_tensor[:, :, 1, 0],sigma)
    stress_tensor_filtered[:, :, 1, 1] = gaussian_with_nans(stress_tensor[:, :, 1, 1],sigma)
    return stress_tensor_filtered

def traction_from_stress(stress_tensor,pixelsize,plot=True,grad_type="diff1",n=1):
    '''
    relation is simple "force balance"
    tx=d(sigma_xx)/dx + d(sigma_xy)/dy
    ty=d(sigma_yy)/dy + d(sigma_yx)/dx
    :param stress_tensor:
    :return:
    '''
    if grad_type=="gradient":
        dsxx_x = np.gradient(stress_tensor[:, :, 0, 0], axis=1)
        dsyy_y = np.gradient(stress_tensor[:, :, 1, 1], axis=0)
        dsxy_x = np.gradient(stress_tensor[:, :, 0, 1], axis=1)
        dsxy_y = np.gradient(stress_tensor[:, :, 0, 1], axis=0)

    # using diff or gradient??
    if grad_type=="diff1":
        dsxx_x = np.diff(stress_tensor[:, :, 0, 0], axis =1 )
        dsyy_y = np.diff(stress_tensor[:, :, 1, 1], axis = 0)
        dsxy_x = np.diff(stress_tensor[:, :, 0, 1], axis = 1)
        dsxy_y = np.diff(stress_tensor[:, :, 0, 1], axis = 0)

        dsxx_x = np.concatenate([dsxx_x,np.zeros((dsxx_x.shape[0],1))],axis=1)
        dsxy_x = np.concatenate([dsxy_x,np.zeros((dsxy_x.shape[0],1))],axis=1)
        dsyy_y = np.concatenate([dsyy_y,np.zeros((1,dsyy_y.shape[1]))],axis=0)
        dsxy_y = np.concatenate([dsxy_y,np.zeros((1,dsxy_y.shape[1]))],axis=0)

    if grad_type=="diffn":
        dsxx_x = stress_tensor[:, n:, 0, 0] - stress_tensor[:, :-n, 0, 0]
        dsyy_y = stress_tensor[n: ,: , 1, 1] - stress_tensor[:-n, :, 1, 1]
        dsxy_x = stress_tensor[:,n:, 0, 1] - stress_tensor[:, :-n, 0, 1]
        dsxy_y = stress_tensor[n:,:, 0, 1] - stress_tensor[:-n, :, 0, 1]

        dsxx_x = np.concatenate([dsxx_x, np.zeros((dsxx_x.shape[0], n))], axis=1)
        dsxy_x = np.concatenate([dsxy_x, np.zeros((dsxy_x.shape[0], n))], axis=1)
        dsyy_y = np.concatenate([dsyy_y,np.zeros((n,dsyy_y.shape[1]))],axis=0)
        dsxy_y = np.concatenate([dsxy_y,np.zeros((n,dsxy_y.shape[1]))],axis=0)
    tx = dsxx_x + dsxy_y
    ty = dsyy_y + dsxy_x

    if plot:
        fig, ax = show_quiver(tx, ty, scale_ratio=0.05, filter=[0, 5], width=0.007, headlength=3, headwidth=4,
                          headaxislength=2, cmap="jet")
        add_title(fig, "traction forces from stress")

    fx, fy = tx*(pixelsize**2), ty*(pixelsize**2)
    return fx,fy

def get_max_values(exp_test=False,mean_normal_list=None):
    max_dict =defaultdict(lambda: None)
    with suppress(NameError): max_dict["def"] = np.max(np.sqrt(u_b**2+v_b**2))
    with suppress(NameError): max_dict["force"] = np.max(np.concatenate([np.sqrt(fx_f**2+fy_f**2),np.sqrt(fx_b**2+fy_b**2)]))
    with suppress(NameError): max_dict["stress"] = np.max(np.concatenate([(stress_tensor_b[:, :, 0, 0] + stress_tensor_b[:, :, 1,1]) / 2,stress_tensor_b[:, :, 1, 0]]))
    if exp_test:
        with suppress(NameError): max_dict["stress"]=np.max(mean_normal_list) if len(mean_normal_list)>0 else None
    return max_dict


def save_def_test(fx_b, fy_b,fx_f,fy_f,out_folder,iteration=""):
    fig1,ax=show_quiver(fx_b, fy_b)
    fig2,ax=show_quiver( fx_f, fy_f )
    fig1.savefig(os.path.join(out_folder,"force_backward%s.png"%iteration))
    fig2.savefig(os.path.join(out_folder, "force_forward%s.png" % iteration))

def def_test():
    for i in np.linspace(0,0.5,5):
        out_folder = "/media/user/GINA1-BK/data_traction_force_microscopy/ev_square_square_edge_no_filter"
        createFolder(out_folder)
        print(i)
        young = 1
        h = 100
        pixelsize = 1
        f_type = "circular"  # display_option
        filter = None
        mask = setup_geometry(im_shape=(100, 100), shape_size=50, shape="rectangle")
        stress_tensor_b = setup_stress_field(mask, distribution="uniform", sigma_n=1, sigma_shear=0, sigma_gf=6)

        fx_b, fy_b = traction_from_stress(stress_tensor_b, pixelsize, plot=False, grad_type="diff1")

        u_b, v_b = finite_thickness_convolution(fx_b, fy_b, pixelsize, h, young, sigma=0.5,
                                                kernel_size=None, force_shift=i)  # somwhat of an approximationn
        fx_f, fy_f = traction_wrapper(u_b, v_b, pixelsize, h, young,mask=mask, filter=filter)  # assuming pixelsize == 1
        save_def_test(fx_b, fy_b,fx_f,fy_f,out_folder,iteration=str(int(i*10)))
from scipy.ndimage import zoom


def expand_mask(mask,i,method="binary_dilation"):
    if method=="binary_dilation":
        zoom_out = binary_dilation(mask, iterations=i) if i > 0 else copy.deepcopy(mask)
    if method == "cv2":
        zoom_out=zoom(mask,(i,i),order=3)
        dx=mask.shape[0]*i-mask.shape[0]
        dy=mask.shape[0]*i-mask.shape[0]
        zoom_out=zoom_out[int(dx/2):-int(dx/2),int(dy/2):-int(dy/2)] # cutting at the edges
        zoom_out[:(mask.shape[0]-1),:(mask.shape[1]-1)]# final cut, might not be necessaryy
    if method == "manual":
        zoom_out = setup_geometry(im_shape=(100, 100), shape_size=i, shape="rectangle")
    return zoom_out


def exp_border(exp_range=[],fx_f=None,fy_f=None,mask=None,young=1,out_folder="",method="binary_dilation"):
    stress_tensors=[]
    mask_exp_list=[]
    mask=mask.astype(int)
    for i in tqdm(exp_range):
        mask_exp2=expand_mask(mask,i,method=method)
        UG_sol, stress_tensor_f = stress_wrapper(mask_exp2.astype(bool), fx_f, fy_f, young, sigma=0.5)
        np.save(os.path.join(out_folder, "stress_tensor_f%s.npy"%str(i)), stress_tensor_f)
        print(os.path.join(out_folder, "stress_tensor_f%s.npy"%str(i)))
        plt.figure();plt.imshow(mask_exp2)
        plt.pause(1)
        stress_tensors.append(stress_tensor_f)
        mask_exp_list.append(mask_exp2)
    mean_normal_list=[(st[:, :, 0, 0] + st[:, :, 1,1])/2 for st in stress_tensors]
    return stress_tensors,mean_normal_list,mask_exp_list

def load_exp_border(exp_range=[],mask=None,out_folder=None,method="binary_dilation"):
    stress_tensors=[]
    mask_exp_list=[]
    for i in tqdm(exp_range):
        mask_exp2=expand_mask(mask,i,method=method)
        stress_tensor_f =np.load(os.path.join(out_folder, "stress_tensor_f%s.npy"%str(i)))
        stress_tensors.append(stress_tensor_f)
        mask_exp_list.append(mask_exp2)
    mean_normal_list=[(st[:, :, 0, 0] + st[:, :, 1,1])/2 for st in stress_tensors]
    return stress_tensors,mean_normal_list,mask_exp_list


def exp_border_real_data():
    out_folder = "/media/user/GINA1-BK/data_traction_force_microscopy/ev_border_expansion_real_data/"
    createFolder(out_folder)
    border_ex_test=(list(range(0,40,2)))
    f_type = "non-circular"
    young = 1
    h = 100
    pixelsize = 1
    #  retrieving clickpoints mask and traction forces
    db = clickpoints.DataFile(
        "/media/user/GINA1-BK/data_traction_force_microscopy/WT_vs_KO_images/KOshift/database.cdb", "r")
    mask = db.getMask(frame=2).data == 1
    db.db.close()
    fx_f, fy_f = try_to_load_traction("/media/user/GINA1-BK/data_traction_force_microscopy/WT_vs_KO_images/KOshift/",
                                      "03", warn=False)
    mask=interpolation(mask, dims=fx_f.shape, min_cell_size=100)
    mask=binary_fill_holes(mask)

    #stress_tensors, mean_normal_list, mask_exp_list = exp_border(exp_range=border_ex_test, fx_f=fx_f, fy_f=fy_f,mask=mask,out_folder=out_folder)
    stress_tensors, mean_normal_list, mask_exp_list = load_exp_border(exp_range=border_ex_test,mask=mask,out_folder=out_folder)
    max_dict = get_max_values(exp_test=True,mean_normal_list=mean_normal_list)
    stress_tensor_b=stress_tensors[0]
    mask_exp = standard_measures(mask=mask,mean_normal_list=mean_normal_list,stress_tensor_b=stress_tensor_b,stress_tensor_f=stress_tensor_f)
    #print(avg_normal_stress_be)
    mask_exp = binary_dilation(mask, iterations=15)
    scalar_comaprisons = full_field_comparision()  # r gives  mostly the spatial distribution
    with suppress(KeyError): del scalar_comaprisons["forces"]
    plot_types = ["test for border expansion"]
    plot_types.extend(["forces_forward", "correlation", "test for border expansion"])
    # plot_types = [ "forces_backward", "full_stress_tensor_backward"]
    general_display(plot_types=plot_types, mask=mask, pixelsize=pixelsize, max_dict=max_dict, f_type=f_type,
                    mean_normal_list=mean_normal_list, mask_exp_list=mask_exp_list,out_folder=out_folder,
                    fx_f=fx_f,fy_f=fy_f,mask_exp=mask_exp,scalar_comaprisons=scalar_comaprisons,border_ex_test=border_ex_test,plot_gt_exp=False)
    plt.close("all")

#exp_border_real_data()

if __name__=="__main__":
    out_folder = "/home/user/Desktop/backup_from_harddrive/data_traction_force_microscopy/ev_square_circ_center_1"
    young = 1
    h = 100
    exp_method = ""
    border_ex_test = []
    pixelsize = 1
    filter="gaussian"
    f_type = "non-circular"  # display_option
    mask = setup_geometry(im_shape=(300, 300), shape_size=150)
    stress_tensor_b = setup_stress_field(mask, distribution="gaussian_flattened_circle", sigma_n=1, sigma_shear=0,
                                         sigma_gf=6)

    fx_b, fy_b = traction_from_stress(stress_tensor_b, pixelsize, plot=False, grad_type="diff1")
    u_b, v_b = finite_thickness_convolution(fx_b, fy_b, pixelsize, h, young, sigma=0.5,
                                            kernel_size=None)  # somwhat of an approximation

    createFolder(out_folder)
    np.save(os.path.join(out_folder,"u_b.npy"),u_b)
    np.save(os.path.join(out_folder,"v_b.npy"),v_b)
    u_b,v_b = np.load(os.path.join(out_folder,"u_b.npy")),np.load(os.path.join(out_folder,"v_b.npy"))
    #compare_scalar_fields(u,u1)
    #u2,v2,z1=finite_thickenss_convolution_exact(fx, fy,pixelsize, h, young, sigma=0.5,kernel_size=None) # best solution possible

    ##forward worklow
    # tractions from deformationa
    fx_f, fy_f = traction_wrapper(u_b, v_b,pixelsize,h,young,mask=mask,filter=filter)# assuming pixelsize == 1

    # stress from tractions
    UG_sol, stress_tensor_f = stress_wrapper(mask, fx_f, fy_f, young, sigma=0.5)
    np.save(os.path.join("/home/user/Desktop/", "stress_tensor_f.npy"), stress_tensor_f)
    np.save(os.path.join(out_folder, "stress_tensor_f.npy"), stress_tensor_f)
    stress_tensor_f=np.load(os.path.join(out_folder, "stress_tensor_f.npy"))
    #stress_tensors, mean_normal_list, mask_exp_list=exp_border(exp_range=[1, 2])###### save this
    # typical stress and force measures


    #stress_tensors, mean_normal_list, mask_exp_list = exp_border(out_folder=out_folder,exp_range=border_ex_test,fx_f=fx_f,fy_f=fy_f,mask=mask,method=exp_method)
    stress_tensors, mean_normal_list, mask_exp_list = load_exp_border(out_folder=out_folder,exp_range=border_ex_test,method=exp_method)
    max_dict = get_max_values(exp_test=len(border_ex_test)>0,mean_normal_list=mean_normal_list)


    #max_dict["force"]=None

    # getting comparison scalar fields
    mask_exp = standard_measures(mask=mask,mean_normal_list=mean_normal_list,stress_tensor_b=stress_tensor_b,stress_tensor_f=stress_tensor_f)
    mask_exp = binary_dilation(mask, iterations=15)
    scalar_comaprisons=full_field_comparision() # r gives  mostly the spatial distribution
    del scalar_comaprisons["forces"]

    #


    plot_types = ["test for border expansion"]
    plot_types.extend(["deformation_backward", "mask", "forces_forward", "forces_backward",
      "shear_forward", "mean_normal_stress_forward",
      "shear_backward", "mean_normal_stress_backward", "full_stress_tensor_forward",
      "full_stress_tensor_backward", "key measures", "deformation_backwards", "correlation"
                                ,"test for border expansion"])
    #plot_types = [ "forces_backward", "full_stress_tensor_backward"]
    general_display(plot_types=plot_types, mask=mask, pixelsize=pixelsize, max_dict=max_dict, f_type=f_type,
                    mean_normal_list=mean_normal_list, mask_exp_list=mask_exp_list, out_folder=out_folder,
                    fx_f=fx_f, fy_f=fy_f, mask_exp=mask_exp, scalar_comaprisons=scalar_comaprisons,
                    border_ex_test=border_ex_test, plot_gt_exp=True,stress_tensor_f=stress_tensor_f,stress_tensor_b=stress_tensor_b)
    plt.close("all")
    a=np.sum(np.sqrt(fx_b[mask.astype(bool)]**2+fy_b[mask.astype(bool)]**2))
    b=np.sum(np.sqrt(fx_f[mask.astype(bool)]**2+fy_f[mask.astype(bool)]**2))

'''
from scipy.interpolate import interp2d
from scipy.optimize import approx_fprime

pixx = np.linspace(0, 1, stress_tensor_b.shape[1])  # coordinate space on wich to perfor the interpolation
pixy = np.linspace(0, 1, stress_tensor_b.shape[0])
# using 2 dimensional interpolation on each component

# this returns a continous function f(x,y)=z, that can be evaluated at any point x,y
sig_xx_inter = interp2d(pixx, pixy, stress_tensor_b[:, :, 0, 0], kind="cubic")
sig_pred = sig_xx_inter(pixx, pixy)
delta = (sig_xx_inter(pixx, pixy) - sig_xx_inter(pixx + 0.001, pixy)) / 0.001
plt.figure();
plt.imshow(delta)
plt.figure();
plt.imshow(sig_pred)


show_quiver(fx_b,fy_b,plot_style="clickpoints")
    f_mag=np.sqrt(fx_b**2+fy_b**2)
    f_angle=fx_b/f_mag # angle towards x axis
    f_mag[f_mag==0]=np.nan
    f_angle[f_angle == 0] = np.nan
    f_mag = gaussian_with_nans(f_mag,sigma=6)
    f_angle = gaussian_with_nans(f_angle,sigma=6)
    f_mag[np.isnan(f_mag)] = 0
    f_angle[np.isnan(f_angle)] = 0

    fx_bn = f_angle * f_mag
    fy_bn = np.sin(np.arccos(f_angle)) * f_mag


    ###

    fx_b, fy_b = traction_from_stress(stress_tensor_b, pixelsize, plot=False, grad_type="diff1")
    f_mag=np.sqrt(fx_b**2+fy_b**2)
    f_angle=np.arctan2(fx_b,fy_b) # angle towards x axis
    f_mag[f_mag == 0] = np.nan
    f_angle[f_angle == 0] = np.nan
    f_mag = gaussian_with_nans(f_mag, sigma=6)
    f_angle = gaussian_with_nans(f_angle, sigma=6)
    f_mag[np.isnan(f_mag)] = 0
    f_angle[np.isnan(f_angle)] = 0
    fx_bn = f_mag * f_angle/(np.sqrt(1+f_angle))
    fy_bn = np.sqrt(f_mag**2 - fx_bn**2)

'''