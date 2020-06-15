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
from pyTFM.plotting import *
from pyTFM.grid_setup_solids_py import grid_setup,interpolation,  correct_torque
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
from plotting_evaluation import *
from pyTFM.grid_setup_solids_py import check_unbalanced_forces
import cv2
from scipy.ndimage.morphology  import distance_transform_edt
from itertools import product


def traction_wrapper(u, v, pixelsize, h, young,mask=None, filter="gaussian",fs=6):
    tx, ty = TFM_tractions(u, v, pixelsize1=pixelsize, pixelsize2=pixelsize,h=h,
                           young=young, sigma=0.49, filter=filter, fs=fs)
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
    mask_exp=binary_dilation(mask,iterations=15) #           addapt to
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
    with suppress(TypeError, NameError):  energy_points_f = strain_energy_points(u_b, v_b, fx_f, fy_f, pixelsize_og, pixelsize_tract)  # contractile energy at any point
    # needs interpolation of mask possibley
    with suppress(TypeError, NameError): cont_energy_f = np.sum(energy_points_f[mask_exp])
    with suppress(TypeError, NameError): contractile_force_f, proj_x, proj_y, center = contractillity(fx_f, fy_f, pixelsize_tract, mask_exp)

    with suppress(TypeError, NameError): shear_b=stress_tensor_b[:, :, 0, 1] / (pixelsize_tract)  # shear component of the stress tensor
    with suppress(TypeError, NameError): mean_normal_stress_b = ((stress_tensor_b[:, :, 0, 0] + stress_tensor_b[:, :, 1, 1]) / 2) / (pixelsize_tract)
    with suppress(TypeError, NameError): mean_normal_stress_b = np.mean(mean_normal_stress_b[mask])
    with suppress(TypeError, NameError): mean_shear_b = np.mean(shear_b[mask])
    # line tension
    # forces
    with suppress(TypeError, NameError):  energy_points_b = strain_energy_points(u_b, v_b, fx_b, fy_b, pixelsize_og, pixelsize_tract)  # contractile energy at any point
    # needs interpolation of mask possibley
    with suppress(TypeError, NameError): cont_energy_b = np.sum(energy_points_b[mask_exp])
    with suppress(TypeError, NameError): contractile_force_b, proj_x, proj_y, center = contractillity(fx_b, fy_b, pixelsize_tract, mask_exp)

    with suppress(TypeError, NameError):
        avg_normal_stress_be=[np.nanmean(ms[mask]) for ms in mean_normal_list]

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
    with suppress(NameError): res_dict["forces_inner"] = compare_scalar_fields(np.sqrt(fx_b[mask_fm] ** 2 + fy_b[mask_fm] ** 2), np.sqrt(fx_f[mask_fm] ** 2 + fy_f[mask_fm] ** 2))
    with suppress(NameError): mean_normal_stress_b = ((stress_tensor_b[:, :, 0, 0] + stress_tensor_b[:, :, 1, 1]) / 2)
    with suppress(NameError): mean_normal_stress_f = ((stress_tensor_f[:, :, 0, 0] + stress_tensor_f[:, :, 1, 1]) / 2)
    with suppress(NameError): res_dict["mean normal stress"] = compare_scalar_fields(mean_normal_stress_b, mean_normal_stress_f)
    with suppress(NameError): res_dict["shear"] = compare_scalar_fields(stress_tensor_b[:, :, 0, 1], stress_tensor_f[:, :, 0, 1])
    return res_dict






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

def get_max_values(u_b=None, v_b=None, fx_f=None, fy_f=None, fx_b=None, fy_b=None, stress_tensor_b=None, exp_test=False, mean_normal_list=None):

    max_dict = defaultdict(lambda: None)
    with suppress(NameError, TypeError): max_dict["def"] = np.max(np.sqrt(u_b**2+v_b**2))
    with suppress(NameError, TypeError): max_dict["force"] = \
        np.max(np.concatenate([np.sqrt(fx_f**2+fy_f**2),np.sqrt(fx_b**2+fy_b**2)]))
    with suppress(NameError, TypeError): max_dict["stress"] = np.max(np.concatenate([(stress_tensor_b[:, :, 0, 0] +
                                                                stress_tensor_b[:, :, 1,1]) / 2,stress_tensor_b[:, :, 1, 0]]))
    if exp_test:
        with suppress(NameError, TypeError): max_dict["stress"]=np.max(mean_normal_list) if len(mean_normal_list)>0 else None
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


def expand_mask(mask,i,m_shape,method="binary_dilation"):
    if method=="binary_dilation":
        zoom_out = binary_dilation(mask, iterations=i) if i > 0 else copy.deepcopy(mask)
    if method == "cv2":
        zoom_out=zoom(mask,(i,i),order=3)
        dx=mask.shape[0]*i-mask.shape[0]
        dy=mask.shape[0]*i-mask.shape[0]
        zoom_out=zoom_out[int(dx/2):-int(dx/2),int(dy/2):-int(dy/2)] # cutting at the edges
        zoom_out[:(mask.shape[0]-1),:(mask.shape[1]-1)]# final cut, might not be necessaryy
    if method == "manual":
        zoom_out = setup_geometry(im_shape=m_shape, shape_size=i, shape="rectangle")
    return zoom_out


def exp_border(exp_range=[],fx_f=None,fy_f=None,mask=None,young=1,out_folder="", method="binary_dilation"):
    tensor_folder = createFolder(os.path.join(out_folder, "stress_tensors"))
    mask_folder = createFolder(os.path.join(out_folder, "masks"))
    stress_tensors=[]
    mask_exp_list=[]
    mask = mask.astype(int)
    for i in tqdm(exp_range):
        mask_exp = expand_mask(mask,i,mask.shape,method=method)
        UG_sol, stress_tensor_f = stress_wrapper(mask_exp.astype(bool), fx_f, fy_f, young, sigma=0.5)
        np.save(os.path.join(tensor_folder, "stress_tensor_f%s.npy"%str(i)), stress_tensor_f)
        np.save(os.path.join(mask_folder, "%s.npy" % str(i)), mask_exp)
        stress_tensors.append(stress_tensor_f)
        mask_exp_list.append(mask_exp)
    mean_normal_list=[(st[:, :, 0, 0] + st[:, :, 1,1])/2 for st in stress_tensors]
    return stress_tensors,mean_normal_list,mask_exp_list

def load_exp_border(exp_range=[],out_folder=None):
    stress_tensors=[]
    mask_exp_list=[]
    for i in tqdm(exp_range):
        stress_tensor_f = np.load(os.path.join(out_folder, "stress_tensors/stress_tensor_f%s.npy"%str(i)))
        mask_exp = np.load(os.path.join(out_folder, "masks/%s.npy" % str(i)))
        stress_tensors.append(stress_tensor_f)
        mask_exp_list.append(mask_exp)
    mean_normal_list=[(st[:, :, 0, 0] + st[:, :, 1,1])/2 for st in stress_tensors]
    return stress_tensors, mean_normal_list, mask_exp_list



def exp_border_real_data():
    out_folder = "/home/user/Desktop/backup_from_harddrive/data_traction_force_microscopy/ev_paper_rd_expansion_fs3"
    createFolder(out_folder)
    border_ex_test=(list(range(0,100,2)))
    f_type = "non-circular"
    young = 1
    h = 100
    pixelsize = 1
    filter="gaussian"
    #  retrieving clickpoints mask and traction forces
    folder="/home/user/Desktop/backup_from_harddrive/data_traction_force_microscopy/WT_vs_KO_images/KOshift/"
    db = clickpoints.DataFile(os.path.join(folder,"database.cdb"), "r")
    mask = db.getMask(frame=2).data == 3
    db.db.close()
    u, v = np.load(os.path.join(out_folder,"u.npy")), np.load(os.path.join(out_folder,"v.npy"))
    fx_f, fy_f = traction_wrapper(u, v, pixelsize, h, young, mask=mask,
                                  filter="gaussian", fs=3) # this filtersize is equal to 3*0.85 ~3.5 µm for real data
    mask = interpolation(mask, dims=fx_f.shape, min_cell_size=100)
    mask = binary_fill_holes(mask)
    np.save(os.path.join(out_folder,"mask.npy"), mask)
    stress_tensors, mean_normal_list, mask_exp_list = exp_border(exp_range=border_ex_test, fx_f=fx_f, fy_f=fy_f, mask=mask, out_folder=out_folder, method="binary_dilation")

    stress_tensor_b=stress_tensors[0]
    max_dict = get_max_values(fx_f=fx_f, fy_f=fy_f, stress_tensor_b=stress_tensor_b,
                             exp_test=len(border_ex_test) > 0, mean_normal_list=mean_normal_list)
    # getting comparison scalar fields
    mask_fm = standard_measures(mask=mask, mean_normal_list=mean_normal_list, stress_tensor_b=stress_tensor_b)
    save_arr = np.array([np.round(np.array(avg_normal_stress_be),5),np.array(border_ex_test)]).T
    np.savetxt(os.path.join(out_folder,"avg_norm_stress.txt"), save_arr, fmt="%.5f", delimiter=",")



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

def cut_arrays(fill,*args):
    ret_list=[]
    for ar in args:
        ret_list.append(ar[fill[0]:fill[1],fill[2]:fill[3]])
    return ret_list

exp_border_real_data()
if __name__=="_main__":
    out_folder = "/home/user/Desktop/backup_from_harddrive/data_traction_force_microscopy/ev_paper_expansion"
    young = 1
    h = 100
    pixelsize = 1
    border_ex_test = list(range(150,230))
    filter = "gaussian"
    exp_method = "manual"
    f_type = "circular"  # display_option
    mask = setup_geometry(im_shape=(300, 300), shape_size=150, shape="rectangle")
    mask_fem = expand_mask(mask, 185, mask.shape, method= "manual")  # highest aggreement (not completely sure
    stress_tensor_b = setup_stress_field(mask, distribution="uniform", sigma_n=1, sigma_shear=0,
                                         sigma_gf=6)
    fx_b, fy_b = traction_from_stress(stress_tensor_b, pixelsize, plot=False, grad_type="diff1")

    u_b, v_b = finite_thickness_convolution(fx_b, fy_b, pixelsize, h, young, sigma=0.5,
                                            kernel_size=None)  # somwhat of an approximation

    createFolder(out_folder)
    np.save(os.path.join(out_folder, "u_b.npy"), u_b)
    np.save(os.path.join(out_folder, "v_b.npy"), v_b)
    u_b, v_b = np.load(os.path.join(out_folder, "u_b.npy")), np.load(os.path.join(out_folder, "v_b.npy"))
    # compare_scalar_fields(u,u1)
    # u2,v2,z1=finite_thickenss_convolution_exact(fx, fy,pixelsize, h, young, sigma=0.5,kernel_size=None) # best solution possible

    ##forward worklow
    # tractions from deformationa
    fx_f, fy_f = traction_wrapper(u_b, v_b, pixelsize, h, young, mask=mask,
                                  filter=filter)  # assuming pixelsize == 1

    # stress from tractions
    UG_sol, stress_tensor_f = stress_wrapper(mask_fem, fx_f, fy_f, young, sigma=0.5)
    np.save(os.path.join("/home/user/Desktop/", "stress_tensor_f.npy"), stress_tensor_f)
    np.save(os.path.join(out_folder, "stress_tensor_f.npy"), stress_tensor_f)
    stress_tensor_f = np.load(os.path.join(out_folder, "stress_tensor_f.npy"))
    # stress_tensors, mean_normal_list, mask_exp_list=exp_border(exp_range=[1, 2])###### save this
    # typical stress and force measures

    stress_tensors, mean_normal_list, mask_exp_list = exp_border(out_folder=out_folder, exp_range=border_ex_test,
                                                                 fx_f=fx_f, fy_f=fy_f, mask=mask, method=exp_method)
    stress_tensors, mean_normal_list, mask_exp_list = load_exp_border(out_folder=out_folder,
                                                                      exp_range=border_ex_test)

    # getting comparison scalar fields
    mask_fm = standard_measures(mask=mask, mean_normal_list=mean_normal_list, stress_tensor_b=stress_tensor_b,
                                stress_tensor_f=stress_tensor_f)
    save_arr = np.array([np.round(np.array(avg_normal_stress_be),5),np.array(border_ex_test)]).T
    np.savetxt(os.path.join(out_folder,"avg_norm_stress.txt"), save_arr, fmt="%.5f", delimiter=",")

    mask_fm = binary_dilation(mask, iterations=25)
    scalar_comaprisons = full_field_comparision()  # r gives  mostly the spatial distribution
    del scalar_comaprisons["forces"]
    key_values = [cont_energy_b, cont_energy_f, contractile_force_b, contractile_force_f,
                  mean_normal_stress_b, mean_normal_stress_f, mean_shear_b, mean_shear_f]

    # cutting to ignore effects close to image edge #### mention this when talking to benn
    fx_b, fy_b, fx_f, fy_f, stress_tensor_b, stress_tensor_f, mask, mask_fem, mask_fm = cut_arrays([20, -20, 20, -20],
                                                                                                   fx_b, fy_b, fx_f,
                                                                                                   fy_f,
                                                                                                   stress_tensor_b,
                                                                                                   stress_tensor_f,
                                                                                                   mask, mask_fem,
                                                                                                   mask_fm)
    max_dict = get_max_values(u_b, v_b, fx_f, fy_f, fx_b, fy_b, stress_tensor_b,
                              exp_test=len(border_ex_test) > 0, mean_normal_list=mean_normal_list)


    plot_types = [""]
    plot_types.extend(["deformation_backward", "mask", "forces_forward", "forces_backward",
                       "shear_forward", "mean_normal_stress_forward",
                       "shear_backward", "mean_normal_stress_backward", "full_stress_tensor_forward",
                       "full_stress_tensor_backward", "key measures", "deformation_backwards", "correlation"])
    # ,"test for border expansion"])
    # plot_types = [ "forces_backward", "full_stress_tensor_backward"]
    plot_types = ["correlation", "key measures", "mean_normal_stress_backward", "mean_normal_stress_forward",
                  "forces_forward", "forces_backward", "mask_outline", "cbars_only"]
    #plot_types = ["correlation", "cbars_only"]
    #max_dict["force"] = 1

    general_display(plot_types=plot_types, mask=mask, pixelsize=pixelsize, max_dict=max_dict, f_type=f_type,
                    out_folder=out_folder, cmap="coolwarm", fx_b=fx_b, fy_b=fy_b,
                    fx_f=fx_f, fy_f=fy_f, mask_fm=mask_fm, mask_fem=mask_fem, scalar_comaprisons=scalar_comaprisons,
                    border_ex_test=border_ex_test, stress_tensor_f=stress_tensor_f,
                    stress_tensor_b=stress_tensor_b, key_values=key_values, plot_gt_exp=False, dm=False, at=False,
                    cb=False)

    plt.close("all")
