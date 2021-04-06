# script concerned with evaluating the general accuracy, effects of geometry, pixelsize and fem grid selection.
from skimage import draw
from skimage import measure
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import copy
from scipy.ndimage import zoom
import clickpoints
from contextlib import suppress
from pyTFM.utilities_TFM import round_flexible, gaussian_with_nans, make_display_mask, createFolder
from pyTFM.plotting import *
from pyTFM.stress_functions import coefficient_of_variation
from pyTFM.grid_setup_solids_py import grid_setup, interpolation, correct_torque
from pyTFM.TFM_functions_for_clickpoints import FEM_simulation, try_to_load_traction
from pyTFM.graph_theory_for_cell_boundaries import mask_to_graph, find_path_circular
import sys
from skimage.measure import regionprops
from collections import defaultdict
from skimage.morphology import binary_erosion
from skimage.draw import circle
from scipy.ndimage import binary_dilation
from scipy.ndimage.morphology import binary_fill_holes
from itertools import chain, product
import os

sys.path.insert(0, '/home/user/Software/pyTFM/analysis_and_testing/')
from simulating_deformation import *
from pyTFM.grid_setup_solids_py import check_unbalanced_forces
import cv2
from scipy.ndimage.morphology import distance_transform_edt
from itertools import product

fields_empty = {key:None for key in["u_b", "v_b", "fx_b", "fy_b", "fx_f", "fy_f", "stress_tensor_f", "stress_tensor_b", "mask_fm", "mask_fem",
         "mask"]}
def load_data(in_folder, pixel_size):
    name_add = "" # "_large

    mask_traction = np.load(os.path.join(in_folder, "traction_mask.npy"))
    mask_traction = binary_fill_holes(mask_traction)
    mask_cell_border = np.load(os.path.join(in_folder, "cell_border_mask.npy"))
    mask_cell_border = binary_fill_holes(mask_cell_border)
    fx = np.load(os.path.join(in_folder, "fx%s.npy"%name_add))
    fy = np.load(os.path.join(in_folder, "fy%s.npy"%name_add))
    tx = fx / ((pixel_size * 10 ** -6) ** 2)
    ty = fy / ((pixel_size * 10 ** -6) ** 2)
    u = np.load(os.path.join(in_folder, "u%s.npy"%name_add))
    v = np.load(os.path.join(in_folder, "v%s.npy"%name_add))

    return mask_traction, mask_cell_border, fx, fy, tx, ty, u, v


def traction_wrapper(u, v, pixelsize, h, young, mask=None, sigma=0.49, filter="gaussian", fs=None, correct_forces=True):
    tx, ty = TFM_tractions(u, v, pixelsize1=pixelsize, pixelsize2=pixelsize, h=h,
                           young=young, sigma=sigma, filter=filter, fs=fs)
    fx, fy = tx * (pixelsize ** 2), ty * (pixelsize ** 2)
    check_unbalanced_forces(fx, fy)
    if correct_forces:
        fx, fy = force_and_torque_correction(fx, fy, mask)  # correct forces

    return fx, fy


def stress_wrapper(mask, fx, fy, young, sigma=0.5,verbose=False):
    nodes, elements, loads, mats = grid_setup(mask, -fx, -fy, young, sigma=0.5, edge_factor=0)  # construct fem grid
    # solve FEM system
    UG_sol, stress_tensor = FEM_simulation(nodes, elements, loads, mats, mask, verbose=verbose,
                                           system_type="colony")  # colony means the pure neumann FEM is applied
    return UG_sol, stress_tensor


def setup_geometry(shape="rectangle", im_shape=(300, 300), shape_size=100):
    if shape == "rectangle":
        start = ((im_shape[0] - shape_size) / 2, (im_shape[1] - shape_size) / 2)
        end = ((im_shape[0] + shape_size) / 2, (im_shape[1] + shape_size) / 2)
        if any([s < 0 for s in start]) or any([end[0] > im_shape[0], end[1] > im_shape[1]]):
            raise Exception("shape is larger the image boundaries")
        x_coords, y_coords = draw.rectangle(start, end, shape=im_shape)
        matrix = np.zeros(im_shape, dtype=int)
        matrix[x_coords.astype(int), y_coords.astype(int)] = 1
    if shape == "circle":
        matrix = np.zeros(im_shape, dtype=int)
        center = regionprops(matrix + 1)[0].centroid
        circ = circle(center[0], center[1], shape_size)
        matrix[circ] = 1

    return matrix


def setup_vector_field(mask, distribution="split_y", x=0, y=0):
    vx = np.zeros(mask.shape)
    vy = np.zeros(mask.shape)
    if distribution == "split_y":
        center = measure.regionprops(mask)[0].centroid
        ind_x, ind_y = np.where(mask)  # all pixels of the mask
        ind_x_above = ind_x[ind_y > int(center[0])]  # all pixels above the centroid
        ind_x_below = ind_x[ind_y < int(center[0])]  # all pixels below the centroid
        ind_y_above = ind_y[ind_y > int(center[0])]  # all pixels above the centroid
        ind_y_below = ind_y[ind_y < int(center[0])]  # all pixels below the centroid
        # loading the forces
        vx[ind_x_above, ind_y_above] = x
        vx[ind_x_below, ind_y_below] = -x
        vy[ind_x_above, ind_y_above] = y
        vy[ind_x_below, ind_y_below] = -y
    return vx, vy


def setup_stress_field(mask, distribution="uniform", sigma_n=1, sigma_shear=0, sigma_gf=4, diameter_gf=4):
    mask = mask.astype(bool)
    sigma_x = np.zeros(mask.shape)
    sigma_y = np.zeros(mask.shape)
    sigma_xy = np.zeros(mask.shape)

    if distribution == "uniform":
        sigma_x[binary_erosion(mask)] = sigma_n  # binary errosion because boundary effects of gradient
        sigma_y[binary_erosion(mask)] = sigma_n
        sigma_xy[binary_erosion(mask)] = sigma_shear
    if distribution == "gaussian_flattened_circle":
        center = regionprops(mask.astype(int))[0].centroid
        shape_length = regionprops(mask.astype(int))[0].equivalent_diameter
        circ = circle(center[0], center[1], radius=shape_length / diameter_gf)
        sigma_x[circ] = 1
        sigma_y[circ] = 1
        sigma_x = gaussian_filter(sigma_x, sigma=sigma_gf)
        sigma_y = gaussian_filter(sigma_y, sigma=sigma_gf)
        sigma_x[~mask] = 0
        sigma_y[~mask] = 0
        # normalizing to get sum on mean of stress of at each pixel to 1
        sigma_x[mask] = (sigma_x[mask] / (np.sum(sigma_x[mask])) * np.sum(mask))
        sigma_y[mask] = (sigma_y[mask] / (np.sum(sigma_y[mask])) * np.sum(mask))
        # sigma_x[binary_erosion(mask)] = sigma_n  # binary errosion because boundary effects of gradient
        # sigma_y[binary_erosion(mask)] = sigma_n

    if distribution == "gaussian_flattened_rectangle":
        sigma_x[binary_erosion(mask, iterations=diameter_gf)] = 1
        sigma_y[binary_erosion(mask, iterations=diameter_gf)] = 1
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
        center = regionprops(mask.astype(int))[0].centroid
        sigma_x[int(center[0]), int(center[1])] = 1
        sigma_y[int(center[0]), int(center[1])] = 1
        sigma_x = gaussian_filter(sigma_x, sigma=sigma_gf)
        sigma_y = gaussian_filter(sigma_y, sigma=sigma_gf)
        mask = mask.astype(bool)
        sigma_x[~mask] = 0
        sigma_y[~mask] = 0
        # normalizing to get sum on mean of stress of at each pixel to 1
        sigma_x[mask] = (sigma_x[mask] / (np.sum(sigma_x[mask])) * np.sum(mask))
        sigma_y[mask] = (sigma_y[mask] / (np.sum(sigma_y[mask])) * np.sum(mask))
        # sigma_x[binary_erosion(mask)] = sigma_n  # binary errosion because boundary effects of gradient
        # sigma_y[binary_erosion(mask)] = sigma_n

    stress_tensor = np.zeros((mask.shape[0], mask.shape[1], 2, 2))
    stress_tensor[:, :, 0, 0] = sigma_x
    stress_tensor[:, :, 0, 1] = sigma_xy
    stress_tensor[:, :, 1, 0] = sigma_xy
    stress_tensor[:, :, 1, 1] = sigma_y

    return stress_tensor


def force_and_torque_correction(f_x, f_y, mask_area):
    f_x_c1 = f_x - np.nanmean(f_x)  # normalizing traction force to sum up to zero (no displacement)
    f_y_c1 = f_y - np.nanmean(f_y)
    f_x_c2, f_y_c2, p = correct_torque(f_x_c1, f_y_c1, mask_area)
    return f_x_c2, f_y_c2,


def standard_measures(mask, pixelsize_tract=1, pixelsize_og=1, mean_normal_list=None, fields=None):
    u_b, v_b, fx_b, fy_b, fx_f, fy_f, mask_fm, stress_tensor_f, stress_tensor_b = [fields[x] for x in
                                                                                   ["u_b", "v_b", "fx_b", "fy_b",
                                                                                    "fx_f", "fy_f", "mask_fm",
                                                                                    "stress_tensor_f",
                                                                                    "stress_tensor_b"]]
    mask = mask.astype(bool)
    mask_fm = mask_fm.astype(bool)

    cv_f, cv_b, cont_energy_b, cont_energy_f, contractile_force_b, contractile_force_f, \
    mean_normal_stress_b, mean_normal_stress_f, mean_shear_b, mean_shear_f, avg_normal_stress, rel_av_norm_stress,max_shear_b,max_shear_f = [
        None for i in range(14)]

    # suppress is equivalent to try and expect...pass
    with suppress(TypeError, NameError): shear_f = stress_tensor_f[:, :, 0, 1] / (
        pixelsize_tract)  # shear component of the stress tensor
    # max_shear_stress
    with suppress(TypeError, NameError):
        max_shear_f = np.sqrt(((stress_tensor_f[:, :, 0, 0] - stress_tensor_f[:, :, 1, 1])/2)**2 + stress_tensor_f[:, :, 0, 1]**2)/pixelsize_tract
        max_shear_f =np.mean(max_shear_f[binary_erosion(mask)])
    with suppress(TypeError, NameError):
        mean_normal_stress_f = ((stress_tensor_f[:, :, 0, 0] + stress_tensor_f[:, :, 1, 1]) / 2) / (
                                                                    pixelsize_tract)
        mean_normal_stress_f = np.mean(mean_normal_stress_f[binary_erosion(mask)]) # mask alone is one tick to big?

    with suppress(TypeError, NameError): mean_shear_f = np.mean(np.abs(shear_f[binary_erosion(mask)]))
    # line tension
    # forces
    with suppress(TypeError, NameError):  energy_points_f = strain_energy_points(u_b, v_b, fx_f, fy_f, pixelsize_og,
                                                                                 pixelsize_tract)  # contractile energy at any point
    # needs interpolation of mask possibley
    with suppress(TypeError, NameError): cont_energy_f = np.sum(energy_points_f[mask_fm])
    with suppress(TypeError, NameError): contractile_force_f, proj_x, proj_y, center = contractillity(fx_f, fy_f,
                                                                                                      pixelsize_tract,
                                                                                                      mask_fm)
    # shear component of the stress tensor
    with suppress(TypeError, NameError): shear_b = stress_tensor_b[:, :, 0, 1] / (
        pixelsize_tract)
    # max_shear_stress
    with suppress(TypeError, NameError):
        max_shear_b = np.sqrt(((stress_tensor_b[:, :, 0, 0] - stress_tensor_b[:, :, 1, 1])/2)**2 + stress_tensor_b[:, :, 0, 1]**2)/pixelsize_tract
        max_shear_b =np.mean(max_shear_b[binary_erosion(mask)])
    with suppress(TypeError, NameError):
        mean_normal_stress_b = ((stress_tensor_b[:, :, 0, 0] + stress_tensor_b[:, :, 1,
                                                    1]) / 2) / (pixelsize_tract)
        mean_normal_stress_b = np.mean(mean_normal_stress_b[binary_erosion(mask)])

    with suppress(TypeError, NameError): mean_shear_b = np.mean(np.abs(shear_b[binary_erosion(mask)]))
    # line tension
    # forces
    with suppress(TypeError, NameError):  energy_points_b = strain_energy_points(u_b, v_b, fx_b, fy_b, pixelsize_og,
                                                                                 pixelsize_tract)  # contractile energy at any point
    # needs interpolation of mask possibley
    with suppress(TypeError, NameError): cont_energy_b = np.sum(energy_points_b[mask_fm])
    with suppress(TypeError, NameError): contractile_force_b, proj_x, proj_y, center = contractillity(fx_b, fy_b,
                                                                                                      pixelsize_tract,
                                                                                                      mask_fm)
    with suppress(TypeError, NameError): cv_f = coefficient_of_variation(binary_erosion(mask), (
            (stress_tensor_f[:, :, 0, 0] + stress_tensor_f[:, :, 1, 1]) / 2) / (pixelsize_tract), 0)
    with suppress(TypeError, NameError): cv_b = coefficient_of_variation(binary_erosion(mask), (
            (stress_tensor_b[:, :, 0, 0] + stress_tensor_b[:, :, 1, 1]) / 2) / (pixelsize_tract), 0)
    with suppress(TypeError, NameError):
        avg_normal_stress_be = [np.nanmean(ms[mask]) for ms in mean_normal_list]
    with suppress(TypeError, NameError):
        rel_av_norm_stress = [x / mean_normal_stress_b for x in avg_normal_stress_be]

    measures = {"mask_fm": mask_fm, "mask": mask, "cv_f": cv_f, "cv_b": cv_b, "cont_energy_b": cont_energy_b,
                "cont_energy_f": cont_energy_f, "contractile_force_b": contractile_force_b,
                "contractile_force_f": contractile_force_f,
                "mean_normal_stress_b": mean_normal_stress_b, "mean_normal_stress_f": mean_normal_stress_f,
                "mean_shear_b": mean_shear_b, "mean_shear_f": mean_shear_f, "max_shear_f":max_shear_f, "max_shear_b":max_shear_b,
                "avg_normal_stress_be": avg_normal_stress_be, "rel_av_norm_stress": rel_av_norm_stress}
    return measures
    # return mean_normal_stress,mean_shear,cont_energy,contractile_force


def contractility_strain_energy_exp(u, v, fx, fy, masks, pixelsize1, pixelsize2):
    # pixelsize1: orgiginal pixel size (deformations are given in pixel with this size)
    # pixelsize2: pixel size of the traction field
    strain_energies = []
    contractilities = []
    for i,m in enumerate(masks):

        m = m.astype(bool)
        ep = strain_energy_points(u, v, fx, fy, pixelsize1, pixelsize2)  # contractile energy at any point
        cenergy = np.sum(ep[m])
        strain_energies.append(cenergy)
        contractile_force, proj_x, proj_y, center = contractillity(fx, fy, pixelsize2, m)
        contractilities.append(contractile_force)


    return np.array(contractilities), np.array(strain_energies)

def simple_height_correction_check(u, v, h_range):  ## implemnent the more advanced one
    '''
    simple function to check the convergence of substrate-height corrected traction force calculated
    :param u:
    :param v:
    :param h_range:
    :return:
    '''

    fig = plt.figure()
    plt.title("traction force convergence")
    ys = []
    hs = []
    for h in tqdm(range(80, 90)):
        tx, ty = ffttc_traction_finite_thickness(u, v, pixelsize1=1, pixelsize2=1, h=h, young=1, sigma=0.49,
                                                 filter="gaussian")
        sum_abs = np.sum(np.sqrt(tx ** 2 + ty ** 2))
        hs.append(h)
        ys.append(sum_abs)
    plt.plot(hs, ys, "o", color="C1")


def compare_scalar_fields(f1, f2, plot_diffrence=False):
    ## correlation coefficient
    r = np.corrcoef(f1.flatten(), f2.flatten())[1, 0]  # r gives  mostly the spatial distribution
    r_squared = r ** 2
    ##mittler absoluter fehler normalirt mit standart abweichung
    # https://de.wikipedia.org/wiki/Mittlerer_absoluter_Fehler
    # https://en.wikipedia.org/wiki/Root-mean-square_deviation#Normalized_root-mean-square_deviation
    abs_error = np.mean(np.abs(f1 - f2))  # absouluter fehler
    xmin = np.nanmin(np.concatenate([f1, f2]))
    xmax = np.nanmax(np.concatenate([f1, f2]))
    xmean = np.nanmean(np.concatenate([f1, f2]))
    abs_error_norm1 = abs_error / (xmax - xmin)
    abs_error_norm2 = abs_error / (xmean)
    ##mean average deviation
    # custom measures tells us: "how much deviates f1 from f2 in % for every pixel
    abs_differences = np.abs(f1 - f2) / f1
    abs_differences[np.isinf(abs_differences)] = np.nan
    avg_deviation = np.nanmean(
        abs_differences)  # not a very good measures because non robust and not exactly zero values are completely ignored

    res = {"r": r, "r_squared": r_squared, "abs_error": abs_error, "abs_error_norm1": abs_error_norm1,
           "abs_error_norm2": abs_error_norm2, "avg_deviation": avg_deviation}
    res = {key: value.astype(float) for key, value in res.items()}
    if plot_diffrence:
        fig, axs = plt.subplots(2, 2)
        axs[0, 0].imshow(f1)
        axs[0, 1].imshow(f2)
        axs[1, 0].imshow(abs_differences)
        axs[1, 1].table(cellText=[[(round_flexible(x, 3))] for x in (res.values())], rowLabels=list(res.keys()),
                        loc='center right', colWidths=[0.3, 0.3])
        axs[1, 1].axis("off")
        plt.show()
    return res


def full_field_comparision(fields):
    fx_b, fy_b, fx_f, fy_f, mask_fm, stress_tensor_f, stress_tensor_b = [fields[x] for x in
                                                                         ["fx_b", "fy_b", "fx_f", "fy_f", "mask_fm",
                                                                          "stress_tensor_f", "stress_tensor_b"]]
    mask_fm = mask_fm.astype(bool)

    res_dict = {}
    with suppress(KeyError): res_dict["forces"] = compare_scalar_fields(np.sqrt(fx_b ** 2 + fy_b ** 2),
                                                                        np.sqrt(fx_f ** 2 + fy_f ** 2))
    with suppress(KeyError): res_dict["forces_inner"] = compare_scalar_fields(
        np.sqrt(fx_b[mask_fm] ** 2 + fy_b[mask_fm] ** 2), np.sqrt(fx_f[mask_fm] ** 2 + fy_f[mask_fm] ** 2))
    with suppress(KeyError): mean_normal_stress_b = ((stress_tensor_b[:, :, 0, 0] + stress_tensor_b[:, :, 1, 1]) / 2)
    with suppress(KeyError): mean_normal_stress_f = ((stress_tensor_f[:, :, 0, 0] + stress_tensor_f[:, :, 1, 1]) / 2)
    with suppress(NameError): res_dict["mean normal stress"] = compare_scalar_fields(mean_normal_stress_b,
                                                                                     mean_normal_stress_f)
    with suppress(KeyError): res_dict["shear"] = compare_scalar_fields(stress_tensor_b[:, :, 0, 1],
                                                                       stress_tensor_f[:, :, 0, 1])
    return res_dict


def gaussian_stress_tensor(stress_tensor, sigma):
    stress_tensor_filtered = np.zeros(stress_tensor.shape)
    stress_tensor_filtered[:, :, 0, 0] = gaussian_with_nans(stress_tensor[:, :, 0, 0], sigma)
    stress_tensor_filtered[:, :, 0, 1] = gaussian_with_nans(stress_tensor[:, :, 0, 1], sigma)
    stress_tensor_filtered[:, :, 1, 0] = gaussian_with_nans(stress_tensor[:, :, 1, 0], sigma)
    stress_tensor_filtered[:, :, 1, 1] = gaussian_with_nans(stress_tensor[:, :, 1, 1], sigma)
    return stress_tensor_filtered


def force_from_stress_field(stress_tensor, pixelsize, plot=True, grad_type="diff1", n=1):
    '''
    relation is simple "force balance"
    tx=d(sigma_xx)/dx + d(sigma_xy)/dy
    ty=d(sigma_yy)/dy + d(sigma_yx)/dx
    :param stress_tensor: 2D stress tensor with 4 elements . ndarray of shape (n,m,2,2)
    :param  pixelsize in m
    :param grad_type: various methods for calculation of the gradient.Gradient uses np.gradient, wich results in a bit
     of interpolation. "diff1" uses diffrence betwen neighbouring pixels and "diffn" uses diffrence between pixels that
     are n steps appart
     :param n distance of pixels in "diffn"
    :return:  fx, fy force field in N (not the traction field!!!!)
    '''
    if grad_type == "gradient":
        dsxx_x = np.gradient(stress_tensor[:, :, 0, 0], axis=1)
        dsyy_y = np.gradient(stress_tensor[:, :, 1, 1], axis=0)
        dsxy_x = np.gradient(stress_tensor[:, :, 0, 1], axis=1)
        dsxy_y = np.gradient(stress_tensor[:, :, 0, 1], axis=0)

    # using diff or gradient??
    if grad_type == "diff1":
        dsxx_x = np.diff(stress_tensor[:, :, 0, 0], axis=1)
        dsyy_y = np.diff(stress_tensor[:, :, 1, 1], axis=0)
        dsxy_x = np.diff(stress_tensor[:, :, 0, 1], axis=1)
        dsxy_y = np.diff(stress_tensor[:, :, 0, 1], axis=0)

        dsxx_x = np.concatenate([dsxx_x, np.zeros((dsxx_x.shape[0], 1))], axis=1)
        dsxy_x = np.concatenate([dsxy_x, np.zeros((dsxy_x.shape[0], 1))], axis=1)
        dsyy_y = np.concatenate([dsyy_y, np.zeros((1, dsyy_y.shape[1]))], axis=0)
        dsxy_y = np.concatenate([dsxy_y, np.zeros((1, dsxy_y.shape[1]))], axis=0)

    if grad_type == "diffn":
        dsxx_x = stress_tensor[:, n:, 0, 0] - stress_tensor[:, :-n, 0, 0]
        dsyy_y = stress_tensor[n:, :, 1, 1] - stress_tensor[:-n, :, 1, 1]
        dsxy_x = stress_tensor[:, n:, 0, 1] - stress_tensor[:, :-n, 0, 1]
        dsxy_y = stress_tensor[n:, :, 0, 1] - stress_tensor[:-n, :, 0, 1]

        dsxx_x = np.concatenate([dsxx_x, np.zeros((dsxx_x.shape[0], n))], axis=1)
        dsxy_x = np.concatenate([dsxy_x, np.zeros((dsxy_x.shape[0], n))], axis=1)
        dsyy_y = np.concatenate([dsyy_y, np.zeros((n, dsyy_y.shape[1]))], axis=0)
        dsxy_y = np.concatenate([dsxy_y, np.zeros((n, dsxy_y.shape[1]))], axis=0)
    tx = dsxx_x + dsxy_y
    ty = dsyy_y + dsxy_x
    fx, fy = tx * (pixelsize ** 2), ty * (pixelsize ** 2)
    return fx, fy


def get_max_values(u_b=None, v_b=None, fx_f=None, fy_f=None, fx_b=None, fy_b=None, stress_tensor_b=None, exp_test=False,
                   mean_normal_list=None):
    max_dict = defaultdict(lambda: None)
    with suppress(NameError, TypeError): max_dict["def"] = np.max(np.sqrt(u_b ** 2 + v_b ** 2))
    with suppress(NameError, TypeError): max_dict["force"] = \
        np.max(np.concatenate([np.sqrt(fx_f ** 2 + fy_f ** 2), np.sqrt(fx_b ** 2 + fy_b ** 2)]))
    with suppress(NameError, TypeError): max_dict["stress"] = np.max(np.concatenate([(stress_tensor_b[:, :, 0, 0] +
                                                                                      stress_tensor_b[:, :, 1, 1]) / 2,
                                                                                     stress_tensor_b[:, :, 1, 0]]))
    if exp_test:
        with suppress(NameError, TypeError): max_dict["stress"] = np.max(mean_normal_list) if len(
            mean_normal_list) > 0 else None
    return max_dict


def save_def_test(fx_b, fy_b, fx_f, fy_f, out_folder, iteration=""):
    fig1, ax = show_quiver(fx_b, fy_b)
    fig2, ax = show_quiver(fx_f, fy_f)
    fig1.savefig(os.path.join(out_folder, "force_backward%s.png" % iteration))
    fig2.savefig(os.path.join(out_folder, "force_forward%s.png" % iteration))


def def_test():
    for i in np.linspace(0, 0.5, 5):
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

        fx_b, fy_b = force_from_stress_field(stress_tensor_b, pixelsize, plot=False, grad_type="diff1")

        u_b, v_b = finite_thickness_convolution(fx_b, fy_b, pixelsize, h, young, sigma=0.5,
                                                kernel_size=None, force_shift=i)  # somwhat of an approximationn
        fx_f, fy_f = traction_wrapper(u_b, v_b, pixelsize, h, young, mask=mask,
                                      filter=filter)  # assuming pixelsize == 1
        save_def_test(fx_b, fy_b, fx_f, fy_f, out_folder, iteration=str(int(i * 10)))




def try_np_load(file):
    try:
        x = np.load(file)
    except FileNotFoundError:
        x = None
    return x


def expand_mask(mask, i, m_shape, method="binary_dilation"):
    if method == "binary_dilation":
        zoom_out = binary_dilation(mask, iterations=i) if i > 0 else copy.deepcopy(mask)
    if method == "cv2":
        zoom_out = zoom(mask, (i, i), order=3)
        dx = mask.shape[0] * i - mask.shape[0]
        dy = mask.shape[0] * i - mask.shape[0]
        zoom_out = zoom_out[int(dx / 2):-int(dx / 2), int(dy / 2):-int(dy / 2)]  # cutting at the edges
        zoom_out = zoom_out[:(mask.shape[0] - 1), :(mask.shape[1] - 1)]  # final cut, might not be necessaryy
    if method == "manual":
        zoom_out = setup_geometry(im_shape=m_shape, shape_size=i, shape="rectangle")
    return zoom_out

def read_header(file):
    with open(file) as f:
        l0 = f.readline()
        l0_split = l0.replace("#", "").strip().split(" ")
        header = {l0_split[0]:float(l0_split[1]),l0_split[2]:float(l0_split[3]), l0_split[4]:l0_split[5]}
    return header

def convert_exp_range(border_ex_range, method, pixelsize):
    # converts list of expansion instruction to µm values based on the pixelsize
    # and the method of expansion
    if method == "binary_dilation":
        range_mu = np.array(border_ex_range) * pixelsize
    if method == "manual":
        # in "manual" a rectangle with the specified edge length is set in the middle of the image
        # we need to correct for zero radius value
        # and we need to devide the values by 2 to get radius
        range_mu = np.array(border_ex_range) * pixelsize/2
        range_mu = range_mu - range_mu[0]
    return range_mu


def exp_border(exp_range=None, fx=None, fy=None, mask=None, young=1, out_folder="", method="binary_dilation", pixelsize =None,verbose=False):
    tensor_folder = os.path.join(out_folder, "stress_tensors")
    mask_folder = os.path.join(out_folder, "masks")
    os.makedirs(tensor_folder, exist_ok=True)
    os.makedirs(mask_folder, exist_ok=True)
    stress_tensors = []
    mask_exp_list = []
    mask = mask.astype(int)
    for i in tqdm(exp_range):
        mask_exp = expand_mask(mask, i, mask.shape, method=method)
        fx_f = fx.copy()
        fy_f = fy.copy()
        fx_f[~mask_exp.astype(bool)] = np.nan  # setting all values outside of mask area to zero
        fy_f[~mask_exp.astype(bool)] = np.nan
        fx_f2 = fx_f - np.nanmean(fx_f)  # normalizing traction force to sum up to zero (no displacement)
        fy_f2 = fy_f - np.nanmean(fy_f)
        fx_f2, fy_f2, p = correct_torque(fx_f2, fy_f2, mask.astype(bool))  # correct forces

        UG_sol, stress_tensor_f = stress_wrapper(mask_exp.astype(bool), fx_f2, fy_f2, young, sigma=0.5,verbose=verbose)
        stress_tensor_f /= pixelsize # conversion to N/m
        np.save(os.path.join(tensor_folder, "stress_tensor_f%d.npy" % i), stress_tensor_f)
        np.save(os.path.join(mask_folder, "%d.npy" % i), mask_exp)
        stress_tensors.append(stress_tensor_f)
        mask_exp_list.append(mask_exp)
    mean_normal_list = [(st[:, :, 0, 0] + st[:, :, 1, 1]) / 2 for st in stress_tensors]

    return stress_tensors, mean_normal_list, mask_exp_list

def load_exp_border(exp_range=[], out_folder=""):
    tensor_folder = os.path.join(out_folder, "stress_tensors")
    mask_folder = os.path.join(out_folder, "masks")
    stress_tensors = []
    mask_exp_list = []
    for i in tqdm(exp_range):
        stress_tensor_f = np.load(os.path.join(tensor_folder, "stress_tensor_f%d.npy" % i))
        mask_exp = np.load(os.path.join(mask_folder, "%d.npy" % i))
        stress_tensors.append(stress_tensor_f)
        mask_exp_list.append(mask_exp)
    mean_normal_list = [(st[:, :, 0, 0] + st[:, :, 1, 1]) / 2 for st in stress_tensors]
    return stress_tensors, mean_normal_list, mask_exp_list


def cut_arrays(fill, arrs, mode="edge", exclude=[]):
    if isinstance(arrs, (list, tuple)):
        ret = []
        if mode == "edge":
            for ar in arrs:
                if isinstance(ar, np.ndarray):
                    ret.append(ar[fill[0]:fill[1], fill[2]:fill[3]])
                else:
                    ret.append(None)
        if mode == "match":
            for ar in arrs:
                if isinstance(ar, np.ndarray):
                    e1 = int((ar.shape[0] - fill[0]) / 2)
                    r1 = (ar.shape[0] - fill[0]) % 2  # potential rest if size diffrence is not divisible by 2
                    e2 = int((ar.shape[1] - fill[1]) / 2)
                    r2 = (ar.shape[1] - fill[1]) % 2  # potential rest if size diffrence is not divisible by 2
                    ret.append(ar[e1:ar.shape[0] - (e1 + r1), e2:ar.shape[1] - (e2 + r2)])
                else:
                    ret.append(None)
    if isinstance(arrs, np.ndarray):
        if mode == "edge":
            ret = arrs[fill[0]:fill[1], fill[2]:fill[3]]

        if mode == "match":
            e1 = int((arrs.shape[0] - fill[0]) / 2)
            r1 = (arrs.shape[0] - fill[0]) % 2  # potential rest if size diffrence is not divisible by 2
            e2 = int((arrs.shape[1] - fill[1]) / 2)
            r2 = (arrs.shape[1] - fill[1]) % 2  # potential rest if size diffrence is not divisible by 2
            ret = arrs[e1:arrs.shape[0] - (e1 + r1), e2:arrs.shape[1] - (e2 + r2)]

    if isinstance(arrs, (dict, defaultdict)):
        ret = {}
        if mode == "edge":
            for key, ar in arrs.items():
                if key not in exclude:
                    if isinstance(ar, np.ndarray):
                        ret[key] = ar[fill[0]:fill[1], fill[2]:fill[3]]
                    else:
                        ret[key] = None
                else:
                    ret[key] = ar

        if mode == "match":
            for key, ar in arrs.items():
                if key not in exclude:
                    if isinstance(ar, np.ndarray):
                        e1 = int((ar.shape[0] - fill[0]) / 2)
                        r1 = (ar.shape[0] - fill[0]) % 2  # potential rest if size diffrence is not divisible by 2
                        e2 = int((ar.shape[1] - fill[0]) / 2)
                        r2 = (ar.shape[1] - fill[0]) % 2  # potential rest if size diffrence is not divisible by 2
                        ret[key] = ar[e1:ar.shape[0] - (e1 + r1), e2:ar.shape[1] - (e2 + r2)]
                    else:
                        ret[key] = None
                else:
                    ret[key] = ar

    return ret


stress_tensor_f = None
