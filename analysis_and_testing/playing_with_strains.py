# script concerned with evaluating the general accuracy, effects of geometry, pixelsize and fem grid selection.
from skimage import draw
from skimage import measure
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import copy
from contextlib import suppress
from pyTFM.utilities_TFM import round_flexible, gaussian_with_nans, make_display_mask, createFolder
from pyTFM.plotting import *
from pyTFM.grid_setup_solids_py import grid_setup
from pyTFM.TFM_functions_for_clickpoints import FEM_simulation
from pyTFM.graph_theory_for_cell_boundaries import mask_to_graph, find_path_circular
import sys
from skimage.morphology import binary_erosion
from skimage.draw import circle
from scipy.ndimage import binary_dilation
from itertools import chain
import os

sys.path.insert(0, '/home/user/Software/tracktion_force_microscopy/tracktion_force_microscopy/analysis_and_testing/')
from analysis_and_testing.simulating_deformation import *
import cv2
from scipy.ndimage.morphology import distance_transform_edt
from itertools import product


def stress_tensor_from_deformation(mask, circ_size=60, type="uniform"):
    if type == "uniform":
        center = measure.regionprops(mask)[0].centroid
        c_x, c_y = np.meshgrid(range(mask.shape[1]), range(mask.shape[0]))  # arrays with all x and y coordinates
        rx = c_x - center[0]  # note maybe its also enough to chose any point as refernece point
        ry = c_y - center[0]

        rx_new = rx * 0.5  # isotropic contraction towards the center by factor 0.5
        ry_new = ry * 0.5
        u = rx - rx_new  # displacement relative to previous
        v = ry - ry_new
        u[~mask.astype(bool)] = 0
        v[~mask.astype(bool)] = 0
    if type == "ring_like_gauss":
        center = measure.regionprops(mask)[0].centroid
        c_x, c_y = np.meshgrid(range(mask.shape[1]), range(mask.shape[0]))  # arrays with all x and y coordinates
        rx = c_x - center[0]  # note maybe its also enough to chose any point as refernece point
        ry = c_y - center[0]

        rx_new = rx * 0.5  # isotropic contraction towards the center by factor 0.5
        ry_new = ry * 0.5
        u = rx - rx_new  # displacement relative to previous
        v = ry - ry_new
        circ = circle(center[0], center[1], circ_size)
        circ_mask = np.zeros(u.shape)
        circ_mask[circ] = 1
        circ_mask = circ_mask.astype(bool)
        u[~circ_mask] = 0
        v[~circ_mask] = 0
        u = gaussian_filter(u, sigma=5)
        v = gaussian_filter(v, sigma=5)
    if type == "ring_like_flat":
        center = measure.regionprops(mask)[0].centroid
        c_x, c_y = np.meshgrid(range(mask.shape[1]), range(mask.shape[0]))  # arrays with all x and y coordinates
        rx = c_x - center[0]  # note maybe its also enough to chose any point as refernece point
        ry = c_y - center[0]

        rx_new = rx * 0.5  # isotropic contraction towards the center by factor 0.5
        ry_new = ry * 0.5
        u = rx - rx_new  # displacement relative to previous
        v = ry - ry_new
        circ = circle(center[0], center[1], circ_size)
        circ_mask = np.zeros(u.shape)
        circ_mask[circ] = 1
        circ_mask = circ_mask.astype(bool)
        disp_lengths = np.sqrt(u ** 2 + v ** 2)
        outer_ring = np.logical_and(~circ_mask, mask.astype(bool))
        u[outer_ring] = (u[outer_ring] / disp_lengths[outer_ring]) * np.max(u[circ_mask])
        v[outer_ring] = (v[outer_ring] / disp_lengths[outer_ring]) * np.max(v[circ_mask])
        # u = gaussian_filter(u, sigma=5)
        # v = gaussian_filter(v, sigma=5)

    # show_quiver(u,v,filter=[0,6])
    # dist=distance_transform_edt(shape1)
    # resized = cv2.resize(dist, (int(im_shape[0]/2),int(im_shape[1]/2)), interpolation=cv2.INTER_CUBIC)
    # center=measure.regionprops(shape1)[0].centroid
    # dist2=np.zeros(im_shape)
    # insert_slices=np.asarray([center[0]-resized.shape[0]/2,center[0]+resized.shape[0]/2,center[1]-resized.shape[1]/2,center[1]+resized.shape[1]/2],dtype=int)
    # dist2[insert_slices[0]:insert_slices[1],insert_slices[2]:insert_slices[3]]=resized

    # dist=(dist*100).astype(np.int32)
    # dist2=(dist2*100).astype(np.int32)
    # u, v, x, y, mask_def, mask_std = calculate_deformation(dist, dist2, window_size=20, overlap=10, std_factor=20)
    return u, v, mask.astype(bool)


def strain_from_deformation(u, v, mask=None):


    du_x = np.gradient(u, axis=1)
    dv_y = np.gradient(v, axis=0)
    dv_x = np.gradient(v, axis=1)
    du_y = np.gradient(u, axis=0)
    e_xx = 0.5 * (du_x + du_x)
    e_xy = 0.5 * (du_y + dv_x)
    e_yx = 0.5 * (dv_x + du_y)
    e_yy = 0.5 * (dv_y + dv_y)
    strain_tensor = np.zeros((e_xx.shape[0], e_xx.shape[1], 2, 2))
    strain_tensor[:, :, 0, 0] = e_xx
    strain_tensor[:, :, 0, 1] = e_yx
    strain_tensor[:, :, 1, 0] = e_xy
    strain_tensor[:, :, 1, 1] = e_yy
    if mask is not None:
        mask = mask.astype(bool)
        strain_tensor[~binary_erosion(mask, iterations=1), :, :] = 0

    return strain_tensor


def stress_from_strain(strain_tensor, E=1, v=0.5):
    # following http://homepages.engineering.auckland.ac.nz/~pkel015/SolidMechanicsBooks/Part_I/BookSM_Part_I/06_LinearElasticity/06_Linear_Elasticity_Complete.pdf
    k_matrix = (E) * np.array([[1 / (1 - v ** 2), v, 0], [v / (1 - v ** 2), 1, 0], [0, 0, 1 / (1 + v)]])
    stress_tensor = np.zeros(strain_tensor.shape)
    for (x, y) in tqdm(product(range(strain_tensor.shape[0]),
                               range(strain_tensor.shape[1]))):  # improve this with array operations
        exx, eyy, exy = strain_tensor[x, y, 0, 0], strain_tensor[x, y, 1, 1], strain_tensor[x, y, 0, 1]
        sxx, syy, sxy = np.matmul(k_matrix, np.array([exx, eyy, exy]))
        stress_tensor[x, y, 0, 0], stress_tensor[x, y, 1, 1], stress_tensor[x, y, 0, 1] = sxx, syy, sxy
    return stress_tensor


def resize_stress_tensor(new_size=(300, 300)):
    new_tensor = np.zeros((new_size[0], new_size[1], 2, 2))
    new_tensor[:, :, 0, 0] = cv2.resize(stress_tensor[:, :, 0, 0], new_size, interpolation=cv2.INTER_CUBIC)
    new_tensor[:, :, 0, 1] = cv2.resize(stress_tensor[:, :, 0, 1], new_size, interpolation=cv2.INTER_CUBIC)
    new_tensor[:, :, 1, 0] = cv2.resize(stress_tensor[:, :, 1, 0], new_size, interpolation=cv2.INTER_CUBIC)
    new_tensor[:, :, 1, 1] = cv2.resize(stress_tensor[:, :, 1, 1], new_size, interpolation=cv2.INTER_CUBIC)
    return new_tensor


if __name__ == "__main__":
    im_shape = (300, 300)
    from plotting_evaluation import display_stress_tensor
    from evaluation_functions import setup_geometry
    mask = setup_geometry(im_shape=im_shape, shape_size=120, shape="circle")

    u, v, mask = stress_tensor_from_deformation(mask, circ_size=60, type="uniform")
    strain_tensor = strain_from_deformation(u, v, mask)
    stress_tensor = stress_from_strain(strain_tensor, E=1, v=0.5)
    stress_tensorRS = resize_stress_tensor(new_size=(300, 300))
    show_quiver(u, v, filter=[0, 8])
    display_stress_tensor(strain_tensor)
    display_stress_tensor(stress_tensor)
    display_stress_tensor(stress_tensorRS)
