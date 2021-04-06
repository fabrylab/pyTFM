# this script calculates the contracility and avg. mena normal stress for diffrent FEM-grid areas
# for a real data example and the standart artififciall sytem

#

import sys
sys.path.append("/home/user/Software/pyTFM/analysis_and_testing")
from evaluation_functions import *
from plotting_evaluation import *
from pyTFM.plotting import show_quiver
import os
from PIL import Image
from pyTFM.frame_shift_correction import correct_stage_drift,croping_after_shift, shift
from pyTFM.parameters_and_strings import default_parameters
from pyTFM.grid_setup_solids_py import correct_forces, interpolation
from pyTFM.stress_functions import all_stress_measures

def shift_arrays(im ,shift_x, shift_y, other_ims=None):
    other_ims = [] if other_ims is None else other_ims
    image_shift = shift(im, shift=(-shift_y, -shift_x), order=5)
    image_shift = croping_after_shift(image_shift, shift_x, shift_y)
    other_ims_shift = [croping_after_shift(im, shift_x, shift_y) for im in other_ims]
    return image_shift, other_ims_shift
# get initial drift

#input parameters
new_analysis =  False
out_folder = "/home/user/Desktop/paper_pyTFM/data_traction_force_microscopy/drift/"
pixel_size = 0.201
std_factor = 15
h = 300
young = 23500
sigma = 0.49
fs = 3

# loading input
im1 = "/home/user/Desktop/paper_pyTFM/data_traction_force_microscopy/data_raw_drift/12_fluo_after.tif"
im2 = "/home/user/Desktop/paper_pyTFM/data_traction_force_microscopy/data_raw_drift/12_fluo_before.tif"
im_cells = "/home/user/Desktop/paper_pyTFM/data_traction_force_microscopy/data_raw_drift/12_BF.tif"
im1 = np.array(Image.open(im1))
im2 = np.array(Image.open(im2))
im_cells = np.array(Image.open(im_cells))
im1 = im1.astype(np.int32)
im2 = im2.astype(np.int32)
mask_border = np.load("/home/user/Desktop/paper_pyTFM/data_traction_force_microscopy/data_raw_drift/cell_border_mask.npy")
mask_border = binary_fill_holes(mask_border)
mask_tractions = np.load("/home/user/Desktop/paper_pyTFM/data_traction_force_microscopy/data_raw_drift/traction_mask.npy")
mask_tractions = binary_fill_holes(mask_tractions)


im1_shift, im2_shift, (mask_border_shift, mask_tractions_shift, im_cells_shift), shift_value = correct_stage_drift(im1, im2,
                    additional_images=[mask_border, mask_tractions, im_cells])
#im2_shift = (np.array(im2_shift))
#im1_shift =  (np.array(im1_shift))
#im_cells_shift = np.array(im_cells_shift)
#mask_border_shift =  (np.array(mask_border_shift)/255).astype(bool)
#mask_tractions_shift =  (np.array(mask_tractions_shift)/255).astype(bool)


## now full with diffrent incomplete shift values along one axis

pixelsize = 0.201
overlap = 16
ws = 20
ws_pix = ws/pixelsize
ov_pix = overlap/pixelsize
young = 49000
height = 300
sigma = 0.49

y_shifts = [55,58.0,58.03,59,63]
shifts = [(shift_value[0], y) for y in y_shifts]
## doesnt work yet
values = []
for ns in tqdm(shifts):
    # shift is applied to im1 and im2, im_cells, mask_border, mask_tractions are cut to the same shape without shift
    im1_shift, (im2_shift,im_cells_shift, mask_border_shift, mask_tractions_shift) = shift_arrays(im1, ns[0],ns[1], [im2, im_cells, mask_border, mask_tractions])
    mask_border_shift = mask_border_shift.astype(bool)
    mask_tractions_shift =  mask_tractions_shift.astype(bool)
    # deformation field
    u, v, mask, mask_std = calculate_deformation(im1.astype(np.int32), im2.astype(np.int32),
                                                 ws_pix, ov_pix,
                                                 std_factor=default_parameters["std_factor"])
    ps_new = pixelsize * np.mean(np.array(im1_shift.shape) / np.array(u.shape))
    # tractions
    mask_tractions_int = interpolation(mask_tractions_shift, u.shape).astype(bool)
    mask_cell_borders_int = interpolation(mask_border_shift, v.shape).astype(bool)

    tx, ty = TFM_tractions(u, v, pixelsize1=pixelsize, pixelsize2=ps_new, h=h, young=young,
                           sigma=sigma, filter=default_parameters["filter_type"], fs=default_parameters["filter_size"])
    # FEM
    fx = tx * (ps_new * 10**-6)**2
    fy = ty * (ps_new * 10**-6)**2
    fx, fy, p = correct_forces(fx, fy, mask_tractions_int)
    UG_sol, stress_tensor = stress_wrapper(mask_tractions_int, fx, fy, young, sigma=0.5,verbose=False)

    # gathering results

    ep = strain_energy_points(u, v, fx, fy, pixelsize, ps_new)  # contractile energy at any point
    energy = np.nansum(ep[mask_tractions_int])
    contractile_force, proj_x, proj_y, center = contractillity(tx, ty, ps_new, mask_tractions_int)
    sigma_max, sigma_min, sigma_max_abs, tau_max, phi_n, phi_shear, sigma_mean = all_stress_measures(stress_tensor,
                                                                                                     px_size=ps_new * 10 ** -6)
    avg_sig_max = np.nanmean(sigma_max[mask_cell_borders_int])
    avg_sig_shear = np.nanmean(tau_max[mask_cell_borders_int])
    avg_sig_mean = np.nanmean(sigma_mean[mask_cell_borders_int])

    values.append([energy, contractile_force, avg_sig_max, avg_sig_mean])

values = list(zip(*values))

plt.figure()
for v in values:
    plt.plot(y_shifts, v)