# this script checks the range of viable PIV window size parameters for a pair of images before and
# after substrate relaxation
## not sure If we want to include this --> opens up whole new discussion on the deformation field algorithms etc..
## but good to know
import numpy as np
import matplotlib.pyplot as plt
from pyTFM.TFM_functions import calculate_deformation, TFM_tractions, strain_energy_points, contractillity
from pyTFM.grid_setup_solids_py import interpolation
from pyTFM.plotting import show_quiver
from PIL import Image#
from tqdm import tqdm
import os
from scipy.ndimage import binary_fill_holes

def load_def(out_folder, ws_range):
    us, vs = [], []
    for ws in ws_range:
        us.append(np.load(os.path.join(out_folder,"%du.npy"%ws)))
        vs.append(np.load(os.path.join(out_folder,"%dv.npy"%ws)))
    return us, vs

#load images
 # µm/pixel
new_analysis =  False
out_folder = "/home/user/Desktop/paper_pyTFM/data_traction_force_microscopy/PIV_ws_check/"
pixel_size = 0.201
std_factor = 15
ws_range = np.arange(5,100,5)
overlaps = ws_range * 0.75
h = 300
young = 23500
sigma = 0.49
fs = 3

im1 = "/home/user/Desktop/paper_pyTFM/data_traction_force_microscopy/data_KO_03/03after_shift.tif"
im2 = "/home/user/Desktop/paper_pyTFM/data_traction_force_microscopy/data_KO_03/03before_shift.tif"
im_cells = "/home/user/Desktop/paper_pyTFM/data_traction_force_microscopy/data_KO_03/03bf_before_shift.tif"
im1 = np.array(Image.open(im1))
im2 = np.array(Image.open(im2))
im_cells = np.array(Image.open(im_cells))
im1 = im1.astype(np.int32)
im2 = im2.astype(np.int32)
mask = np.load("/home/user/Desktop/paper_pyTFM/data_traction_force_microscopy/data_KO_03/traction_mask.npy")
mask = binary_fill_holes(mask)


if new_analysis:
    os.makedirs(out_folder,exist_ok=True)
    for ws, ov in tqdm(zip(ws_range, overlaps)):
        window_size_pix = ws / pixel_size
        overlap_pix = 0.75 *  window_size_pix
        u, v, mask, mask_std = calculate_deformation(im1, im2, window_size=window_size_pix, overlap=overlap_pix, std_factor=std_factor)
        np.save(os.path.join(out_folder,"%du.npy"%ws), u)
        np.save(os.path.join(out_folder,"%dv.npy"%ws), v)
    show_quiver(u, v, filter=[0,15])

us, vs =  load_def(out_folder, ws_range)
strain_energies = []
contractile_forces = []
for u, v in zip(us,vs):
    #show_quiver(u, v, filter=[0, 7])
    pixelsize1 = pixel_size
    pixelsize2 =  pixel_size * np.mean(np.array(im1.shape) / np.array(u.shape))
    mask_interp = interpolation(mask, u.shape)
    tx, ty = TFM_tractions(u, v, pixelsize1=pixelsize1 * 10**-6, pixelsize2=pixelsize2 * 10**-6, h=h,
                           young=young, sigma=sigma, filter="gaussian", fs=fs * 10**-6)

    ep = strain_energy_points(u, v, tx, ty, pixelsize1, pixelsize2)
    cenergy = np.sum(ep[mask_interp])

    contractile_force, proj_x, proj_y, center = contractillity(tx, ty, pixelsize2, mask_interp)

    strain_energies.append(cenergy)
    contractile_forces.append(contractile_force)

    #show_quiver(tx,ty,filter=[0,3])



plt.figure()
plt.plot(ws_range,strain_energies)
plt.ylim(0, np.max(strain_energies) * 1.2)

