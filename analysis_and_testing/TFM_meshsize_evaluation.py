# this script evaluates the effect of diffrent gridsize ( size of an individual FEA element) on the stress calcualtion
# a high resolution traction field is downsized (by 3rd order spline interpolation) and used to calculate stress tensors

import numpy as np
import matplotlib.pyplot as plt
from pyTFM.TFM_functions import calculate_deformation
from pyTFM.plotting import show_quiver
from PIL import Image#
from tqdm import tqdm
import os


from scipy.ndimage import binary_fill_holes, zoom
from pyTFM.TFM_functions import calculate_deformation, TFM_tractions, strain_energy_points, contractillity
from pyTFM.grid_setup_solids_py import interpolation, grid_setup,FEM_simulation
from pyTFM.stress_functions import all_stress_measures
import re
import sys

sys.path.insert(0, "/home/user/Software/pyTFM/analysis_and_testing")
from evaluation_functions import *
from plotting_evaluation import *

def stress_wrapper(mask, fx, fy, young=1, sigma=0.5, verbose=False):
    nodes, elements, loads, mats = grid_setup(mask, -fx, -fy, young, sigma=sigma, edge_factor=0)  # construct fem grid
    # solve FEM system
    UG_sol, stress_tensor = FEM_simulation(nodes, elements, loads, mats, mask, verbose=verbose,
                                           system_type="colony")  # colony means that the pure neumann FEM is applied
    return UG_sol, stress_tensor

def zoom_fields(tx,ty,u,v, l,pixelsize,mask,order=3):
    zl = init_length / l
    tx_int = zoom(tx, zl, order=order)
    ty_int = zoom(ty, zl, order=order)

    # check if correct interpolation range
    check = pixel_size * np.mean(np.array(mask.shape) / np.array(tx_int.shape))
    print("target element size: ", l, "actual element size", check)

    # check how much the interpolation influences the strain energy and the contractility
    pixelsize1 = pixelsize
    pixelsize2 = pixel_size * np.mean(np.array(mask.shape) / np.array(tx_int.shape))
    area2 = (pixelsize2 * 10 ** - 6) ** 2
    u_int = zoom(u, zl, order=order)
    v_int = zoom(v, zl, order=order)

    mask_interp = interpolation(mask, ty_int.shape)

    return tx_int, ty_int, u_int, v_int, mask_interp, pixelsize1, pixelsize2, area2


def try_grid_sizes(lengths, tx, ty, u, v,  out_folder, pixelsize, verbose=False):
    strain_energies = []
    contractile_forces = []
    order = 3  # doesnt seem to matter for fluctuations of contractillity and strain energy
    for l in tqdm(lengths):
        tx_int, ty_int, u_int, v_int, mask_interp, pixelsize1, pixelsize2, area2= zoom_fields(tx, ty, u, v, l, pixelsize, mask,order=3)

        ep = strain_energy_points(u_int, v_int, tx_int, ty_int, pixelsize1, pixelsize2)
        cenergy = np.sum(ep[mask_interp])

        contractile_force, proj_x, proj_y, center = contractillity(tx_int, ty_int, pixelsize2, mask_interp)

        strain_energies.append(cenergy)
        contractile_forces.append(contractile_force)

        # calcualte stresses with the interpolated force field
        UG_sol, stress_tensor = stress_wrapper(mask_interp, tx_int * area2, ty_int * area2, young=1, sigma=0.5,
                                               verbose=verbose)
        l_label = str(np.round(l, 2)).replace(".", "_")
        np.save(os.path.join(out_folder, "stress_tensor%s.npy" % l_label), stress_tensor)

        # monitor contractility /strain energy change

    strain_energies = np.array(strain_energies)
    contractilities = np.array(contractile_forces)
    out = np.array([lengths, strain_energies, contractilities]).T
    np.savetxt(os.path.join(out_folder, "res.txt"), out, fmt="%.4e", delimiter=' ',
               header="grid_element_length strain_energy contractility")


def read_lenghs_from_file_names(out_folder):
    names = os.listdir(out_folder)
    numbers1 = [re.search("(\d.*\d)", x).group(1) for x in names if re.search("(\d.*\d)", x)]

    numbers2 = [re.search("[a-zA-Z](\d)\.", x).group(1) for x in names if re.search("[a-zA-Z](\d)\.", x)]
    numbers = numbers1 + numbers2
    numbers = [float(n.replace("_", ".")) for n in numbers]
    return sorted(numbers)

def load_grid_size_test(out_folder, in_folder):
    #res = np.loadtxt(os.path.join(out_folder, "res.txt"), delimiter=' ', skiprows=1)
    mask, mask_cell_border, fx, fy, tx, ty, u, v = load_data(in_folder, pixel_size)



    contractilities = []
    strain_energies = []
    lengths = read_lenghs_from_file_names(out_folder)
    stress_tensors = []

    for l in tqdm(lengths):

        tx_int, ty_int, u_int, v_int, mask_interp, pixelsize1, pixelsize2, area2 = zoom_fields(tx, ty, u, v, l, pixel_size,
                                                                                        mask,
                                                                                        order=3)
        ep = strain_energy_points(u_int, v_int, tx_int, ty_int, pixelsize1, pixelsize2)
        cenergy = np.sum(ep[mask_interp])
        strain_energies.append(cenergy)
        contractile_force, proj_x, proj_y, center = contractillity(tx_int, ty_int, pixelsize2, mask_interp)
        contractilities.append(contractile_force)

        l = int(l) if int(l) ==  l else l
        l_label = str(np.round(l, 2)).replace(".", "_")
        print(l_label)
        stress_tensors.append(np.load(os.path.join(out_folder, "stress_tensor%s.npy" % l_label)))
    return lengths, stress_tensors, strain_energies ,contractilities


pixel_size = 0.201
new_analysis = False
# loading the traction field

in_folder = "/home/user/Desktop/paper_pyTFM/data_traction_force_microscopy/data_KO_03/"
out_folder = "/home/user/Desktop/paper_pyTFM/data_traction_force_microscopy/meshsize_evaluation/"

mask, mask_cell_border, fx, fy ,tx, ty, u, v = load_data(in_folder, pixel_size)
os.makedirs(out_folder,exist_ok=True)


init_length = pixel_size * np.mean(np.array(mask.shape)/np.array(u.shape))
mask_int = interpolation(mask,u.shape)#initial  grid length
#tx, ty = TFM_tractions(u, v, pixelsize1=pixel_size, pixelsize2=init_length, h=300,
#                           young=18500, sigma=0.49, filter="gaussian", fs=3)
#np.save("/home/user/Desktop/paper_pyTFM/data_traction_force_microscopy/data_KO_03/fy_large", ty*init_length**2)
lengths = [init_length] + list(np.arange(1,20)) +[1.3, 1.7, 2.5]# [init_length] + [0.4, 0.6, 0.8] + list(np.arange(1,20)) +[1.3, 1.7, 2.5]
lengths = sorted(lengths)

out_folder = "/home/user/Desktop/paper_pyTFM/data_traction_force_microscopy/meshsize_evaluation/"
os.makedirs(out_folder,exist_ok=True)
# iterate over FEA element size

if new_analysis:
    try_grid_sizes(lengths, tx, ty, u, v, out_folder, pixel_size, verbose=False)
lengths, stress_tensors, strain_energies ,contractilities = load_grid_size_test(out_folder, in_folder)

sig_means, sig_maxs, sig_shears = [], [], []
pixel_sizes = []
for st in stress_tensors:
    pixel_size2 = pixel_size * np.mean(np.array(mask.shape) / np.array(st.shape[:2]))
    pixel_sizes.append(pixel_size2)
    mask_interp = interpolation(mask,st.shape[:2])
    sigma_max, sigma_min, sigma_max_abs, tau_max, phi_n, phi_shear, sigma_mean = all_stress_measures(st, px_size=pixel_size2 * 10 ** -6)
    sig_means.append(np.nanmean(sigma_mean[mask_interp]))
    sig_maxs.append(np.nanmean(sigma_max[mask_interp]))
    sig_shears.append(np.nanmean(tau_max[mask_interp]))



# strain energies and contractilities to confirm that interpolation doesnt introduce large errors
plt.figure()
plt.plot(lengths, strain_energies, "o-", label="strain energy")
plt.ylim(0, np.max(strain_energies) * 1.2)
plt.legend()
plt.xlabel("grid element edge length [µm]")
plt.savefig(os.path.join(out_folder, "strain_energies.png"))
plt.figure()
plt.plot(lengths,contractilities, "o-",  label="contractillity")
plt.ylim(0, np.max(contractilities) * 1.2)
plt.legend()
plt.xlabel("grid element edge length [µm]")
plt.savefig(os.path.join(out_folder, "contractilities.png"))

# stresses
plt.figure()
plt.plot(lengths, sig_means, "o-", label="avg. mean stress", lw=3, ms=7)
plt.plot(lengths, sig_maxs, "o-", label="avg. max stress", lw=3, ms=7)
plt.plot(lengths, sig_shears, "o-", label="avg. shear stress", lw=3, ms=7)
plt.xlabel("grid element edge length [µm]")
plt.ylabel("stresses [N/m]")


plt.gca().tick_params(axis="both", which="both", color="black", length=4, width=2, labelsize=20,
                      labelcolor="black")
set_axis_attribute(plt.gca(), "set_color", "black")
set_axis_attribute(plt.gca(), "set_linewidth", 2)
plt.tight_layout()
plt.savefig(os.path.join(out_folder, "streses_nol.png"))
plt.legend()
plt.savefig(os.path.join(out_folder, "streses.png"))
plt.savefig(os.path.join(out_folder, "streses.svg"))


# stresses
plt.figure()
plt.plot(lengths, sig_means, "o-", label="avg. mean stress")
plt.plot(lengths, sig_maxs, "o-", label="avg. max stress")
plt.plot(lengths, sig_shears, "o-", label="avg. shear stress")
plt.xlabel("grid element edge length [µm]")
plt.ylabel("stresses [N/m]")
plt.legend()
plt.savefig(os.path.join(out_folder, "streses.png"))



# TODO: could think about including the line tension
# TODO: use smaller grid sizes ....
# why is the contractility and stuff decreasing?
0.201 * (20-19)