# simualtes the the foramtion of two forces, and predicst the force from the deformations with TFM
# plot the output strain energy vs the input strain energy for spread out forces



from scipy.ndimage import binary_erosion
import sys
sys.path.insert(0,'/home/user/Software/pyTFM/analysis_and_testing/')
from simulating_deformation import *
from playing_with_strains import *
from plotting_evaluation import *
from general_evaluation import *
from pyTFM.TFM_functions import strain_energy_points
from pyTFM.utilities_TFM import createFolder
import time
import os
import matplotlib.pyplot as plt
import numpy as np
from skimage.filters import gaussian
from skimage.draw import circle
from tqdm import tqdm


figsize = (7, 7)
arrow_scale = 0.13  # 0.1
arrow_width = 0.004
headlength = 4
headwidth = 4  # 6
headaxislength = 3
paras = {"figsize": figsize,
      "arrow_scale": arrow_scale,
      "arrow_width": arrow_width,
      "headlength": headlength,
      "headwidth": headwidth,
      "headaxislength": headaxislength,
         "cmap":"coolwarm",
         "vmax":None, "filter":[0,6]
         }


def show_forces_forward(x=None, y=None, figsize=None, arrow_scale=None, arrow_width=None, filter=None
                        , headlength=None, headwidth=None, headaxislength=None, cmap=None, vmax=None, cb=None):
    fig, ax = show_quiver(x, y, figsize=figsize, scale_ratio=arrow_scale - 0.03, filter=filter,
                          width=arrow_width,
                          headlength=headlength, headwidth=headwidth,
                          headaxislength=headaxislength, cbar_tick_label_size=30, cmap=cmap,
                          cbar_style="not-clickpoints",
                          vmin=0, vmax=vmax, plot_cbar=cb)
    return fig,ax


pixelsize=1
h=1000
young=1
new_calculation = False

output_folder = "/home/user/Desktop/backup_from_harddrive/data_traction_force_microscopy/gaussfilter_strain_energy"
createFolder(output_folder)
#fx_b, fy_b = np.zeros((200,200)), np.zeros((200,200))
#l1 = np.array([np.arange(35,45),np.arange(35,45)])
#l2 = np.array([np.arange(55,65),np.arange(55,65)])



out_folder2=createFolder(os.path.join(output_folder,"force_plots"))
cont_input = []
cont_output_filtered = []
cont_output_unfiltered = []
rs=list(range(1,12,2))
ar_size=600
mask=np.ones((ar_size,ar_size))
mask = binary_erosion(np.ones((ar_size, ar_size)), iterations=20).astype(bool)

filter_size=3
plt.ioff()

name_add = ""
plotting=True
rs = list(range(1,18,2))
if new_calculation:
    for r in tqdm(rs):
        fx_b, fy_b = np.zeros((ar_size, ar_size)), np.zeros((ar_size, ar_size))
        c1 = circle(r=ar_size / 2, c=200, radius=r)
        fx_b[c1[0], c1[1]] = -1
        c2 = circle(r=ar_size / 2, c=400, radius=r)
        fx_b[c2[0], c2[1]] = 1
        fx_b = fx_b / np.sum(np.abs(fx_b))

        u_b, v_b = deformation_by_upsampling(fx_b, fy_b, factor=10, kernel_size=(3000, 3000),
                                             method="fftconvolve")  # 3000x3000 is maximum
        fig, a = show_quiver(u_b, v_b)
        fig.savefig(os.path.join(out_folder2, str(r) + "def.svg"))
        fx_f_nf, fy_f_nf = traction_wrapper(u_b, v_b, pixelsize, h, young, mask=mask, sigma=0.5,
                                            filter=None, fs=None)  # assuming pixelsize == 1
        # show_quiver(fx_f_nf, fy_f_nf)
        fx_f, fy_f = traction_wrapper(u_b, v_b, pixelsize, h, young, mask=mask, sigma=0.5,
                                      filter="gaussian", fs=filter_size)  # assuming pixelsize == 1
        if plotting:
            paras["filter"] = [0, 1]
            f, a = show_forces_forward(x=fx_b, y=fy_b, **paras)
            f.savefig(os.path.join(out_folder2, str(r) + "input.svg"))
            paras["filter"] = [0, 1] if r < 15 else [0, 4]
            f, a = show_forces_forward(x=fx_f_nf, y=fy_f_nf, **paras)
            f.savefig(os.path.join(out_folder2, str(r) + "output_unfiiltered.svg"))
            f, a = show_forces_forward(x=fx_f, y=fy_f, **paras)
            f.savefig(os.path.join(out_folder2, str(r) + "output_filtered.svg"))
        plt.close("all")
        energy_points = strain_energy_points(u_b, v_b, fx_b, fy_b, 1, 1)  # contractile energy at any point
        cont_input.append(np.sum(energy_points[mask]))
        energy_points = strain_energy_points(u_b, v_b, fx_f_nf, fy_f_nf, 1, 1)  # contractile energy at any point
        cont_output_unfiltered.append(np.sum(energy_points[mask]))
        energy_points = strain_energy_points(u_b, v_b, fx_f, fy_f, 1, 1)  # contractile energy at any point
        cont_output_filtered.append(np.sum(energy_points[mask]))
    out_arr = np.stack([np.array(cont_input), np.array(cont_output_unfiltered), np.array(cont_output_filtered)]).T
    header = "input,unfiltered,filtered"
    np.savetxt(os.path.join(output_folder, "out.txt"), out_arr, header=header, delimiter=",")

la = np.loadtxt(os.path.join(output_folder, "out.txt"),delimiter=",")
cont_input = la[:,0]
cont_output_unfiltered = la[:,1]
cont_output_filtered = la[:,2]

cont_input_norm=np.array(cont_input)/np.array(cont_input)
cont_out_nf_norm=np.array(cont_output_unfiltered)/np.array(cont_input)
cont_out_f_norm=np.array(cont_output_filtered)/np.array(cont_input)


plt.ion()
fig=plt.figure()
plt.plot(rs,cont_input_norm,label="strain energy of input forces", color="C1", lw=3)
plt.plot(rs,cont_out_f_norm,label="strain energy of output forces", color="C2", lw=3)
plt.gca().tick_params(axis="both", which="both", color="black", length=4, width=2, labelsize=20, labelcolor="black")
plt.legend()
set_axis_attribute(plt.gca(), "set_color", "black")
set_axis_attribute(plt.gca(), "set_linewidth", 2)

fig.savefig(os.path.join(output_folder,name_add+"obj_size_dependance1.svg"))



fig=plt.figure()
plt.plot(rs,cont_input_norm,label="strain energy of input forces", color="C1", lw=3)
plt.plot(rs,cont_out_nf_norm,label="strain energy of output forces unfiltered", color="C2", lw=3)
plt.plot(rs,cont_out_f_norm,label="strain energy of output forces filtered", color="C3", lw=3)
plt.legend()
plt.gca().tick_params(axis="both", which="both", color="black", length=4, width=2, labelsize=20, labelcolor="black")
set_axis_attribute(plt.gca(), "set_color", "black")
set_axis_attribute(plt.gca(), "set_linewidth", 2)
fig.savefig(os.path.join(output_folder,name_add+"obj_size_dependance2.svg"))

