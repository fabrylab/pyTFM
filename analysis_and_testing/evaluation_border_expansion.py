import sys
import numpy as np
import os
import re
import matplotlib.pyplot as plt
sys.path.append("/home/user/Software/pyTFM/analysis_and_testing")
from evaluation_functions import *
from plotting_evaluation import *

stress_tensor_f, u_b, v_b = None, None, None


out_folder = "/home/user/Desktop/backup_from_harddrive/data_traction_force_microscopy/ev_paper_rd_expansion"


display_bd = None
end_shape = (380,390) # shape of real data/ and that looks nice in artfiical system, need to adjust for every experiment
#end_shape =  (260,260) # for artificial system: (260,260) , for real data: (390,400)
end_shape_ex = (380,390)  #shape fior the border exapnsion plots
border_ex_factor =1 # 0.5 for artificial system, 1 for real data; this factor compensates for the fact that realdata border exapnsion is done by binary delation
# and the artificial system by increasing the shape diameter. 1 dilation step adds 2 pixel to the diameter
young = 1
h = 100
pixelsize = 0.85 # 0.85 for "real data"
border_ex_test = np.loadtxt(os.path.join(out_folder,"avg_norm_stress.txt"), delimiter=",")[:,1].astype(int)
# read factor for filename conversion:
with open(os.path.join(out_folder,"avg_norm_stress.txt")) as f:
    l = f.readline()
    if "conversion_to_filenames=" in l:
        filename_factor = int(re.search("#conversion_to_filenames=(\d{1,4}).*", l).group(1))
    else:
        filename_factor = 1


filter = "gaussian"
exp_method = "manual"
f_type = "circular"  # display_option


stress_tensors, mean_normal_list, mask_exp_list = load_exp_border(out_folder=out_folder, exp_range=border_ex_test,filename_factor=filename_factor)
try:
    mask = np.load(os.path.join(out_folder, "mask.npy"))
except FileNotFoundError:
    mask = mask_exp_list[0]

stress_tensor_b = setup_stress_field(mask, distribution="uniform", sigma_n=1, sigma_shear=0)

fx_b, fy_b = force_from_stress_field(stress_tensor_b, pixelsize, plot=False, grad_type="diff1")
try:
    u_b, v_b = np.load(os.path.join(out_folder, "u_b.npy")), np.load(os.path.join(out_folder, "v_b.npy"))
except FileNotFoundError:
    try:
        u_b, v_b = np.load(os.path.join(out_folder, "u.npy")), np.load(os.path.join(out_folder, "v.npy"))
    except FileNotFoundError:
        pass

try:
    fx_f, fy_f = np.load(os.path.join(out_folder, "fx_f.npy")), np.load(os.path.join(out_folder, "fy_f.npy"))
except:

    fx_f, fy_f = traction_wrapper(u_b, v_b, pixelsize, h, young, mask=mask,
                                  filter=filter, fs=3 * pixelsize)

try:
    stress_tensor_f = np.load(os.path.join(out_folder, "stress_tensor_f.npy"))
except FileNotFoundError:
    pass

# calculating a list of contractilities and Strain energies on all FEM masks
contractilities, strain_energies = contractility_strain_energy_exp(u_b,v_b,fx_f,fy_f, mask_exp_list)
contractilities *= 10**12

avg_normal_stress_be = [np.nanmean(ms[mask.astype(bool)]) for ms in mean_normal_list]

border_ex_test = border_ex_test - np.min(border_ex_test) # setting to zero (only usefull for these sepecific cases)
border_ex_test = border_ex_test * border_ex_factor
# picking mask of fem grid to be displayed
if not isinstance(display_bd, int):
    display_id = np.argmax(avg_normal_stress_be)
else:
    display_id = np.argmin(np.abs(np.array(border_ex_test) - display_bd))
stress_tensor_f = stress_tensors[display_id]
mask_fem = mask_exp_list[display_id]


mask_fm = mask_fem #

fx_f_exp = cut_arrays(end_shape_ex, fx_f, mode="match")
fy_f_exp = cut_arrays(end_shape_ex, fy_f, mode="match")
mask_exp = cut_arrays(end_shape_ex, mask, mode="match")

fields = {"u_b":u_b, "v_b":v_b, "fx_b":fx_b, "fy_b":fy_b,"fx_f":fx_f, "fy_f":fy_f, "stress_tensor_b":stress_tensor_b, "stress_tensor_f":stress_tensor_f,
        "mask": mask, "mask_fem":mask_fem, "mask_fm": mask_fm, "fx_f_exp":fx_f_exp, "fy_f_exp":fy_f_exp,"mask_exp": mask_exp}#

measures = standard_measures(mask=mask, mean_normal_list=mean_normal_list, fields=fields)
measures["contractile_force_b"] *= 10**12
measures["contractile_force_f"] *= 10**12
scalar_comaprisons = full_field_comparision(fields)

# getting comparison scalar fields

# cutting to ignore effects close to image edge #### mention this when talking to benn

# used to be [20, -20, 20, -20] with mode="edge
fields = cut_arrays(end_shape, fields, mode="match", exclude =["fx_f_exp","fy_f_exp","mask_exp"])
mask_exp_list = cut_arrays(end_shape_ex, mask_exp_list, mode="match")



#max_force = np.max(np.sqrt(fx_f**2 + fy_f**2)) #0.06
max_force = 1
max_stress = np.max(avg_normal_stress_be)
max_dict = {"stress": max_stress, "force": max_force*1.1}


plot_types = ["cbars_only"]
plot_types = ["correlation", "key measures", "mean_normal_stress_backward", "mean_normal_stress_forward",
              "forces_forward", "forces_backward", "mask_outline", "cbars_only", "test for border expansion", "be5","be2", "be2A"]
#plot_types = ["key measures"]
#plot_types = ["be5", "be2A",  "be2"]
#plot_types = ["correlation", "cbars_only"]
# max_dict["force"] = 1
#plt.ioff()
plot_types = ["be5"]
general_display(plot_types=plot_types, pixelsize=pixelsize, max_dict=max_dict, f_type=f_type,
                out_folder=out_folder, cmap="coolwarm", scalar_comaprisons=scalar_comaprisons,
                border_ex_test=border_ex_test, be_avm_list=avg_normal_stress_be, mask_exp_list=mask_exp_list,
                mean_normal_list=mean_normal_list, strain_energies=strain_energies, contractilities=contractilities,
                fields=fields, key_values=measures, plot_gt_exp=False, dm=True, at=False,
                    cb=False)
#plt.ion()
plt.close("all")
'''
# plotting a single colorbar
fig = plt.figure(figsize=(3.2, 4.75))
plt.gca().set_axis_off()
cbar = add_colorbar(vmin=0, vmax=max_dict["force"]*0.1, aspect=8, shrink=1, cbar_axes_fraction=1.2, cmap="coolwarm",
                    cbar_style="not-clickpoints")
set_axis_attribute(cbar.ax, "set_color", "black")
cbar.ax.tick_params(axis="both", which="both", color="black", length=4, width=2, labelsize=20, labelcolor="black")
fig.savefig(os.path.join(out_folder, "cbars_only" + "force"  + ".svg"))




plt.close("all")
plt.show()


'''

