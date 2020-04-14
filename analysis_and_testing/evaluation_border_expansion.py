import sys
sys.path.append("/home/user/Software/pyTFM/analysis_and_testing")
from general_evaluation import *
from plotting_evaluation import *


out_folder = "/home/user/Desktop/backup_from_harddrive/data_traction_force_microscopy/ev_paper_expansion4"

young = 1
h = 100
pixelsize = 1
border_ex_test = np.loadtxt(os.path.join(out_folder,"avg_norm_stress.txt"), delimiter=",")[:,1].astype(int)

filter = "gaussian"
exp_method = "manual"
f_type = "non-circular"  # display_option


stress_tensors, mean_normal_list, mask_exp_list = load_exp_border(out_folder=out_folder, exp_range=border_ex_test)
try:
    mask = np.load(os.path.join(out_folder, "mask.npy"))
except FileNotFoundError:
    mask = mask_exp_list[0]

stress_tensor_b = setup_stress_field(mask, distribution="uniform", sigma_n=1, sigma_shear=0,
                                      sigma_gf=6)
fx_b, fy_b = traction_from_stress(stress_tensor_b, pixelsize, plot=False, grad_type="diff1")
try:
    u_b, v_b = np.load(os.path.join(out_folder, "u_b.npy")), np.load(os.path.join(out_folder, "v_b.npy"))
except FileNotFoundError:
    try:
        u_b, v_b = np.load(os.path.join(out_folder, "u.npy")), np.load(os.path.join(out_folder, "v.npy"))
    except FileNotFoundError:
        u_b, v_b = None, None

try:
    fx_f, fy_f = traction_wrapper(u_b, v_b, pixelsize, h, young, mask=mask,
                                  filter=filter)  # assuming pixelsize == 1
except:
    fx_f, fy_f =  np.load(os.path.join(out_folder, "fx_f.npy")), np.load(os.path.join(out_folder, "fy_f.npy"))

try:
    stress_tensor_f = np.load(os.path.join(out_folder, "stress_tensor_f.npy"))
except FileNotFoundError:
    stress_tensor_f = None




avg_normal_stress_be = [np.nanmean(ms[mask.astype(bool)]) for ms in mean_normal_list]
avg_normal_stress_be = np.array(avg_normal_stress_be)/np.max(avg_normal_stress_be) # normalizing
border_ex_test=border_ex_test-np.min(border_ex_test) # setting to zero (only usefull for these sepecific cases)

# getting comparison scalar fields

# cutting to ignore effects close to image edge #### mention this when talking to benn
fx_f, fy_f, stress_tensor_b, stress_tensor_f,\
                 mask = cut_arrays([20, -20, 20, -20], [fx_f, fy_f, stress_tensor_b,
                                stress_tensor_f, mask])
mask_exp_list = cut_arrays([20, -20, 20, -20], mask_exp_list)

max_force = np.max(np.sqrt(fx_f**2 + fy_f**2)) #0.06
max_stress = np.max(avg_normal_stress_be)
max_dict = {"stress": max_stress, "force": max_force*1.1}


plot_types = ["cbars_only"]
plot_types = ["be5", "be3", "be2"]
plot_types = ["be5"]
#plot_types = ["correlation", "cbars_only"]
# max_dict["force"] = 1
plt.ioff()
general_display(plot_types=plot_types, mask=mask, pixelsize=pixelsize, max_dict=max_dict, f_type=f_type,
                    out_folder=out_folder, cmap="coolwarm", fx_b=fx_b, fy_b=fy_b,
                    fx_f=fx_f, fy_f=fy_f,
                    border_ex_test=border_ex_test, be_avm_list=avg_normal_stress_be, mask_exp_list=mask_exp_list,
                    mean_normal_list=mean_normal_list, stress_tensor_f=stress_tensor_f,
                    stress_tensor_b=stress_tensor_b, plot_gt_exp=False, dm=False, at=False,
                    cb=False)
plt.ion()
plt.close("all")


fig = plt.figure(figsize=(3.2, 4.75))
plt.gca().set_axis_off()
cbar = add_colorbar(vmin=0, vmax=max_dict["force"]*0.1, aspect=8, shrink=1, cbar_axes_fraction=1.2, cmap="coolwarm",
                    cbar_style="not-clickpoints")
set_axis_attribute(cbar.ax, "set_color", "black")
cbar.ax.tick_params(axis="both", which="both", color="black", length=4, width=2, labelsize=20, labelcolor="black")
fig.savefig(os.path.join(out_folder, "cbars_only" + "force"  + ".svg"))