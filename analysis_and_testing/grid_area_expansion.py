# this script calculates the contracility and avg. mena normal stress for diffrent FEM-grid areas
# for a real data example and the standart artififciall sytem

#

import sys
sys.path.append("/home/user/Software/pyTFM/analysis_and_testing")
from evaluation_functions import *
from plotting_evaluation import *
from pyTFM.plotting import show_quiver
import os






def exp_border_real_data(out_folder, u, v,tx, ty, mask, border_ex_range, pixelsize):
    # TODO: unit seems to be inconsistent
    border_ex_range =  sorted(border_ex_range)
    exp_method = "binary_dilation"
    # convert traction to forces
    pixelsize1 = pixelsize * 10**-6
    pixelsize2 = pixelsize1 * np.mean(np.array(mask.shape) / np.array(tx.shape))
    fx, fy = tx * (pixelsize2 ** 2), ty * (pixelsize2 ** 2)

    np.save(os.path.join(out_folder, "fx.npy"), fx)
    np.save(os.path.join(out_folder, "fy.npy"), fy)
    np.save(os.path.join(out_folder, "u.npy"), u)
    np.save(os.path.join(out_folder, "v.npy"), v)
    # interpolation of the mask
    mask_int = interpolation(mask, dims=fx.shape)
    # border expansion test
    stress_tensors, mean_normal_list, mask_exp_list = exp_border(exp_range=border_ex_range, fx=fx, fy=fy,
                                mask=mask_int, out_folder=out_folder, method=exp_method, pixelsize=pixelsize2, verbose=False)
    # saving the a bit of output
    avg_normal_stresses = [np.nanmean(ms[mask_int.astype(bool)]) for ms in mean_normal_list]
    exp_mu = convert_exp_range(border_ex_range, exp_method, pixelsize2 / (10**-6))
    save_arr = np.array([np.round(np.array(avg_normal_stresses), 5), np.array(border_ex_range), exp_mu]).T
    h = 300
    header = "pixelsize %.4f h %d exp_method %s" % (pixelsize, h, exp_method)
    np.savetxt(os.path.join(out_folder, "avg_norm_stress.txt"), save_arr, fmt="%.7f", delimiter=",", header=header)


def exp_border_artificial_data(out_folder, border_ex_range, pixelsize=1, young=1, h=100):

    border_ex_range =  sorted(border_ex_range)
    # defining the size, shape and
    im_shape = (400, 400)  # (400, 400)
    stress_filed_size = border_ex_range[0]  # 150
    stress_filed_shape = "rectangle"
    stress_field_distribution = "uniform"
    sigma_n = 1
    sigma_shear = 0
    filter = "gaussian"
    exp_method = "manual" # defines weather grid is drawn new (yields rectangular shape, or whether grid is
    # expanded with binary expansion (yields form with 8 edges, stop sign like

    # inintial mask defining the "cell area"
    mask = setup_geometry(im_shape=im_shape, shape_size=stress_filed_size, shape=stress_filed_shape)
    # corresponding stress field
    stress_tensor_b = setup_stress_field(mask, distribution=stress_field_distribution, sigma_n=sigma_n,
                                         sigma_shear=sigma_shear)
    # forces from derivation of the stress field
    fx_b, fy_b = force_from_stress_field(stress_tensor_b, pixelsize, plot=False, grad_type="diff1")
    # deformation field based on fourrier transformed, height corrected boussinesq solution
    u, v = fourrier_deformation(fx_b, fy_b, pixelsize, young, sigma=0.49, h=h)
    # reconstructed froces with height corrected TFM algorithm
    fx_f, fy_f = traction_wrapper(u, v, pixelsize, h, young, mask=mask, filter=filter, fs=2 * pixelsize)

    np.save(os.path.join(out_folder, "mask.npy"), mask)
    np.save(os.path.join(out_folder, "fx_b.npy"), fx_b)
    np.save(os.path.join(out_folder, "fy_b.npy"), fy_b)
    np.save(os.path.join(out_folder, "u.npy"), u)
    np.save(os.path.join(out_folder, "v.npy"), v)
    np.save(os.path.join(out_folder, "stress_tensor_b.npy"), stress_tensor_b)
    np.save(os.path.join(out_folder, "fx.npy"), fx_f)
    np.save(os.path.join(out_folder, "fy.npy"), fy_f)

    stress_tensors, mean_normal_list, mask_exp_list = exp_border(exp_range=border_ex_range,
                                                                 fx=fx_f, fy=fy_f, mask=mask,out_folder=out_folder,
                                                                 method=exp_method, pixelsize=pixelsize)

    avg_normal_stresses = [np.nanmean(ms[mask.astype(bool)]) for ms in mean_normal_list]
    exp_mu = convert_exp_range(border_ex_range, exp_method, pixelsize)
    save_arr = np.array([np.round(np.array(avg_normal_stresses), 5), np.array(border_ex_range), exp_mu]).T
    header = "pixelsize %.4f h %d exp_method %s" % (pixelsize, h, exp_method)
    np.savetxt(os.path.join(out_folder, "avg_norm_stress.txt"), save_arr, fmt="%.7f", delimiter=",", header=header)

def read_txt_file(file):
    header = read_header(file)
    res = np.loadtxt(file, delimiter=",", skiprows=0)
    avg_normal_stresses = res[:, 0]
    border_ex_range = res[:, 1].astype(int)
    border_ex_range_mu = res[:, 2]
    return header, avg_normal_stresses, border_ex_range, border_ex_range_mu


def evaluate_sizes(in_folder):
    max_values = []
    max_indices = []
    zero_exp_values = []
    sizes = []
    max_contractility_recov = []
    for folder in os.listdir(in_folder):
        sub_folder = os.path.join(in_folder, folder)
        file = os.path.join(sub_folder, "avg_norm_stress.txt")
        if os.path.exists(file):
            size = int(folder)
            u_s, v_s, fx_s, fy_s, stress_tensors_s, avg_normal_stresses_s, mask_exp_list_s, header_s \
                , border_ex_range_s, border_ex_range_mu_s = load_exp_data(sub_folder)
            fx_b, fy_b = np.load(os.path.join(sub_folder,"fx_b.npy")), np.load(os.path.join(sub_folder,"fy_b.npy"))
            contractilities_b, strain_energies_b = contractility_strain_energy_exp(u_s, v_s, fx_b, fy_b,
                                                                                   mask_exp_list_s, 1, 1)
            contractilities_s, strain_energies_s = contractility_strain_energy_exp(u_s, v_s, fx_s, fy_s,
                                                                                   mask_exp_list_s, 1, 1)
            cont_recov = np.array(contractilities_s) / np.array(contractilities_b)
            # plt.ioff()
            # general_display(plot_types=["exp_test2_single"], pixelsize=header_s["pixelsize"],
            #                 out_folder=sub_folder, cmap="coolwarm",
            #                 avm_list=avg_normal_stresses_s, exp_list_mu=border_ex_range_mu_s,
            #                 contractilities=contractilities_s,
            #                 strain_energies=strain_energies_s,
            #                 plot_gt_exp=False, at=False)
            # plt.ion()
            #plt.close("all")
            zero_exp_values.append(avg_normal_stresses_s[0])
            max_values.append(np.max(avg_normal_stresses_s))
            max_indices.append(border_ex_range_mu_s[np.argmax(avg_normal_stresses_s)])
            sizes.append(size)
            max_contractility_recov.append(np.max(cont_recov))

    sort_indices = np.argsort(sizes)
    sizes = np.array(sizes)[sort_indices]
    max_values = np.array(max_values)[sort_indices]
    max_indices = np.array(max_indices)[sort_indices]
    zero_exp_values = np.array(zero_exp_values)[sort_indices]
    max_contractility_recov =  np.array(max_contractility_recov)[sort_indices]

    out_arr = np.array([zero_exp_values, max_values, max_indices, sizes]).T
    header = "avg_normal_stress_no_expansion[N/m] avg_normal_stress_best_expansion[N/m] size_best_expansion[µm] object_size[µm]"
    np.savetxt(os.path.join(out_folder_art_sizes,"results.txt"), out_arr, header=header)

    general_display(plot_types=["exp_size_stress", "exp_size_index", "exp_size_force"],
                    out_folder=out_folder_art_sizes, cmap="coolwarm",
                    max_stresses=max_values, sizes=sizes, max_index=max_indices, zero_exp_values=zero_exp_values, cont_rec=max_contractility_recov,
                    plot_gt_exp=False, at=False)


    general_display(plot_types=["exp_size_stress", "exp_size_index", "exp_size_force"],
                        out_folder=out_folder_art_sizes, cmap="coolwarm",
                        max_stresses=max_values, sizes=sizes, max_index=max_indices, zero_exp_values=zero_exp_values, cont_rec=max_contractility_recov,
                        plot_gt_exp=False, at=False, ext=".png")

    #plt.plot(sizes, max_values)
    #plt.plot(sizes, zero_exp_values)

def load_exp_data(folder):

    # loading the deformation and force field field
    u = try_np_load(os.path.join(folder, "u.npy"))
    v = try_np_load(os.path.join(folder, "v.npy"))
    fx = try_np_load(os.path.join(folder, "fx.npy"))
    fy = try_np_load(os.path.join(folder, "fy.npy"))

    # loading list of expansion values and avg. normal stresses
    header, avg_normal_stresses, border_ex_range, border_ex_range_mu = read_txt_file(os.path.join(folder, "avg_norm_stress.txt"))

    # loading the masks and stress tensors
    tensor_folder = os.path.join(folder, "stress_tensors")
    mask_folder = os.path.join(folder, "masks")
    stress_tensors = []
    mask_exp_list = []
    for i in tqdm(border_ex_range):
        stress_tensor_f = np.load(os.path.join(tensor_folder, "stress_tensor_f%d.npy" % i))
        mask_exp = np.load(os.path.join(mask_folder, "%d.npy" % i))
        stress_tensors.append(stress_tensor_f)
        mask_exp_list.append(mask_exp)

    return u, v, fx, fy, stress_tensors, avg_normal_stresses, mask_exp_list, header, border_ex_range, border_ex_range_mu


def fill_fields(mask, mask_exp, fx, fy):
    fields = fields_empty.copy()
    fields["mask"] = mask
    fields["mask_exp"] =  mask_exp
    fields["fx_f_exp"] = fx
    fields["fy_f_exp"] = fy
    return fields


new_analysis_artificial = False
new_analysis_real_data = False
new_analysis_artificial_size = False
out_folder_art =  "/home/user/Desktop/paper_pyTFM/data_traction_force_microscopy/grid_area_expansion_artificial_data"
out_folder_art_sizes =  "/home/user/Desktop/paper_pyTFM/data_traction_force_microscopy/grid_area_expansion_artificial_data/size"
out_folder_rl = "/home/user/Desktop/paper_pyTFM/data_traction_force_microscopy/grid_area_expansion_real_data"
in_folder = "/home/user/Desktop/paper_pyTFM/data_traction_force_microscopy/data_KO_03/"


# todo: repeat this analysis
pixel_size_rl = 0.201
mask_traction, mask_cell_border, fx, fy, tx, ty, u, v = load_data(in_folder, pixel_size_rl)
if new_analysis_real_data:
    # expansion with real data real data
    os.makedirs(out_folder_rl, exist_ok=True)
    border_ex_range = list(range(0, 150, 4))
    #border_ex_range = list(range(0, 2, 2))
    exp_border_real_data(out_folder_rl, u, v,tx, ty, mask_cell_border, border_ex_range, pixel_size_rl)


if new_analysis_artificial_size:
    # TODO: stopped at 120!! 08.12.2020
    for ex_start in tqdm([2,5,10,15, 20,25, 30, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150]):
        out_folder = os.path.join(out_folder_art_sizes, str(ex_start))
        os.makedirs(out_folder, exist_ok=True)
        border_ex_range = list(range(ex_start, ex_start + 8, 1)) + list(range(ex_start+8, ex_start + 30, 4))
        exp_border_artificial_data(out_folder, border_ex_range)

## new calculations
if new_analysis_artificial:
    os.makedirs(out_folder_art, exist_ok=True)
    border_ex_range = sorted(np.unique(np.array(list(range(150, 260, 4)) + list(range(180, 190, 1)))))
    #border_ex_range = [50,60]
    exp_border_artificial_data(out_folder_art, border_ex_range)



#die

### evaluation of diffrent object sizes
evaluate_sizes(out_folder_art_sizes)
## this looks pretty good..-> obviously wont work that to well for single cells
#key notes: no added noise on piv field wich can widely influence the result
# 20 µm seems to very small already
# this is worst case, where forces are gnerated strictly on the cell edge--> really may be much better
# wich could easily be remidied with higher resolution PIV and TFM



## plotting and stuff
u_art, v_art, fx_art, fy_art, stress_tensors_art, avg_normal_stresses_art, mask_exp_list_art, header_art\
    , border_ex_range_art, border_ex_range_mu_art= load_exp_data(out_folder_art)
u_rl, v_rl, fx_rl, fy_rl, stress_tensors_rl, avg_normal_stresses_rl, mask_exp_list_rl, header_rl\
    , border_ex_range_rl, border_ex_range_mu_rl= load_exp_data(out_folder_rl)


contractilities_art, strain_energies_art = contractility_strain_energy_exp(u_art,v_art,fx_art,fy_art, mask_exp_list_art, 1,1)
pixelsize2 = pixel_size_rl * np.mean(np.array(mask_cell_border.shape) / np.array(fx_rl.shape))
contractilities_rl, strain_energies_rl =  contractility_strain_energy_exp(u_rl, v_rl, fx_rl, fy_rl, mask_exp_list_rl, pixel_size_rl, pixelsize2)

contractilities_art_norm = contractilities_art/contractilities_art.max()
contractilities_rl_norm = contractilities_rl/contractilities_rl.max()
avg_normal_stresses_art_norm = avg_normal_stresses_art/1
avg_normal_stresses_rl_norm = avg_normal_stresses_rl/avg_normal_stresses_rl.max()




out_array_art = np.array([contractilities_art, contractilities_art_norm,  avg_normal_stresses_art, avg_normal_stresses_art_norm ,  border_ex_range_mu_art]).T
out_array_rl = np.array([contractilities_rl, contractilities_rl_norm,  avg_normal_stresses_rl, avg_normal_stresses_rl_norm,  border_ex_range_mu_rl]).T
header = "contracility[N] contractility_normalized avg_normal_stress[N/m] acg_normal_stress_normalized expansion[µm]"

np.savetxt(os.path.join(out_folder_art,"plot_numbers.txt"), X=out_array_art, header=header)
np.savetxt(os.path.join(out_folder_rl,"plot_numbers.txt"), X=out_array_rl, header=header)





## maybe read that from the load_exp data // also immediate conversion to µm in mask_exp_list
f_type = "circular"

display_bd = None
end_shape_rl = (380,390) # shape of real data/ and that looks nice in artfiical system, need to adjust for every experiment
end_shape_art = (260,260) # for artificial system: (260,260) , for real data: (390,400)
end_shape_ex = (380,390)  #shape for the border exapnsion plots

disp_index = -1


fields_art = fill_fields(mask_exp_list_art[0], mask_exp_list_art[0], fx_art, fy_art)
fields_rl = fill_fields(mask_cell_border, mask_exp_list_rl[0], fx_rl, fy_rl)
max_dict = {}
max_dict["stress"] = np.max(avg_normal_stresses_rl)
plot_types = ["border_ex2", "border_ex5"]
general_display(plot_types=plot_types, pixelsize=header_art["pixelsize"], max_dict=max_dict, f_type=f_type,
                out_folder=out_folder_art, cmap="coolwarm",
                avm_list_rl=avg_normal_stresses_rl, contractilities_rl=contractilities_rl, strain_energies_rl=strain_energies_rl, exp_list_mu_rl=border_ex_range_mu_rl,
                avm_list_art=avg_normal_stresses_art,contractilities_art=contractilities_art,strain_energies_art=strain_energies_art, exp_list_mu_art=border_ex_range_mu_art,
                fields_rl=fields_rl,fields_art=fields_art,border_ex_test_rl=border_ex_range_rl,border_ex_test_art=border_ex_range_art, exp_masks_rl=mask_exp_list_rl, exp_masks_art=mask_exp_list_art,
                fields=None, key_values=None, plot_gt_exp=False, dm=True, at=False,
                    cb=False)




plt.close("all")