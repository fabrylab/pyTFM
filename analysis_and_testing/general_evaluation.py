import sys
sys.path.append("/home/user/Software/pyTFM/analysis_and_testing")
from evaluation_functions import *
from plotting_evaluation import *
from pyTFM.plotting import show_quiver
import os


#exp_border_real_data()
if __name__ == "__main__":
    new_calculation = False
    single_tensor = True  # set true if you dont want use only one FEM grid size
    FEM_grid_size_single = 185  # only if single_tensor
    FEM_grid_method_single = "manual"  # defines weather grid is drawn new (yields rectangular shape, or whether grid is
    # expanded with binary expansion (yields form with 8 edges, stop sign like
    display_bd = 185  # which border expanision index to display/ if none the one with maximum mean normal stress is chosen
    # also ony if single_tensor == False
    out_folder = "/home/user/Desktop/paper_pyTFM/data_traction_force_microscopy/ev_paper_ft_400"
    os.makedirs(out_folder, exist_ok=True)
    young = 1
    force_factor = 1
    h = 100
    pixelsize = 1

    im_shape = (400, 400)  # (400, 400)
    stress_filed_size = 150  # 150
    stress_filed_shape = "rectangle"
    stress_field_distribution = "uniform"
    sigma_n = 1
    sigma_shear = 0
    # border_ex_test = list(range(150,350,4))
    border_ex_test = sorted(np.unique(np.array(list(range(150, 325, 4)) + list(range(180, 190, 1)))))
    border_ex_test = []
    # border_ex_test = list(range(150, 158, 4))
    filter = "gaussian"
    exp_method = "manual"
    f_type = "circular"  # display_option

    mask = setup_geometry(im_shape=im_shape, shape_size=stress_filed_size, shape=stress_filed_shape)
    mask_fm = binary_dilation(mask, iterations=25)
    np.save(os.path.join(out_folder, "mask.npy"), mask)

    stress_tensor_b = setup_stress_field(mask, distribution=stress_field_distribution, sigma_n=sigma_n,
                                         sigma_shear=sigma_shear)

    # mask = binary_dilation(mask)
    masks = [binary_dilation(mask, iterations=i).astype(bool) for i in range(1, 15, 1)] + [
        binary_dilation(mask, iterations=i).astype(bool) for i in range(15, 50, 2)]
    masks_id = list(range(1,15,1))  +  list(range(15, 50, 2))
    # masks = [binary_dilation(mask, iterations=i).astype(bool) for i in range(1, 50, 2)]
    if new_calculation:
        np.save(os.path.join(out_folder, "stress_tensor_b.npy"), stress_tensor_b)
        fx_b, fy_b = force_from_stress_field(stress_tensor_b, pixelsize, plot=False, grad_type="diff1")
        np.save(os.path.join(out_folder, "fx_b.npy"), fx_b)
        np.save(os.path.join(out_folder, "fy_b.npy"), fy_b)

        u_b, v_b = fourrier_deformation(fx_b, fy_b, pixelsize, young, sigma=0.49)

        np.save(os.path.join(out_folder, "u_b.npy"),  u_b)
        np.save(os.path.join(out_folder, "v_b.npy"),  v_b)
        # forward workflow
        # tractions from deformation
        # assuming pixelsize == 1
        fx_f, fy_f = traction_wrapper(u_b, v_b, pixelsize, h, young, mask=mask,
                                      filter=filter, fs=2 * pixelsize)

        np.save(os.path.join(out_folder, "fx_f.npy"), fx_f)
        np.save(os.path.join(out_folder, "fy_f.npy"), fy_f)
        # stress from tractions
        if single_tensor:
            mask_fem = expand_mask(mask, FEM_grid_size_single, mask.shape,
                                   method=FEM_grid_method_single)  # highest aggreement (not completely sure
            UG_sol, stress_tensor_f = stress_wrapper(mask_fem, fx_f, fy_f, young, sigma=0.5)
            np.save(os.path.join(out_folder, "stress_tensor_f.npy"), stress_tensor_f)

        # stress_tensors, mean_normal_list, mask_exp_list=exp_border(exp_range=[1, 2])###### save this
        # typical stress and force measures
        stress_tensors, mean_normal_list, mask_exp_list = exp_border(out_folder=out_folder,
                                                                     exp_range=border_ex_test,
                                                                     fx_f=fx_f, fy_f=fy_f, mask=mask,
                                                                     method=exp_method)
        # getting comparison scalar fields
        # side note:   mask_exp=binary_dilation(mask,iterations=15)
        avg_normal_stress_be = [np.nanmean(ms[mask.astype(bool)]) for ms in mean_normal_list]
        save_arr = np.array([np.round(np.array(avg_normal_stress_be), 5), np.array(border_ex_test)]).T
        np.savetxt(os.path.join(out_folder, "avg_norm_stress.txt"), save_arr, fmt="%.5f", delimiter=",")
    else:
        avg_normal_stress_be = None
    fx_b, fy_b = np.load(os.path.join(out_folder, "fx_b.npy")), np.load(os.path.join(out_folder, "fy_b.npy"))
    u_b, v_b = np.load(os.path.join(out_folder, "u_b.npy")), np.load(
        os.path.join(out_folder, "v_b.npy"))
    fx_f, fy_f = np.load(os.path.join(out_folder, "fx_f.npy")), np.load(os.path.join(out_folder, "fy_f.npy"))

    mask = np.load(os.path.join(out_folder, "mask.npy"))
    stress_tensors, mean_normal_list, mask_exp_list = load_exp_border(out_folder=out_folder,
                                                                      exp_range=border_ex_test)

    stress_tensor_b = np.load(os.path.join(out_folder, "stress_tensor_b.npy"))
    if single_tensor:
        mask_fem = expand_mask(mask, FEM_grid_size_single, mask.shape,
                               method=FEM_grid_method_single)  # highest aggreement (not completely sure
        stress_tensor_f = np.load(os.path.join(out_folder, "stress_tensor_f.npy"))
    else:
        if not isinstance(display_bd, int):
            display_id = np.argmax([np.nanmean(ms[mask.astype(bool)]) for ms in mean_normal_list])
        else:
            display_id = np.argmin(np.abs(np.array(border_ex_test) - display_bd))
        stress_tensor_f = stress_tensors[display_id]
        mask_fem = mask_exp_list[display_id]

    # save_arr = np.array([np.round(np.array(avg_normal_stress_be), 5), np.array(border_ex_test)]).T
    # np.savetxt(os.path.join(out_folder, "avg_norm_stress.txt"), save_arr, fmt="%.5f", delimiter=",")

    fields = {"u_b": u_b, "v_b": v_b, "fx_b": fx_b, "fy_b": fy_b, "fx_f": fx_f, "fy_f": fy_f,
              "stress_tensor_b": stress_tensor_b,
              "stress_tensor_f": stress_tensor_f,
              "mask": mask, "mask_fem": mask_fem, "mask_fm": mask_fm}
    measures = standard_measures(mask=mask, mean_normal_list=mean_normal_list, fields=fields)
    scalar_comparisons = full_field_comparision(fields)  # r gives  mostly the spatial distribution
    del scalar_comparisons["forces"]

    # cutting to ignore effects close to image edge #### mention this when talking to benn
    fields = cut_arrays([20, -20, 20, -20], fields, mode="edge")
    mask_exp_list = cut_arrays([20, -20, 20, -20], mask_exp_list, mode="edge")
    max_dict = get_max_values(u_b, v_b, fx_f, fy_f, fx_b, fy_b, stress_tensor_b,
                              exp_test=len(border_ex_test) > 0, mean_normal_list=mean_normal_list)
    max_dict["stress"] = 1


    # plot_types = [ "forces_backward", "full_stress_tensor_backward"]
    measures["contractile_force_f"] = 5.22e-10
    measures ["mean_normal_stress_f"] = 0.93
    plot_types = [ "key measures", "cbars_only"]
              #    "forces_forward", "forces_backward", "mask_outline", "cbars_only", "test for border expansion",
            #      "be5", "be2"]
    # plot_types = ["key measures"]
    # max_dict["force"] = 1

    general_display(plot_types=plot_types, pixelsize=pixelsize, max_dict=max_dict, f_type=f_type,
                    out_folder=out_folder, cmap="coolwarm", scalar_comaprisons=scalar_comparisons,
                    border_ex_test=border_ex_test, be_avm_list=avg_normal_stress_be, mask_exp_list=mask_exp_list,
                    mean_normal_list=mean_normal_list, strain_energies=None, contractilities=None,
                    fields=fields, key_values=measures, plot_gt_exp=False, dm=False, at=False,
                    cb=False)
    plt.close("all")












