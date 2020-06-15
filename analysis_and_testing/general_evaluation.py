from evaluation_functions import *
from plotting_evaluation import *
def exp_border_real_data():
    out_folder = "/home/user/Desktop/backup_from_harddrive/data_traction_force_microscopy/ev_paper_rd_expansion"
    createFolder(out_folder)
    border_ex_test=(list(range(0,80,1)))
    f_type = "non-circular"
    young = 1
    h = 100
    pixelsize = 1
    #  retrieving clickpoints mask and traction forces
    folder="/home/user/Desktop/backup_from_harddrive/data_traction_force_microscopy/WT_vs_KO_images/KOshift/"
    db = clickpoints.DataFile(os.path.join(folder,"database.cdb"), "r")
    mask = db.getMask(frame=2).data == 3
    db.db.close()
    fx_f, fy_f =np.load(os.path.join(out_folder,"u.npy")), np.load(os.path.join(out_folder,"v.npy"))
    mask = interpolation(mask, dims=fx_f.shape, min_cell_size=100)
    mask = binary_fill_holes(mask)
    stress_tensors, mean_normal_list, mask_exp_list = exp_border(exp_range=border_ex_test, fx_f=fx_f, fy_f=fy_f, mask=mask, out_folder=out_folder, method="binary_dilation")
    stress_tensors, mean_normal_list, mask_exp_list = load_exp_border(exp_range=border_ex_test,out_folder=out_folder)


    stress_tensor_b=stress_tensors[0]
    max_dict = get_max_values(fx_f=fx_f, fy_f=fy_f, stress_tensor_b=stress_tensor_b,
                             exp_test=len(border_ex_test) > 0, mean_normal_list=mean_normal_list)

    # getting comparison scalar fields
    mask_fm = standard_measures(mask=mask, mean_normal_list=mean_normal_list, stress_tensor_b=stress_tensor_b)
    save_arr = np.array([np.round(np.array(avg_normal_stress_be),5),np.array(border_ex_test)]).T
    save_arr[:,1]=save_arr[:,1]/2 # to get "exapansion radius"
    np.savetxt(os.path.join(out_folder,"avg_norm_stress.txt"), save_arr, fmt="%.5f", delimiter=",")



    mask_exp = binary_dilation(mask, iterations=15)
    scalar_comaprisons = full_field_comparision()  # r gives  mostly the spatial distribution
    with suppress(KeyError): del scalar_comaprisons["forces"]
    plot_types = ["test for border expansion"]
    plot_types.extend(["forces_forward", "correlation", "test for border expansion"])
    # plot_types = [ "forces_backward", "full_stress_tensor_backward"]
    general_display(plot_types=plot_types, mask=mask, pixelsize=pixelsize, max_dict=max_dict, f_type=f_type,
                    mean_normal_list=mean_normal_list, mask_exp_list=mask_exp_list,out_folder=out_folder,
                    fx_f=fx_f,fy_f=fy_f,mask_exp=mask_exp,scalar_comaprisons=scalar_comaprisons,border_ex_test=border_ex_test,plot_gt_exp=False)
    plt.close("all")


#exp_border_real_data()
if __name__ == "__main__":

    new_calcualtion = False
    single_tensor = False # set true if you dont want use only one FEM grid size
    FEM_grid_size_single = None # only if single_tensor
    FEM_grid_method_single = "manual" # defines weather grid is drawn new (yields rectangular shape, or whether grid is
    # expanded with binary expansion (yields form with 8 edges, stop sign like
    display_bd = 185 # which border expanision index to display/ if none the one with maximum mean normal stress is chosen
    # also ony if single_tensor == False
    out_folder = "/home/user/Desktop/backup_from_harddrive/data_traction_force_microscopy/ev_paper_expansion"
    createFolder(out_folder)
    young = 1
    h = 100
    pixelsize = 1
    ks = (2500, 2500)
    def_factor = 10
    im_shape = (650, 650)
    stress_filed_size = 150
    stress_filed_shape = "rectangle"
    stress_field_distribution = "uniform"
    sigma_n = 1
    sigma_shear = 0
    # border_ex_test = list(range(150,350,4))
    border_ex_test = sorted(np.unique(np.array(list(range(150, 325, 4)) + list(range(180, 190, 1)))))
    # border_ex_test = list(range(150, 158, 4))
    filter = "gaussian"
    exp_method = "manual"
    f_type = "circular"  # display_option


    mask = setup_geometry(im_shape=im_shape, shape_size=stress_filed_size, shape=stress_filed_shape)
    np.save(os.path.join(out_folder, "mask.npy"), mask)

    stress_tensor_b = setup_stress_field(mask, distribution=stress_field_distribution, sigma_n=sigma_n, sigma_shear=sigma_shear)
    if new_calcualtion:
        np.save(os.path.join(out_folder, "stress_tensor_b.npy"), stress_tensor_b)

        fx_b, fy_b = traction_from_stress(stress_tensor_b, pixelsize, plot=False, grad_type="diff1")
        np.save(os.path.join(out_folder, "fx_b.npy"), fx_b)
        np.save(os.path.join(out_folder, "fy_b.npy"), fy_b)

        u_b, v_b = deformation_by_upsampling(fx_b, fy_b, factor=10, pixelsize=pixelsize, sigma=0.49, young=young, h=h,
                                             kernel_size=ks, method="fftconvolve")

        # u_b, v_b = finite_thickness_convolution(fx_b, fy_b, pixelsize, h, young, sigma=0.49,
        #                                        kernel_size=None)  # somwhat of an approximation

        # createFolder(out_folder)
        np.save(os.path.join(out_folder, "u_b.npy"), u_b)
        np.save(os.path.join(out_folder, "v_b.npy"), v_b)

        # compare_scalar_fields(u,u1)
        # u2,v2,z1=finite_thickenss_convolution_exact(fx, fy,pixelsize, h, young, sigma=0.5,kernel_size=None) # best solution possible

        ##forward worklow
        # tractions from deformation
        fx_f, fy_f = traction_wrapper(u_b, v_b, pixelsize, h, young, mask=mask,
                                      filter=filter, fs=3)  # assuming pixelsize == 1
        np.save(os.path.join(out_folder, "fx_f.npy"), fx_f)
        np.save(os.path.join(out_folder, "fy_f.npy"), fy_f)
        # stress from tractions
        if single_tensor:
            mask_fem = expand_mask(mask, FEM_grid_size_single, mask.shape, method=FEM_grid_method_single)  # highest aggreement (not completely sure
            UG_sol, stress_tensor_f = stress_wrapper(mask_fem, fx_f, fy_f, young, sigma=0.5)
            np.save(os.path.join(out_folder, "stress_tensor_f.npy"), stress_tensor_f)


        # stress_tensors, mean_normal_list, mask_exp_list=exp_border(exp_range=[1, 2])###### save this
        # typical stress and force measures

        stress_tensors, mean_normal_list, mask_exp_list = exp_border(out_folder=out_folder, exp_range=border_ex_test,
                                                                     fx_f=fx_f, fy_f=fy_f, mask=mask, method=exp_method)
        # getting comparison scalar fields
        ### does this work???
        # side note:   mask_exp=binary_dilation(mask,iterations=15)

        mask_fm = standard_measures(mask=mask, mean_normal_list=mean_normal_list, stress_tensor_b=stress_tensor_b,
                                    stress_tensor_f=stress_tensor_f)
        save_arr = np.array([np.round(np.array(avg_normal_stress_be), 5), np.array(border_ex_test)]).T
        np.savetxt(os.path.join(out_folder, "avg_norm_stress.txt"), save_arr, fmt="%.5f", delimiter=",")



    fx_b, fy_b = np.load(os.path.join(out_folder, "fx_b.npy")), np.load(os.path.join(out_folder, "fy_b.npy"))
    u_b, v_b = np.load(os.path.join(out_folder, "u_b.npy")), np.load(
        os.path.join(out_folder, "v_b.npy"))
    fx_f, fy_f = np.load(os.path.join(out_folder, "fx_f.npy")), np.load(os.path.join(out_folder, "fy_f.npy"))

    mask = np.load(os.path.join(out_folder, "mask.npy"))
    stress_tensors, mean_normal_list, mask_exp_list = load_exp_border(out_folder=out_folder,
                                                                      exp_range=border_ex_test)

    stress_tensor_b = np.load(os.path.join(out_folder, "stress_tensor_b.npy"))
    if single_tensor:
        mask_fem = expand_mask(mask, FEM_grid_size_single, mask.shape, method=FEM_grid_method_single)  # highest aggreement (not completely sure
        stress_tensor_f = np.load(os.path.join(out_folder, "stress_tensor_f.npy"))
    else:
        if not isinstance(display_bd, int):
            display_id = np.argmax([np.nanmean(ms[mask.astype(bool)]) for ms in mean_normal_list])
        else:
            display_id = np.argmin(np.abs(np.array(border_ex_test) - display_bd))
        stress_tensor_f = stress_tensors[display_id]
        mask_fem = mask_exp_list[display_id]




    # getting comparison scalar fields


    ####### is mask_fm defined?? ######
    measures = standard_measures(mask=mask, mean_normal_list=mean_normal_list, stress_tensor_b=stress_tensor_b,
                                stress_tensor_f=stress_tensor_f)
    #save_arr = np.array([np.round(np.array(avg_normal_stress_be), 5), np.array(border_ex_test)]).T
    #np.savetxt(os.path.join(out_folder, "avg_norm_stress.txt"), save_arr, fmt="%.5f", delimiter=",")

    fields = {"u_b":u_b, "v_b": v_b, "fx_b": fx_b, "fy_b": fy_b, "fx_f": fx_f, "fy_f": fy_f, "stress_tensor_b": stress_tensor_b,
              "stress_tensor_f": stress_tensor_f,
              "mask": mask, "mask_fem": mask_fem, "mask_fm": measures["mask_fm"]}

    scalar_comaprisons = full_field_comparision(fields)  # r gives  mostly the spatial distribution
    del scalar_comaprisons["forces"]
    key_values = [cont_energy_b, cont_energy_f, contractile_force_b, contractile_force_f,
                  mean_normal_stress_b, mean_normal_stress_f, mean_shear_b, mean_shear_f]

    # cutting to ignore effects close to image edge #### mention this when talking to benn
    fields = cut_arrays([20, -20, 20, -20], fields, mode="edge")
    mask_exp_list = cut_arrays([20, -20, 20, -20], mask_exp_list, mode="edge")
    max_dict = get_max_values(u_b, v_b, fx_f, fy_f, fx_b, fy_b, stress_tensor_b,
                              exp_test=len(border_ex_test) > 0, mean_normal_list=mean_normal_list)
    max_dict["stress"] = 1

    plot_types = [""]
    plot_types.extend(["deformation_backward", "mask", "forces_forward", "forces_backward",
                       "shear_forward", "mean_normal_stress_forward",
                       "shear_backward", "mean_normal_stress_backward", "full_stress_tensor_forward",
                       "full_stress_tensor_backward", "key measures", "deformation_backwards", "correlation"])
    # ,"test for border expansion"])
    # plot_types = [ "forces_backward", "full_stress_tensor_backward"]
    plot_types = ["correlation", "key measures", "mean_normal_stress_backward", "mean_normal_stress_forward",
                  "forces_forward", "forces_backward", "mask_outline", "cbars_only", "test for border expansion", "be5"]
    # plot_types = ["key measures"]
    # max_dict["force"] = 1

    general_display(plot_types=plot_types, pixelsize=pixelsize, max_dict=max_dict, f_type=f_type, fields=fields,
                    out_folder=out_folder, cmap="coolwarm", scalar_comaprisons=scalar_comaprisons,
                    border_ex_test=border_ex_test, be_avm_list=rel_av_norm_stress, mask_exp_list=mask_exp_list,
                    mean_normal_list=mean_normal_list, key_values=measures,  plot_gt_exp=False, dm=True,
                    at=False,
                    cb=False)

    plt.close("all")








