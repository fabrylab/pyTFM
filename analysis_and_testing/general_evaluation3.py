import sys
sys.path.append("/home/user/Software/pyTFM/analysis_and_testing")
from evaluation_functions import *
from plotting_evaluation import *
from pyTFM.plotting import show_quiver

def exp_border_real_data():
    out_folder = "/home/user/Desktop/backup_from_harddrive/data_traction_force_microscopy/ev_paper_expansion_test"
    createFolder(out_folder)
    border_ex_test=(list(range(0,80,1)))
    f_type = "non-circular"
    young = 1
    h = 100
    pixelsize = 1
    #  retrieving clickpoints mask and traction forces
    folder="/home/user/Desktop/biophysDS/abauer/backup_from_harddrive/data_traction_force_microscopy/WT_vs_KO_images/KOshift/"
    db = clickpoints.DataFile(os.path.join(folder,"database.cdb"), "r")
    mask = db.getMask(frame=2).data == 3
    db.db.close()
    fx_f, fy_f = np.load(os.path.join(out_folder,"fx.npy")), np.load(os.path.join(out_folder,"fy.npy"))
    mask = interpolation(mask, dims=fx_f.shape, min_cell_size=100)
    mask = binary_fill_holes(mask)
    stress_tensors, mean_normal_list, mask_exp_list = exp_border(exp_range=border_ex_test, fx_f=fx_f, fy_f=fy_f, mask=mask, out_folder=out_folder, method="binary_dilation")
    stress_tensors, mean_normal_list, mask_exp_list = load_exp_border(exp_range=border_ex_test,out_folder=out_folder)


    stress_tensor_b=stress_tensors[0]
    max_dict = get_max_values(fx_f=fx_f, fy_f=fy_f, stress_tensor_b=stress_tensor_b,
                             exp_test=len(border_ex_test) > 0, mean_normal_list=mean_normal_list)


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
    plt.ioff()
    new_calculation = True
    single_tensor = False  # set true if you dont want use only one FEM grid size
    FEM_grid_size_single = 160  # only if single_tensor
    FEM_grid_method_single = "manual"  # defines weather grid is drawn new (yields rectangular shape, or whether grid is
    # expanded with binary expansion (yields form with 8 edges, stop sign like
    display_bd = 160  # which border expanision index to display/ if none the one with maximum mean normal stress is chosen
    # also ony if single_tensor == False
    out_folder = "/home/user/Desktop/backup_from_harddrive/data_traction_force_microscopy/ev_paper_expansion_ft_400"
    #createFolder(out_folder)
    young = 1
    h = 100
    sigma = 0.49
    pixelsize = 1
    ks = (2500, 2500)  # (3000,3000)

    #folders = ((400,50),(400,20),(400,10),(650,10),(2000,10))
    #folders = [(400, 10)]

    # out_folder =os.path.join(out_folder_og, "ev_paper_expansion_ft_"+str(f[0])+"_"+str(f[1]))
    print(out_folder)
    createFolder(out_folder)
    im_shape = (400, 400)  # (650, 650)
    stress_filed_size = 150  # 150
    stress_filed_shape = "rectangle"
    stress_field_distribution = "uniform"
    sigma_n = 1
    sigma_shear = 0
    # border_ex_test = list(range(150,350,4))
    # border_ex_test = sorted(np.unique(np.array(list(range(150, 325, 4)) + list(range(180, 190, 1)))))
    # border_ex_test = [150, 153, 155, 157, 160, 165, 170, 180, 190, 200, 300]
    # border_ex_test = [150, 153, 155, 157, 160, 165, 170, 180, 190, 200, 300]
    border_ex_test = [150] + list(range(155, 170)) + list(range(170, 220, 10)) + [250, 300]
    border_ex_test=border_ex_test[18:]
    # border_ex_test = list(range(150, 158, 4))
    filter = "gaussian"
    exp_method = "manual"
    f_type = "circular"  # display_option

    mask = setup_geometry(im_shape=im_shape, shape_size=stress_filed_size, shape=stress_filed_shape)

    mask_fem = expand_mask(mask, FEM_grid_size_single, mask.shape, method=FEM_grid_method_single)
    mask_fm = mask_fem
    np.save(os.path.join(out_folder, "mask.npy"), mask)

    stress_tensor_b = setup_stress_field(mask, distribution=stress_field_distribution, sigma_n=sigma_n,
                                         sigma_shear=sigma_shear)

    # mask = binary_dilation(mask)
    masks = [binary_dilation(mask, iterations=i).astype(bool) for i in range(1, 15, 1)] + [
        binary_dilation(mask, iterations=i).astype(bool) for i in range(15, 50, 2)]
    # masks = [binary_dilation(mask, iterations=i).astype(bool) for i in range(1, 50, 2)]
    if new_calculation:
        np.save(os.path.join(out_folder, "stress_tensor_b.npy"), stress_tensor_b)
        fx_b, fy_b = force_from_stress_field(stress_tensor_b, pixelsize, plot=False, grad_type="diff1")
        np.save(os.path.join(out_folder, "fx_b.npy"), fx_b)
        np.save(os.path.join(out_folder, "fy_b.npy"), fy_b)

        # u_b, v_b = deformation_by_upsampling_sparse(fx_b, fy_b, factor=10, pixelsize=pixelsize, sigma=0.49,
        #                                            young=young,
        #                                            h=h, use_numba=False, use_exact=False)

        u_b, v_b = fourrier_deformation(fx_b, fy_b, pixelsize, young, sigma=sigma, h=h, use_finite_height=True)

        # u_b, v_b = deformation_by_upsampling(fx_b, fy_b, factor=10, pixelsize=pixelsize, sigma=0.49, young=young,
        #                                 h=h, kernel_size=ks, method="fftconvolve")
        # show_quiver(u_b,v_b)
        # plt.figure();plt.imshow(np.sqrt(u_b**2+v_b**2))
        # createFolder(out_folder)

        #u_b = np.load(os.path.join(out_folder, "u_b.npy"))
        #v_b = np.load(os.path.join(out_folder, "v_b.npy"))
        np.save(os.path.join(out_folder, "u_b.npy"), u_b)
        np.save(os.path.join(out_folder, "v_b.npy"), v_b)

        # forward workflow
        # tractions from deformation
        # assuming pixelsize == 1
        fx_f, fy_f = traction_wrapper(u_b, v_b, pixelsize, h, young, mask=mask, sigma=sigma,
                                      filter=filter, fs= 3 * pixelsize)
        #show_quiver(fx_f, fy_f)
        #contractile_force1, proj_x1, proj_y1, center1 = contractillity(fx_f, fy_f, 1, mask.astype(bool))
        #contractile_force2, proj_x2, proj_y2, center2 = contractillity(fx_b, fy_b, 1, mask.astype(bool))
        #np.sum(np.sqrt(fx_f**2+ fy_f**2))
        #np.sum(np.sqrt(fx_b ** 2 + fy_b ** 2))
        #plt.figure();
        #plt.imshow(np.sqrt(fx_f ** 2 + fy_f ** 2) - np.sqrt(fx_b ** 2 + fy_b ** 2)); plt.colorbar()

        c, e = contractility_strain_energy_exp(u_b, v_b, fx_f, fy_f, masks)


        # contractile_force, proj_x, proj_y, center = contractillity(fx_f, fy_f, 1, binary_dilation(mask,iterations=40))
        # contractile_force, proj_x, proj_y, center = contractillity(fx_b, fy_b, 1, binary_dilation(mask,iterations=40))
        c_in, e_in = contractility_strain_energy_exp(u_b, v_b, fx_b, fy_b, masks)

        plt.figure()
        plt.plot(u_b[int(u_b.shape[0] / 2), :])
        plt.savefig(os.path.join(out_folder, "def_slice.png"))
        plt.figure()
        plt.plot(fx_f[int(fx_f.shape[0] / 2), :])
        plt.savefig(os.path.join(out_folder, "traction_slice.png"))
        # plt.figure()
        # plt.plot(c/c.max(),label="force")

        # plt.plot(e/e.max(),label="energy")
        # plt.plot(e_in / e_in.max(), label="input energy")
        # plt.legend()

        plt.figure()
        plt.plot(c, label="output force")
        plt.plot(c_in, label="input force")
        plt.legend()
        plt.ylim((0, c_in.max() * 1.2))
        plt.savefig(os.path.join(out_folder, "contractility_radial.png"))
        print(c.max() / c_in.max())

        fig, ax = show_quiver(fx_f, fy_f, filter=[0, 4])
        fig.savefig(os.path.join(out_folder, "tractions_output.png"))
        # ax.imshow(masks[4], alpha=0.2)
        fig, ax = show_quiver(u_b, v_b, filter=[0, 12])
        fig.savefig(os.path.join(out_folder, "deformations.png"))

        plt.figure()
        plt.plot(e, label="output strain energy")
        plt.plot(e_in, label="input strain energy")
        plt.legend()
        plt.ylim((0, e_in.max() * 1.2))
        plt.savefig(os.path.join(out_folder, "stain_energy_radial.png"))

        # ax.imshow(masks[4], alpha=0.2)
        # plt.figure();
        # fx_f, fy_f = fx_f/fx_f.max(), fy_f /fx_f.max()
        # plt.plot(fx_f[int(fx_f.shape[0] / 2), :])

        np.save(os.path.join(out_folder, "fx_f.npy"), fx_f)
        np.save(os.path.join(out_folder, "fy_f.npy"), fy_f)
        # stress from tractions
        if single_tensor:
            UG_sol, stress_tensor_f = stress_wrapper(mask_fem, fx_f, fy_f, young, sigma=0.5)
            np.save(os.path.join(out_folder, "stress_tensor_f.npy"), stress_tensor_f)

        # stress_tensors, mean_normal_list, mask_exp_list=exp_border(exp_range=[1, 2])###### save this
        # typical stress and force measures

        stress_tensors, mean_normal_list, mask_exp_list = exp_border(out_folder=out_folder,
                                                                     exp_range=border_ex_test,
                                                                     fx=fx_f, fy=fy_f, mask=mask,
                                                                     method=exp_method)
        # getting comparison scalar fields
        ### does this work???
        # side note:   mask_exp=binary_dilation(mask,iterations=15)
        avg_normal_stress_be = [np.nanmean(ms[mask.astype(bool)]) for ms in mean_normal_list]
        save_arr = np.array([np.round(np.array(avg_normal_stress_be), 5), np.array(border_ex_test)]).T
        np.savetxt(os.path.join(out_folder, "avg_norm_stress.txt"), save_arr, fmt="%.5f", delimiter=",")

    try:

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

        plot_types = [""]
        plot_types.extend(["deformation_backward", "mask", "forces_forward", "forces_backward",
                           "shear_forward", "mean_normal_stress_forward",
                           "shear_backward", "mean_normal_stress_backward", "full_stress_tensor_forward",
                           "full_stress_tensor_backward", "key measures", "deformation_backwards", "correlation"])
        # ,"test for border expansion"])
        # plot_types = [ "forces_backward", "full_stress_tensor_backward"]
        plot_types = ["correlation", "key measures", "mean_normal_stress_backward", "mean_normal_stress_forward",
                      "forces_forward", "forces_backward", "mask_outline", "cbars_only",
                      "test for border expansion",
                      "be5", "be2"]
        # plot_types = ["key measures"]
        # max_dict["force"] = 1
        try:
            general_display(plot_types=plot_types, pixelsize=pixelsize, max_dict=max_dict, f_type=f_type,
                            out_folder=out_folder, cmap="coolwarm", scalar_comaprisons=scalar_comparisons,
                            border_ex_test=border_ex_test, be_avm_list=avg_normal_stress_be,
                            mask_exp_list=mask_exp_list,
                            mean_normal_list=mean_normal_list, strain_energies=None, contractilities=None,
                            fields=fields, key_values=measures, plot_gt_exp=False, dm=False, at=False,
                            cb=False)

        except:
            pass
    except:
        pass

    plt.close("all")








