### function integrating Traktion force microscopy into a clcikpoints database

from andreas_TFM_package.grid_setup_solids_py import *
from andreas_TFM_package.functions_for_cell_colonie import *
from andreas_TFM_package.solids_py_stress_functions import *
from andreas_TFM_package.utilities_TFM import *
from andreas_TFM_package.TFM_functions import *
import solidspy.postprocesor as pos
import solidspy.assemutil as ass
from peewee import IntegrityError
import clickpoints
import os
import re
import warnings
from tqdm import tqdm


# preset names of layers
layer_name_dict=defaultdict(list)
layer_name_dict["deformation"]="def_plots"
layer_name_dict["tracktion_force"]="traction_plots"
layer_name_dict["FEM_analysis"]="FEM_output"

# some message to be printed

calculation_messages_multiple={"deformation":"calculating deformation on frames %s",
                    "tracktion_force":"calculating tracktion_forces on frames %s",
                    "FEM_analysis":"FEM_analysis on frames %s",
                    "get_contractility_contractile_energy": "contractility/contractile energy on frames %s"
                                    }
calculation_messages_single={"deformation":"calculating deformation on current frame",
                    "tracktion_force":"calculating tracktion_forces on current frame",
                    "FEM_analysis":"FEM_analysis on current frame",
                     "get_contractility_contractile_energy": "contractility/contractile energy on current frame"
                                      }
# units of the returned stress and energy measures:
units={	"avarage line stress":"N/m",
	"avarage cell force":"N",
	"avarage cell pressure"	:"N/m",
	"avarage cell shear":"N/m",
	"std line stress":"N/m",
	"std cell force":"N",
	"std cell pressure":"N/m",
	"std cell shear":"N/m",
"contractility":"N",
"contractile energy":"J",
"contractility per cell":"N",
"contractile energy per cell":"J",
"contractility per area":"N/m2",
"contractile energy per area":"J/m2",
"avarage normal stress":"N/m",
"avarage shear stress":"N/m",
"area":"m2",
"sum deformations":"pixels",
"sum traction forces":"N/m2"
}



def default_parameters():
    sigma = 0.49  # poison ratio
    young = 49000  # youngsmodulus
    pixelsize = 0.201  # µm/pixel
    window_size = 100  # in pixel , keep at 10 µm
    overlapp = 85  # window_size/2
    std_factor = 15
    #pixelsize2 = pixelsize * (window_size - overlapp)  # pixelsize of the deformation image
    h = 300
    parameter_dict = make_paramters_dict_tfm(sigma=sigma, young=young, pixelsize=pixelsize, window_size=window_size,
                                             overlapp=overlapp,
                                             std_factor=std_factor, h=h, pixel_factor=window_size - overlapp)
    return parameter_dict

class Mask_Error(Exception):
    pass




def write_output_file_with_nice_prefix(values,value_type, file_path):

    if value_type=="parameters":
        with open(file_path, "w+") as f:
            f.write("analysis_paramters\n")
            for parameter, value in values.items():
                f.write(parameter + "\t" + str(value)+ "\n")

    if value_type == "results":
        # list of available frames sorted
        frames=list(values.keys())
        frames.sort(key=lambda x:int(x))

        with open(file_path, "a+") as f:
            for frame in frames:
                res_part=values[frame]
                for name,res in res_part.items():
                    res_round,pref=find_prefix(res)
                    f.write(frame+"\t"+name+"\t"+str(res_round)+" "+pref+units[name]+"\n")


def write_output_file(values,value_type, file_path):

    if value_type=="parameters":
        with open(file_path, "w+") as f:
            f.write("analysis_paramters\n")
            for parameter, value in values.items():
                f.write(parameter + "\t" + str(value)+ "\n")
    if value_type == "results":
        # list of available frames sorted
        frames=list(values.keys())
        frames.sort(key=lambda x:int(x))

        with open(file_path, "w+") as f:
            for frame in frames:
                res_part=values[frame]
                for name,res in res_part.items():
                    res_unpack,warn=unpack_list(res)
                    warn_empty=(warn=="")
                    f.write(frame+"\t"+name+"\t"+str(round_flexible(res_unpack))+"\t"+units[name]+"\t"*warn_empty+warn+"\n")



def except_error(func, error, **kwargs):  # take functio and qkwarks
    '''
    wraper to handle errors and return false if the exception is encountered
    :param func:
    :param error:
    :param kwargs:
    :param return_values:
    :return:
    '''
    # if function doesnt return any values

    try:
        values = func(**kwargs)
    except error:
        return False
    return values


def try_to_load_mask(db,frame,mtype="membrane"):

    try:
        mask = db.getMask(frame=frame,layer=1).data
        # extract only one type of mask
        if mtype=="membrane": #
            mask=mask==1
        if mtype=="contractility_select":
            mask = mask == 2
    except AttributeError: # checks if mask object exists
        raise Mask_Error("no mask of the cell membrane found for frame " + str(frame))

    if np.sum(mask)==0: # checks if mask is empty
        raise Mask_Error("no mask of the cell membrane found for frame " + str(frame))
    else:
        return mask

def try_to_load_deformation(path, frame, warn=False):
    '''loading the deformations fro a given frame. If deformations are not found either raises an error
    or a warning (warn=True) and returns None for u and v

    :param path:
    :param frame:
    :param warn:
    :return:
    '''
    try:
        u = np.load(os.path.join(path, frame + "u.npy"))
        v = np.load(os.path.join(path, frame + "v.npy"))
    except FileNotFoundError:
        if warn:
            warnings.warn("no deformations found for frame " + frame)
        else:
            raise FileNotFoundError("no deformations found for frame " + frame)
        return (None, None)
    return (u, v)


def try_to_load_traction(path, frame, warn=False):
    '''loading the tractions from a given frame. If tractions are not found either raises an error
    or a warning (warn=True) and returns None for u and v

    :param path:
    :param frame:
    :param warn:
    :return:
    '''
    try:
        t_x = np.load(os.path.join(path, frame + "tx.npy"))
        t_y = np.load(os.path.join(path, frame + "ty.npy"))
    except FileNotFoundError as e:
        if warn:
            warnings.warn("no traction forces found for frame " + frame)
        else:
            raise FileNotFoundError("no traction forces found for frame " + frame)
        return (None, None)
    return (t_x, t_y)


def check_small_mask(mask, warn_thresh=300):
    '''
    warns when mask is very small and raises error if amsk is too small
    '''
    mask_sum=np.sum(mask)
    if mask_sum==0:
        raise Mask_Error("mask empty after filtering/interpolation")
        return "error"
    if mask_sum<warn_thresh:
        warnings.warn("mask is very small, consider increasing the resoultion of deformation and traction field")
        return "mask for cellcolny is very small"
    return ""

def setup_database_for_tfm(folder, name, return_db=False):

    db = clickpoints.DataFile(os.path.join(folder, "database3.cdb"), "w")
    images = [x for x in os.listdir(folder) if ".tif" in x]
    frames = [get_group(re.search('(\d{1,4})', x), 0) for x in images]
    sort_id_list=make_rank_list(frames,dtype=int)# list of sort indexes (frames) of images in the database
  
    layer_list = ["images_after", "images_before","membranes","def_plots","traction_plots","FEM_output"]
    base_layer = db.getLayer(layer_list[0], create=True, id=0)
    for l in layer_list[1:]:
        db.getLayer(l, base_layer=base_layer, create=True)
    db.setMaskType("membrane", color="#99EA44",index=1)
    db.setMaskType("contractility_select", color="#ff0000",index=2)
    path = db.setPath(folder, 1)  # for setting the path to the images

    for id, (sort_index_id,frame, im) in enumerate(zip(sort_id_list,frames, images)):
        if "after" in im:
            db.setImage(id=id, filename=im, sort_index=sort_index_id
                        , layer="images_after", path=1)
            db.setAnnotation(filename=im, comment=frame + "after__"+str(sort_index_id)+"sid")
            # setting new layer and sorting in at the same time
        if "before" in im and not "bf_" in im:
            db.setImage(id=id, filename=im
                        , sort_index=sort_index_id, layer="images_before", path=1)
            db.setAnnotation(filename=im, comment=frame + "before__"+str(sort_index_id)+"sid")
        if "bf_before" in im:  # change names
            db.setImage(id=id, filename=im, sort_index=sort_index_id,
                        layer="membranes", path=1)
            db.setAnnotation(filename=im, comment=frame+"bf__"+str(sort_index_id)+"sid")
    # delete_empty_layers(db) # not necessary
    if return_db:
        return frames, db
    else:
        db.db.close()
        return frames


def assert_and_set_default_TFM_mode(parameter_dict):
    valid_values = ["finite_thikness", "infinite_thikness"]
    if not "TFM_mode" in parameter_dict.keys(): # setting default value for tfm mode
        parameter_dict["tfm_mode"]="finite_thikness"
        print("using finite thikness TFM")
    else:
        v_v=["'%s'"%s for s in valid_values]
        separator=" or " if len(valid_values)==2 else ", "
        assert (parameter_dict["TFM_mode"] in valid_values), "invalid choice for TFM mode. Use "+separator.join(v_v)

def create_layers_on_demand(db, layer_list):

    '''
    :param db: clickpointsw database
    :param layer_list: list of layer names that should be created
    :return:
    '''
    layer_list=make_iterable(layer_list) 
    base_layer = db.getLayer(id=1)
    layer_names = [l.name for l in db.getLayers()]  # all existing layer

    for pl in layer_list:
        if pl not in layer_names:
            db.getLayer(pl, base_layer=base_layer, create=True)


def get_file_order_and_frames(db):

    # string to search for frame: max for numbers at the beginning
    s_frame='^(\d{1,4})'  
    #string to search for sort index (frame in the cdb database)
    # max for numbers after__ 
    s_sid='__(\d{1,4})'

    file_order = defaultdict(list) # name of the image (contains before and after): id in cdb database
    frames_ref_dict={} #frame of image: frame in cdb database
    
    for an in db.getAnnotations():
        img_frame,name,cdb_frame=get_group(re.search(s_frame+'(\w{1,})'+s_sid, an.comment), "all")
        frames_ref_dict[img_frame]=int(cdb_frame)
        file_order[img_frame+name].append(an.image_id)
    unique_frames = np.unique(list(frames_ref_dict.keys()))
    
    return unique_frames, file_order,frames_ref_dict

def get_frame_from_annotation(db,cdb_frame):
    '''

    :param db:
    :param cdb_frame: number of the frame (sort index in the cdb database
    :return: str_frame: frame as a string in the typical convention
    '''
    annotations=[]
    for an in db.getAnnotations(frame=cdb_frame):
        annotations.append(an.comment)
    str_frame = np.unique([get_group(re.search('(\d{1,4})', x), 0) for x in annotations])
    return str_frame[0]




def deformation(frame, parameter_dict,res_dict, db, frames_ref_dict=None, file_order=None, path=None,single=True,**kwargs):
    if single:
        frames, file_order,frames_ref_dict = get_file_order_and_frames(db)
        path = db.getPath(id=1).path
        layer_list = ["def_plots"]
        create_layers_on_demand(db, layer_list)
        print(calculation_messages_single["deformation"])


    # deformation for 1 frame
    im1 = db.getImage(id=file_order[frame + "after"]).data  ## thats very slow though
    im2 = db.getImage(id=file_order[frame + "before"]).data
    u, v, x, y, mask, mask_std = calculate_deformation(im1.astype(np.int32), im2.astype(np.int32),
                                                       parameter_dict["window_size"]
                                                       , parameter_dict["overlapp"],
                                                       std_factor=parameter_dict["std_factor"])
    res_dict[frame]["sum deformations"] = np.sum(np.sqrt(u ** 2 + v** 2))
    # plotting
    plt.ioff()
    dpi = 200
    fig1 = show_quiver_clickpoints(u, v, filter=[0, int(np.ceil(u.shape[0] / 40))], scale_ratio=0.2,
                                   headwidth=3, headlength=3, width=0.002,
                                   figsize=(im1.shape[1] / dpi, im1.shape[0] / dpi), cbar_str="deformation\n[pixel]")
    # saving plots
    fig1.savefig(os.path.join(path, frame + "deformation.png"), dpi=200)
    # closing figure objects
    plt.close(fig1)
    # adding plots to data base
    # if entry already exist dont change anything (image is overwritten, but database entry stays the same)
    except_error(db.setImage, IntegrityError, filename=frame + "deformation.png", layer="def_plots", path=1,
                 sort_index=frames_ref_dict[frame])
    # saving raw files
    np.save(os.path.join(path, frame + "u.npy"), u)
    np.save(os.path.join(path, frame + "v.npy"), v)
    return None, frame

def plot_deformations(folder):
    files=os.listdir(folder)
    files_dict=defaultdict(dict)
    for f in files:
        if f[2]=="u":
            files_dict[f[:2]]["u"]=f
        if f[2] == "v":
            files_dict[f[:2]]["v"]=f
    for frame, files in files_dict.items():
        u=np.load(os.path.join(folder,files["u"]))
        v = np.load(os.path.join(folder,files["v"]))
        dpi = 200
        fig1 = show_quiver_clickpoints(u, v, filter=[0, int(np.ceil(u.shape[0] / 40))], scale_ratio=0.2,
                                       headwidth=3, headlength=3, width=0.002,
                                       figsize=(2022 / dpi, 2011 / dpi),
                                       cbar_str="deformation\n[pixel]")
        fig1.savefig(os.path.join(folder, frame + "deformation.png"), dpi=200)




def get_contractility_contractile_energy(frame, parameter_dict,res_dict, db,frames_ref_dict=None, path=None,
                                         single=True,per_cell=True,per_area=True,**kwargs):

    if single:
        frames, file_order,frames_ref_dict = get_file_order_and_frames(db)
        path = db.getPath(id=1).path
        layer_list = ["FEM_output"]
        create_layers_on_demand(db, layer_list)
        print(calculation_messages_single["get_contractility_contractile_energy"])

    mask=try_to_load_mask(db,frames_ref_dict[frame],mtype="contractility_select")

    u,v=try_to_load_deformation(path, frame, warn=True)
    t_x,t_y=try_to_load_traction(path, frame, warn=False)

    # filling holes
    mask=binary_fill_holes(mask)
    warn="selected area is very small" if np.sum(mask) < 1000 else "" # set a warning written in the outputfile

    ps_new = parameter_dict["pixelsize_beads_image"] * np.mean(  # should be equivalent to "pixelsize_def_image"
        np.array(mask.shape) / np.array(t_x.shape))
    # interpolation to size of traction force array
    mask_int = interpolation(mask, t_x.shape)


    contractile_force, proj_x, proj_y,center=contractility(t_x, t_y, ps_new, mask_int)
    # calculate contractile energy if deformations are provided
    if isinstance(u,np.ndarray):
        contr_energy=contractile_energy(u,v,t_x,t_y,parameter_dict["pixelsize_beads_image"], ps_new,mask_int)
    else:
        contr_energy=None

    # writing results to dictionray

    res_dict[frame]["contractility"]=[contractile_force,warn]
    res_dict[frame]["contractile energy"]= [contr_energy,warn]



    # normalize with number of cells, these are exactly the cells recognized in FEM analysis
    # could be problematic for large pixelsize of deformation image
    # maybe improve this a bit
    if per_cell:
        mask_membrane = try_to_load_mask(db, frames_ref_dict[frame], mtype="membrane")
        mask_membrane = remove_small_holes(mask_membrane, 100)
        mask_membrane = remove_small_objects(label(mask_membrane), 1000) > 0  # removing other small bits
        mask_int_mem = interpolation(mask_membrane, t_x.shape)
        mask_area, mask_boundaries, borders = prepare_mask(mask_int_mem,min_cell_size=5)
        n_cells=len(borders.cell_ids)
        res_dict[frame]["contractility per cell"] = [contractile_force/n_cells,warn]
        res_dict[frame]["contractile energy per cell"] = [contr_energy/n_cells,warn]
    if per_area:
        mask_membrane = try_to_load_mask(db, frames_ref_dict[frame], mtype="membrane")
        mask_membrane = binary_fill_holes(mask_membrane)
        area=np.sum(mask_membrane)*((ps_new*10**-6)**2)
        res_dict[frame]["contractility per area"] = [contractile_force / area,warn]
        res_dict[frame]["contractile energy per area"] = [contr_energy / area,warn]
        res_dict[frame]["area"] =  area
    return (contractile_force, contr_energy), frame





def tracktion_force(frame, parameter_dict,res_dict, db, path=None, frames_ref_dict=None, single=True,im_shape=0,**kwargs):

    assert_and_set_default_TFM_mode(parameter_dict)
    if single:
        frames, file_order,frames_ref_dict = get_file_order_and_frames(db)
        im_shape=db.getImage(0).data.shape
        path = db.getPath(id=1).path
        layer_list = ["traction_plots"]
        create_layers_on_demand(db, layer_list)
        print(calculation_messages_single["tracktion_force"])

    # trying to laod deformation
    u,v=try_to_load_deformation(path, frame, warn=False)
    ps_new = parameter_dict["pixelsize_beads_image"] * np.mean(  # should be equivalent to "pixelsize_def_image"
        np.array(u.shape) / np.array(u.shape))  # pixelsize of fem grid in µm
    # using tfm with or without finite thikness correction
    if parameter_dict["TFM_mode"] == "finite_thikness":
        tx, ty = ffttc_traction_finite_thickness(u, v, pixelsize1=parameter_dict["pixelsize_beads_image"],
                                                     pixelsize2=ps_new,
                                                     h=parameter_dict["h"], young=parameter_dict["young"],
                                                     sigma=parameter_dict["sigma"],
                                                     filter="gaussian")
        # unit is N/m**2
        if np.isnan(tx).all():
            warnings.warn("falling back to inifinte thikness assumption due to nummerical issues")
            parameter_dict["TFM_mode"] = "infinite_thikness"
    if parameter_dict["TFM_mode"] == "infinite_thikness":
        tx, ty = ffttc_traction(u, v, pixelsize1=parameter_dict["pixelsize_beads_image"],
                                                 pixelsize2=ps_new,
                                                 young=parameter_dict["young"],
                                                 sigma=parameter_dict["sigma"],
                                                 filter="gaussian")
    res_dict[frame]["sum traction forces"]=np.sum(np.sqrt(tx**2+ty**2))
    # plotting
    plt.ioff()
    dpi = 200
    fig2 = show_quiver_clickpoints(tx, ty, filter=[0, int(np.ceil(u.shape[0] / 40))], scale_ratio=0.2,
                                   headwidth=3, headlength=3, width=0.002,
                                   figsize=(im_shape[1] / dpi, im_shape[0] / dpi), cbar_str="tracktion\n[Pa]")
    # saving plots
    fig2.savefig(os.path.join(path, frame + "traction.png"), dpi=200)
    # closing figure objects
    plt.close(fig2)
    # adding plots to data base
    # if entry already exist dont change anything (image is overwritten, but database entry stays the same)
    except_error(db.setImage, IntegrityError, filename=frame + "traction.png", layer="traction_plots", path=1,
                 sort_index=frames_ref_dict[frame])
    # saving raw files
    np.save(os.path.join(path, frame + "tx.npy"), tx)
    np.save(os.path.join(path, frame + "ty.npy"), ty)
    return None, frame




def FEM_analysis(frame,parameter_dict,res_dict,db,path=None, frames_ref_dict=None,single=True,**kwargs):
    # performing MSM/finite elements analysis
    if single:
        frames, file_order,frames_ref_dict = get_file_order_and_frames(db)
        path = db.getPath(id=1).path
        layer_list = ["FEM_output"]
        create_layers_on_demand(db, layer_list)
        print(calculation_messages_single["FEM_analysis"])


    # trying to load mask, skip if not found
    mask = try_to_load_mask(db, frames_ref_dict[frame],mtype="membrane")
    im_shape=mask.shape
    # trying to traction forces, skip if not found
    t_x, t_y = try_to_load_traction(path, frame, warn=False)
    # some pre clean up of the mask
    mask = remove_small_holes(mask, 100)
    mask = remove_small_objects(label(mask), 1000) > 0  # removing other small bits
    # interpolation to size of traction force array
    mask_int = interpolation(mask, t_x.shape)
    warn=check_small_mask(mask_int)  ## imprve to get full error catching when mask is e.g. dicontinous
    # further preparation of mask data
    mask_area, mask_boundaries, borders = prepare_mask(mask_int, min_cell_size=5) # min_cell_size at least 2 or none
    # plt.figure();plt.imshow(mask_area)
    ps_new = parameter_dict["pixelsize_beads_image"] * np.mean(  # should be equivalent to "pixelsize_def_image"
        np.array(mask.shape) / np.array(mask_area.shape))  # pixelsize of fem grid in µm (?? is this ok??)

    # preparing force mao
    f_x = t_x * ((ps_new * (10 ** -6)) ** 2)  # point force for each node from tractions
    f_y = t_y * ((ps_new * (
                10 ** -6)) ** 2)  ## this factor is just for numerical reasons.... ## alsow try to mae this more efficient
    f_x[~mask_area] = np.nan  # setting all values outside of maske area to zero
    f_y[~mask_area] = np.nan
    f_x_c1 = f_x - np.nanmean(f_x)  # normalizing traction force to sum up to zero (no displacement)
    f_y_c1 = f_y - np.nanmean(f_y)
    f_x_c2, f_y_c2, p = correct_torque(f_x_c1, f_y_c1, mask_area)
    # get_torque1(f_y,f_x,mask_area)

    # setup of the grid
    nodes, elements, loads, mats = grid_setup(mask_area, f_x_c2, f_y_c2, 1, sigma=0.5)
    DME, IBC, neq = ass.DME(nodes, elements)  # boundary conditions asembly??
    print("Number of elements: {}".format(elements.shape[0]))
    print("Number of equations: {}".format(neq))

    # System assembly
    KG = ass.assembler(elements, mats, nodes, neq, DME, sparse=True)
    RHSG = ass.loadasem(loads, IBC, neq)

    # System solution with custom conditions
    UG_sol, rx = custom_solver(KG, RHSG, mask_area,
                               verbose=True)  # solver with constratinst to zero translation and zero rotation

    if not (np.allclose(KG.dot(UG_sol) / KG.max(), RHSG / KG.max())):
        print("The system is not in equilibrium!")
    UC = pos.complete_disp(IBC, nodes, UG_sol)  # uc are x and y displacements
    E_nodes, S_nodes = pos.strain_nodes(nodes, elements, mats, UC)  # stresses and strains
    stress_tensor = calculate_stress_tensor(S_nodes, nodes, dims=mask_area.shape)  # assembling the stress tensor

    # average shear and normal stress on the colony area
    avg_shear, avg_normal_stress = calculate_mean_stress_measure(mask_area, stress_tensor, ps_new)



    res_dict[frame]["avarage normal stress"] = [avg_normal_stress,warn]
    res_dict[frame]["avarage shear stress"] = [avg_shear,warn]
    ### other possible stress measures, just for a nice picture
    #sigma_max, sigma_min, tau_max, phi_n, phi_shear, sigma_avg = all_stress_measures(S_nodes, nodes,
     #                                                                                dims=mask_area.shape)
    #sigma_max_abs = np.maximum(np.abs(sigma_min), np.abs(sigma_max))  ### highest possible norm of the stress tensor



    # retrievieving spline representation of bourders
    lines_spline_points = borders.lines_spline_points
    lines_splines = borders.lines_splines
    lines_points = borders.lines_points
    # plot linestresses over border as continous curves:
    lines_interpol, min_v, max_v = interpolation_for_stress_and_normal_vector(lines_splines, lines_points,
                                                                              stress_tensor, pixel_length=ps_new,
                                                                              interpol_factor=6)
    # calcualting measures for stress on cell borders
    avg_line_stress=mean_stress_vector_norm(lines_interpol, borders, norm_level="points", vtype="t_vecs")
    avg_cell_force=mean_stress_vector_norm(lines_interpol, borders, norm_level="cells", vtype="t_vecs")
    avg_cell_pressure=mean_stress_vector_norm(lines_interpol, borders, norm_level="cells", vtype="tn")
    avg_cell_shear=mean_stress_vector_norm(lines_interpol, borders, norm_level="cells", vtype="ts")

    res_dict[frame]["avarage line stress"]=[avg_line_stress[1],warn]
    res_dict[frame]["avarage cell force"] =[avg_cell_force[1],warn]
    res_dict[frame]["avarage cell pressure"] =[avg_cell_pressure[1],warn]
    res_dict[frame]["avarage cell shear"] =[avg_cell_shear[1],warn]
    res_dict[frame]["std line stress"] = avg_line_stress[2]
    res_dict[frame]["std cell force"] = avg_cell_force[2]
    res_dict[frame]["std cell pressure"] = avg_cell_pressure[2]
    res_dict[frame]["std cell shear"] = avg_cell_shear[2]

    #plot of line stresses to the database
    plt.ioff()
    dpi = 200

    fig=plot_continous_boundary_stresses(mask_area.shape, borders.edge_lines, lines_interpol, min_v,
                                     max_v, figsize=(im_shape[1] / dpi, im_shape[0] / dpi),
                                     cbar_str="line stress\n[N/µm]")

    # saving the the plot
    fig.savefig(os.path.join(path, frame + "border_stress_img.png"), dpi=200)
    plt.close(fig)
    # adding figure to database
    except_error(db.setImage, IntegrityError, filename=os.path.join(path, frame + "border_stress_img.png"),
                 layer="FEM_output", path=1, sort_index=frames_ref_dict[frame])

    #saving the fem solution
    np.save(os.path.join(path, frame + "MSM_sol_defo.npy"), UG_sol)

    return None, frame





def get_db_info_for_analysis(db):

    all_frames, file_order, frames_ref_dict = get_file_order_and_frames(db)
    path = db.getPath(id=1).path
    im_shape = db.getImage(0).data.shape
    db_info={ "file_order":file_order,
            "frames_ref_dict":frames_ref_dict,
            "path":path,
            "im_shape":im_shape
            }
    return db_info,all_frames


def apply_to_frames(db, parameter_dict, analysis_function,res_dict,frames=[],db_info=None):
    '''
    wrapper to apply analysis function on all frames
    :param db: clcikpoints database
    :param parameter_dict: parameters for piv deforamtion calcualtion: (windowsize, overlapp), sigma and youngs modulus
    of the gel (or of the cell sheet when applying FEM), hight of the gel
    :param func: function that is analyzed
    :param frames: list f frames (of the cdb database) to be analyze e.g [0,1,2]
    :param db_info: dicitionary with the keys "path","frames_ref_dict","im_shape","file_order" as constructed from
    get db_info_for_analysis
    :param res_dict: ictionary of results to be filled up adn appended
    :return:
    '''


    if not isinstance(db_info,dict):
        db_info, all_frames=get_db_info_for_analysis(db)
    frames=make_iterable(frames) # if only one frame as string
    print(calculation_messages_multiple[analysis_function.__name__] % str(frames))
    create_layers_on_demand(db, layer_name_dict[analysis_function.__name__])
    for frame in tqdm(frames,total=len(frames)):
        try:
            analysis_function(frame, parameter_dict,res_dict, db=db,**db_info, single=False)
        except (Mask_Error,FileNotFoundError,IndexError): ## IndexError is only supposed to be there for a short time!!!!!
            pass
    return res_dict

if __name__=="__main__":
    ## setting upnecessary paramteres
    db=clickpoints.DataFile("/media/user/GINA1-BK/data_traktion_force_microscopy/WT_vs_KO_images_10_09_2019/wt_vs_ko_images_Analyzed/KOshift/database3.cdb","r")
    parameter_dict = default_parameters()
    analysis_function=FEM_analysis
    res_dict=defaultdict(dict)
    db_info, all_frames = get_db_info_for_analysis(db)
    apply_to_frames(db, parameter_dict, get_contractility_contractile_energy,res_dict, frames=all_frames)
    write_output_file(res_dict, "results", "/media/user/GINA1-BK/data_traktion_force_microscopy/WT_vs_KO_images_10_09_2019/wt_vs_ko_images_Analyzed/WTshift/out_test.txt")
    # calculating the deformation field and adding to data base
'''
def add_deforamtion_and_traction_force(db, parameter_dict):
    layer_list = ["def_plots", "traction_plots", "test_layer"]
    # check if layer already exist and create otherwise

    # correct ordering of files from image anotations

    # plotting
    dpi = 200
    fig1 = show_quiver_clickpoints(u, v, filter=[0, int(np.ceil(u.shape[0] / 40))], scale_ratio=0.2,
                                   headwidth=3, headlength=3, width=0.002,
                                   figsize=(im1.shape[1] / dpi, im1.shape[0] / dpi), cbar_str="deformation\n[pixel]")
    fig2 = show_quiver_clickpoints(tx_h, ty_h, filter=[0, int(np.ceil(u.shape[0] / 40))], scale_ratio=0.2,
                                   headwidth=3, headlength=3, width=0.002,
                                   figsize=(im1.shape[1] / dpi, im1.shape[0] / dpi), cbar_str="tracktion\n[Pa]")
    fig3 = show_quiver_clickpoints(tx, ty, filter=[0, int(np.ceil(u.shape[0] / 40))], scale_ratio=0.2,
                                   headwidth=3, headlength=3, width=0.002,
                                   figsize=(im1.shape[1] / dpi, im1.shape[0] / dpi), cbar_str="tracktion\n[Pa]")
    # saving plots
    fig1.savefig(os.path.join(path, frame + "deformation.png"), dpi=200)
    fig2.savefig(os.path.join(path, frame + "traction.png"), dpi=200)
    fig3.savefig(os.path.join(path, frame + "traction_test.png"), dpi=200)
    # closing figure objects
    plt.close(fig1)
    plt.close(fig2)
    plt.close(fig3)
    # adding plots to data base
    # if entry already exist dont change anything (image is overwritten, but database entry stays the same)
    except_error(db.setImage, IntegrityError, filename=frame + "deformation.png", layer="def_plots", path=1,
                 sort_index=int(frame) - 1)
    except_error(db.setImage, IntegrityError, filename=frame + "traction.png", layer="traction_plots", path=1,
                 sort_index=int(frame) - 1)
    except_error(db.setImage, IntegrityError, filename=frame + "traction_test.png", layer="test_layer", path=1,
                 sort_index=int(frame) - 1)

    # saving raw files
    np.save(os.path.join(path, frame + "u.npy"), u)
    np.save(os.path.join(path, frame + "v.npy"), v)
    np.save(os.path.join(path, frame + "tx_h.npy"), tx_h)
    np.save(os.path.join(path, frame + "ty_h.npy"), ty_h)
'''
