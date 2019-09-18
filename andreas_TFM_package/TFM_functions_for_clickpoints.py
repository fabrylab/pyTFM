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


# default parameters for calculation
default_parameters={"sigma":0.49, # poison ratio
"young":49000, # youngs modulus
"pixelsize":0.201, #pixelsize of the image with beads in  µm/pixel
"window_size":20,  # in µm
"overlapp":100, # set bigger then window_size/2
"std_factor":85,  # additional filter for extreme values in deformation field
"h":300, # hight of the substrate in µm
"TFM_mode":"finite_thikness",  # mode of traction force microscopy ("finite_thikness" or "infinite_thikness")
}

default_fig_parameters={
"cbar_str": {"deformation":"deformation\n[pixel]","tracktion":"tracktion\n[Pa]","FEM":"line stress\n[N/µm]"},  # label of the colorbar
"cmap": "rainbow",  # colormap for displaying magnitudes
"vmin": {"deformation":None,"tracktion":None,"FEM":None},  # minimal value displayed in the colormap
"vmax": {"deformation":None,"tracktion":None,"FEM":None},  # maximal value displayed in the colormap
"cbar_width": "2%",  # width of the color bar in % of the main image. Must be string with % at the end.
"cbar_height": "50%",  # height of the color bar in % of the main image. Must be string with % at the end.
"cbar_borderpad": 4,  # distance between the edge of the image and the colorbar (in pixels???)
"scale_ratio": 0.2,  # scale arrows so that the longest arrow is "maximum image dimension"*"scale ratio" long
"headwidth": 3,  # width of the arrow heads in pixels?
"headlength": 3,  # length of the arrow heads in pixels?
"width": 0.002,  # width of the arrow shaft /what unit?
"plot_t_vecs":{"FEM":False}, # plotting the stress vectors on the cell borders
"plot_n_arrows":{"FEM":False}, # plotting normal vectors on cell borders
"linewidth":{"FEM":4}, # linewidth when plotting the cell border stresses
"cm_cmap":{"FEM":cm.jet}, # color map for plotting the cell border stresses. Needs a color maps object
"border_arrow_filter":{"FEM":1}, #plot only every n'th arrow for  border stresses
#"filter":[0,4],
#"figsize":(10,10),
}

# default parameters for plotting
def set_fig_parameters(shape, fig_shape,dpi, default_fig_parameters,figtype):
    fig_parameters = {
        # filtering: 1.minimal length of arrow, 2. draw only every n'th arrow (in x and y direction)
        "filter": [0, int(np.ceil(shape[0] / 50))],
        # figsize, so that saving the figure with dpi=dpi, gives an image of the shape fig_shape[0]
        # used to match the other images in the database
        "figsize": (fig_shape[1] / dpi, fig_shape[0] / dpi),
    }
    for key,value in default_fig_parameters.items(): #adding, or potentially overwriting other paramters
        if isinstance(value,dict): # check if values for multiple plot types exist
            if figtype in value.keys():
                fig_parameters[key] = value[figtype]
        else:
            fig_parameters[key]=value


    return fig_parameters





# preset names of layers
layer_name_dict=defaultdict(list)
layer_name_dict["deformation"]="def_plots"
layer_name_dict["tracktion_force"]="traction_plots"
layer_name_dict["FEM_analysis"]="FEM_output"

# some message to be printed
calculation_messages=defaultdict(lambda:"%s")
calculation_messages["deformation"]="calculating deformation on frames %s"
calculation_messages["tracktion_force"]="calculating tracktion_forces on frames %s"
calculation_messages["FEM_analysis"]="FEM_analysis on frames %s"
calculation_messages["get_contractility_contractile_energy"]="contractility/contractile energy on frames %s"
calculation_messages["general_properites"]="getting colony properties on frames %s"


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
"avarage normal stress":"N/m",
"avarage shear stress":"N/m",
"area":"m2",
"n_cells":"",
"sum deformations":"pixels",
"sum traction forces":"N/m2",
"sum deformations on cell colony":"pixels"
}








class Mask_Error(Exception):
    pass


def write_output_file(values,value_type, file_path,with_prefix=False):

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
                    res_unpack, warn = unpack_list(res)
                    warn_empty = (warn == "")
                    f.write(frame + "\t" + name + "\t" + str(round_flexible(res_unpack)) + "\t" + units[
                        name] + "\t" * warn_empty + warn + "\n")




def except_error(func, error, **kwargs):  # take functino and qkwarks
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


def try_to_load_mask(db,frame,mtype="cell colony",warn_thresh=300):
    '''
    loading the mask from clickpoints and checking if the size is reasonable. Returns warnings or errors.

    :param db:
    :param frame:
    :param mtype:
    :return:
    '''
    # trying to load the mask from clickpoints
    try:
        mask = db.getMask(frame=frame,layer=1).data
        # extract only one type of mask
        if mtype=="cell colony": #
            mask=mask==1
        if mtype=="contractility selection":
            mask = mask == 2
    # raising error if no mask object in clickpoints exist
    except AttributeError:
        raise Mask_Error("no mask of the cell membrane found for frame " + str(frame))

    # checking if mask is empty
    if np.sum(mask)==0:
        raise Mask_Error("mask_empty for frame " + str(frame))
    # checking if mask is susiciously small
    elif np.sum(binary_fill_holes(mask))<warn_thresh:
        print("mask for %s is very small"%mtype)
        warn= "selected area is very small"
    # if everything is alright warn is empty string
    else: warn=""

    return mask, warn


def try_to_load_deformation(path, frame, warn=False):
    '''
    loading the deformations fro a given frame. If deformations are not found either raises an error
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
    '''
    loading the tracktions from a given frame. If tracktions are not found either raises an error
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



def setup_database_for_tfm(folder, name, return_db=False):

    db = clickpoints.DataFile(os.path.join(folder,name), "w")
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

    # string to search for frame: allows max 4 numbers at the beginning
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

def get_db_info_for_analysis(db,unpack=False):

    all_frames, file_order, frames_ref_dict = get_file_order_and_frames(db)
    path = db.getPath(id=1).path
    im_shape = db.getImage(0).data.shape
    db_info = {"file_order": file_order,
               "frames_ref_dict": frames_ref_dict,
               "path": path,
               "im_shape": im_shape
               }

    return db_info,all_frames


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



def general_properties(frame, parameter_dict,res_dict, db,db_info=None,
                                         single=True,**kwargs):
    '''
    Number of cells, area of the cell colony...
    :param frame:
    :param parameter_dict:
    :param res_dict:
    :param db:
    :return:
    '''
    if single:
        db_info, all_frames = get_db_info_for_analysis(db)
        print(calculation_messages["general_properites"]%frame)

    t_x, t_y = try_to_load_traction(db_info["path"], frame, warn=False) # just to get correct shape
    ps_new = parameter_dict["pixelsize"] * np.mean(np.array(db_info["im_shape"]) / np.array(t_x.shape))

    mask_membrane, warn = try_to_load_mask(db, db_info["frames_ref_dict"][frame], mtype="cell colony",
                                           warn_thresh=1500/(ps_new/parameter_dict["pixelsize"]))
    area = np.sum(binary_fill_holes(mask_membrane)) * ((parameter_dict["pixelsize"] * 10 ** -6) ** 2)
    mask_int_mem = interpolation(mask_membrane, t_x.shape, min_cell_size=100)
    mask_area, mask_boundaries, borders = prepare_mask(mask_int_mem,min_cell_size=5)
    n_cells=len(borders.cell_ids)
    res_dict[frame]["area"]=[area,warn]
    res_dict[frame]["n_cells"]=[n_cells,warn]





def deformation(frame, parameter_dict,res_dict, db,db_info=None,single=True,**kwargs):

    if single:
        db_info,all_frames=get_db_info_for_analysis()
        create_layers_on_demand(db, ["def_plots"])
        print(calculation_messages["deformation"]%frame)


    # deformation for 1 frame
    im1 = db.getImage(id=db_info["file_order"][frame + "after"]).data  ## thats very slow though
    im2 = db.getImage(id=db_info["file_order"][frame + "before"]).data
    u, v, x, y, mask, mask_std = calculate_deformation(im1.astype(np.int32), im2.astype(np.int32),
                                                       parameter_dict["window_size"]
                                                       , parameter_dict["overlapp"],
                                                       std_factor=parameter_dict["std_factor"])
    res_dict[frame]["sum deformations"] = np.sum(np.sqrt(u ** 2 + v** 2))


    # plotting
    plt.ioff()
    dpi = 200

    fig_parameters = set_fig_parameters(u.shape, db_info["im_shape"], dpi,default_fig_parameters,figtype="deformation")
    fig1 = show_quiver_clickpoints(u, v,**fig_parameters)
    # saving plots
    fig1.savefig(os.path.join(db_info["path"], frame + "deformation.png"), dpi=200)
    # closing figure objects
    plt.close(fig1)
    # adding plots to data base
    # if entry already exist dont change anything (image is overwritten, but database entry stays the same)
    except_error(db.setImage, IntegrityError, filename=frame + "deformation.png", layer="def_plots", path=1,
                 sort_index=db_info["frames_ref_dict"][frame])
    # saving raw files
    np.save(os.path.join(db_info["path"], frame + "u.npy"), u)
    np.save(os.path.join(db_info["path"], frame + "v.npy"), v)

    # sum of deformations on the cellcolony area
    ps_new = parameter_dict["pixelsize"] * np.mean(np.array(db_info["im_shape"]) / np.array(u.shape))
    mask_membrane, warn = try_to_load_mask(db, db_info["frames_ref_dict"][frame], mtype="cell colony",
                                           warn_thresh=1500 / (ps_new / parameter_dict["pixelsize"]))
    mask_int = interpolation(binary_fill_holes(mask_membrane), dims=u.shape, min_cell_size=100)
    res_dict[frame]["sum deformations on cell colony"] = np.sum(np.sqrt(u[mask_int] ** 2 + v[mask_int] ** 2))
    return None, frame





def get_contractility_contractile_energy(frame, parameter_dict,res_dict, db,db_info=None,
                                         single=True,**kwargs):

    if single:
        db_info, all_frames = get_db_info_for_analysis(db)
        print(calculation_messages["get_contractility_contractile_energy"]%frame)

    mask,warn=try_to_load_mask(db,db_info["frames_ref_dict"][frame],mtype="contractility selection",warn_thresh=1000)
    u,v=try_to_load_deformation(db_info["path"], frame, warn=True)
    t_x,t_y=try_to_load_traction(db_info["path"], frame, warn=False)

    # filling holes
    mask=binary_fill_holes(mask)
    ps_new = parameter_dict["pixelsize"] * np.mean(np.array(mask.shape) / np.array(t_x.shape))
    # interpolation to size of traction force array
    mask_int = interpolation(mask, t_x.shape)
    contractile_force, proj_x, proj_y,center=contractility(t_x, t_y, ps_new, mask_int)
    # calculate contractile energy if deformations are provided
    if isinstance(u,np.ndarray):
        contr_energy=contractile_energy(u,v,t_x,t_y,parameter_dict["pixelsize"], ps_new,mask_int)
    else:
        contr_energy=None
    print("contractile energy=",round_flexible(contr_energy),"contractillity=",round_flexible(contractile_force))

    # writing results to dictionray

    res_dict[frame]["contractility"]=[contractile_force,warn]
    res_dict[frame]["contractile energy"]= [contr_energy,warn]

    return (contractile_force, contr_energy), frame




def tracktion_force(frame, parameter_dict,res_dict, db, db_info=None, single=True,**kwargs):

    if single:
        db_info, all_frames = get_db_info_for_analysis(db)
        create_layers_on_demand(db, ["traction_plots"])
        print(calculation_messages["tracktion_force"]%frame)

    # trying to laod deformation
    u,v=try_to_load_deformation(db_info["path"], frame, warn=False)
    ps_new = parameter_dict["pixelsize"] * np.mean(  # should be equivalent to "pixelsize_def_image"
        np.array(u.shape) / np.array(u.shape))  # pixelsize of fem grid in µm
    # using tfm with or without finite thikness correction
    if parameter_dict["TFM_mode"] == "finite_thikness":
        tx, ty = ffttc_traction_finite_thickness(u, v, pixelsize1=parameter_dict["pixelsize"],
                                                     pixelsize2=ps_new,
                                                     h=parameter_dict["h"], young=parameter_dict["young"],
                                                     sigma=parameter_dict["sigma"],
                                                     filter="gaussian")
        # unit is N/m**2
        if np.isnan(tx).all():
            warnings.warn("falling back to inifinte thikness assumption due to nummerical issues")
            parameter_dict["TFM_mode"] = "infinite_thikness"
    if parameter_dict["TFM_mode"] == "infinite_thikness":
        tx, ty = ffttc_traction(u, v, pixelsize1=parameter_dict["pixelsize"],
                                                 pixelsize2=ps_new,
                                                 young=parameter_dict["young"],
                                                 sigma=parameter_dict["sigma"],
                                                 filter="gaussian")
    res_dict[frame]["sum traction forces"]=np.sum(np.sqrt(tx**2+ty**2))
    # plotting
    plt.ioff()
    dpi = 200
    fig_parameters = set_fig_parameters(u.shape, db_info["im_shape"], dpi,default_fig_parameters,figtype="tracktion")
    fig2 = show_quiver_clickpoints(tx, ty, **fig_parameters)
    # saving plots
    fig2.savefig(os.path.join(db_info["path"], frame + "traction.png"), dpi=dpi)
    plt.close(fig2) # closing figure objects
    # adding plots to data base
    # if entry already exist don't change anything (image is overwritten, but database entry stays the same)
    except_error(db.setImage, IntegrityError, filename=frame + "traction.png", layer="traction_plots", path=1,
                 sort_index=db_info["frames_ref_dict"][frame])
    # saving raw files
    np.save(os.path.join(db_info["path"], frame + "tx.npy"), tx)
    np.save(os.path.join(db_info["path"], frame + "ty.npy"), ty)
    return None, frame




def FEM_analysis(frame,parameter_dict,res_dict,db,single=True,db_info=None,**kwargs):

    # performing MSM/finite elements analysis
    if single:
        db_info, all_frames = get_db_info_for_analysis(db)
        create_layers_on_demand(db, ["FEM_output"])
        print(calculation_messages["FEM_analysis"]%frame)

    # trying to traction forces, raise error if not found
    t_x, t_y = try_to_load_traction(db_info["path"], frame, warn=False)
    ps_new = parameter_dict["pixelsize"] * np.mean(
        np.array(db_info["im_shape"]) / np.array(t_x.shape))  # pixelsize of fem grid in µm
    # trying to load mask, raise error if not found
    mask, warn = try_to_load_mask(db, db_info["frames_ref_dict"][frame], mtype="cell colony",
                                           warn_thresh=1500 / (ps_new/parameter_dict["pixelsize"]))

    # interpolation to size of traction force array
    mask_int = interpolation(mask, t_x.shape,min_cell_size=100) # this cell size in pixels of beads image

    # further preparation of mask data
    mask_area, mask_boundaries, borders = prepare_mask(mask_int, min_cell_size=5) # min_cell_size at least 2 or none
    # plt.figure();plt.imshow(mask_area)

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
    fig_parameters = set_fig_parameters(mask_area.shape, db_info["im_shape"], dpi,default_fig_parameters,figtype="FEM")
    fig=plot_continous_boundary_stresses(mask_area.shape, borders.edge_lines, lines_interpol, min_v,
                                     max_v,**fig_parameters)

    # saving the the plot
    fig.savefig(os.path.join(db_info["path"], frame + "border_stress_img.png"), dpi=200)
    plt.close(fig)
    # adding figure to database
    except_error(db.setImage, IntegrityError, filename=os.path.join(db_info["path"], frame + "border_stress_img.png"),
                 layer="FEM_output", path=1, sort_index=db_info["frames_ref_dict"][frame])

    #saving the fem solution (deformations on the cell sheet)
    np.save(os.path.join(db_info["path"], frame + "MSM_sol_defo.npy"), UG_sol)

    return None, frame



def apply_to_frames(db, parameter_dict, analysis_function,res_dict,frames=[],db_info=None):
    '''
    wrapper to apply analysis function on all frames
    :param db: clickpoints database
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
    print(calculation_messages[analysis_function.__name__] % str(frames))
    create_layers_on_demand(db, layer_name_dict[analysis_function.__name__])
    for frame in tqdm(frames,total=len(frames)):
        try:
            analysis_function(frame, parameter_dict,res_dict, db=db,db_info=db_info, single=False)
        except (Mask_Error,FileNotFoundError,IndexError): ## IndexError is only supposed to be there for a short time!!!!!
            pass
    return res_dict


### code to work on clickpoint outside of the addon
if __name__=="__main__":
    ## setting upnecessary paramteres
    db=clickpoints.DataFile("/media/user/GINA1-BK/data_traktion_force_microscopy/WT_vs_KO_images_10_09_2019/wt_vs_ko_images_Analyzed/WTshift/database3.cdb","r")
    parameter_dict = default_parameters
    res_dict=defaultdict(dict)
    db_info, all_frames = get_db_info_for_analysis(db)
    apply_to_frames(db, parameter_dict, FEM_analysis,res_dict, frames="01",db_info=db_info)
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
