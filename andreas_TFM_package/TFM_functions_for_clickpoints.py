### function integrating Traktion force microscopy into a clcikpoints database

from andreas_TFM_package.grid_setup_solids_py import *
from andreas_TFM_package.functions_for_cell_colonie import *
from andreas_TFM_package.solids_py_stress_functions import *
from andreas_TFM_package.utilities_TFM import *
from andreas_TFM_package.TFM_functions import *
from andreas_TFM_package.parameters_and_strings import *
import solidspy.postprocesor as pos
import solidspy.assemutil as ass
import solidspy.solutil as sol
from peewee import IntegrityError
import clickpoints
import os
import re
import warnings
from tqdm import tqdm
import itertools





# preset names of layers
layer_name_dict=defaultdict(list)
layer_name_dict["deformation"]="def_plots"
layer_name_dict["traction_force"]="traction_plots"
layer_name_dict["FEM_analysis"]="FEM_output"

class Mask_Error(Exception):
    pass


def write_output_file(values,value_type, file_path,with_prefix=False):

    if value_type=="parameters":
        with open(file_path, "w+") as f:
            f.write("analysis_paramters\n")
            for parameter, value in values.items():
                if parameter!="mask_labels":
                    f.write(parameter + "\t" + str(value)+ "\n")
    if value_type == "results":
        # list of available frames sorted
        frames=list(values.keys())
        frames.sort(key=lambda x:int(x))

        with open(file_path, "a+") as f:
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

    try:
        values = func(**kwargs)
    except error as e:
        print(e)
        return False
    return values





def try_to_load_mask(db,frame,raise_error=True,mtype="cell colony",warn_thresh=300):
    '''
    loading the mask from clickpoints and checking if the size is reasonable. Returns warnings or errors.

    :param db:
    :param frame:
    :param raise_error raises an error if true, else returns None
    :param mtype:
    :return:
    '''
    warn=""

    if mtype not in ["cell type1","cel type2","membrane","contractillity_colony"]:
        raise Mask_Error("unsupported mask type")

    # trying to load the mask from clickpoints
    if mtype not in [m.name for m in db.getMaskTypes()]:
        if raise_error:
            raise Mask_Error("no mask for mask type %s found" % (mtype))
        else:
            return None, ""

    try:
        mask = db.getMask(frame=frame,layer=1).data
        # extract only one type of mask
        if mtype in ["membrane","cell type1"] : #
            mask = mask == 1
        if mtype in ["contractillity_colony","cell type2"]:
            mask = mask == 2

    # raising error if no mask object in clickpoints exist

    except AttributeError:
        if raise_error:
            raise Mask_Error("no mask found in frame %s " % (str(frame)))
        else:
            return None,""

    # checking if mask is empty
    if np.sum(mask)==0:
        if raise_error:
            raise Mask_Error("mask empty for mask type %s found in frame %s " % (mtype,str(frame)))
    # checking if mask is suspiciously small
    elif isinstance(warn_thresh,(int,float)):
        if np.sum(binary_fill_holes(mask))<warn_thresh:
            print("mask for %s is very small"%mtype)
            warn= "selected area is very small"
    # if everything is alright warn is empty string

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
    loading the tractions from a given frame. If tractions are not found either raises an error
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


def warn_incorrect_files(frames):
    '''
    throws a waring when it  more or less then three images per frame are found.
    :param frames:
    :return:
    '''
    frames = np.array(frames)
    unique_frames, counts = np.unique(frames, return_counts=True)
    problems=np.where(counts!=3)[0]
    if len(problems)>0:
        warn="There seems to be a problem with the your images."
        for p_id in problems:
            warn+="Found %s files for frame %s. " %(counts[p_id],unique_frames[p_id])
        warnings.warn(warn+"Excpeted only three files per frame.")

def warn_small_FEM_area(mask_area,threshold):
    warn=""
    area=np.sum(mask_area)
    if area<threshold:
        warnings.warn("FEM grid is very small (%d pixel). Consider increasing resolution of deformation and traction field."%area)
        warn="small FEM grid"
    return warn



def setup_masks(db,db_info,parameter_dict):
   #mask_type=[m.name for m in db.getMaskTypes()]
    if parameter_dict["FEM_mode"]=="colony":
        db.deleteMaskTypes("cell type1")
        db.deleteMaskTypes("cell type2")
        if "membrane" not in db_info["mask_types"]:
            db.setMaskType("membrane", color="#99EA44",index=1)
        if "contractillity_colony" not in db_info["mask_types"]:
            db.setMaskType("contractillity_colony", color="#ff0000",index=2)
    if parameter_dict["FEM_mode"] == "cell layer":
        db.deleteMaskTypes("membrane", color="#99EA44", index=1)
        db.deleteMaskTypes("contractillity_colony", color="#ff0000", index=2)
        if "cell type1"  not in db_info["mask_types"]:
            db.setMaskType("cell type1", color="#1322ff", index=1)
        if "cell type2" not in db_info["mask_types"]:
            db.setMaskType("cell type2",color="#ebff05",index=2)

def guess_TFM_mode(db_info,parameter_dict):

    cl_cond1= len(db_info["mask_types"])==2 and  "cell type1" in db_info["mask_types"] and  "cell type2" in db_info["mask_types"]
    cl_cond2= len(db_info["mask_types"])==1 and  ("cell type1" in db_info["mask_types"] or  "cell type2" in db_info["mask_types"])

    co_cond1 = len(db_info["mask_types"]) == 2 and "membrane" in db_info["mask_types"] and "contractillity_colony" in db_info[
        "mask_types"]
    co_cond2 = len(db_info["mask_types"]) == 1 and (
                "membrane" in db_info["mask_types"] or "contractillity_colony" in db_info["mask_types"])

    cond_empty=len(db_info["mask_types"]) == 0

    undertermined=False
    if cl_cond1 or cl_cond2:
        mode = "cell layer"
    elif co_cond1 or  co_cond2:
        mode = "colony"
    elif cond_empty:
        mode=parameter_dict["FEM_mode"]
    else:
        warnings.warn("failed to guess analysis mode. Try to select it manually")
        mode=parameter_dict["FEM_mode"]
        undertermined=True

    return mode,undertermined



def setup_database_for_tfm(folder, name, return_db=False):

    '''
    Sorting images into a clickpoints database. Frames are identified by leading numbers. Layers are identified by
    the file name.
    :param folder: Folder where images are searched.
    :param name: Name of the database. Needs to end with .cdb.
    :param return_db: Choose weather function returns the database object, or weather the connection to the
    database is closed
    :return:
    '''

    # creating a new cdb database, will override an existing one.
    db = clickpoints.DataFile(os.path.join(folder,name), "w")
    # regex patterns to sort files into layers. If any of these matches, the file will  be sorted into a layer.
    # keys: name of the layer, values: list of regex patterns
    file_endings = "(png|jpg|tif|swg)" # all allowed file endings
    layer_search = {"images_after": [re.compile("\d{1,4}after.*" + file_endings)],
                    "images_before": [re.compile("\d{1,4}before.*" + file_endings)],
                    "membranes": [re.compile("\d{1,4}mebrane.*" + file_endings),
                                          re.compile("\d{1,4}bf_before.*" + file_endings)]
                            }
    # filtering all files in the folder
    all_patterns=list(itertools.chain(*layer_search.values()))
    images = [x for x in os.listdir(folder) if any([pat.match(x) for pat in all_patterns])]
    # identifying frames by evaluating the leading number.
    frames = [get_group(re.search('(\d{1,4})', x), 0) for x in images] # extracting frame
    # generating a list of sort_ids for the clickpoints database (allows you to miss some frames)
    sort_id_list=make_rank_list(frames,dtype=int)# list of sort indexes (frames) of images in the database
    warn_incorrect_files(frames) # checking if there where more or less then three images per frame

    # initializing layer in the database
    layer_list = ["images_after", "images_before","membranes","def_plots","traction_plots","FEM_output"]  # change this..
    base_layer = db.getLayer(layer_list[0], create=True, id=0)
    for l in layer_list[1:]:
        db.getLayer(l, base_layer=base_layer, create=True)
    path = db.setPath(folder, 1)  # setting the path to the images

    # sorting images into layers
    for id, (sort_index_id,frame, im) in enumerate(zip(sort_id_list,frames, images)):
        if any([pat.match(im) for pat in layer_search["images_after"]]):
            db.setImage(id=id, filename=im, sort_index=sort_index_id
                        , layer="images_after", path=1)
            db.setAnnotation(filename=im, comment=frame + "after__"+str(sort_index_id)+"sid")
            # setting new layer and sorting in at the same time
        if any([pat.match(im) for pat in layer_search["images_before"]]):
            db.setImage(id=id, filename=im
                        , sort_index=sort_index_id, layer="images_before", path=1)
            db.setAnnotation(filename=im, comment=frame + "before__"+str(sort_index_id)+"sid")
        if any([pat.match(im) for pat in layer_search["membranes"]]):
            db.setImage(id=id, filename=im, sort_index=sort_index_id,
                        layer="membranes", path=1)
            db.setAnnotation(filename=im, comment=frame+"bf__"+str(sort_index_id)+"sid")
    # delete_empty_layers(db) # not necessary

    # return database connection or close
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



def get_db_info_for_analysis(db):

    all_frames, file_order, frames_ref_dict = get_file_order_and_frames(db)
    path = db.getPath(id=1).path
    im_shapes={} #exact list of image shapes
    for frame in all_frames:
        im_shapes[frame]=db.getImage(frame=frames_ref_dict[frame],layer="membranes").data.shape

    mask_types=[m.name for m in db.getMaskTypes()] # list mask types
    db_info = {"file_order": file_order,
               "frames_ref_dict": frames_ref_dict,
               "path": path,
               "im_shape": im_shapes,
               "mask_types":mask_types
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



def fill_patches_for_cell_layer(frame, parameter_dict,res_dict, db,db_info=None,single=True,**kwargs):
    # trying to load the mask from clickpoints
    print(db_info["frames_ref_dict"][frame])
    try:
        image=db.getImage(frame=db_info["frames_ref_dict"][frame]) #this is a work around, dont know why i cant use getMask directly
        mask = db.getMask(image).data
    # raising error if no mask object in clickpoints exist
    except AttributeError:
            raise Mask_Error("no mask of the cell membrane found for frame " + str(frame))
    mask_part = mask == 1  # type of "cell type1"

    ### lets not talk about this...
    labels=label(mask_part,background=-1)

    edge_labels=np.hstack([labels[0, :], labels[-1, :],labels[:, 0], labels[:, -1]]).astype(int)
    edge_mask_values=np.hstack([mask_part[0, :], mask_part[-1, :],mask_part[:, 0], mask_part[:, -1]]).astype(int)
    background_patches_edge=np.unique(edge_labels[edge_mask_values==0])
    fore_ground_patches_edge=np.unique(edge_labels[edge_mask_values==1])
    found_labels=np.concatenate([background_patches_edge,fore_ground_patches_edge])
    neigbouring=[]
    for la in fore_ground_patches_edge:
        expansion=np.logical_and(binary_dilation(labels == la), ~(labels == la))
        neigbouring.append(np.unique(labels[expansion]))
    neigbouring=np.concatenate(neigbouring)
    neigbouring=np.array([n for n in neigbouring if n not in found_labels])

    neigbouring2=[]
    for la in neigbouring:
        expansion = np.logical_and(binary_dilation(labels == la), ~(labels == la))
        neigbouring2.append(np.unique(labels[expansion]))
    neigbouring2 = np.concatenate(neigbouring2)
    neigbouring2 = np.array([n for n in neigbouring2 if n not in neigbouring and n not in found_labels])

    all_labels=np.concatenate([neigbouring,neigbouring2,fore_ground_patches_edge])
    mask_new = np.zeros_like(labels)
    for la in all_labels:
        mask_new[labels==la]=1
    mask_new=mask_new.astype(bool)
    mask[mask_new]=1
    mask[~mask_new] = 2
    mask=mask.astype(np.uint8)
    db.setMask(image=image,data=mask) # udapting mask in clickpoints



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

    mtypes = [m for m in ["cell type1", "cell type2"] if m in parameter_dict["area masks"] and m in db_info["mask_types"]]
    sum_on_area(mtypes,frame,res_dict,parameter_dict, db,db_info, label="area", sumtype="area")

    u, v = try_to_load_deformation(db_info["path"], frame, warn=False)  # just to get correct shape
    ps_new = parameter_dict["pixelsize"] * np.mean(np.array(db_info["im_shape"][frame]) / np.array(u.shape))


    # loads an calculates exact area and cell count for mask cell colony if it exists
    mask_membrane, warn = try_to_load_mask(db, db_info["frames_ref_dict"][frame], mtype="membrane",
                                           warn_thresh=1500/(ps_new/parameter_dict["pixelsize"]))

    area = np.sum(binary_fill_holes(mask_membrane)) * ((parameter_dict["pixelsize"] * 10 ** -6) ** 2)
    res_dict[frame]["area on colony"] = [area, warn]

    mask_area, borders = prepare_mask(mask_membrane,u.shape,min_cell_size=500)
    n_cells=len(borders.cell_ids)
    res_dict[frame]["colony n_cells"]=[n_cells,warn]



def sum_on_area(masks,frame,res_dict,parameter_dict, db,db_info,label,x=None,y=None,sumtype="abs",):

    mtypes=make_iterable(masks)
    for mtype in mtypes:
        mask_membrane, warn = try_to_load_mask(db, db_info["frames_ref_dict"][frame], mtype=mtype,
                                               warn_thresh=1500)
        mask_membrane=cut_mask_from_edge(mask_membrane,parameter_dict["edge_padding"])
        label2=default_parameters["mask_labels"][mtype]
        if sumtype=="abs":
            mask_int = interpolation(binary_fill_holes(mask_membrane), dims=x.shape, min_cell_size=100)
            res_dict[frame]["%s on %s"%(label,label2)]= np.sum(np.sqrt(x[mask_int] ** 2 + y[mask_int] ** 2))
        if sumtype=="mean":
            mask_int = interpolation(binary_fill_holes(mask_membrane), dims=x.shape, min_cell_size=100)
            res_dict[frame]["%s on %s" % (label, label2)] = np.mean(x[mask_int])
        if sumtype=="area": # area of original mask, without interpolation
            area = np.sum(binary_fill_holes(mask_membrane)) * ((parameter_dict["pixelsize"] * 10 ** -6) ** 2)
            res_dict[frame]["%s of %s" % (label, label2)] = area



def deformation(frame, parameter_dict,res_dict, db,db_info=None,single=True,**kwargs):

    if single:
        db_info,all_frames=get_db_info_for_analysis()
        create_layers_on_demand(db, ["def_plots"])
        print(calculation_messages["deformation"]%frame)


    # deformation for 1 frame
    im1 = db.getImage(id=db_info["file_order"][frame + "after"]).data  ## thats very slow though
    im2 = db.getImage(id=db_info["file_order"][frame + "before"]).data

    # overlapp and windowsize in pixels
    window_size_pix=int(np.ceil(parameter_dict["window_size"] / parameter_dict["pixelsize"]))
    overlapp_pix=int(np.ceil(parameter_dict["overlapp"] / parameter_dict["pixelsize"]))
    u, v, x, y, mask, mask_std = calculate_deformation(im1.astype(np.int32), im2.astype(np.int32),
                                                      window_size_pix, overlapp_pix,
                                                       std_factor=parameter_dict["std_factor"])
    db_info["defo_shape"]=u.shape
    res_dict[frame]["sum deformations"] = np.sum(np.sqrt(u ** 2 + v** 2))


    # plotting
    plt.ioff()
    dpi = 200
    fig_parameters = set_fig_parameters(u.shape, db_info["im_shape"][frame], dpi,default_fig_parameters,figtype="deformation")
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
    # summing deformation over certain areas
    mtypes=[m for m in db_info["mask_types"] if m in parameter_dict["area masks"]]
    sum_on_area(mtypes,frame,res_dict,parameter_dict, db,db_info,label="sum of deformations",x=u,y=v,sumtype="abs")

    return None, frame


def make_contractile_energy_plot(u,v,t_x,t_y,frame, parameter_dict,res_dict, db,db_info):

    ps_new = parameter_dict["pixelsize"] * np.mean(np.array(db_info["im_shape"][frame]) / np.array(db_info["defo_shape"]))
    pixelsize1 = parameter_dict["pixelsize"] * 10 ** -6  # conversion to m
    pixelsize2 = ps_new * 10 ** -6
    energy_points = 0.5 * (pixelsize2 ** 2) * (np.sqrt((t_x * u * pixelsize1) ** 2 + (
            t_y * v * pixelsize1) ** 2))

    bg = np.percentile(energy_points, 30)  # value of a background point
    energy_points -= bg
    plt.ioff()
    dpi = 200
    fig_parameters = set_fig_parameters(db_info["defo_shape"], db_info["im_shape"][frame], dpi,
                                        default_fig_parameters,
                                        figtype="energy_points")
    fig = show_map_clickpoints(energy_points, **fig_parameters)

    # saving the the plot
    fig.savefig(os.path.join(db_info["path"], frame + "energy_distribution.png"), dpi=200)
    plt.close(fig)
    # adding figure to database
    except_error(db.setImage, IntegrityError,
                 filename=os.path.join(db_info["path"], frame + "energy_distribution.png"),
                 layer="FEM_output", path=1, sort_index=db_info["frames_ref_dict"][frame])


def cut_mask_from_edge(mask,cut_factor):
    dims=mask.shape
    inds=[int(dims[0]*cut_factor),int(dims[0]-(dims[0]*cut_factor)),int(dims[1]*cut_factor),int(dims[1]-(dims[1]*cut_factor))]
    mask[:inds[0], :] = 0
    mask[inds[1]:, :] = 0
    mask[:, inds[2]] = 0
    mask[:, inds[3]:] = 0
    return mask


def get_contractillity_contractile_energy(frame, parameter_dict,res_dict, db,db_info=None,
                                         single=True, **kwargs):

    if single:
        db_info, all_frames = get_db_info_for_analysis(db)
        print(calculation_messages["get_contractillity_contractile_energy"]%frame)


    u, v = try_to_load_deformation(db_info["path"], frame, warn=True)
    t_x, t_y = try_to_load_traction(db_info["path"], frame, warn=False)
    db_info["defo_shape"]=t_x.shape

    # select mask and plotting according to FEM type
    if parameter_dict["FEM_mode"]=="cell layer":
        mtypes = ["cell type1", "cell type2"]
        #total contractile energy as a map
        make_contractile_energy_plot(u, v, t_x, t_y, frame, parameter_dict, res_dict, db, db_info)

    else:
        mtypes=["contractillity_colony"]

    # iterating though mask that are selected for summation
    for mtype in mtypes:
        contractile_force=None
        contr_energy = None
        mask,warn=try_to_load_mask(db,db_info["frames_ref_dict"][frame],mtype=mtype,warn_thresh=1000)
        ps_new = parameter_dict["pixelsize"] * np.mean(np.array(mask.shape) / np.array(t_x.shape))
        # removing some fraction of the mask close to the border
        mask=cut_mask_from_edge(mask,parameter_dict["edge_padding"])
            # filling holes
        mask = binary_fill_holes(mask)
        # interpolation to size of traction force array
        mask_int = interpolation(mask, t_x.shape)
        # calculate contractillity only in "colony" mode
        if mtype=="contractillity_colony":
            contractile_force, proj_x, proj_y,center=contractillity(t_x, t_y, ps_new, mask_int)
            res_dict[frame]["contractillity on " + default_parameters["mask_labels"][mtype]] = [contractile_force, warn]
        # calculate contractile energy if deformations are provided
        if isinstance(u,np.ndarray):
            contr_energy=contractile_energy(u,v,t_x,t_y,parameter_dict["pixelsize"], ps_new,mask_int)
            res_dict[frame]["contractile energy on " + default_parameters["mask_labels"][mtype]] = [contr_energy, warn]
        print("contractile energy=",round_flexible(contr_energy),"contractillity=",round_flexible(contractile_force))


    return (contractile_force, contr_energy), frame




def traction_force(frame, parameter_dict,res_dict, db, db_info=None, single=True,**kwargs):

    if single:
        db_info, all_frames = get_db_info_for_analysis(db)
        create_layers_on_demand(db, ["traction_plots"])
        print(calculation_messages["traction_force"]%frame)

    # trying to laod deformation
    u,v=try_to_load_deformation(db_info["path"], frame, warn=False)
    db_info["defo_shape"] = u.shape
    ps_new = parameter_dict["pixelsize"] * np.mean(  # should be equivalent to "pixelsize_def_image"
        np.array(u.shape) / np.array(u.shape))  # pixelsize of fem grid in µm
    # using tfm with or without finite thickness correction
    if parameter_dict["TFM_mode"] == "finite_thickness":
        tx, ty = ffttc_traction_finite_thickness(u, v, pixelsize1=parameter_dict["pixelsize"],
                                                     pixelsize2=ps_new,
                                                     h=parameter_dict["h"], young=parameter_dict["young"],
                                                     sigma=parameter_dict["sigma"],
                                                     filter="gaussian")
        # unit is N/m**2
        if np.isnan(tx).all():
            warnings.warn("falling back to inifinte thickness assumption due to nummerical issues")
            parameter_dict["TFM_mode"] = "infinite_thickness"
    if parameter_dict["TFM_mode"] == "infinite_thickness":
        tx, ty = ffttc_traction(u, v, pixelsize1=parameter_dict["pixelsize"],
                                                 pixelsize2=ps_new,
                                                 young=parameter_dict["young"],
                                                 sigma=parameter_dict["sigma"],
                                                 filter="gaussian")

    # plotting
    plt.ioff()
    dpi = 200
    fig_parameters = set_fig_parameters(u.shape, db_info["im_shape"][frame], dpi,default_fig_parameters,figtype="traction")
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

    mtypes = [m for m in db_info["mask_types"] if m in parameter_dict["area masks"]]
    sum_on_area(mtypes,frame,res_dict,parameter_dict, db,db_info, label="sum of traction forces", x=tx, y=ty, sumtype="abs")
    return None, frame


def FEM_setup_cell_layer(frame,parameter_dict,db,db_info=None,**kwargs):
    '''

    :param frame:
    :param parameter_dict:
    :param db:
    :param db_info:
    :param kwargs:
    :return:
    '''
    t_x, t_y = try_to_load_traction(db_info["path"], frame, warn=False)
    db_info["defo_shape"] = t_x.shape
    ps_new = parameter_dict["pixelsize"] * np.mean(
        np.array(db_info["im_shape"][frame]) / np.array(t_x.shape))  # pixelsize of fem grid in µm

    # preparing forces/ dont apply any further correction here
    f_x = t_x * ((ps_new * (10 ** -6)) ** 2)  # point force for each node from tractions
    f_y = t_y * ((ps_new * (10 ** -6)) ** 2)
    # use the full field of few for fem grid
    mask_area = np.ones(t_x.shape).astype(bool)
    # grid setup
    nodes, elements, loads, mats = grid_setup(mask_area, -f_x, -f_y, 1, sigma=0.5)  # note the negative signe

    # see if we can find any borders##
    #mask, warn = try_to_load_mask(db, db_info["frames_ref_dict"][frame], raise_error=False, mtype="cell colony",
    #                              warn_thresh=1500 / (ps_new / parameter_dict["pixelsize"]))
    #
    #
    warn,borders=None,None

    return nodes, elements, loads, mats,mask_area, ps_new, warn, borders


def FEM_setup_colony(frame,parameter_dict,db,db_info=None,**kwargs):

    # trying to traction forces, raise error if not found
    t_x, t_y = try_to_load_traction(db_info["path"], frame, warn=False)
    db_info["defo_shape"]=t_x.shape
    ps_new = parameter_dict["pixelsize"] * np.mean(
        np.array(db_info["im_shape"][frame]) / np.array(t_x.shape))  # pixelsize of fem grid in µm

    # preparing forces
    f_x = t_x * ((ps_new * (10 ** -6)) ** 2)  # point force for each node from tractions
    f_y = t_y * ((ps_new * (10 ** -6)) ** 2)

    # trying to load cell colony mask, raise error if not found
    mask, warn = try_to_load_mask(db, db_info["frames_ref_dict"][frame],raise_error=True, mtype="membrane",
                                           warn_thresh=1500 / (ps_new/parameter_dict["pixelsize"]))

    # preparation of mask data
    mask_area, borders = prepare_mask(mask,t_x.shape, min_cell_size=500) # min_cell_size at least 2 or none
    warn=warn_small_FEM_area(mask_area,threshold=1000)

   # coorecting force for torque and net force
    f_x[~mask_area] = np.nan  # setting all values outside of maske area to zero
    f_y[~mask_area] = np.nan
    f_x_c1 = f_x - np.nanmean(f_x)  # normalizing traction force to sum up to zero (no displacement)
    f_y_c1 = f_y - np.nanmean(f_y)
    f_x_c2, f_y_c2, p = correct_torque(f_x_c1, f_y_c1, mask_area)
    # get_torque1(f_y,f_x,mask_area)

    # setup of the grid
    nodes, elements, loads, mats = grid_setup(mask_area, -f_x_c2, -f_y_c2, 1, sigma=0.5) # note the negative signe

    return nodes, elements, loads, mats, mask_area, ps_new, warn, borders

def FEM_simulation(nodes, elements, loads, mats, mask_area, parameter_dict,**kwargs):

    DME, IBC, neq = ass.DME(nodes, elements)  # boundary conditions asembly??
    print("Number of elements: {}".format(elements.shape[0]))
    print("Number of equations: {}".format(neq))

    # System assembly
    KG = ass.assembler(elements, mats, nodes, neq, DME, sparse=True)
    RHSG = ass.loadasem(loads, IBC, neq)

    # System solution with custom conditions
    if parameter_dict["FEM_mode"]=="colony":
        # solver with constraints to zero translation and zero rotation
        UG_sol, rx = custom_solver(KG, RHSG, mask_area, verbose=False)

    # System solution with default solver
    if parameter_dict["FEM_mode"] == "cell layer":
        UG_sol = sol.static_sol(KG, RHSG)  # automatically detect sparce matrix

    if not (np.allclose(KG.dot(UG_sol) / KG.max(), RHSG / KG.max())):
        print("The system is not in equilibrium!")
    # average shear and normal stress on the colony area
    UC = pos.complete_disp(IBC, nodes, UG_sol)  # uc are x and y displacements
    E_nodes, S_nodes = pos.strain_nodes(nodes, elements, mats, UC)  # stresses and strains
    stress_tensor = calculate_stress_tensor(S_nodes, nodes, dims=mask_area.shape)  # assembling the stress tensor

    return  UG_sol,stress_tensor

def FEM_analysis_average_stresses(frame,res_dict,parameter_dict, db,db_info,stress_tensor,ps_new,borders=None,**kwargs):

    # analyzing the FEM results with average stresses
    shear=stress_tensor[:,:,0,1] # shear component of the stress tensor
    mean_normal_stress =(stress_tensor[:,:,0,0]+stress_tensor[:,:,1,1])/2 # mena normal component of the stress tensor
    shear=shear*ps_new*10**6# conversion to N/m
    mean_normal_stress=mean_normal_stress*ps_new*10**6# conversion to N/m

    if parameter_dict["FEM_mode"]=="cell layer":
        # plotting mean normal stress as a map
        #plt.ioff()
        #dpi = 200
        #fig_parameters = set_fig_parameters(db_info["defo_shape"],db_info["im_shape"][frame], dpi, default_fig_parameters,
        #                                   figtype="FEM_cell_layer")
        #fig = show_map_clickpoints(mean_normal_stress, **fig_parameters)

        # saving the the plot
        #fig.savefig(os.path.join(db_info["path"], frame + "avg_normal.png"), dpi=200)
        #plt.close(fig)
        # adding figure to database
        #except_error(db.setImage, IntegrityError,
        #             filename=os.path.join(db_info["path"], frame + "avg_normal.png"),
        #             layer="FEM_output", path=1, sort_index=db_info["frames_ref_dict"][frame])

        mtypes=[m for m in db_info["mask_types"] if m in parameter_dict["area masks"]]
        sum_on_area(mtypes,frame,res_dict,parameter_dict, db,db_info, "average normal stress", x=mean_normal_stress, sumtype="mean")
        sum_on_area(mtypes,frame,res_dict,parameter_dict, db,db_info, "average shear stress", x=shear, y=None, sumtype="mean")


    if parameter_dict["FEM_mode"]=="colony":
        res_dict[frame]["average normal stress colony"] = np.mean(mean_normal_stress[borders.mask_area])
        res_dict[frame]["average shear stress colony"] = np.mean(shear[borders.mask_area])

    ### other possible stress measures, just for a nice picture
    #sigma_max, sigma_min, tau_max, phi_n, phi_shear, sigma_avg = all_stress_measures(S_nodes, nodes,
     #                                                                                dims=mask_area.shape)
    #sigma_max_abs = np.maximum(np.abs(sigma_min), np.abs(sigma_max))  ### highest possible norm of the stress tensor

def FEM_analysis_borders(frame, res_dict, db,db_info, stress_tensor, ps_new, warn, mask_area, borders=None,
                             **kwargs):

    # retrieving spline representation of bourders

    lines_splines = borders.lines_splines
    lines_points = borders.lines_points
    # plot lines tresses over border as continous curves:
    lines_interpol, min_v, max_v = interpolation_for_stress_and_normal_vector(lines_splines, lines_points,
                                                                              stress_tensor, pixel_length=ps_new,
                                                                              interpol_factor=1)
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
    fig_parameters = set_fig_parameters(mask_area.shape, db_info["im_shape"][frame], dpi,default_fig_parameters,figtype="FEM")
    fig=plot_continous_boundary_stresses(mask_area.shape, borders.edge_lines, lines_interpol, min_v,
                                     max_v,**fig_parameters)

    # saving the the plot
    fig.savefig(os.path.join(db_info["path"], frame + "border_stress_img.png"), dpi=200)
    plt.close(fig)
    # adding figure to database
    except_error(db.setImage, IntegrityError, filename=os.path.join(db_info["path"], frame + "border_stress_img.png"),
                 layer="FEM_output", path=1, sort_index=db_info["frames_ref_dict"][frame])

    return None, frame

def FEM_full_analysis(frame, parameter_dict,res_dict, db, single=True, db_info=None, **kwargs):
    # performing full MSM/finite elements analysis
    if single:
        db_info, all_frames = get_db_info_for_analysis(db)
        create_layers_on_demand(db, ["FEM_output"])
        print(calculation_messages["FEM_analysis"] % frame)

    #wrapper to flexibly perform FEM analysis




    if parameter_dict["FEM_mode"]=="colony":
        nodes, elements, loads, mats, mask_area, ps_new, warn, borders = FEM_setup_colony(frame, parameter_dict, db,
                                                                                          db_info=db_info, **kwargs)
        UG_sol, stress_tensor = FEM_simulation(nodes, elements, loads, mats, mask_area,parameter_dict)
        np.save(os.path.join(db_info["path"], frame + "fem_sol.npy"), UG_sol)

        FEM_analysis_average_stresses(frame,res_dict,parameter_dict, db,db_info,stress_tensor,ps_new,borders)
        FEM_analysis_borders(frame, res_dict, db,db_info, stress_tensor, ps_new, warn, mask_area, borders=borders,
                             **kwargs)

    if parameter_dict["FEM_mode"]=="cell layer":
        nodes, elements, loads, mats, mask_area, ps_new, warn, borders = FEM_setup_cell_layer(frame, parameter_dict, db,
                                                                               db_info=db_info, **kwargs)
        UG_sol, stress_tensor = FEM_simulation(nodes, elements, loads, mats, mask_area, parameter_dict,frame=frame)
        np.save(os.path.join(db_info["path"], frame + "fem_sol.npy"), UG_sol)
        #if isinstance(borders,Cells_and_Lines):
        #    FEM_analysis_borders(frame, res_dict, db,db_info, stress_tensor, ps_new, warn, mask_area, borders=borders,
        #                          **kwargs)
        FEM_analysis_average_stresses(frame,res_dict,parameter_dict, db,db_info,stress_tensor,ps_new)





def apply_to_frames(db, parameter_dict, analysis_function,res_dict,frames=[],db_info=None,**kwargs):
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
            analysis_function(frame, parameter_dict,res_dict, db=db,db_info=db_info, single=False,**kwargs)
        except Exception as e:
            if type(e) in (Mask_Error,FileNotFoundError,FindingBorderError):
                print(e)
            else:
                raise(e)

    return res_dict


### code to work on clickpoint outside of the addon
if __name__=="__main__":
    ## setting upnecessary paramteres
    #db=clickpoints.DataFile("/home/user/Desktop/Monolayers_new_images/monolayers_new_images/KO_DC1_tomatoshift/database.cdb","r")
    db = clickpoints.DataFile(
        "/media/user/GINA1-BK/data_traktion_force_microscopy/WT_vs_KO_images_10_09_2019/wt_vs_ko_images_Analyzed/KOshift/database3.cdb", "r")

    parameter_dict = default_parameters
    res_dict=defaultdict(dict)
    db_info, all_frames = get_db_info_for_analysis(db)
    parameter_dict["overlapp"]=10
    parameter_dict["FEM_mode"] = "colony"
    #apply_to_frames(db, parameter_dict, deformation, res_dict, frames=all_frames, db_info=db_info)
    #apply_to_frames(db, parameter_dict, traction_force, res_dict, frames=all_frames, db_info=db_info)
    apply_to_frames(db, parameter_dict, deformation,res_dict, frames="01",db_info=db_info)
    apply_to_frames(db, parameter_dict, traction_force, res_dict, frames="01", db_info=db_info)
    apply_to_frames(db, parameter_dict, general_properties, res_dict, frames="01", db_info=db_info)
    apply_to_frames(db, parameter_dict, FEM_full_analysis, res_dict, frames="01", db_info=db_info)
    apply_to_frames(db, parameter_dict, get_contractillity_contractile_energy, res_dict, frames="01", db_info=db_info)

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
                                   figsize=(im1.shape[1] / dpi, im1.shape[0] / dpi), cbar_str="traction\n[Pa]")
    fig3 = show_quiver_clickpoints(tx, ty, filter=[0, int(np.ceil(u.shape[0] / 40))], scale_ratio=0.2,
                                   headwidth=3, headlength=3, width=0.002,
                                   figsize=(im1.shape[1] / dpi, im1.shape[0] / dpi), cbar_str="traction\n[Pa]")
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
