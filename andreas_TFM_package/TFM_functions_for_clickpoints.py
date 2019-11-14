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


class Mask_Error(Exception):
    pass
class ShapeMismatchError(Exception):
    pass


def write_output_file(values,value_type, file_path,with_prefix=False,new_file=False):
    if new_file:
        if os.path.exists(file_path): # try some other out names when one already exists
            for i in range(100000):
                file_path=os.path.join(os.path.split(file_path)[0],"out"+str(i)+".txt")
                if not os.path.exists(file_path):
                    break
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
                    warn_empty = (warn != "")
                    f.write(frame + "\t" + name + "\t" + str(round_flexible(res_unpack)) + "\t" + units[
                        name] + "\t" * warn_empty + warn + "\n")
    return file_path

def except_error(func, error,print_error=True, **kwargs):  # take functino and qkwarks
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
        if print_error:
            print(e)
        return False
    return values

def check_shape(x,y):
    s1=np.array(x.shape)
    s2=np.array(y.shape)
    if not all(s1==s2):
        raise ShapeMismatchError("shape of input arrays is unequal. Try recalculating the corresponding arrays.")


def try_mask_load(db,frame,raise_error=True,mtype="cell colony"):
    try:
        mask = db.getMask(frame=frame, layer=1).data
        # extract only one type of mask
        if mtype in ["membrane", "cell type1"]:  #
            mask = mask == 1
        if mtype in ["contractillity_colony", "cell type2"]:
            mask = mask == 2

    # raising error if no mask object in clickpoints exist

    except AttributeError:
        if raise_error:
            raise Mask_Error("no mask found in frame %s for type %s" % (str(frame), mtype))
        else:
            return None
    return mask


def warn_small_FEM_area(mask_area,threshold):
    warn=""
    area=np.sum(mask_area)
    if area<threshold:
        warnings.warn("FEM grid is very small (%d pixel). Consider increasing resolution of deformation and traction field."%area)
        warn="small FEM grid"
    return warn


def check_small_or_empty_mask(mask,frame, mtype,warn_thresh,raise_error):
    # checking if mask is empty
    warn=""
    if np.sum(mask)==0:
        if raise_error:
            raise Mask_Error("mask empty for mask type %s found in frame %s " % (mtype,str(frame)))
    # checking if mask is suspiciously small
    elif isinstance(warn_thresh,(int,float)):
        if np.sum(binary_fill_holes(mask))<warn_thresh:
            print("mask for %s is very small"%mtype)
            warn= "selected area is very small"
    return warn


def load_mask(db,frame,raise_error=True,mtype="cell colony",warn_thresh=300,fill_holes=False):
    '''
    loading the mask from clickpoints and checking if the size is reasonable. Returns warnings or errors.

    :param db:
    :param frame:
    :param raise_error raises an error if true, else returns None
    :param mtype:
    :return:
    '''
    warn=""
    # check if amsk type is supported
    if mtype not in ["cell type1","cell type2","membrane","contractillity_colony"]:
        raise Mask_Error("unsupported mask type")

    # trying to load the mask from clickpoints
    mask=try_mask_load(db,frame,raise_error=True,mtype=mtype)

    #check if mask is completey empty or suspiciously small
    warn = check_small_or_empty_mask(mask,frame, mtype,warn_thresh,raise_error)
    # if everything is alright warn is empty string

    # optional filling holes in the mask, only makes sense for some mask type
    if mtype in ["membrane","contractillity_colony"] and fill_holes:
        mask=binary_fill_holes(mask)
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





def create_layers_on_demand(db,db_info, layer_list):

    '''
    :param db: clickpointsw database
    :param layer_list: list of layer names that should be created
    :return:
    '''
    layer_list=make_iterable(layer_list)
    if any([l not in db_info["layers"] for l in layer_list]):
        base_layer = db.getLayer(id=1)
        for pl in layer_list:
            if pl not in db_info["layers"]:
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
    layers=[l.name for l in db.getLayers()]
    try:
        path=db.getOption("folder")
    except:
        path = db.getPath(id=1).path
    if path==".": # if empty path object in clickpoints use the path where clickpoints is saved
        path=os.path.split(db._database_filename)[0]
    im_shapes={} #exact list of image shapes
    for frame in all_frames:
        im_shapes[frame]=db.getImage(frame=frames_ref_dict[frame],layer="membranes").data.shape

    mask_types=[m.name for m in db.getMaskTypes()] # list mask types
    db_info = {"file_order": file_order,
               "frames_ref_dict": frames_ref_dict,
               "path": path,
               "im_shape": im_shapes,
               "mask_types":mask_types,
               "layers":layers
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
    str_frame = np.unique([get_group(re.search('(\d{1,4})', x), 1) for x in annotations])
    return str_frame[0]

def cut_mask_from_edge(mask,cut_factor):
    sum_mask1=np.sum(mask)
    dims=mask.shape
    inds=[int(dims[0]*cut_factor),int(dims[0]-(dims[0]*cut_factor)),int(dims[1]*cut_factor),int(dims[1]-(dims[1]*cut_factor))]
    mask[:inds[0], :] = 0
    mask[inds[1]:, :] = 0
    mask[:, :inds[2]] = 0
    mask[:, inds[3]:] = 0

    sum_mask2 = np.sum(mask)
    if sum_mask2<sum_mask1:
        warn="mask was cut becuase close to image edge"
    else:
        warn=""
    return mask, warn



def add_plot(plot_type, values,plot_function,frame,db_info,default_fig_parameters,parameter_dict,db):
    #values: values (args) that are needed as input for the plotting functions

    values=make_iterable_args(values) # returns list if input is not a list or other iterable

    if plot_type in default_fig_parameters["plots"][parameter_dict["FEM_mode"]]:  # checking if this should be plotted
        layer = default_fig_parameters["plots_layers"][plot_type]
        file_name=default_fig_parameters["file_names"][plot_type]

        create_layers_on_demand(db, db_info, [layer])
        plt.ioff()
        dpi = 200
        fig_parameters = set_fig_parameters(db_info["defo_shape"], db_info["im_shape"][frame], dpi,
                                            default_fig_parameters,
                                            figtype=plot_type)
        fig = plot_function(*values, **fig_parameters)

        # saving the the plot
        fig.savefig(os.path.join(db_info["path"], frame + file_name), dpi=200)
        plt.close(fig)

        # adding the plot to the database
        except_error(db.setImage, IntegrityError, print_error=True, filename=frame + file_name,
                     layer=layer, path=1,sort_index=db_info["frames_ref_dict"][frame])

def general_properties(frame, parameter_dict,res_dict, db,db_info=None,
                                         **kwargs):
    '''
    Number of cells, area of the cell colony...
    :param frame:
    :param parameter_dict:
    :param res_dict:
    :param db:
    :return:
    '''
    mtypes = [m for m in parameter_dict["area masks"] if m in db_info["mask_types"]]
    cut = True if parameter_dict["FEM_mode"] == "cell layer" else False  # cut close to image edges
    sum_on_area(mtypes,frame,res_dict,parameter_dict, db,db_info, label="area", sumtype="area",cut=cut)

    if not "defo_shape" in db_info.keys():
        u, v = try_to_load_deformation(db_info["path"], frame, warn=False)  # just to get correct shape
        db_info["defo_shape"] = u.shape

    ps_new = parameter_dict["pixelsize"] * np.mean(np.array(db_info["im_shape"][frame]) / np.array(db_info["defo_shape"]))


    # loads an calculates exact area and cellcount for mask cell colony if it exists
    # loading the mask, raising errors if not present or empty, warning if to small, also filling holes for some mask types
    if "membrane" in db_info["mask_types"]:
        mask_membrane, warn = load_mask(db, db_info["frames_ref_dict"][frame], mtype="membrane",
                                               warn_thresh=1500/(ps_new/parameter_dict["pixelsize"]),fill_holes=False)

        mask_area, borders = prepare_mask(mask_membrane,db_info["defo_shape"],min_cell_size=500)
        n_cells=len(borders.cell_ids)
        res_dict[frame]["colony n_cells"]=[n_cells,warn]


def sum_on_area(masks,frame,res_dict,parameter_dict, db,db_info,label,x=None,y=None,sumtype="abs",cut=True):

    mtypes=make_iterable(masks)
    for mtype in mtypes:
        mask_membrane, warn = load_mask(db, db_info["frames_ref_dict"][frame], mtype=mtype,
                                               warn_thresh=1500,fill_holes=True)
        if cut:
            mask_membrane,warn_edge=cut_mask_from_edge(mask_membrane,parameter_dict["edge_padding"])
        label2=default_parameters["mask_labels"][mtype]
        if sumtype=="abs":
            mask_int = interpolation(mask_membrane, dims=x.shape, min_cell_size=100)
            res_dict[frame]["%s on %s"%(label,label2)]= np.sum(np.sqrt(x[mask_int] ** 2 + y[mask_int] ** 2))
        if sumtype=="mean":
            mask_int = interpolation(mask_membrane, dims=x.shape, min_cell_size=100)
            res_dict[frame]["%s on %s" % (label, label2)] = np.mean(x[mask_int])
        if sumtype=="area": # area of original mask, without interpolation
            area = np.sum(mask_membrane) * ((parameter_dict["pixelsize"] * 10 ** -6) ** 2)
            res_dict[frame]["%s of %s" % (label, label2)] = area



def deformation(frame, parameter_dict,res_dict, db,db_info=None,**kwargs):


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

    # adding plot of derformation field to the database
    add_plot("deformation", (u,v),show_quiver_clickpoints,frame,db_info,default_fig_parameters,parameter_dict,db)

    # saving raw files
    np.save(os.path.join(db_info["path"], frame + "u.npy"), u)
    np.save(os.path.join(db_info["path"], frame + "v.npy"), v)
    # summing deformation over certain areas
    mtypes=[m for m in db_info["mask_types"] if m in parameter_dict["area masks"]]
    cut = True if parameter_dict["FEM_mode"]=="cell layer" else False # cut close to image edges
    sum_on_area(mtypes,frame,res_dict,parameter_dict, db,db_info,label="sum deformations",x=u,y=v,sumtype="abs",cut=cut)

    return None, frame




def get_contractillity_contractile_energy(frame, parameter_dict,res_dict, db,db_info=None,
                                          **kwargs):

    u, v = try_to_load_deformation(db_info["path"], frame, warn=True)
    t_x, t_y = try_to_load_traction(db_info["path"], frame, warn=False)
    db_info["defo_shape"]=t_x.shape
    ps_new = parameter_dict["pixelsize"] * np.mean(np.array(db_info["im_shape"][frame]) / np.array(t_x.shape))

    # select mask
    mtypes=[m for m in db_info["mask_types"] if m in ["cell type1","cell type2","contractillity_colony"]]

    if isinstance(u, np.ndarray):
        energy_points = contractile_energy_points(u, v, t_x, t_y, parameter_dict["pixelsize"], ps_new)  # contractile energy at any point
        # plotting contractile energy (only happens if enable in default_fig_parameters
        add_plot("energy_points",energy_points,show_map_clickpoints,frame,db_info,default_fig_parameters,parameter_dict,db)

        # iterating though mask that are selected for summation
    for mtype in mtypes:
        contractile_force=None
        contr_energy = None
        mask,warn=load_mask(db,db_info["frames_ref_dict"][frame],mtype=mtype,warn_thresh=1000,fill_holes=True)

        # removing some fraction of the mask close to the border
        # interpolation to size of traction force array
        mask_int = interpolation(mask, t_x.shape)

        mask_int,warn_edge=cut_mask_from_edge(mask_int,parameter_dict["edge_padding"])
        # calculate contractillity only in "colony" mode
        if mtype=="contractillity_colony":
            warn+=" "*(len(warn_edge)>0) +warn_edge
            contractile_force, proj_x, proj_y,center=contractillity(t_x, t_y, ps_new, mask_int)
            res_dict[frame]["contractillity on " + default_parameters["mask_labels"][mtype]] = [contractile_force, warn]
        # calculate contractile energy if deformations are provided
        if isinstance(u,np.ndarray):
            check_shape(u, t_x)
            contr_energy = np.sum(energy_points[mask_int.astype(bool)])  # sum of contractile energy on on mask

            res_dict[frame]["contractile energy on " + default_parameters["mask_labels"][mtype]] = [contr_energy, warn]
        print("contractile energy=",round_flexible(contr_energy),"contractillity=",round_flexible(contractile_force))


    return (contractile_force, contr_energy), frame




def traction_force(frame, parameter_dict,res_dict, db, db_info=None,**kwargs):

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

    # add a plot of the trackitoon filed to the database
    add_plot("traction", (tx,ty),show_quiver_clickpoints,frame,db_info,default_fig_parameters,parameter_dict,db)

    # saving raw files
    np.save(os.path.join(db_info["path"], frame + "tx.npy"), tx)
    np.save(os.path.join(db_info["path"], frame + "ty.npy"), ty)

    mtypes = [m for m in db_info["mask_types"] if m in parameter_dict["area masks"]]
    cut = True if parameter_dict["FEM_mode"] == "cell layer" else False  # cut close to image edges
    sum_on_area(mtypes,frame,res_dict,parameter_dict, db,db_info, label="sum traction forces", x=tx, y=ty, sumtype="abs",cut=cut)
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
    #mask, warn = load_mask(db, db_info["frames_ref_dict"][frame], raise_error=False, mtype="cell colony",
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
    mask, warn = load_mask(db, db_info["frames_ref_dict"][frame],raise_error=True, mtype="membrane",
                                           warn_thresh=1500 / (ps_new/parameter_dict["pixelsize"]),fill_holes=False)

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

    plot_fields(nodes, fields=[S_nodes[:, 0], S_nodes[:, 1], S_nodes[:, 2]], dims=mask_area.shape,
                titles=["x_stress", "y_stress", "xy_stress"], cbar_str="stress in N/pixel", origin="upper",
                mask=mask_area)  # ,mask_overlay=mask_int)
    return  UG_sol,stress_tensor





def FEM_analysis_average_stresses(frame,res_dict,parameter_dict, db,db_info,stress_tensor,ps_new,borders=None,**kwargs):

    # analyzing the FEM results with average stresses
    shear=stress_tensor[:,:,0,1] # shear component of the stress tensor
    mean_normal_stress =(stress_tensor[:,:,0,0]+stress_tensor[:,:,1,1])/2 # mean normal component of the stress tensor
    shear=shear/(ps_new*10**-6)# conversion to N/m
    mean_normal_stress=mean_normal_stress/(ps_new*10**-6)# conversion to N/m
    add_plot("stress_map", mean_normal_stress, show_map_clickpoints, frame, db_info, default_fig_parameters,
             parameter_dict,db)

    if parameter_dict["FEM_mode"]=="cell layer":
        mtypes=[m for m in db_info["mask_types"] if m in parameter_dict["area masks"]]
        sum_on_area(mtypes,frame,res_dict,parameter_dict, db,db_info, "average normal stress", x=mean_normal_stress, sumtype="mean",cut=True)
        sum_on_area(mtypes,frame,res_dict,parameter_dict, db,db_info, "average shear stress", x=shear, y=None, sumtype="mean",cut=True)


    if parameter_dict["FEM_mode"]=="colony":
        res_dict[frame]["average normal stress colony"] = np.mean(mean_normal_stress[borders.mask_area])
        res_dict[frame]["average shear stress colony"] = np.mean(shear[borders.mask_area])

    ### other possible stress measures, just for a nice picture
    #sigma_max, sigma_min, tau_max, phi_n, phi_shear, sigma_avg = all_stress_measures(S_nodes, nodes,
     #                                                                                dims=mask_area.shape)
    #sigma_max_abs = np.maximum(np.abs(sigma_min), np.abs(sigma_max))  ### highest possible norm of the stress tensor

def FEM_analysis_borders(frame, res_dict, db,db_info,parameter_dict, stress_tensor, ps_new, warn, borders=None,
                             **kwargs):

    # retrieving spline representation of bourders

    lines_splines = borders.lines_splines
    lines_points = borders.lines_points
    # plot lines tresses over border as continous curves:
    lines_interpol, min_v, max_v = interpolation_for_stress_and_normal_vector(lines_splines, lines_points,
                                                                              stress_tensor, pixel_length=ps_new,
                                                                              interpol_factor=1)
    # plotting the stress at cell borders
    add_plot("FEM_borders", (borders.mask_area.shape,borders.edge_lines, lines_interpol, min_v, max_v), plot_continous_boundary_stresses, frame,
             db_info, default_fig_parameters, parameter_dict, db)

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



    return None, frame

def FEM_full_analysis(frame, parameter_dict,res_dict, db, db_info=None, **kwargs):
    # performing full MSM/finite elements analysis
    #wrapper to flexibly perform FEM analysis



    if parameter_dict["FEM_mode"]=="colony":
        nodes, elements, loads, mats, mask_area, ps_new, warn, borders = FEM_setup_colony(frame, parameter_dict, db,
                                                                                          db_info=db_info, **kwargs)
        UG_sol, stress_tensor = FEM_simulation(nodes, elements, loads, mats, mask_area,parameter_dict)
        np.save(os.path.join(db_info["path"], frame + "stress_tensor.npy"), stress_tensor)

        FEM_analysis_average_stresses(frame,res_dict,parameter_dict, db,db_info,stress_tensor,ps_new,borders)
        FEM_analysis_borders(frame, res_dict, db,db_info,parameter_dict, stress_tensor, ps_new, warn, borders=borders,
                             **kwargs)

    if parameter_dict["FEM_mode"]=="cell layer":
        nodes, elements, loads, mats, mask_area, ps_new, warn, borders = FEM_setup_cell_layer(frame, parameter_dict, db,
                                                                               db_info=db_info, **kwargs)
        UG_sol, stress_tensor = FEM_simulation(nodes, elements, loads, mats, mask_area, parameter_dict,frame=frame)
        np.save(os.path.join(db_info["path"], frame + "stress_tensor.npy"), stress_tensor)
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
    for frame in tqdm(frames,total=len(frames)):
        try:
            analysis_function(frame, parameter_dict,res_dict, db=db,db_info=db_info,**kwargs)
        except Exception as e:
            if type(e) in (Mask_Error,FileNotFoundError,FindingBorderError,ShapeMismatchError):
                print(e)
            else:
                raise(e)

    return res_dict


### code to work on clickpoint outside of the addon
if __name__=="__main__":
    ## setting up necessary paramteres
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
    #apply_to_frames(db, parameter_dict, deformation,res_dict, frames="04",db_info=db_info)
    #apply_to_frames(db, parameter_dict, traction_force, res_dict, frames="04", db_info=db_info)
    #apply_to_frames(db, parameter_dict, FEM_full_analysis, res_dict, frames="01", db_info=db_info)
    #apply_to_frames(db, parameter_dict, traction_force, res_dict, frames="02", db_info=db_info)
    #apply_to_frames(db, parameter_dict, FEM_full_analysis, res_dict, frames="01", db_info=db_info)
    apply_to_frames(db, parameter_dict, get_contractillity_contractile_energy, res_dict, frames=all_frames, db_info=db_info)

    #write_output_file(res_dict, "results", "/media/user/GINA1-BK/data_traktion_force_microscopy/WT_vs_KO_images_10_09_2019/wt_vs_ko_images_Analyzed/WTshift/out_test.txt")
    # calculating the deformation field and adding to data base
