### function integrating Traktion force microscopy into a clcikpoints database

import os
import warnings

import clickpoints
from peewee import DoesNotExist
from peewee import IntegrityError
from pyTFM.TFM_functions import *
from pyTFM.frame_shift_correction import *
from pyTFM.grid_setup_solids_py import *
from pyTFM.parameters_and_strings import *
from pyTFM.plotting import *
from pyTFM.stress_functions import *
from pyTFM.utilities_TFM import *
from skimage.morphology import label
from tqdm import tqdm



class Mask_Error(Exception):
    pass


class ShapeMismatchError(Exception):
    pass


class cells_masks():
    def __init__(self, frames, db, db_info, parameter_dict):
        self.frames = []
        self.db = db
        self.db_info = db_info
        self.parameter_dict = parameter_dict
        self.indices = {m.name: m.index for m in self.db.getMaskTypes()}
        self._masks_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: None)))
        self._warns_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: "")))
        self.warn_thresh = 1000
        self.add_frames(frames, min_size=parameter_dict["min_obj_size"])

    # dictionaries properties because iterating through these dictionaries often adds empty entries,due to the
    # default dict structure. This will for example generate empty cell (cell colonies in a frame when trying to load a
    # mask that is not existing
    @property
    def masks_dict(self):
        return copy.deepcopy(self._masks_dict)

    @property
    def warns_dict(self):
        return copy.deepcopy(self._warns_dict)

    def add_frames(self, frames, min_size=1000):
        print("loading masks: ")
        frames = make_iterable(frames)
        for frame in tqdm(frames):
            # loading mask
            mask = try_mask_load(self.db, self.db_info["frames_ref_dict"][frame], raise_error=False, mtype="all")
            if not isinstance(mask, np.ndarray):  # checking if mask instance could be loaded
                for mask_name, index in self.indices.items():
                    self._masks_dict[frame][0][mask_name] = None  # writing to dict
                    self._warns_dict[frame][0][mask_name] = "no mask found"
                continue
            # optional cutting close to image edge
            mask_full = mask > 0  # mask as one block
            mask_full = binary_fill_holes(mask_full)  # filling the whole area
            mask_cut, warn_edge = cut_mask_from_edge_wrapper(self.parameter_dict["edge_padding"], mask,
                                                             self.parameter_dict, cut=True)  # just for the warning here
            mask_full = remove_small_objects(mask_full, min_size=min_size)  # cleanup
            mask_full = label(mask_full)  # finding individual objects
            regions = regionprops(mask_full)

            if len(regions) == 0:  # checking if mask is completely empty
                for mask_name, index in self.indices.items():
                    self._masks_dict[frame][0][mask_name] = None  # writing to dict
                    self._warns_dict[frame][0][mask_name] = "no mask found"
                    self._masks_dict[frame][0]["com"] = None
                continue

            for obj_id, r in enumerate(regions):  # iterating through all objects
                for mask_name, index in self.indices.items():  # iterating through all max types
                    shape = mask.shape
                    new_mask = np.zeros(shape)
                    new_mask[r.coords[:, 0], r.coords[:, 1]] = 1
                    coords_final = np.where(np.logical_and(mask == index, new_mask))
                    # extracting only one mask type in only one object
                    warn = check_mask_size(mask, self.warn_thresh, print_out=True, mask_name=mask_name, frame=frame)
                    self._masks_dict[frame][obj_id][mask_name] = (
                        coords_final, shape)  # writing to dict only the coordinates of true values
                    self._warns_dict[frame][obj_id][mask_name] = warn + " " * (warn != "") + warn_edge
                    self._masks_dict[frame][obj_id]["com"] = r.centroid

    def get_com_frame(self, frame):
        ret = []
        for cell_id in self.masks_dict[frame].keys():
            com = self.masks_dict[frame][cell_id]["com"]  # do i need to use the property here??
            ret.append((cell_id, com))
        return ret

    def reconstruct_mask(self, frame, cell_id, mtype, raise_error=True, fill_holes=False, cut_close_to_edge=True):
        # takes about 0.2 seconds for one type of mask
        # note: small objects have already been removed
        mask = None
        ind_shape = self.masks_dict[frame][cell_id][mtype]  # do i need to use the property here??
        if isinstance(ind_shape, tuple):
            indices = ind_shape[0]
            shape = ind_shape[1]
            if isinstance(indices, tuple) and isinstance(shape, tuple):
                if len(indices[0]) > 0:  # can also be empty tuple sometimes
                    mask = np.zeros(shape)
                    mask[indices] = 1
                    if fill_holes:
                        mask = binary_fill_holes(mask)
                    if cut_close_to_edge:
                        mask, warn_edge = cut_mask_from_edge_wrapper(self.parameter_dict["edge_padding"], mask,
                                                                     self.parameter_dict, cut=True)
                    if np.sum(mask) > 0:
                        return mask.astype(bool)
        else:
            if raise_error:
                raise Mask_Error("no mask found in frame %s for type %s" % (str(frame), mtype))
        return mask

    def reconstruct_masks_frame(self, frame, mtype="mask", raise_error=True, fill_holes=False, obj_ids=[], indices=None):
        # retrieves the a list of all masks and cell ids in one frame
        ret = []
        mtypes = make_iterable(mtype)
        for cell_id in self.masks_dict[frame].keys():
            if cell_id in obj_ids or len(obj_ids) == 0:  # only select certain objects
                for mtype in mtypes:
                    mask = self.reconstruct_mask(frame, cell_id, mtype, raise_error=raise_error, fill_holes=fill_holes)
                    warn = self.get_warning(frame, cell_id, mtype)
                    ret.append((cell_id, mask, mtype, warn))
        return ret

    def reconstruct_masks_frame_add(self, frame, mtypes, raise_error=True, fill_holes=False, obj_ids=[]):
        ret = self.reconstruct_masks_frame(frame, mtypes, raise_error=raise_error, fill_holes=fill_holes,
                                           obj_ids=obj_ids)
        shape = [mask.shape for cell_id, mask, mtype, warn in ret if isinstance(mask, np.ndarray)]
        if len(shape) == 0:
            return ret
        ret_dict = defaultdict(lambda: np.zeros(shape[0]).astype(bool))
        warn_dict = defaultdict(lambda: "")
        if len(ret) > 1:
            for cell_id, mask, mtype, warn in ret:
                if not isinstance(mask, np.ndarray):
                    continue
                ret_dict[cell_id] = np.logical_or(ret_dict[cell_id], mask)
                if warn not in warn_dict[cell_id]:
                    warn_dict[cell_id] += ", " + warn
        ret_final = [(cell_id, mask, "SUM", warn) for (cell_id, mask), warn in
                     zip(ret_dict.items(), warn_dict.values())]
        ret_final = ret if len(ret_final) == 0 else ret_final
        return ret_final

    def get_warning(self, frame, cell_id, mtype):
        return self.warns_dict[frame][cell_id][mtype]


def check_mask_size(mask, warn_tresh, print_out=True, mask_name="", frame=""):
    warn = "" if np.sum(mask.astype(bool)) > warn_tresh else "small mask"
    if print_out and warn != "":
        print(warn + " of type %s in %s" % (mask_name, frame))
    return warn


def write_output_file(values, value_type, file_path, new_file=False):
    if new_file:
        if os.path.exists(file_path):  # try some other out names when one already exists
            for i in range(100000):
                file_path = os.path.join(os.path.split(file_path)[0], "out" + str(i) + ".txt")
                if not os.path.exists(file_path):
                    break
    if value_type == "parameters":
        with open(file_path, "w+") as f:
            f.write("analysis_paramters\n")
            for parameter, value in values.items():
                if parameter not in ["mask_properties", "FEM_mode_id", "fig_parameters", "cv_pad"]:
                    f.write(parameter + "\t" + str(value) + "\n")
    if value_type == "results":
        # list of available frames sorted
        frames = list(values.keys())
        frames.sort(key=lambda x: int(x))  # extra sorting step
        with open(file_path, "a+") as f:
            for frame in frames:
                for name, res_list in values[frame].items():
                    for res_single in res_list:
                        cell_id, res, warn = res_single
                        f.write(
                            frame + "\t" + str(cell_id) + "\t" + name + "\t" + str(round_flexible(res)) + "\t" + units[
                                name] + "\t" * (warn != "") + warn + "\n")
    return file_path


def check_shape(x, y):
    s1 = np.array(x.shape)
    s2 = np.array(y.shape)
    if not all(s1 == s2):
        raise ShapeMismatchError("shape of input arrays is unequal. Try recalculating the corresponding arrays.")


def try_mask_load(db, frame, raise_error=True, mtype="cell colony", ret_type="None"):
    try:
        mask = db.getMask(frame=frame, layer=1).data
        # extract only one type of mask
        if not mtype == "all":
            index = db.getMaskType(mtype).index
            mask = mask == index
    # raising error if no mask object in clickpoints exist
    except AttributeError:
        if raise_error:
            raise Mask_Error("no mask found in frame %s for type %s" % (str(frame), mtype))
        else:
            if ret_type == "zeros":
                return np.zeros(db.getImage(frame=frame).data.shape)
            else:
                return None
    return mask


def warn_small_FEM_area(mask_area, threshold):
    warn = ""
    area = np.sum(mask_area)
    if area < threshold:
        warnings.warn(
            "FEM grid is very small (%d pixel). Consider increasing resolution of deformation and traction field." % area)
        warn = "small FEM grid"
    return warn




def check_small_or_empty_mask(mask, frame, mtype, warn_thresh=None, raise_error=True, add_str_error="",
                              add_str_warn=""):
    # checking if mask is empty
    warn = ""
    if np.sum(mask) == 0:
        if raise_error:
            raise Mask_Error("mask empty for mask type %s in frame %s" % (mtype, str(frame)) + " " + add_str_error)
    # checking if mask is suspiciously small
    elif isinstance(warn_thresh, (int, float)):
        if np.sum(binary_fill_holes(mask)) < warn_thresh:
            print("mask for %s is very small" % mtype + add_str_error)
            warn = "selected area is very small" + " " + add_str_error
    return warn


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


def create_layers_on_demand(db, db_info, layer_list):
    '''
    :param db: clickpointsw database
    :param layer_list: list of layer names that should be created
    :return:
    '''
    layer_list = make_iterable(layer_list)
    if any([l not in db_info["layers"] for l in layer_list]):
        base_layer = db.getLayer(id=1)
        for pl in layer_list:
            if pl not in db_info["layers"]:
                db.getLayer(pl, base_layer=base_layer, create=True)


def split_dict_str(string, convert_key=False, convert_value=False, **kwargs):
    string = string.strip("{").strip("}")
    string_list = string.split(",")
    dict_obj = {}
    for e in string_list:
        key, value = [sub_str.strip(" ") for sub_str in e.split(":")]
        if convert_key:
            key = try_int_strip(key)
        else:
            key = key.strip("'")
        if convert_value:
            value = try_int_strip(value)
        else:
            value = value.strip("'")
        dict_obj[key] = value

    return dict_obj


def split_list_str(string, **kwargs):
    # for example when original input was list
    string = string.strip("[").strip("]")
    string_list = string.split(", ")
    list_obj = [v.strip("'") for v in string_list]
    return list_obj


def get_option_wrapper(db, key, unpack_funct=None, empty_return=list, **kwargs):
    try:
        if unpack_funct:
            return unpack_funct(db.table_option.select().where(db.table_option.key == key).get().value, **kwargs)
        else:
            return db.table_option.select().where(db.table_option.key == key).get().value
    except:
        return empty_return()


def get_db_info_for_analysis(db):
    # list of all frames, order according to their sort index
    unique_frames = get_option_wrapper(db, "unique_frames", split_list_str)
    file_order = get_option_wrapper(db, "file_order", split_dict_str, convert_value=True)
    frames_ref_dict = get_option_wrapper(db, "frames_ref_dict", split_dict_str, empty_return=dict, convert_value=True)
    id_frame_dict = get_option_wrapper(db, "id_frame_dict", split_dict_str, empty_return=dict, convert_key=True)

    # inverse of frames_ref_dict
    cbd_frames_ref_dict = {value: key for key, value in frames_ref_dict.items()}

    layers = [l.name for l in db.getLayers()]
    try:
        path = get_option_wrapper(db, "folder", None)
    except:
        path = db.getPath(id=1).path
    # if empty path object in clickpoints use the path where clickpoints is saved
    if path == ".":
        path = os.path.split(db._database_filename)[0]

    # TODO: think about optimizing that// this esentially loads one image per frame each time it is called
    im_shapes = {}  # exact list of image shapes
    for frame in unique_frames:
        images = db.getImages(frame=frames_ref_dict[frame])
        if len(images) > 0:
            im_shapes[frame] = db.getImages(frame=frames_ref_dict[frame])[0].data.shape
        else:
            im_shapes[frame] = ()
    # list of exiting mask types
    mask_types = [m.name for m in db.getMaskTypes()]
    db_info = {"file_order": file_order,
               "frames_ref_dict": frames_ref_dict,
               "path": path,
               "im_shape": im_shapes,
               "mask_types": mask_types,
               "layers": layers,
               "unique_frames": unique_frames,
               "id_frame_dict": id_frame_dict,
               "cbd_frames_ref_dict": cbd_frames_ref_dict
               }
    return db_info, unique_frames


def write_image(db, layer, sort_index, filename):
    try:  # deleting all entries that might have been in this sot_index/layer slot
        prev_img = db.getImages(layer=layer, frame=sort_index)
        for im in prev_img:
            db.deleteImages(id=im.id)
    except DoesNotExist:
        pass
    # setting a new image
    except_error(db.setImage, IntegrityError, print_error=True, filename=filename,
                 layer=layer, path=1, sort_index=sort_index)


def add_plot(plot_type, values, plot_function, frame, db_info, parameter_dict, db):
    # values: values (args) that are needed as input for the plotting functions
    fig_parameters = parameter_dict["fig_parameters"]
    if plot_type in fig_parameters["plots"][parameter_dict["FEM_mode"]]:  # checking if this should be plotted
        layer = fig_parameters["plots_layers"][plot_type]
        file_name = fig_parameters["file_names"][plot_type]

        create_layers_on_demand(db, db_info, [layer])
        plt.ioff()
        fig_parameters_use = set_fig_parameters(db_info["defo_shape"], db_info["im_shape"][frame], fig_parameters,
                                                figtype=plot_type)
        fig, ax = plot_function(*values, **fig_parameters_use)

        # saving the the plot
        print("saving to " + os.path.join(db_info["path"], frame + file_name))
        fig.savefig(os.path.join(db_info["path"], frame + file_name), facecolor=fig.get_facecolor(),
                    edgecolor=fig.get_facecolor(), dpi=fig_parameters_use["resolution"])
        plt.close(fig)
        # adding the plot to the database
        write_image(db, layer=layer, sort_index=db_info["frames_ref_dict"][frame], filename=frame + file_name)

def check_empty_mask(mask, mtype="---", frame="---", obj_id="---", add_str_error="", raise_error=False):
    if not isinstance(mask, np.ndarray):
        wrong=True
    elif np.sum(mask) == 0:
        wrong = True
    else:
        wrong=False
    if wrong:
        string = "couldn't identify mask %s in frame %s patch %s" % (str(mtype), str(frame), str(obj_id))
        if raise_error:
            raise Mask_Error(string +  " " + add_str_error)
        else:
            print(string)
    return wrong

def calc_on_area(frame, res_dict, parameter_dict, label, masks, mask_types=None, obj_ids=[], x=None, y=None,
                 sumtype="abs", add_cut_factor=None, db_info=None, label2=None):
    fill_holes = True if parameter_dict["FEM_mode"] == "colony" else False
    mask_iter = masks.reconstruct_masks_frame(frame, mask_types, obj_ids=obj_ids, raise_error=False,
                                              fill_holes=fill_holes)
    for obj_id, mask, mtype, warn in mask_iter:
        if check_empty_mask(mask, mtype, frame, obj_id):
            continue
        if isinstance(add_cut_factor, float):  # additional cutting when in layer-mode
            mask, w = cut_mask_from_edge(mask, add_cut_factor, parameter_dict["TFM_mode"] == "colony")
        if check_empty_mask(mask, mtype, frame, obj_id):
            continue
        if label2 == None:
            _label2 = default_parameters["mask_properties"][mtype]["label"]
        else:
            _label2 = label2
        _label2 = _label2.strip()
        if sumtype == "abs":
            mask_int = interpolation(mask, dims=x.shape, min_cell_size=100)
            res_dict[frame]["%s %s" % (label, _label2)].append(
                [obj_id, np.sum(np.sqrt(x[mask_int] ** 2 + y[mask_int] ** 2)), warn])
        if sumtype == "mean":
            mask_int = interpolation(mask, dims=x.shape, min_cell_size=100)
            res_dict[frame]["%s %s" % (label, _label2)].append([obj_id, np.mean(x[mask_int]), warn])
        if sumtype == "mean_abs":
            mask_int = interpolation(mask, dims=x.shape, min_cell_size=100)
            res_dict[frame]["%s %s" % (label, _label2)].append([obj_id, np.mean(np.abs(x[mask_int])), warn])
        if sumtype == "area":  # area of original mask, without interpolation
            area = np.sum(mask) * ((parameter_dict["pixelsize"] * 10 ** -6) ** 2)
            res_dict[frame]["%s %s" % (label, _label2)].append([obj_id, area, warn])
        if sumtype == "cv":
            mask_int = interpolation(mask, dims=x.shape, min_cell_size=100)
            ps_new = parameter_dict["pixelsize"] * np.mean(np.array(db_info["im_shape"][frame]) / np.array(x.shape))
            border_pad = int(np.round(parameter_dict["cv_pad"] / ps_new))
            res_dict[frame]["%s %s" % (label, _label2)].append(
                [obj_id, coefficient_of_variation(mask_int, x, border_pad), warn])



def general_properties(frame, parameter_dict, res_dict, db, db_info=None, masks=None, **kwargs):

    '''
    Number of cells, area of the cell colony...
    :param frame:
    :param parameter_dict:
    :param res_dict:
    :param db:
    :return:
    '''

    # area of the cells/cell colony
    # retrieve which masks should be used to calculate the area
    use_type = "area_colony" if parameter_dict["FEM_mode"] == "colony" else "area_layer"
    mtypes = [m for m in db_info["mask_types"] if m in get_masks_by_key(default_parameters, "use", use_type)]
    # place holder for shape if not defined in defo-shape, needed for counting cells
    int_shape = db_info["defo_shape"] if "defo_shape" in db_info.keys() else (int(db_info["im_shape"][frame][0] * 0.2),
                                                                              int(db_info["im_shape"][frame][1] * 0.2))
    # calculate the area of each mask
    calc_on_area(frame, res_dict, parameter_dict, "area", masks, mask_types=mtypes, sumtype="area")
    # calculating the cell count for each colony
    mask_iter = masks.reconstruct_masks_frame(frame, "Cell Boundary", raise_error=True, fill_holes=False)
    for obj_id, mask_membrane, mtype, warn in mask_iter:
        if check_empty_mask(mask_membrane, mtype, frame, obj_id):
            continue
        borders = find_borders(mask_membrane, int_shape, raise_error=False, type=parameter_dict["FEM_mode"],
                               min_length=parameter_dict["min_line_length"])
        if not isinstance(borders, Cells_and_Lines):
            print("couldn't identify cell borders in frame %s patch %s" % (str(frame), str(obj_id)))
            continue
        n_cells = borders.n_cells
        res_dict[frame]["cell number"].append([obj_id, n_cells, warn])
    # center of mass (centroid) of each object// a simple check if objects have been identified correctly
    for obj_id, com in masks.get_com_frame(frame):
        res_dict[frame]["center of object"].append([obj_id, str(np.round(com, 2)), ""])


def find_full_im_path(cdb_image, base_folder):
    rel_path = cdb_image.path.path
    filename = cdb_image.filename
    return os.path.join(os.path.join(base_folder, rel_path), filename)


def simple_shift_correction(frame, parameter_dict, res_dict, db, db_info=None, **kwargs):
    # load images from database
    im_a = db.getImage(id=db_info["file_order"][frame + "images_after"])
    im_b = db.getImage(id=db_info["file_order"][frame + "images_before"])

    image_after = im_a.data
    image_before = im_b.data
    try:
        im_m = db.getImage(id=db_info["file_order"][frame + "membranes"])
        image_membrane = im_m.data
        im_m_path = find_full_im_path(im_m, db_info["path"])
        n_frames = 3
    except KeyError:
        n_frames = 2

    # get paths for saving later
    im_a_path = find_full_im_path(im_a, db_info["path"])
    im_b_path = find_full_im_path(im_b, db_info["path"])

    # performing drift correction
    if n_frames == 2:
        b_save, a_save, [], drift = correct_stage_drift(image_before, image_after,
                                                        additional_images=[])
    if n_frames == 3:
        b_save, a_save, [m_save], drift = correct_stage_drift(image_before, image_after,
                                                              additional_images=[image_membrane])

    print("\nframe %s: found drift of %s" % (frame, str(np.round(drift, 3))))
    b_save.save(im_b_path)
    a_save.save(im_a_path)

    if n_frames == 3:
        m_save.save(im_m_path)


def simple_segmentation(frame, parameter_dict, res_dict, db, db_info=None, masks=None, seg_threshold=0,
                        seg_type="cell_area", im_filter=None, **kwargs):
    # only for cell layer mode
    if seg_type == "cell_area":
        mtypes = get_masks_by_key(default_parameters, "use", "stress_layer")
        index_mem = parameter_dict["mask_properties"]["Cell Boundary"]["index"]
        index_ct1 = parameter_dict["mask_properties"][mtypes[0]]["index"]
        index_ct2 = parameter_dict["mask_properties"][mtypes[1]]["index"]
        if all([m in db_info["mask_types"] for m in mtypes]):
            im = db.getImage(id=db_info["file_order"][frame + "membranes"]).data
            if not isinstance(im_filter, np.ndarray):
                im_filter = gaussian_filter(im, sigma=10)
            mask = im_filter > seg_threshold
            og_mask = try_mask_load(db, db_info["frames_ref_dict"][frame], raise_error=False, mtype="all",
                                    ret_type="zeros")
            memb_mask = og_mask == index_mem
            og_mask[mask] = index_ct1
            og_mask[~mask] = index_ct2
            og_mask[memb_mask] = index_mem
            db.setMask(frame=db_info["frames_ref_dict"][frame], data=og_mask.astype("uint8"))
    # only for colony mode
    if seg_type == "membrane":
        mtypes = get_masks_by_key(default_parameters, "use", "borders")
        index_mem = parameter_dict["mask_properties"]["Cell Boundary"]["index"]
        if all([m in db_info["mask_types"] for m in mtypes]):
            im = db.getImage(id=db_info["file_order"][frame + "membranes"]).data.astype(float)
            if not isinstance(im_filter, np.ndarray):
                im_filter = gaussian_filter(im, sigma=2) - gaussian_filter(im, sigma=10)
            mask = im_filter > seg_threshold
            mask = remove_small_objects(mask, min_size=1000)  # same as when actually preparing borders
            og_mask = try_mask_load(db, db_info["frames_ref_dict"][frame], raise_error=False, mtype="all",
                                    ret_type="zeros")
            og_mask[og_mask == index_mem] = 0  # deleting old membrane mask
            og_mask[mask] = index_mem
            db.setMask(frame=db_info["frames_ref_dict"][frame], data=og_mask.astype("uint8"))
    return im_filter


def cover_entire_image_with_mask(frame, parameter_dict, res_dict, db, db_info=None, masks=None,
                                 mask_type="Tractions", mode="fill out image", **kwargs):
    mtype = db.getMaskTypes(name=mask_type)[:]
    if len(mtype) == 0:
        print("mask " + mask_type + " not found")
        return
    index = mtype[0].index

    # filling the entire image with the index
    if mode == "fill out image":
        data = np.ones(db_info["im_shape"][frame]) * index
        db.setMask(frame=db_info["frames_ref_dict"][frame], data=data.astype(np.uint8))
    # marking the outer perimeter of the image with the mask
    # perimeter directly at the edge doesnt work for binary fill holes method later one
    # perimeter is set index-pixels distant to the image edge, so that you can set multiple masks with this method
    if mode == "encircle image":
        data = try_mask_load(db, db_info["frames_ref_dict"][frame], raise_error=False, mtype="all", ret_type="zeros")
        data[(index, -(index+1)), index:-index] = index
        data[index:-index, (index, -(index+1))] = index
        db.setMask(frame=db_info["frames_ref_dict"][frame], data=data.astype(np.uint8))

    return None, frame


def deformation(frame, parameter_dict, res_dict, db, db_info=None, masks=None, **kwargs):
    # deformation for 1 frame
    im1 = db.getImage(id=db_info["file_order"][frame + "images_after"]).data
    im2 = db.getImage(id=db_info["file_order"][frame + "images_before"]).data

    # overlap and windowsize in pixels
    window_size_pix = int(np.ceil(parameter_dict["window_size"] / parameter_dict["pixelsize"]))
    overlap_pix = int(np.ceil(parameter_dict["overlap"] / parameter_dict["pixelsize"]))
    u, v, mask, mask_std = calculate_deformation(im1.astype(np.int32), im2.astype(np.int32),
                                                 window_size_pix, overlap_pix,
                                                 std_factor=parameter_dict["std_factor"])
    db_info["defo_shape"] = u.shape
    res_dict[frame]["sum deformations image"].append(["image", np.sum(np.sqrt(u ** 2 + v ** 2)), ""])  # propably remove that

    # adding plot of deformation field to the database
    add_plot("deformation", (u, v), show_quiver, frame, db_info, parameter_dict, db)

    # saving raw files
    np.save(os.path.join(db_info["path"], frame + "u.npy"), u)
    np.save(os.path.join(db_info["path"], frame + "v.npy"), v)
    # summing deformation over certain areas
    mtypes = [m for m in db_info["mask_types"] if m in get_masks_by_key(default_parameters, "use", "defo")]
    calc_on_area(frame, res_dict, parameter_dict, "sum deformations", masks, mask_types=mtypes, x=u, y=v, sumtype="abs")
    return None, frame


def get_contractillity_contractile_energy(frame, parameter_dict, res_dict, db, db_info=None, masks=None, **kwargs):
    u, v = try_to_load_deformation(db_info["path"], frame, warn=True)
    t_x, t_y = try_to_load_traction(db_info["path"], frame, warn=False)
    if u.shape != t_x.shape:
        raise ShapeMismatchError("Traction field and deformation field have different shapes. "
                                 "Try to repeat the calculation of the traction field.")
    db_info["defo_shape"] = t_x.shape
    ps_new = parameter_dict["pixelsize"] * np.mean(np.array(db_info["im_shape"][frame]) / np.array(t_x.shape))
    # select masks
    mtypes = [m for m in db_info["mask_types"] if m in get_masks_by_key(default_parameters, "use", "forces")]
    if isinstance(u, np.ndarray):
        energy_points = strain_energy_points(u, v, t_x, t_y, parameter_dict["pixelsize"],
                                             ps_new)  # contractile energy at any point
        # plotting contractile energy (only happens if enable in default_fig_parameters)
        add_plot("energy_points", [energy_points/((ps_new*10**-6)**2)], show_map_clickpoints, frame, db_info, parameter_dict, db)

    # iterating though mask that are selected for summation
    mask_iter = masks.reconstruct_masks_frame(frame, mtypes, raise_error=False, fill_holes=True)
    contractile_force = None
    contr_energy = None
    for obj_id, mask, mtype, warn in mask_iter:
        if check_empty_mask(mask, mtype, frame, obj_id):
            continue
        # interpolation to size of traction force array
        mask_area = np.sum(mask) * (parameter_dict["pixelsize"] ** 2)
        mask_int = interpolation(mask, t_x.shape)
        # calculate contractillity only in "colony" mode
        if parameter_dict["FEM_mode"] == "colony":
            contractile_force, proj_x, proj_y, center = contractillity(t_x, t_y, ps_new, mask_int)
            #  + default_parameters["mask_properties"][mtype]["label"]
            res_dict[frame]["contractility"].append([obj_id, contractile_force, warn])
        # calculating the area of this mask
        label = "area" + parameter_dict["mask_properties"][mtype]["label"]
        res_dict[frame][label].append([obj_id, mask_area, warn])

        # calculate contractile energy if deformations are provided
        if isinstance(u, np.ndarray):
            check_shape(u, t_x)
            contr_energy = np.sum(energy_points[mask_int])  # sum of contractile energy on on mask
            #+ default_parameters["mask_properties"][mtype]["label"]
            res_dict[frame]["strain energy"].append(
                [obj_id, contr_energy, warn])
        print("strain energy=", round_flexible(contr_energy), "contractility=", round_flexible(contractile_force))

    return (contractile_force, contr_energy), frame


def traction_force(frame, parameter_dict, res_dict, db, db_info=None, masks=None, **kwargs):
    # trying to laod deformation
    u, v = try_to_load_deformation(db_info["path"], frame, warn=False)
    db_info["defo_shape"] = u.shape
    ps_new = parameter_dict["pixelsize"] * np.mean(np.array(db_info["im_shape"][frame]) / np.array(u.shape))

    # using tfm with or without finite thickness correction
    if parameter_dict["TFM_mode"] == "finite_thickness":
        tx, ty = TFM_tractions(u, v, pixelsize1=parameter_dict["pixelsize"],
                               pixelsize2=ps_new,
                               h=parameter_dict["h"], young=parameter_dict["young"],
                               sigma=parameter_dict["sigma"],
                               filter=parameter_dict["filter_type"], fs=parameter_dict["filter_size"])

    if parameter_dict["TFM_mode"] == "infinite_thickness":
        tx, ty = ffttc_traction(u, v, pixelsize1=parameter_dict["pixelsize"],
                                pixelsize2=ps_new,
                                young=parameter_dict["young"],
                                sigma=parameter_dict["sigma"],
                                filter=parameter_dict["filter_type"], fs=parameter_dict["filter_size"])

    # add a plot of the traction filed to the database
    add_plot("traction", (tx, ty), show_quiver, frame, db_info, parameter_dict, db)

    # saving raw files
    np.save(os.path.join(db_info["path"], frame + "tx.npy"), tx)
    np.save(os.path.join(db_info["path"], frame + "ty.npy"), ty)

    mtypes = [m for m in db_info["mask_types"] if m in get_masks_by_key(default_parameters, "use", "forces")]
    calc_on_area(frame, res_dict, parameter_dict, "sum traction forces", masks, mask_types=mtypes, x=tx, y=ty,
                 sumtype="abs")
    return None, frame


def FEM_grid_setup(frame, parameter_dict, mask_grid, db_info=None, warn="", **kwargs):
    '''
    :param frame:
    :param parameter_dict:
    :param db:
    :param db_info:
    :param kwargs:
    :return:
    '''

    # loading forces and update shape info
    t_x, t_y = try_to_load_traction(db_info["path"], frame, warn=False)
    db_info["defo_shape"] = t_x.shape  # pixelsize of fem grid in µm
    ps_new = parameter_dict["pixelsize"] * np.mean(np.array(db_info["im_shape"][frame]) / np.array(t_x.shape))
    # preparig the mask
    mask_area = prepare_mask_FEM(mask_grid, t_x.shape)  # area for FEM analysis
    warn_grid = warn_small_FEM_area(mask_area, threshold=1000)
    warn = warn + " " + warn_grid

    # FEM grid setupplt.
    # preparing forces
    f_x = t_x * ((ps_new * (10 ** -6)) ** 2)  # point force for each node from tractions
    f_y = t_y * ((ps_new * (10 ** -6)) ** 2)
    if parameter_dict["FEM_mode"] == "colony":
        # using mask for grid setup
        # trying to load cell colony mask, raise error if not found
        # coorecting force for torque and net force
        f_x_c2, f_y_c2, p = correct_forces(f_x, f_y, mask_area)
        # get_torque1(f_y,f_x,mask_area)
        nodes, elements, loads, mats = grid_setup(mask_area, -f_x_c2, -f_y_c2, 1, sigma=parameter_dict["sigma_cells"],
                                                  edge_factor=parameter_dict["edge_padding"])  # note the negative signe

    if parameter_dict["FEM_mode"] == "cell layer":
        nodes, elements, loads, mats = grid_setup(mask_area, -f_x, -f_y, 1, sigma=parameter_dict["sigma_cells"],
                                                  edge_factor=parameter_dict["edge_padding"])  # note the negative signe

    return nodes, elements, loads, mats, mask_area, warn, ps_new


def FEM_analysis_average_stresses(frame, res_dict, parameter_dict, db, db_info, stress_tensor, ps_new, masks, obj_id,
                                  **kwargs):
    # analyzing the FEM results with average stresses
    sigma_max, sigma_min, sigma_max_abs, tau_max, phi_n, phi_shear, sigma_mean = all_stress_measures(stress_tensor,
                                                                                                     px_size=ps_new * 10 ** -6)

    if parameter_dict["FEM_mode"] == "cell layer":
        use_type = "stress_layer"
        add_cut_factor = parameter_dict["edge_padding"] + parameter_dict["padding_cell_layer"]
    else:
        use_type = "stress_colony"
        add_cut_factor = None

    mtypes = [m for m in db_info["mask_types"] if m in get_masks_by_key(default_parameters, "use", use_type)]

    calc_on_area(frame, res_dict, parameter_dict, "mean normal stress", masks, mask_types=mtypes, obj_ids=[obj_id],
                 x=sigma_mean, sumtype="mean_abs", add_cut_factor=add_cut_factor)
    calc_on_area(frame, res_dict, parameter_dict, "max normal stress", masks, mask_types=mtypes, obj_ids=[obj_id],
                 x=sigma_max_abs, sumtype="mean_abs", add_cut_factor=add_cut_factor)
    calc_on_area(frame, res_dict, parameter_dict, "max shear stress", masks, mask_types=mtypes, obj_ids=[obj_id],
                 x=tau_max,
                 sumtype="mean_abs", add_cut_factor=add_cut_factor)

    calc_on_area(frame, res_dict, parameter_dict, "cv mean normal stress", masks, mask_types=mtypes, obj_ids=[obj_id],
                 x=sigma_mean, sumtype="cv", add_cut_factor=add_cut_factor, db_info=db_info)
    calc_on_area(frame, res_dict, parameter_dict, "cv max normal stress", masks, mask_types=mtypes, obj_ids=[obj_id],
                 x=sigma_max_abs, sumtype="cv", add_cut_factor=add_cut_factor, db_info=db_info)
    calc_on_area(frame, res_dict, parameter_dict, "cv max shear stress", masks, mask_types=mtypes, obj_ids=[obj_id],
                 x=tau_max, sumtype="cv", add_cut_factor=add_cut_factor, db_info=db_info)

    return sigma_mean


def FEM_analysis_borders(frame, res_dict, db, db_info, parameter_dict, stress_tensor, ps_new, borders, obj_id, warn,
                         **kwargs):
    # retrieving spline representation of borders
    lines_splines = borders.lines_splines
    line_lengths = borders.line_lengths
    # plot lines tresses over border as continuous curves:
    lines_interpol, min_v, max_v = lineTension(lines_splines, line_lengths,
                                               stress_tensor, pixel_length=ps_new,
                                               interpol_factor=1)
    plot_values = [borders.inter_shape, borders.edge_lines, lines_interpol, min_v, max_v]

    # norm of the line tension vector
    line_tension_norm = mean_stress_vector_norm(lines_interpol, borders, norm_level="points", vtype="t_vecs",
                                                exclude_colony_edge=True)
    # normal component of the line tension vector
    line_tension_n = mean_stress_vector_norm(lines_interpol, borders, norm_level="points", vtype="t_normal",
                                             exclude_colony_edge=True)
    # shear component of the line tension vector
    line_tension_sh = mean_stress_vector_norm(lines_interpol, borders, norm_level="points", vtype="t_shear",
                                              exclude_colony_edge=True)

    res_dict[frame]["average magnitude line tension"].append([obj_id, line_tension_norm[1], warn])
    res_dict[frame]["std magnitude line tension"].append([obj_id, line_tension_norm[2], ""])
    res_dict[frame]["average normal line tension"].append([obj_id, line_tension_n[1], warn])
    res_dict[frame]["std normal line tension"].append([obj_id, line_tension_n[2], ""])
    res_dict[frame]["average shear line tension"].append([obj_id, line_tension_sh[1], warn])
    res_dict[frame]["std shear line tension"].append([obj_id, line_tension_sh[2], ""])

    # currently cells are not detected in cell layer mode// could ne implemented though...
    if parameter_dict["FEM_mode"] == "colony":
        avg_cell_force = mean_stress_vector_norm(lines_interpol, borders, norm_level="cells", vtype="t_vecs",
                                                 exclude_colony_edge=True)
        avg_cell_pressure = mean_stress_vector_norm(lines_interpol, borders, norm_level="cells", vtype="t_normal",
                                                    exclude_colony_edge=True)
        avg_cell_shear = mean_stress_vector_norm(lines_interpol, borders, norm_level="cells", vtype="t_shear",
                                                 exclude_colony_edge=True)
        res_dict[frame]["average cell force"].append([obj_id, avg_cell_force[1], warn])
        res_dict[frame]["average cell pressure"].append([obj_id, avg_cell_pressure[1], warn])
        res_dict[frame]["average cell shear"].append([obj_id, avg_cell_shear[1], warn])
        res_dict[frame]["std cell force"].append([obj_id, avg_cell_force[2], ""])
        res_dict[frame]["std cell pressure"].append([obj_id, avg_cell_pressure[2], ""])
        res_dict[frame]["std cell shear"].append([obj_id, avg_cell_shear[2], ""])

    return None, frame, plot_values


def FEM_full_analysis(frame, parameter_dict, res_dict, db, db_info=None, masks=None, **kwargs):
    # performing full MSM/finite elements analysis
    # wrapper to flexibly perform FEM analysis
    # masks for FEM grid
    plot_values = []
    m_stresses = []
    FEM_type = "FEM_layer" if parameter_dict["FEM_mode"] == "cell layer" else "FEM_colony"
    # masks that make up the FEM_are
    FEM_area_masks = get_masks_by_key(parameter_dict, "use", FEM_type)
    # add relevant masks, usefull in cell layer mode
    if FEM_type == "cell layer":
        mask_iter_grid = masks.reconstruct_masks_frame_add(frame, FEM_area_masks, raise_error=False, fill_holes=True)
    else:
        mask_iter_grid = masks.reconstruct_masks_frame(frame, FEM_area_masks, raise_error=False, fill_holes=True)
    # masks for line tension along cell-cell borders
    mask_iter_borders = masks.reconstruct_masks_frame(frame, "Cell Boundary", raise_error=False, fill_holes=False)
    for obj_id, mask_grid, mtype, warn in mask_iter_grid:
        if check_empty_mask(mask_grid, mtype, frame, obj_id):
            continue
        nodes, elements, loads, mats, mask_area, warn, ps_new = FEM_grid_setup(frame, parameter_dict, mask_grid,
                                                                               db_info=db_info, **kwargs)
        # FEM solution
        UG_sol, stress_tensor = FEM_simulation(nodes, elements, loads, mats, mask_area)
        np.save(os.path.join(db_info["path"], frame + "stress_tensor.npy"), stress_tensor / (ps_new * 10**-6))

        # analyzing stresses and stress distribution
        mean_normal_stress = FEM_analysis_average_stresses(frame, res_dict, parameter_dict, db, db_info, stress_tensor,
                                                           ps_new, masks, obj_id)
        m_stresses.append(mean_normal_stress)

        # finding cell borders
        mask_borders, warn_borders = next(((v[1], v[3]) for v in mask_iter_borders if v[0] == obj_id), (None, ""))

        if check_empty_mask(mask_borders, "Cell Boundary", frame, obj_id):
            continue
        # additional cutting when in layer mode
        if parameter_dict["FEM_mode"] == "layer":
            mask_borders, w = cut_mask_from_edge(mask_borders,
                                                 parameter_dict["edge_padding"] + parameter_dict[
                                                     "padding_cell_layer"],
                                                 parameter_dict["TFM_mode"] == "colony")
        warn += " " + warn_borders
        borders = find_borders(mask_borders, mask_area.shape, raise_error=False, type=parameter_dict["FEM_mode"],
                               min_length=parameter_dict["min_line_length"])
        if not isinstance(borders, Cells_and_Lines):
            print("couldn't identify cell borders in frame %s patch %s" % (str(frame), str(obj_id)))
            continue
        # analyzing line tension
        k, f, pv = FEM_analysis_borders(frame, res_dict, db, db_info, parameter_dict, stress_tensor, ps_new,
                                        borders, obj_id, warn,
                                        **kwargs)
        plot_values.append(pv)
    # plotting the stress at cell borders
    if len(plot_values) > 0:
        add_plot("FEM_borders", [plot_values], plot_continuous_boundary_stresses, frame, db_info, parameter_dict, db)
    # plotting the stress on the colony area

    if len(m_stresses) > 0:
        m_stresses = np.sum(m_stresses, axis=0)
        add_plot("stress_map", [m_stresses], show_map_clickpoints, frame, db_info, parameter_dict, db)


def provide_basic_objects(db, frames, parameter_dict, db_info, masks, res_dict):
    frames = make_iterable(frames)
    if not isinstance(db_info, dict):
        db_info, all_frames = get_db_info_for_analysis(db)
    if not isinstance(masks, cells_masks):
        masks = cells_masks(frames, db, db_info, parameter_dict)  # create new masks object if necessary
    else:  # add frames that are missing
        for frame in frames:
            if frame not in masks.masks_dict.keys():
                masks.add_frames(frame, parameter_dict["min_obj_size"])

    if isinstance(res_dict, defaultdict):
        if isinstance(res_dict.default_factory(), defaultdict):
            if isinstance(res_dict.default_factory().default_factory(), list):
                return db_info, masks, res_dict
    res_dict = defaultdict(lambda: defaultdict(list))
    return db_info, masks, res_dict


def apply_to_frames(db, parameter_dict, analysis_function, leave_basics=False, res_dict=None, frames=[], db_info=None,
                    masks=None, **kwargs):
    '''
    wrapper to apply analysis function on all frames
    :param db: clickpoints database
    :param parameter_dict: parameters for piv deformation calculation: (window size, overlap), sigma and Young's modulus
    of the gel (or of the cell sheet when applying FEM), hight of the gel
    :param func: function that is analyzed
    :param frames: list f frames (of the cdb database) to be analyze e.g [0,1,2]
    :param db_info: dictionary with the keys "path","frames_ref_dict","im_shape","file_order" as constructed from
    get db_info_for_analysis
    :param res_dict: dictionary of results to be filled up adn appended
    :return:
    '''

    frames = make_iterable(frames)
    if not leave_basics:
        db_info, masks, res_dict = provide_basic_objects(db, frames, parameter_dict, db_info, masks, res_dict)
    print(calculation_messages[analysis_function.__name__] % str(frames))
    for frame in tqdm(frames, total=len(frames)):
        try:
            analysis_function(frame, parameter_dict, res_dict, db=db, db_info=db_info, masks=masks, **kwargs)
        except Exception as e:
            if type(e) in (Mask_Error, FileNotFoundError, FindingBorderError, ShapeMismatchError):
                print(e)
            else:
                raise (e)
    return db_info, masks, res_dict


### code to work on clickpoint outside of the addon
if __name__ == "__main__":

    ## setting up necessary paramteres
    # db=clickpoints.DataFile("/home/user/Desktop/Monolayers_new_images/monolayers_new_images/KO_DC1_tomatoshift/database.cdb","r")
    db = clickpoints.DataFile("/home/andy/test_data_pyTFM/KOshift/database.cdb", "r")
    parameter_dict = default_parameters
    res_dict = defaultdict(lambda: defaultdict(list))
    db_info, all_frames = get_db_info_for_analysis(db)

    # db_info, masks, res_dict = apply_to_frames(db, parameter_dict, deformation, res_dict, frames="12",
    #                                            db_info=db_info, masks=None)
    #db_info, masks, res_dict = apply_to_frames(db, parameter_dict, general_properties, res_dict=res_dict,
   #                                            frames=all_frames,
    #                                           db_info=db_info, masks=None)
    db_info, masks, res_dict = apply_to_frames(db, parameter_dict, FEM_full_analysis, res_dict=res_dict, frames="01",
                                               db_info=db_info, masks=None)

    # parameter_dict["filter_type"]=None
# default_fig_parameters["file_names"]["traction"] = "t_none.png"
# db_info, masks, res_dict = apply_to_frames(db, parameter_dict, traction_force, res_dict, frames="12",
#                                       db_info=db_info, masks=None)

'''
    ###### problem: produces empty entries when try to access non-exisitng str
    masks = cells_masks(all_frames, db, db_info, parameter_dict)
    db_info, masks, res_dict = apply_to_frames(db, parameter_dict, general_properties, res_dict, frames="01",
                                               db_info=db_info, masks=masks)
    #
   # db_info, masks, res_dict = apply_to_frames(db, parameter_dict, simple_segmentation, res_dict, frames="1",
   #                                            db_info=db_info, masks=masks,seg_threshold=0,seg_type="cell_area",)
    db_info, masks, res_dict = apply_to_frames(db, parameter_dict, FEM_full_analysis, res_dict, frames="01", db_info=db_info,masks=masks)
    # # db_info, masks, res_dict = apply_to_frames(db, parameter_dict, traction_force, res_dict, frames="1", db_info=db_info, masks=masks)

    #db_info, masks, res_dict = apply_to_frames(db, parameter_dict, traction_force, res_dict, frames="01", db_info=db_info,masks=masks)
    db_info,masks,res_dict=apply_to_frames(db, parameter_dict, get_contractillity_contractile_energy, res_dict, frames="01", db_info=db_info,masks=masks)
   # print(res_dict)
        #apply_to_frames(db, parameter_dict, FEM_full_analysis, res_dict, frames="12", db_info=db_info)

        #apply_to_frames(db, parameter_dict, FEM_full_analysis, res_dict, frames=all_frames, db_info=db_info)

    write_output_file(res_dict, "results", "/home/user/Desktop/backup_from_harddrive/data_traction_force_microscopy/WT_vs_KO_images/WTshift/out_test.txt",new_file=True)
        # calculating the deformation field and adding to data base
'''
#TDOD: deal with oroblems when database is moved around/ folders are renamed --> make option with this directory
#TODO insert tables as images for read the docs