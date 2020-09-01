# function to build and manipulate the clickpoints database
import itertools

import clickpoints
from pyTFM.parameters_and_strings import *
from pyTFM.utilities_TFM import *
from skimage.morphology import label, binary_dilation


class Mask_Error(Exception):
    pass


def guess_TFM_mode(db_info, parameter_dict):
    # enter type in ass add_option
    cl_cond1 = len(db_info["mask_types"]) == 2 and "cell type1" in db_info["mask_types"] and "cell type2" in db_info[
        "mask_types"]
    cl_cond2 = len(db_info["mask_types"]) == 1 and (
                "cell type1" in db_info["mask_types"] or "cell type2" in db_info["mask_types"])

    co_cond1 = len(db_info["mask_types"]) == 2 and "Cell Boundary" in db_info["mask_types"] and "Tractions" in \
               db_info[
                   "mask_types"]
    co_cond2 = len(db_info["mask_types"]) == 1 and (
            "Cell Boundary" in db_info["mask_types"] or "Tractions" in db_info["mask_types"])

    cond_empty = len(db_info["mask_types"]) == 0

    undertermined = False
    if cl_cond1 or cl_cond2:
        mode = "cell layer"
    elif co_cond1 or co_cond2:
        mode = "colony"
    elif cond_empty:
        mode = parameter_dict["FEM_mode"]
    else:
        warnings.warn("failed to guess analysis mode. Setting to 'FEM_mode'")
        mode = parameter_dict["FEM_mode"]
        undertermined = True

    return mode, undertermined


def setup_database_for_tfm(folder, name):
    '''
    Sorting images into a clickpoints database. Frames are identified by leading numbers. Layers are identified by
    the file name.
    :param folder: Folder where images are searched.
    :param name: Name of the database. Needs to end with .cdb.
    :param return_db: Choose weather function returns the database object, or weather the connection to the
    database is closed
    :param key1,key2,key3: regular expression that define how to sort images. Can be single string
    or a list. If any of the regex is matched for one key, the image will be classified accordingly.
    Don't include the file ending. Typical image endings (.png,.tif ... ) are added automatically.
    key1: image after bead removal, key2: image before bead removal, key3: image of the
    cells.
    :param frame_key: reguar expression that defines how the frame number is searched. You must
    mark the group that contains the frame with parenthesis "()".
    :return:
    '''

    # creating a new cdb database, will override an existing one.
    db = clickpoints.DataFile(os.path.join(folder, name), "w")
    folders = {"folder1_txt": os.getcwd(),
               "folder2_txt": os.getcwd(),
               "folder3_txt": os.getcwd(),
               "folder_out_txt": os.getcwd()}
    search_keys = {"after": "\d{1,4}after", "before": "\d{1,4}before",
                   "cells": "\d{1,4}bf_before",
                   "frames": "(\d{1,4})"}
    setup_database_internal(db, search_keys, folders)


def setup_database_internal(db, keys_dict, folders_dict):
    '''
    Sorting images into a clickpoints database. Frames are identified by leading numbers. Layers are identified by
    the file name.
    :param folder: Folder where images are searched.
    :param name: Name of the database. Needs to end with .cdb.
    :param return_db: Choose weather function returns the database object, or weather the connection to the
    database is closed
    :param key1,key2,key3: regular expression that define how to sort images. Can be single string
    or a list. If any of the regex is matched for one key, the image will be classified accordingly.
    Don't include the file ending. Typical image endings (.png,.tif ... ) are added automatically.
    key1: image after bead removal, key2: image before bead removal, key3: image of the
    cells.
    :param frame_key: regular expression that defines how the frame number is searched. You must
    mark the group that contains the frame with parenthesis "()".
    :return:
    '''

    # regex patterns to sort files into layers. If any of these matches, the file will  be sorted into a layer.
    # keys: name of the layer, values: list of regex patterns

    key1 = keys_dict["after"]
    key2 = keys_dict["before"]
    key3 = keys_dict["cells"]
    key3 = None if key3 in ["None","none"] else key3
    key_frame = keys_dict["frames"]

    key1 = make_iterable(key1)
    key2 = make_iterable(key2)
    key3 = make_iterable(key3)


    file_endings = "(.*\.png|.*\.jpg|.*\.tif|.*\.swg)"  # all allowed file endings
    layer_search = {"images_after": {"folder": folders_dict["folder_after"],
                                     "file_key": [re.compile(k + file_endings) for k in key1]},
                    "images_before": {"folder": folders_dict["folder_before"],
                                      "file_key": [re.compile(k + file_endings) for k in key2]},
                    "membranes": {"folder": folders_dict["folder_cells"],
                                  "file_key": [re.compile(k + file_endings) if k is not None else None for k in key3]}
                    }

    # filtering all files in the folder
    # all_patterns = list(itertools.chain(*layer_search.values()))
    images = {}
    for layer in layer_search.keys():
        folder = layer_search[layer]["folder"]
        skey = layer_search[layer]["file_key"]
        if skey[0] is not None and folder is not None:
           images[layer] = [os.path.join(folder, x) for x in os.listdir(folder) if any([pat.search(x) for pat in skey]) and re.search(key_frame, x)]

    # number of layers either 3 or 2 if no image of the cells is provided
    expected_layers = len(images.keys())

    # breaking if no images at all are found
    if all([len(ims) == 0 for ims in images.values()]):
        warnings.warn("no images found")
        return

    # finding the frames of all images
    all_frames = []
    image_frames = {}
    for layer, ims in images.items():
        frs = [get_group(re.search(key_frame, os.path.split(x)[1]), 1) for x in ims]
        all_frames.extend(frs)
        image_frames[layer] = {i:f for i, f in zip(ims, frs)}

    # defining which frame belongs to which sort index
    all_s_id, frames_ids = make_rank_list(all_frames)

    # merge to single dictionary that associates one file with a sort index
    # and layer in the clickpoints database
    image_sort_id= defaultdict(dict)
    for layer, im_frame in image_frames.items():
        for im, fr in im_frame.items():
            image_sort_id[frames_ids[fr]][layer] = im

    # checking if there where more or less then three images per frame
    image_sort_id, frames_ids = filter_incorrect_files(image_sort_id, frames_ids, expected_layers)
    ids_frames = {s_id: frame for frame, s_id in frames_ids.items()}

    # creating the layers
    layer_list = default_layer_names[:expected_layers]
    base_layer = db.getLayer(layer_list[0], create=True, id=0)
    for l in layer_list[1:]:
        db.getLayer(l, base_layer=base_layer, create=True)
    # inserting path to output text file
    db.setPath(folders_dict["folder_out"], id=1)

    # sorting images into layers
    file_order = {}
    id_frame_dict = {}

    for sort_index in list(sorted(image_sort_id.keys())):
        for layer in layer_list:
            try:
                frame = ids_frames[sort_index]
                im = image_sort_id[sort_index][layer]
                print("file:", im, "frame:", frame, "layer:", layer)
                image_object = db.setImage(filename=im, sort_index=sort_index,
                                           layer=layer)
                file_order[frame + layer] = image_object.id
                id_frame_dict[image_object.id] = frame
            except Exception as e:
                print("Someting went wrong when setting images in frame %s:"%ids_frames[sort_index])
                print(e)

    # writing meta information to the database:
    # dictionaries mapping the layer and the frame to the image id (file_order),
    # image to the frame (id_frame_dict), and the frame to the sort_index (frames_ref_dict)
    # and a list of sorted (by the associated sort indices) frames
    unique_frames = sorted(frames_ids.items() , key=lambda x: x[1])
    unique_frames = [x[0] for x in unique_frames]

    db._AddOption(key="frames_ref_dict", value=frames_ids)
    db.setOption(key="frames_ref_dict", value=frames_ids)
    db._AddOption(key="file_order", value=file_order)
    db.setOption(key="file_order", value=file_order)
    db._AddOption(key="unique_frames", value=unique_frames)
    db.setOption(key="unique_frames", value=unique_frames)
    db._AddOption(key="id_frame_dict", value=id_frame_dict)
    db.setOption(key="id_frame_dict", value=id_frame_dict)


def check_existing_masks(db, parameter_dict):
    current_mask_types = [m.name for m in db.getMaskTypes()]
    FEM_mode = parameter_dict["FEM_mode"]
    other_masks = [m for m in current_mask_types if m not in get_masks_by_key(default_parameters, "FEM_mode", FEM_mode)]
    return other_masks


def setup_masks(db, db_info, parameter_dict, delete_all=False, delete_specific=[]):
    if delete_all:  # delete every existing mask
        db.deleteMaskTypes()
    for m_type in delete_specific:  # delete a specific set of mask, usually provided by "check_exsisting_masks"
        db.deleteMaskTypes(m_type)
    # setting new masks
    FEM_mode = parameter_dict["FEM_mode"]
    mtypes = get_masks_by_key(default_parameters, "FEM_mode", FEM_mode)
    for mask_name, color, index in zip(*get_properties_masks(default_parameters, mtypes, ["name", "color", "index"])):
        db.setMaskType(mask_name, color=color, index=index)
    # update db info
    db_info["mask_types"] = [m.name for m in db.getMaskTypes()]


def fill_patches_for_cell_layer(frame, parameter_dict, res_dict, db, db_info=None, **kwargs):
    # trying to load the mask from clickpoints
    print(db_info["frames_ref_dict"][frame])
    try:
        image = db.getImage(
            frame=db_info["frames_ref_dict"][frame])  # this is a work around, dont know why i cant use getMask directly
        mask = db.getMask(image).data
    # raising error if no mask object in clickpoints exist
    except AttributeError:
        raise Mask_Error("no mask of the cell membrane found for frame " + str(frame))
    mask_part = mask == 1  # type of "cell type1"

    ### lets not talk about this...
    labels = label(mask_part, background=-1)

    edge_labels = np.hstack([labels[0, :], labels[-1, :], labels[:, 0], labels[:, -1]]).astype(int)
    edge_mask_values = np.hstack([mask_part[0, :], mask_part[-1, :], mask_part[:, 0], mask_part[:, -1]]).astype(int)
    background_patches_edge = np.unique(edge_labels[edge_mask_values == 0])
    fore_ground_patches_edge = np.unique(edge_labels[edge_mask_values == 1])
    found_labels = np.concatenate([background_patches_edge, fore_ground_patches_edge])
    neigbouring = []
    for la in fore_ground_patches_edge:
        expansion = np.logical_and(binary_dilation(labels == la), ~(labels == la))
        neigbouring.append(np.unique(labels[expansion]))
    neigbouring = np.concatenate(neigbouring)
    neigbouring = np.array([n for n in neigbouring if n not in found_labels])

    neigbouring2 = []
    for la in neigbouring:
        expansion = np.logical_and(binary_dilation(labels == la), ~(labels == la))
        neigbouring2.append(np.unique(labels[expansion]))
    neigbouring2 = np.concatenate(neigbouring2)
    neigbouring2 = np.array([n for n in neigbouring2 if n not in neigbouring and n not in found_labels])

    all_labels = np.concatenate([neigbouring, neigbouring2, fore_ground_patches_edge])
    mask_new = np.zeros_like(labels)
    for la in all_labels:
        mask_new[labels == la] = 1
    mask_new = mask_new.astype(bool)
    mask[mask_new] = 1
    mask[~mask_new] = 2
    mask = mask.astype(np.uint8)
    db.setMask(image=image, data=mask)  # udapting mask in clickpoints


def filter_incorrect_files(frames, frame_sort_index, expected=3):

    '''
    throws a waring when there more or less then three images per frame are found.
    :param frames:
    :return:
    '''
    frames_cp = frames.copy()
    deleted_sort_indices = []
    for sort_index, layers in frames.items():
        if len(layers) != expected:
            frame = [f for f, s_id in frame_sort_index.items() if s_id == sort_index]
            frame = frame[0] if len(frame)==1 else "unidentified"
            warn = "Found %s files for frame %s. " % (str(len(layers)), frame)
            warnings.warn(warn + "Expected %s files per frame." % str(expected))
            # delete all images from the uncorrect sort index
            deleted_sort_indices.append(sort_index)

    # how to reorganize sort indices after some are removed
    old_new_si = {s_id: s_id for s_id in frame_sort_index.values()}
    # iterating through all deleted sort indices
    # at every iteration the deleted sort_index is dropped and all new sort indices, that refer to an old
    # sort index that is larger then then deleted one are subtracted by -1
    for s_id in deleted_sort_indices:
        old_new_si =  {si_old: si_new - (si_old > s_id) for si_old, si_new in old_new_si.items() if si_old != s_id}

    # updating the dictionaries accordingly
    frames_cp = {old_new_si[s_id]: layer_im for s_id, layer_im in frames_cp.items() if s_id not in deleted_sort_indices}
    frame_sort_index_cp = {frame: old_new_si[s_id] for frame, s_id in frame_sort_index.items() if s_id not in deleted_sort_indices}

    return frames_cp, frame_sort_index_cp