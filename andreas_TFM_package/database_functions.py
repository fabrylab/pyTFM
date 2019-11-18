# function to build and manipulate the clickpoints database
import re
import warnings
import clickpoints
import os
from andreas_TFM_package.utilities_TFM import *
from skimage.morphology import label,binary_dilation
import itertools


class Mask_Error(Exception):
    pass






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



def setup_database_for_tfm(folder, name, return_db=False,key1="\d{1,4}after",
                           key2="\d{1,4}before",key3=["\d{1,4}mebrane","\d{1,4}bf_before"],frame_key='(\d{1,4})'):

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
    db = clickpoints.DataFile(os.path.join(folder,name), "w")
    # regex patterns to sort files into layers. If any of these matches, the file will  be sorted into a layer.
    # keys: name of the layer, values: list of regex patterns
    key1 = make_iterable(key1)
    key2 = make_iterable(key2)
    key3 = make_iterable(key3)

    file_endings = "(.*\.png|.*\.jpg|.*\.tif|.*\.swg)" # all allowed file endings
    layer_search = {"images_after": [re.compile(k + file_endings) for k in key1],
                    "images_before": [re.compile(k + file_endings) for k in key2],
                    "membranes": [re.compile(k + file_endings) for k in key3]
                            }
    # filtering all files in the folder
    all_patterns=list(itertools.chain(*layer_search.values()))
    images = [x for x in os.listdir(folder) if any([pat.match(x) for pat in all_patterns])]
    # identifying frames by evaluating the leading number.
    frames = [get_group(re.search(frame_key, x), 1) for x in images] # extracting frame
    # generating a list of sort_ids for the clickpoints database (allows you to miss some frames)
    sort_id_list=make_rank_list(frames,dtype=int)# list of sort indexes (frames) of images in the database
    warn_incorrect_files(frames) # checking if there where more or less then three images per frame

    # initializing layer in the database
    layer_list = ["images_after", "images_before","membranes"]
    base_layer = db.getLayer(layer_list[0], create=True, id=0)
    for l in layer_list[1:]:
        db.getLayer(l, base_layer=base_layer, create=True)
    path = db.setPath(folder, 1)  # setting the path to the images

    # sorting images into layers
    for id, (sort_index_id,frame, im) in enumerate(zip(sort_id_list,frames, images)):

        if any([pat.match(im) for pat in layer_search["images_after"]]):
            layer="images_after"
            comment=frame + "after__"+str(sort_index_id)+"sid"
        if any([pat.match(im) for pat in layer_search["images_before"]]):
            layer = "images_before"
            comment = frame + "before__"+str(sort_index_id)+"sid"
        if any([pat.match(im) for pat in layer_search["membranes"]]):
            layer = "membranes"
            comment = frame+"bf__"+str(sort_index_id)+"sid"
        print("file:", im, "frame:", frame, "layer", "layer:", layer)
        db.setImage(id=id, filename=im, sort_index=sort_index_id,
                    layer=layer, path=1)
        db.setAnnotation(filename=im, comment=comment)
    # return database connection or close
    if return_db:
        return frames, db
    else:
        db.db.close()
        return frames




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
    db_info["mask_types"]=[m.name for m in db.getMaskTypes()] # update db info

def fill_patches_for_cell_layer(frame, parameter_dict,res_dict, db,db_info=None,**kwargs):
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
