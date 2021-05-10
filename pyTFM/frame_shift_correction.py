# correction for frame shift in the images of beads

import copy
import os
import re
from collections import defaultdict

import numpy as np
from PIL import Image
from pyTFM.utilities_TFM import createFolder, get_group
from scipy.ndimage.interpolation import shift
from skimage.feature import register_translation
from skimage.registration import phase_cross_correlation


def normalizing(img):
    img = img - np.percentile(img, 1)  # 1 Percentile
    img = img / np.percentile(img, 99.99)  # norm to 99 Percentile
    img[img < 0] = 0.0
    img[img > 1] = 1.0
    return img


def croping_after_shift(image, shift_x, shift_y):
    if shift_x <= 0:
        image = image[:, int(np.ceil(-shift_x)):]
    else:
        image = image[:, :-int(np.ceil(shift_x))]
    if shift_y <= 0:
        image = image[int(np.ceil(-shift_y)):, :]
    else:
        image = image[:-int(np.ceil(shift_y)), :]
    return np.array(image, dtype=np.float)


def check_subdirs(after_identifier, before_identifier, dirs):
    after_subdir = None
    before_subdir = None
    for sd in dirs:
        if re.match(after_identifier, sd):
            after_subdir = sd
        if re.match(before_identifier, sd):
            before_subdir = sd
    return after_subdir, before_subdir


def check_files_dict(files_dict):
    # warning if not before and after images where found and deleting corresponding entry
    require_keys = ["before", "after"]
    files_dict_cp = copy.deepcopy(files_dict)
    for exp, frames_dict in files_dict_cp.items():
        for frame, img_files in frames_dict.items():
            for key in require_keys:
                if key not in img_files.keys():
                    print("couldn't find %s image for frame %s in experiment %s" % (key, frame, exp))
                    try:
                        del files_dict[exp][frame]
                    except KeyError:
                        pass


def find_files_for_shifting(folder, identifier_list):
    after_folder_identifier, before_folder_identifier, after_file_identifier, before_file_identifier, bf_file_identifier = identifier_list
    files_dict = {}
    for subdir, dirs, files in os.walk(folder):
        after_subdir, before_subdir = check_subdirs(after_folder_identifier, before_folder_identifier, dirs)
        if after_subdir is None or before_subdir is None:
            continue  # making sure the folder contains both "before" and "after" images
        else:
            experiment = os.path.split(subdir)[1]  # selecting experiment name: the current directory
            files_dict[experiment] = defaultdict(dict)
            after_subdir_f = os.path.join(subdir, after_subdir)
            before_subdir_f = os.path.join(subdir, before_subdir)

            for file in os.listdir(after_subdir_f):
                frame_after = get_group(re.search(after_file_identifier, file), 1)
                frame_bf = get_group(re.search(bf_file_identifier, file), 1)
                if frame_after:  # finding fluo images of beads
                    files_dict[experiment][frame_after]["after"] = os.path.join(after_subdir_f, file)
                if frame_bf:
                    files_dict[experiment][frame_bf]["bf"] = os.path.join(after_subdir_f, file)

            for file in os.listdir(before_subdir_f):
                frame_before = get_group(re.search(before_file_identifier, file), 1)
                frame_bf = get_group(re.search(bf_file_identifier, file), 1)
                if frame_before:  # finding fluo images of beads
                    files_dict[experiment][frame_before]["before"] = os.path.join(before_subdir_f, file)
                if frame_bf:
                    files_dict[experiment][frame_bf]["bf"] = os.path.join(before_subdir_f, file)
    check_files_dict(files_dict)
    return files_dict


def cut_images(folder, files_dict, names=["after_shift.tif", "before_shift.tif", "bf_before_shift.tif"]):
    for exp, frames_dict in files_dict.items():
        new_folder = os.path.join(folder, exp + "_shift")
        createFolder(new_folder)
        for frame, img_files in frames_dict.items():
            print("experiment", exp)
            print("frame", frame)
            print("image_file", img_files)

            img_b = np.asarray(Image.open(img_files["before"]))
            img_a = np.asarray(Image.open(img_files["after"]))
            imb_b_BF = np.asarray(Image.open(img_files["bf"]))

            b, a [bf], (shift_x, shift_y) = correct_stage_drift(img_b, img_a, additional_images=[imb_b_BF])

            b_save.save(os.path.join(new_folder, frame + "after_shift.tif"))
            a_save.save(os.path.join(new_folder, frame + "before_shift.tif"))
            bf_save.save(os.path.join(new_folder, frame + "bf_before_shift.tif"))


# correcting frame shift between images of beads before and after cell removal.

# the correction is done by finding the shift between two images using image registration. Then the images are cropped
# to the common field of view. If this script finds further images of the cells, it wil also cropp them to this field
# of view. The output is saved to the input folder. For each "experiment" a new folder is created. An experiment is
# identified as a directory that contains one folder for the images before cell removal and one folder with images after
# the cell removal.

def correct_stage_drift(image1, image2, additional_images=[]):

    # find shift with image registration
    shift_values = phase_cross_correlation(image1, image2, upsample_factor=100)

    shift_y = shift_values[0][0]
    shift_x = shift_values[0][1]

    # using interpolation to shift subpixel precision, image2 is the reference
    image1_shift = shift(image1, shift=(-shift_y, -shift_x), order=5)

    # normalizing and converting to image format
    b = normalizing(croping_after_shift(image1_shift, shift_x, shift_y))
    a = normalizing(croping_after_shift(image2, shift_x, shift_y))
    b_save = Image.fromarray(b * 255)
    a_save = Image.fromarray(a * 255)

    # doing the same with additional images
    additional_images_save = []
    for add_image in additional_images:
        add_image_shift = shift(add_image, shift=(-shift_y, -shift_x), order=5)
        add_image_norm = normalizing(croping_after_shift(add_image_shift, shift_x, shift_y))
        add_image_save = Image.fromarray(add_image_norm * 255)
        additional_images_save.append(add_image_save)

    return b_save, a_save, additional_images_save, (shift_x, shift_y)
