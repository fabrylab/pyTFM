# correcting frame shift between images of beads before and after cell removal
# also searching for the corresponding bright field image and cutting this image to the same field of view.



from andreas_TFM_package.frame_shift_correction import *
import re

folder = "/media/user/GINA1-BK/test_tfm_structure" # enter a folder
# output names for the frame shift corrected images. THe order is always: [after bead removal, before bead removal,
#third image showing the cells]
names=["after_shift.tif","before_shift.tif","bf_before_shift.tif"]
# the script expects two sub folder: one containing images after cell removal and one containing images
# before cell removal. Bright field or other images of the cells are searched in both of these folders
# identifier for the folder with images after cell removal:
after_folder_identifier = re.compile("after", re.IGNORECASE)

# identifier for the folder with images before cell removal:
before_folder_identifier = re.compile("before", re.IGNORECASE)

# the script searches for image files in the folders above and only if both folders can be found in the same directory
# "after" images are only searched for in the "after" folder, same with "before". "brightfield" is searched in both
# folders.

# identifiers for files: th frame number needs to be marked as a group with "()"
after_file_identifier = re.compile("(\d{0,3})_{0,1}fluo", re.IGNORECASE)  # needs frame identifier
before_file_identifier = re.compile("(\d{0,3})_{0,1}fluo", re.IGNORECASE)
bf_file_identifier = re.compile("(\d{0,3})_{0,1}BF_before", re.IGNORECASE)
# putting identifiers in one list
identifier_list=[after_folder_identifier, before_folder_identifier, after_file_identifier, before_file_identifier, bf_file_identifier]

# finding files and sorting  for frames and experiments
# this functions iterates through a directory tree. If it identifies subdirectories containing images
# before and after cell removal (specified by the "after_folder_identifier" and "before_folder_identifier")
# it notes down the name of the directory as the "experiment". Then the function enters both subdirectories and searches
# for the "before", "after" and "brightfield" images. It identifies the frame of the images. Any frame for which one of
# the images is missing is discarded entirely

files_dict=find_files_for_shifting(folder, identifier_list)

# frame shift correction an cutting to a common field of view
# using image registration
cut_images(folder,files_dict,names)













