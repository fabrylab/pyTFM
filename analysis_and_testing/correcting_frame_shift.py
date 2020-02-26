# correcting frame shift between images of beads before and after cell removal.

# the correction is done by finding the shift between two images using image registration. Then the images are cropped
# to the common field of view. If this script finds further images of the cells, it wil also cropp them to this field
# of view. The output is saved to the input folder. For each "experiment" a new folder is created. An experiment is
# identified as a directory that contains one folder for the images before cell removal and one folder with images after
# the cell removal.


# importing functions for drift correction and finding and grouping images
from pyTFM.frame_shift_correction import *
import re

# First set an input folder. The script will search the entire tree of this folder.
# Its better to enter a directory with r"". This is needed to interpret backslashes correctly.
#folder = r"/media/user/GINA1-BK/test_tfm_structure" # enter a custom folder
folder=os.getcwd() # or use the current directory

# The script expects two sub folders: one containing images after cell removal and one containing images
# before cell removals.
# you can set the identifier for these two folders like this:

# The first argument of re.compile must be a regular expression
# regular expressions allow you to search stings. Go to https://docs.python.org/3/library/re.html for a documentation
# identifier for the folder with images after cell removal:
after_folder_identifier = re.compile("after", re.IGNORECASE)
# identifier for the folder with images before cell removal:
before_folder_identifier = re.compile("before", re.IGNORECASE)

# The script searches for image files in the folders above only if both folders can be found in the same directory.
# "after" images are only searched for in the "after" folder, and "before" only in the "before" folder.
# "brightfield" is searched in both folders. While searching images the script also identifies the frame number
# of the image. You can set the identifiers for files like this:
# The frame number needs to be marked as a group with "()". Not including a group will cause an error.
after_file_identifier = re.compile("(\d{0,3})_{0,1}fluo", re.IGNORECASE)
before_file_identifier = re.compile("(\d{0,3})_{0,1}fluo", re.IGNORECASE)
bf_file_identifier = re.compile("(\d{0,3})_{0,1}BF_before", re.IGNORECASE)
# putting identifiers in one list
identifier_list=[after_folder_identifier, before_folder_identifier, after_file_identifier, before_file_identifier, bf_file_identifier]

# You can also define the output names for the frame shift corrected images. The order is:
# [after bead removal, before bead removal,third image showing the cells]
names=["after_shift.tif","before_shift.tif","bf_before_shift.tif"]

# finding files and sorting for frames and experiments
# this functions iterates through a directory tree. If it identifies subdirectories containing images
# before and after cell removal (specified by the "after_folder_identifier" and "before_folder_identifier")
# it notes down the name of the directory as the "experiment". Then the function enters both subdirectories and searches
# for the "before", "after" and "brightfield" images. It identifies the frame of the images. Any frame for which one of
# the images is missing is discarded entirely
files_dict=find_files_for_shifting(folder, identifier_list)

# frame shift correction an cutting to a common field of view
# using image registration
cut_images(folder,files_dict,names)













