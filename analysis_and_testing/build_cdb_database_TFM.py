from pyTFM.database_functions import *  # must be on top becauseof some matplotlib backend issues
import os
folder=os.getcwd()
#folder="/home/user/Desktop/Andreas-Python/tracktion_force_microscopy/exmaple_analysis/KO_shift_part/"
#folder="/home/user/Software/tracktion_force_microscopy/tracktion_force_microscopy/test_images_database_setup/"
setup_database_for_tfm(folder,"database.cdb",return_db=False)




