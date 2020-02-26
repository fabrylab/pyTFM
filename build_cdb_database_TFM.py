from pyTFM.database_functions import *
import os
folder=os.getcwd()
#folder="/home/user/Desktop/Andreas-Python/tracktion_force_microscopy/exmaple_analysis/KO_shift_part/"
setup_database_for_tfm(folder,"database.cdb",return_db=False)
