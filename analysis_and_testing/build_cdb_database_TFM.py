from andreas_TFM_package.database_functions import *  # must be on top becauseof some matplotlib backend issues
import os
folder=os.getcwd()
#folder="/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/out_clickpoints"
frames=setup_database_for_tfm(folder,"database.cdb",return_db=False)
