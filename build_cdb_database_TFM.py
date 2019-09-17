# script to set up a clickpoints database suitable for using the TFM addon

from andreas_TFM_package.TFM_functions_for_clickpoints import *
#folder="/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/out_clickpoints"
folder=os.getcwd()
frames=setup_database_for_tfm(folder,"database.cdb",return_db=False)
