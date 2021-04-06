
import clickpoints
from pyTFM.database_functions import *
from pyTFM.TFM_functions_for_clickpoints import *

if __name__ == "__main__":
    import clickpoints.launch
    print(clickpoints.__file__)
    clickpoints.launch.main("/home/user/Desktop/paper_pyTFM/data_traction_force_microscopy/WT_vs_KO_images/KOshift/database.cdb")

