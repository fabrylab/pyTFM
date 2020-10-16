
import clickpoints
from pyTFM.database_functions import *
from pyTFM.TFM_functions_for_clickpoints import *

if __name__ == "__main__":
    import clickpoints.launch
    print(clickpoints.__file__)
    clickpoints.launch.main("/home/user/Software/pyTFM/test_dataset/KO/04after.tif")
