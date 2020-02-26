
import clickpoints
from pyTFM.database_functions import *

db=clickpoints.DataFile("/home/user/Software/tracktion_force_microscopy"
                        "/tracktion_force_microscopy/test_images_database_setup/database.cdb")

folders = {"folder1_txt": "/home/user/Software/tracktion_force_microscopy/tracktion_force_microscopy/test_images_database_setup/",
                "folder2_txt": "/home/user/Software/tracktion_force_microscopy/tracktion_force_microscopy/test_images_database_setup/",
                "folder3_txt": "/home/user/Software/tracktion_force_microscopy/tracktion_force_microscopy/test_images_database_setup/",
                "folder_out_txt": "/home/user/Software/tracktion_force_microscopy/tracktion_force_microscopy/test_images_database_setup/"}
search_keys = {"after": "\d{1,4}after", "before": "\d{1,4}before",
                    "cells": "\d{1,4}bf_before",
                    "frames": "(\d{1,4})"}
db.getPath(id=1)
db._AddOption(key="key",value="vlaue")