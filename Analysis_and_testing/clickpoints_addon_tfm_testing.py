from TFM_functions_for_clickpoints import *  # must be on top becauseof some matplotlib backend issues
import clickpoints
import os
import re
import numpy as np
from utilities import  get_group,createFolder
from PIL import Image
import copy
from collections import defaultdict
from tqdm import tqdm
from solidspy import solids_GUI
from collections import defaultdict
from datetime import datetime as dt


# def delete_empty_layers(db): ## improve this??
#     for l in db.getLayers():
#         if len(l.images)==0: # this is slow?
#             db.deleteLayers(id=l.id)

# def remove_png_images(db):
#     del_ids=[]
#     for i in db.getImages():
#         filename = i.filename
#         id = i.id
#         if filename.endswith(".png"):
#             del_ids.append(id)
#     db.deleteImages(id=del_ids)

#folder="/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/out_clickpoints/"
#images=[x for x in os.listdir(folder) if ".tif" in x]
#frames=[int(get_group(re.search("(\d{1,4})*",x),0))-1 for x in images ]
#images1=[os.path.join(folder,x) for x in images]
#zero_pad_images(images1,os.path.join(folder,"padded"))




folder="/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/out_clickpoints"
db=clickpoints.DataFile(os.path.join(folder,"database3.cdb"),"r") # genarating new clickpoints file

# calculating deformation with piv

sigma=0.49#poison ratio
young=25536# youngsmodulus
pixelsize = 4.09 / 40 #µm/pixel
window_size=100   # in pixel , keep at 10 µm
overlapp=60   # window_size/2
std_factor=15
pixelsize2=pixelsize*(window_size - overlapp)# pixelsize of the deformation image

h=300

parameter_dict=make_paramters_dict_tfm( sigma=sigma, young=young, pixelsize=pixelsize, window_size=window_size, overlapp=overlapp,
                         std_factor=std_factor, h=h, pixel_factor = window_size - overlapp)
# calculating the deformation field and adding to data base


apply_to_all_frames(db,parameter_dict,analysis_function=deformation) # calculation of deformations on all frames
apply_to_all_frames(db,parameter_dict,analysis_function=tracktion_force) # calculation of deformations on all frames


sigma=0.49#poison ratio
pixelsize = 4.09 / 40 #µm/pixel
parameter_dict = {"sigma":0.49,"pixelsize":4.09 / 40}
folder="/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/out_clickpoints"

apply_to_all_frames(db,parameter_dict,analysis_function=FEM_analysis) # calculation of deformations on all frames
db.db.close()
### keep as external scrip..

#db.
#
# to do. override images if already existing