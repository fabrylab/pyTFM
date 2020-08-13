
import clickpoints
from pyTFM.database_functions import *
from pyTFM.TFM_functions_for_clickpoints import *

# some illustration of line tension and stress tensors
tensor = np.array([[1,1.5],[1.5,2]])
def lt(tensor,angle):
    vx = np.cos(angle)
    vy = np.sin(angle)
    v = np.array([vx,vy])
    l = np.matmul(tensor, v)
    norm = np.linalg.norm(l)
    normal_component = np.matmul(l, v)
    shear = norm**2 - normal_component**2
    return norm, normal_component, shear
angles = np.linspace(0,1*np.pi,10000)
res = [(lt(tensor,a)) for a in angles ]
norms = [r[0] for r in res]
norm_components = [r[1] for r in res]
shear_components = [r[2] for r in res]
plt.figure()
plt.plot(angles,norms, label="norm")
plt.plot(angles,norm_components, label="normal component")
plt.plot(angles,np.abs(norm_components), label="abs normal component")
plt.plot(angles,shear_components, label="shear component")
plt.legend()

np.mean(norms)
np.mean(norm_components)

## testing filling the whole area with masks
db = clickpoints.DataFile("/home/andy/Desktop/KOshift/database.cdb","r")
parameter_dict = default_parameters
db_info, all_frames = get_db_info_for_analysis(db)
db_info, masks, res_dict = apply_to_frames(db, parameter_dict, FEM_full_analysis, frames=all_frames[0],
                                         db_info=db_info, masks=None)

db_info, masks, res_dict = apply_to_frames(db, parameter_dict, general_properties, frames=all_frames,
                                         db_info=db_info, masks=None)
db_info, masks, res_dict = apply_to_frames(db, parameter_dict, get_contractillity_contractile_energy, frames=all_frames,
                                         db_info=db_info, masks=None)










## testing filling the whole area with masks
db = clickpoints.DataFile("/home/andy/Desktop/KOshift/database5.cdb","r")
db.deleteMasks()
db.setMaskType(name="mask",color="FFFFFF")
parameter_dict = default_parameters
db_info, all_frames = get_db_info_for_analysis(db)

# db_info, masks, res_dict = apply_to_frames(db, parameter_dict, deformation, res_dict, frames="12",
#                                            db_info=db_info, masks=None)
db_info, masks, res_dict = apply_to_frames(db, parameter_dict, cover_entire_image_with_mask, res_dict=None,
                                           frames=all_frames, leave_basics=True,
                                           db_info=db_info, masks=None, mask_type="mask", mode="encircle image")

masks = cells_masks(all_frames[0], db, db_info, parameter_dict)
mask_list = masks.reconstruct_masks_frame(all_frames[0], "mask", obj_ids=[0], raise_error=False,
                                              fill_holes=True)
plt.figure();plt.imshow(db.getMasks()[0].data)
plt.figure();plt.imshow(mask_list[0][1])

# testing the image selection
db = clickpoints.DataFile("/home/andy/Desktop/KOshift/database5.cdb","w")

folders = {"folder_after": "/home/andy/Desktop/KOshift",
                "folder_before": "/home/andy/Desktop/KOshift",
                "folder_cells": "/home/andy/Desktop/KOshift",
                "folder_out": "/home/andy/Desktop/KOshift"}
search_keys = {"after": "\d{1,4}after", "before": "\d{1,4}before",
                    "cells": "\d{1,4}bf_before",
                    "frames": "(\d{1,4})"}



# testing meta info collection
setup_database_internal(db,search_keys,folders)
db_info, all_frames = get_db_info_for_analysis(db)



# renaming new convention:
import os
import copy
import shutil
rename_dict = {"x-1_y-1":0,"x-1_y0":1,"x0_y-1":2,"x0_y0":3,"x0_y1":4,"x1_y0":5,"x1_y1":6,"x-1_y1":7,"x1_y-1":8}
folder = "/home/user/Desktop/plate1/relaxed/"
files = os.listdir(folder)
for f in files:
    match = [k for k in rename_dict.keys() if k in f]
    if len(match)==1:
        f_new = f.replace(match[0],str(rename_dict[match[0]]))
        shutil.move(os.path.join(folder,f),os.path.join(folder,f_new))




def test_n(phi, st):
    n = np.array([np.cos(phi), np.sin(phi)])
    n_f = np.abs(np.dot(np.matmul(st, n), n))
    s_f = np.sqrt(np.linalg.norm(np.matmul(st, n))**2-n_f**2)
    return n_f, s_f


st = np.array([[5000,0],[0,0.5]])
r = np.linspace(0,2*np.pi,1000)
nvs = [test_n(i,st)[0] for i in r]
svs = [test_n(i,st)[1] for i in r]
print(np.mean(nvs))
print(np.mean(svs))

plt.figure()
plt.plot(r,nvs)
plt.figure()
plt.plot(r,svs)

