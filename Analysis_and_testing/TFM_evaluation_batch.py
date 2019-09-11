import numpy as np
import matplotlib.pyplot as plt
import os

import openpiv.tools
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from scipy.ndimage import zoom
from skimage.filters import rank
from skimage.morphology import cube
from scipy.ndimage.filters import uniform_filter,median_filter,gaussian_filter
import openpiv.tools
import openpiv.process
import openpiv.scaling
import openpiv.validation
import openpiv.filters
import matplotlib.pyplot as plt
import numpy as np
import copy
import os
import imageio
from tqdm import tqdm
import re
from skimage.filters import gaussian

def get_xy_for_quiver(u):
    xs=np.zeros(np.shape(u))
    for i in range(np.shape(u)[0]):
        xs[i,:]=np.arange(0,np.shape(u)[1],1)
    ys = np.zeros(np.shape(u))
    for j in range(np.shape(u)[1]):  ## is inverted in other skript...
        ys[:, j] = np.arange(0,np.shape(u)[0], 1)
    return xs, ys


def createFolder(directory):
    '''
    function to create directories, if they dont already exist
    '''
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)


from TFM_functions import *
in_folder=r"/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/"
out_folder=r"/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/" ## output folder


min_scale_value_tract=0
max_scale_value_tract=70
min_scale_value_disp=0
max_scale_value_disp=40
std_factor=15
sigma=0.49#poison ratio
young=25536# youngsmodulus
pixelsize = 6.25 / 40 #µm/pixel
window_size=64   # in pixel , keep at 10 µm
overlapp=32    # window_size/2
pixel_factor = window_size - overlapp


# listing files
files_dict={}
file_list=os.listdir(in_folder)
file_list=[x for x in file_list if "shift.tif" in x]
pos_list = [re.search("(\d{1,3})\D.*", x).group(1) for x in file_list]  # just first two positions (maybe use three
pos_list=np.unique(pos_list)


for pos in pos_list:
    files_dict[pos] = [os.path.join(in_folder, file) for file in file_list if re.match(pos+"\D.*",file)]


for pos,files in tqdm(files_dict.items(),total=len(files_dict.items())):
    if not pos =="11":
        continue
    print(pos)
    file2=[x for x in files if "after" in os.path.split(x)[1]][0] # cahnged order
    file1=[x for x in files if "before" in os.path.split(x)[1] and not "bf" in os.path.split(x)[1]][0]
    file=[x for x in files if "bf_before" in os.path.split(x)[1]][0]

    print("calculating deformation")
    u, v, x, y, mask,mask_std = calculate_deformation(file1, file2, window_size, overlapp,  std_factor=20)
    print("calculating traktion")
    tx, ty = ffttc_traction(u, v, young, pixelsize, bf_image=False, filter="mean")
    u_no_std, v_no_std, x, y, mask, mask_no_std = calculate_deformation(file1, file2, window_size, overlapp, std_factor=200)
    #saving output and images
    out_folder2=os.path.join(in_folder,pos)
    createFolder( out_folder2)
    np.save(os.path.join(out_folder2,"u.npy"),u)
    np.save(os.path.join(out_folder2, "v.npy"), v)
    np.save(os.path.join(out_folder2, "x.npy"),x)
    np.save(os.path.join(out_folder2, "y.npy"),y)
    np.save(os.path.join(out_folder2, "mask.npy"),mask)
    np.save(os.path.join(out_folder2, "tx.npy"),tx)
    np.save(os.path.join(out_folder2, "ty.npy"),ty)

    #defo image:
    filter=np.sqrt((u**2+v**2))<1
    u_filter=copy.deepcopy(u)
    u_filter[filter]=0
    v_filter = copy.deepcopy(v)
    v_filter[filter] = 0
    fig=plt.figure()
    plt.imshow(np.sqrt(u**2+v**2),vmin=min_scale_value_disp,vmax=max_scale_value_disp)
    cbar=plt.colorbar()
    cbar.set_label("displacement in pixel")
    x1,y1=get_xy_for_quiver(u)
    scale=0.5#
    plt.quiver(x1,y1,u_filter/scale,v_filter/scale,scale=1,headwidth=3, scale_units='xy',width=0.002)
    plt.arrow(np.shape(u)[0] * 1.2, 0 - np.shape(u)[1] * 0.2, 5 / scale, 0, head_width=0, facecolor="black",
              edgecolor="black",
              clip_on=False, head_starts_at_zero=False, overhang=2)
    plt.text(np.shape(u)[0] * 1.25, 0 - np.shape(u)[1] * 0.15, "%d pixel displacement" % 5, fontsize=10,
             color="black")

    fig.savefig(os.path.join(out_folder2,pos+"deformation.png"))
    plt.close()


    #defo image_no_std:
    filter=np.sqrt((u_no_std**2+v_no_std**2))<1
    u_filter=copy.deepcopy(u_no_std)
    u_filter[filter]=0
    v_filter = copy.deepcopy(v_no_std)
    v_filter[filter] = 0
    fig=plt.figure()
    plt.imshow(np.sqrt(u_no_std**2+v_no_std**2),vmin=min_scale_value_disp,vmax=max_scale_value_disp)
    cbar=plt.colorbar()
    cbar.set_label("displacement in pixel")
    x1,y1=get_xy_for_quiver(u_no_std)
    scale=0.5#
    plt.quiver(x1,y1,u_filter/scale,v_filter/scale,scale=1,headwidth=3, scale_units='xy',width=0.002)

    plt.arrow(np.shape(u)[0] * 1.2, 0 - np.shape(u)[1] * 0.2, 5 / scale, 0, head_width=0, facecolor="black",
              edgecolor="black",
              clip_on=False, head_starts_at_zero=False, overhang=2)
    plt.text(np.shape(u)[0] * 1.25, 0 - np.shape(u)[1] * 0.15, "%d pixel displacement" % 5, fontsize=10,
             color="black")


    fig.savefig(os.path.join(out_folder2,pos+"deformation_no_std.png"))
    plt.close()


    # traction image:

    fig=plt.figure()
    plt.imshow(np.sqrt((tx/1000) ** 2 + (ty/1000) ** 2),vmin=min_scale_value_tract, vmax=max_scale_value_tract)
    cbar = plt.colorbar()
    cbar.set_label("traktion forces in kPa")
    x1, y1 = get_xy_for_quiver(tx)
    scale=2# #
    plt.quiver(x1, y1, (tx/1000) / scale, (ty/1000) / scale, scale=1, scale_units='xy',width=0.002)
    plt.arrow(np.shape(tx)[0] * 1.2, 0 - np.shape(tx)[1] * 0.2,  30 / scale, 0, head_width=0, facecolor="black",
              edgecolor="black",
              clip_on=False, head_starts_at_zero=False, overhang=2)
    plt.text(np.shape(tx)[0] * 1.25, 0 - np.shape(tx)[1] * 0.15, "%d pixel displacement" %  30, fontsize=10,
             color="black")
    fig.savefig(os.path.join(out_folder2, pos + "traktion.png"))
    plt.close()

    #saving gif
    images = [plt.imread(file1),plt.imread(file2)]
    imageio.mimsave(os.path.join(out_folder2, pos + ".gif"), images)


    # std filtering image
    fig=plt.figure()
    plt.imshow(mask_std)
    fig.savefig(os.path.join(out_folder2, pos + "std_filtering.png"))
    plt.close()

### energy calculation and gui









