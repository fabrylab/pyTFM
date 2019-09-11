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
from mpl_toolkits.axes_grid1 import make_axes_locatable
import imageio
from tqdm import tqdm
import re
from skimage.filters import gaussian

from TFM_functions import *
try:
    from functions_for_cell_colonie import show_quiver,filter_values
except:
    pass


def get_scale_bar_length(fig, im, scale, bar_length):
    '''
    accesory funcion to get the length of the scale bar in terms of
    :param fig: figure object
    :param im: object from imshow
    :param scale: scale used to rescale the quiver arrows
    :param bar_length: length of the scale bar in figure coordinates
    :return:
    '''
    # transform scale bar length in figure coordinates (0.1) to axis coordinates
    f = fig.transFigure.transform([bar_length, bar_length]) - fig.transFigure.transform([0, 0])
    k = im.axes.transData.inverted().transform(f) - im.axes.transData.inverted().transform([0, 0])
    # apply scaling factor that is used for quiver arrows
    l = k[0] / scale
    return np.round(l, 1)

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


in_folder=r"/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/"
out_folder=r"/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/out_clickpoints" ## output folder


min_scale_value_tract=0
max_scale_value_tract=2
min_scale_value_disp=0
max_scale_value_disp=20

sigma=0.49#poison ratio
young=25536# youngsmodulus
pixelsize = 4.09 / 40 #µm/pixel
window_size=100   # in pixel , keep at 10 µm
overlapp=85   # window_size/2
std_factor=15
pixel_factor = window_size - overlapp
outputmode=2 # 1 makeing individual folders, 2 throwing everything in one folder
h=300


# listing files
files_dict={}
file_list=os.listdir(in_folder)
file_list=[x for x in file_list if "shift.tif" in x]
pos_list = [re.search("(\d{1,3})\D.*", x).group(1) for x in file_list]  # just first two positions (maybe use three
pos_list=np.unique(pos_list)


for pos in pos_list:
    files_dict[pos] = [os.path.join(in_folder, file) for file in file_list if re.match(pos+"\D.*",file)]

plt.ioff()
for pos,files in tqdm(files_dict.items(),total=len(files_dict.items())):
    #if not pos =="11":
    #   continue
    print(pos)
    file1=[x for x in files if "after" in os.path.split(x)[1]][0] # chnged order
    file2=[x for x in files if "before" in os.path.split(x)[1] and not "bf" in os.path.split(x)[1]][0]
    file=[x for x in files if "bf_before" in os.path.split(x)[1]][0]
    print("calculating deformation")
    u, v, x, y, mask,mask_std = calculate_deformation(file1, file2, window_size, overlapp,  std_factor=std_factor)

    tx, ty = ffttc_traction_finite_thickness(u, v,pixelsize,pixelsize*pixel_factor,h=h,young=young,sigma=0.49, filter="gaussian") ## try this at one point....

    #saving output and images
    #seting folder structure
    if outputmode==1:
        out_folder2=os.path.join(in_folder,pos)
        createFolder( out_folder2)
        pos_name=""

    else:
        out_folder2=out_folder
        createFolder( out_folder2)
        pos_name=str(pos)


    # saving output data
    np.save(os.path.join(out_folder2, pos_name+"u.npy"),u)
    np.save(os.path.join(out_folder2, pos_name+"v.npy"), v)
    np.save(os.path.join(out_folder2, pos_name+"x.npy"),x)
    np.save(os.path.join(out_folder2, pos_name+"y.npy"),y)
    np.save(os.path.join(out_folder2, pos_name+"mask.npy"),mask)
    #np.save(os.path.join(out_folder2, pos_name+"tx_h.npy"),tx_h)
    #np.save(os.path.join(out_folder2, pos_name+"ty_h.npy"),ty_h)
    np.save(os.path.join(out_folder2, pos_name+"tx_h.npy"),tx)
    np.save(os.path.join(out_folder2, pos_name+"ty_h.npy"),ty)


    pixx = np.arange(np.shape(tx)[0])
    pixy = np.arange(np.shape(tx)[1])
    xv, yv = np.meshgrid(pixy, pixx)


    # defo image:

    # displaying only every 9th arrow  as mean of sourrounding arrows

    def_abs =np.sqrt((u**2+v**2))
    u_show = gaussian(u, sigma=3)  # u_show = copy.deepcopy(u)  for unfiltered arrows
    v_show = gaussian(v, sigma=3)  # v_show = copy.deepcopy(v) for unfiltered arrows
    select_x = ((xv - 1) % 1) == 0
    select_y = ((yv - 1) % 1) == 0
    selct_size = def_abs > 1
    select = select_x * select_y * selct_size
    u_show[~select] = 0
    v_show[~select] = 0

    fig=plt.figure()
    im=plt.imshow(def_abs,vmin=min_scale_value_disp,vmax=max_scale_value_disp,cmap="rainbow")
    plt.axis("off")
    x1,y1=get_xy_for_quiver(u)
    ratio = 0.2  # ratio of length of biggest arrow to max axis lenght
    scale = ratio * np.max(np.shape(u)) / np.max(np.sqrt((u_show**2+v_show**2)))  # automatic sacleing in dependace of the image size
    plt.quiver(x1,y1,u_show*scale,v_show*scale,scale=1,headwidth=3, scale_units='xy',width=0.002, angles='xy')

    # arrow now with fixed position and lenght with respect ot figure
    plt.arrow(0.75,0.96, 0.2 , 0, head_width=0, facecolor="black",
              edgecolor="black",
              clip_on=False, head_starts_at_zero=False, overhang=2, transform=fig.transFigure)
    plt.text(np.shape(u)[0] * 1, 0 - np.shape(u)[1] * 0.05, "%.1f pixel displacement" % get_scale_bar_length(fig,im,scale,0.2), fontsize=10,
             color="black")

    ## making colorbar and stuff
    ax = plt.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(mappable=im, cax=cax)
    cbar.set_label("displacement in pixel")

    fig.savefig(os.path.join(out_folder2,pos+"deformation.png"))
    plt.close()



    # traction image:

    # displaying only every 9th arrow with as mean of sourrounding arrows

    t_abs=np.sqrt((tx/1000) ** 2 + (ty/1000) ** 2)
    tx_show = gaussian(tx/1000, sigma=3) # tx_show = copy.deepcopy(tx/1000)  for unfiltered arrows
    ty_show = gaussian(ty/1000, sigma=3) # ty_show = copy.deepcopy(ty/1000) for unfiltered arrows
    select_y = ((yv - 1) % 1) == 0
    select_x = ((xv - 1) % 1) == 0


    selct_size = t_abs > 0.1
    select = select_x * select_y * selct_size
    tx_show[~select] = 0
    ty_show[~select] = 0

    fig=plt.figure()
    im=plt.imshow(t_abs,vmin=min_scale_value_tract, vmax=max_scale_value_tract,cmap="rainbow")
    plt.axis("off")
    x1, y1 = get_xy_for_quiver(tx)
    ratio = 0.2  # ratio of length of biggest arrow to max axis lenght
    scale = ratio * np.max(np.shape(tx_show)) / np.max(np.sqrt((tx_show) ** 2 + (ty_show) ** 2))  # automatic sacleing in dependace of the image size
    plt.quiver(x1, y1, tx_show * scale, ty_show * scale, scale=1, headwidth=3, scale_units='xy', width=0.002, angles='xy')

    # arrow now with fixed position and lenght with respect ot figure
    plt.arrow(0.75, 0.96, 0.2, 0, head_width=0, facecolor="black",
              edgecolor="black",
              clip_on=False, head_starts_at_zero=False, overhang=2, transform=fig.transFigure)
    plt.text(np.shape(u)[0] * 1, 0 - np.shape(u)[1] * 0.05,
             "%.1f kPa traction force" % get_scale_bar_length(fig, im, scale, 0.2), fontsize=10,
             color="black")

    ## making colorbar and stuff
    ax = plt.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(mappable=im, cax=cax)
    cbar.set_label("traction force in kPa")
    fig.savefig(os.path.join(out_folder2, pos + "traktion.png"))
    plt.close()

    #saving gif
    images = [plt.imread(file1),plt.imread(file2)]
    imageio.mimsave(os.path.join(out_folder2, pos + ".gif"), images,fps=2)

### energy calculation and gui







