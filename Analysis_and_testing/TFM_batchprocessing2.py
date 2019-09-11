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
in_folder=r"E:\TFM\20190222_I\40x\WTshift"
out_folder=r"E:\TFM\20190222_I\40x\WTshift" ## output folder


min_scale_value_tract=0
max_scale_value_tract=70
min_scale_value_disp=0
max_scale_value_disp=20
sigma=0.49#poison ratio
young=25536# youngsmodulus
pixelsize = 6.25 / 40
window_size=64
overlapp=32
pixel_factor = window_size - overlapp
#file1="/media/user/GINA1-BK/traktion_force_microscopy/new_data_from_magdalena26_02_19/WTshift/12after_shift.tif"
#file2="/media/user/GINA1-BK/traktion_force_microscopy/new_data_from_magdalena26_02_19/WTshift/12before_shift.tif"
#file3="/media/user/GINA1-BK/traktion_force_microscopy/new_data_from_magdalena26_02_19/WTshift/12bf_before_shift.tif"

# listing files
files_dict={}
file_list=os.listdir(in_folder)
file_list=[x for x in file_list if "shift.tif" in x]
pos_list = [re.search("(\d{1,3})\D.*", x).group(1) for x in file_list]  # just first two positions (maybe use three
pos_list=np.unique(pos_list)


for pos in pos_list:
    files_dict[pos] = [os.path.join(in_folder, file) for file in file_list if re.match(pos+"\D.*",file)]




for pos,files in tqdm(files_dict.items(),total=len(files_dict.items())):
    #if not pos ==12:
    #    continue
    print(pos)
    file2=[x for x in files if "after" in os.path.split(x)[1]][0] # cahnged order
    file1=[x for x in files if "before" in os.path.split(x)[1]][0]
    file=[x for x in files if "bf_before" in os.path.split(x)[1]][0]
    print("calculating deformation")
    u, v, x, y, mask = calculate_deformation(file1, file2, window_size, overlapp)
    print("calculating traktion")
    tx, ty = ffttc_traction(u, v, young, pixelsize, bf_image=False, filter="mean")

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
    plt.arrow(75, -14, 10/scale, 0, head_width=0, facecolor="black", edgecolor="black",
              clip_on=False, head_starts_at_zero=False, overhang=2)
    plt.text(75, -10, "10 pixel displacement", color="black")
    fig.savefig(os.path.join(out_folder2,pos+"deformation.png"))
    plt.close()

    # traction image:


    fig=plt.figure()
    plt.imshow(np.sqrt((tx/1000) ** 2 + (ty/1000) ** 2),vmin=min_scale_value_tract, vmax=max_scale_value_tract)
    cbar = plt.colorbar()
    cbar.set_label("traktion forces in kPa")
    x1, y1 = get_xy_for_quiver(tx)
    scale=2# #
    plt.quiver(x1, y1, (tx/1000) / scale, (ty/1000) / scale, scale=1, scale_units='xy',width=2)
    plt.arrow(75, -14, 30 / scale, 0, head_width=0, facecolor="black", edgecolor="black",
              clip_on=False, head_starts_at_zero=False, overhang=2)
    plt.text(75, -10, "30 kPa Force", color="black")
    fig.savefig(os.path.join(out_folder2, pos + "traktion.png"))
    plt.close()

    #saving gif
    images = [plt.imread(file1),plt.imread(file2)]
    imageio.mimsave(os.path.join(out_folder2, pos + ".gif"), images)

### energy calculation and gui





### energy calculation



'''

date = datetime.datetime.now()  # this constructsa string containing the date
date_str = "_%d_%d_%d_%.2d%.2d%.2d" % (date.day, date.month, date.year, date.hour, date.minute, date.second)
output_file = os.path.join(out_folder, date_str + "output.txt")
with open((output_file), "w") as f:
    f.write("filename\tstrain_energy\n")


tx=tx/1000
ty=ty/1000
u_shift = u - np.mean(u)
v_shift = v - np.mean(v)

# gui


plt.close("all")
fig = plt.figure(figsize=(12,4))
ax1 = plt.axes([0.45, 0.25, 0.3, 0.6]) # ax image bf
ax2 = plt.axes([0.05, 0.25, 0.3, 0.6]) # ax img traction forces
ax_col=plt.axes([0.65, 0.25, 0.1, 0.6]) # ax for colorbar
ax_col.set_axis_off()
ax_text=plt.axes([0.45, 0.1, 0.1, 0.1]) # ax for colorbar
ax_text.set_axis_off()
im = ax1.imshow(np.sqrt(tx ** 2 + ty ** 2))
fig.colorbar(mappable=im,ax=ax_col)

x_small, y_small = get_xy_for_quiver(tx)
ax1.quiver(x_small, y_small, tx, ty)


ax2.imshow(bf_image)


pixx = np.arange(np.shape(u_shift)[0])
pixy = np.arange(np.shape(u_shift)[1])
xv, yv = np.meshgrid(pixy, pixx)
pix = np.vstack((xv.flatten(), yv.flatten())).T
ind = []
energy = 0
plt_str = 0


def onselect(verts):
    global array, pix, ind, energy, plt_str,ax_text,energy_points,mask,verts_out,p
    p = path.Path(verts)
    verts_out=verts
    ind = p.contains_points(pix, radius=1)
    print(ind)  # boolean mask
    mask = np.reshape(np.array(ind), (np.shape(u_shift)))

    energy_points = 0.5 * pixelsize * pixelsize * (np.sqrt((tx * (u_shift / pixel_factor) * pixelsize) ** 2 + (
            ty * (v_shift / pixel_factor) * pixelsize) ** 2)) / 1000
    bg = np.percentile(energy_points, 50)
    energy = np.sum(energy_points[mask]) - bg * np.sum(mask)
    print(energy)


    # writing energy text
    if plt_str != 0:
        plt_str.remove()
    plt_str = ax_text.text(0, 0, "strain energy = " + str(np.round(energy, 2))+" pJ")
    fig.canvas.draw_idle()


lasso = LassoSelector(ax1, onselect)


# closing button
def close(event):
    plt.close("all")


button_ax_c = plt.axes([0.6, 0.1, 0.1, 0.075])
button_ax_c.xaxis.set_ticks_position('none')
button_ax_c.yaxis.set_ticks_position('none')
button_ax_c.set_xticks([])
button_ax_c.set_yticks([])
b_close = Button(button_ax_c, 'close')
b_close.on_clicked(close)


# saving button

def save(event):
    global energy, output_file
    with open((output_file), "a") as f:
        f.write(file1 + "\t" + str(np.round(energy, 2)) + "\n")
    print("saved line " + file1 + "\t" + str(np.round(energy, 2)) + "\n")


button_ax_s = plt.axes([0.8, 0.1, 0.1, 0.075])
button_ax_s.xaxis.set_ticks_position('none')
button_ax_s.yaxis.set_ticks_position('none')
button_ax_s.set_xticks([])
button_ax_s.set_yticks([])
b_save = Button(button_ax_s, 'save')
b_save.on_clicked(save)

'''






