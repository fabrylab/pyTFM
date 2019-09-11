import openpiv.tools
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import numpy as np
import copy
import os
from scipy.ndimage import zoom
from skimage.filters import rank
from skimage.morphology import cube
from scipy.ndimage.filters import uniform_filter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import LassoSelector,Button
from matplotlib import path
import os
import time
import datetime




def get_xy_for_quiver(u):
    xs=np.zeros(np.shape(u))
    for i in range(np.shape(u)[0]):
        xs[i,:]=np.arange(0,np.shape(u)[1],1)
    ys = np.zeros(np.shape(u))
    for j in range(np.shape(u)[1]):  ## is inverted in other skript...
        ys[:, j] = np.arange(0,np.shape(u)[0], 1)
    return xs, ys


pixelsize=6.45/40
window_size=64
overlapp=32
pixel_factor=window_size-overlapp
pixelsize=pixelsize*pixel_factor  ### u is given in terms of original pixe


file1="/media/user/GINA1-BK/traktion_force_microscopy/new_data_from_magdalena26_02_19/WTshift/11after_shift.tif"
file1="/media/user/GINA1-BK/traktion_force_microscopy/new_data_from_magdalena26_02_19/WTshift/12after_shift.tif"
folder1="/media/user/GINA1-BK/traktion_force_microscopy/new_data_from_magdalena26_02_19/WTshift/11"
folder2="/media/user/GINA1-BK/traktion_force_microscopy/new_data_from_magdalena26_02_19/WTshift/12"


# makeing output file

date = datetime.datetime.now()  # this constructsa string containing the date
date_str = "_%d_%d_%d_%.2d%.2d%.2d" % (date.day, date.month, date.year, date.hour, date.minute, date.second)
output_file=os.path.join(folder1,date_str+"output.txt")
with open((output_file),"w") as f:
    f.write("filename\tstrain_energy\n")




tx1_filter=np.load(os.path.join(folder1,"tx.npy"))/1000
ty1_filter=np.load(os.path.join(folder1,"ty.npy"))/1000
u1=np.load(os.path.join(folder1,"u.npy"))   ## dx
v1=np.load(os.path.join(folder1,"v.npy"))    ##dy

tx2_filter=np.load(os.path.join(folder2,"tx.npy"))/1000
ty2_filter=np.load(os.path.join(folder2,"ty.npy"))/1000
u2=np.load(os.path.join(folder2,"u.npy"))   ## dx
v2=np.load(os.path.join(folder2,"v.npy"))    ##dy




## one could further filter these arrows with a maens filter


u1_shift=u1-np.mean(u1)
v1_shift=v1-np.mean(v1)

# selecting mask for eneergy calculation
plt.close("all")
fig = plt.figure()
ax1=plt.axes([0,0.3,0.8,0.6])
im=ax1.imshow(np.sqrt(tx1_filter**2+ty1_filter**2))
fig.colorbar(im)
x_small,y_small=get_xy_for_quiver(tx1_filter)
ax1.quiver(x_small,y_small,tx1_filter,-ty1_filter)

pixx = np.arange(np.shape(u1_shift)[0])
pixy=np.arange(np.shape(u1_shift)[1])
xv, yv = np.meshgrid(pixx ,pixy)
pix = np.vstack( (xv.flatten(), yv.flatten()) ).T
ind=[]
energy=0
plt_str=0
index=0

def onselect(verts):
    global array, pix, ind, energy, plt_str
    p = path.Path(verts)
    ind = p.contains_points(pix, radius=1)
    print(ind)   # boolean mask
    mask = np.reshape(np.array(ind), np.shape(u1_shift))

    energy_points = 0.5 * pixelsize * pixelsize * (np.sqrt((tx1_filter * (u1_shift / pixel_factor) * pixelsize) ** 2 + (
                ty1_filter * (v1_shift / pixel_factor) * pixelsize) ** 2)) / 1000
    bg = np.percentile(energy_points, 50)
    energy = np.sum(energy_points[mask]) - bg * np.sum(mask)
    print(energy)

    # writing energy text
    if plt_str!=0:
        plt_str.remove()
    plt_str=plt.text(-7,-5.5,"strain energy = "+str(np.round(energy,2)))
    fig.canvas.draw_idle()

lasso = LassoSelector(ax1, onselect)


#closing button
def close(event):
    plt.close("all")

button_ax_c=plt.axes([0.8, 0.8, 0.1, 0.075])
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
        f.write(file1+"\t"+str(np.round(energy,2))+"\n")
    print("saved line "+ file1+"\t"+str(np.round(energy,2))+"\n")

button_ax_s = plt.axes([0.8, 0.6, 0.1, 0.075])
button_ax_s.xaxis.set_ticks_position('none')
button_ax_s.yaxis.set_ticks_position('none')
button_ax_s.set_xticks([])
button_ax_s.set_yticks([])
b_save = Button(button_ax_s, 'save')
b_save.on_clicked(save)

# next button

def next(event):
    global energy, output_file,ax1,index
    index+=1
    ax1.clear()
    ax1 = plt.axes([0, 0.3, 0.8, 0.6])

    im = ax1.imshow(np.sqrt(tx2_filter ** 2 + ty2_filter ** 2))
    fig.colorbar(im)
    x_small, y_small = get_xy_for_quiver(tx1_filter)
    ax1.quiver(x_small, y_small, tx2_filter, -ty2_filter)
button_ax_n = plt.axes([0.8, 0.4, 0.1, 0.075])
button_ax_n.xaxis.set_ticks_position('none')
button_ax_n.yaxis.set_ticks_position('none')
button_ax_n.set_xticks([])
button_ax_n.set_yticks([])
b_next = Button(button_ax_n, 'next')
b_next.on_clicked(next)


