import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os

import openpiv.tools
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import numpy as np
import copy
import os
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
import time
import cv2 as cv
from imageio import imread
from skimage.filters import gaussian


from TFM_functions import *
out_file=r"E:\TFM\20190222_I\40x\IFBmshift.txt"    ## output folder

file1=r"E:\TFM\20190222_I\40x\IFBmshift\03before_shift.tif" # after image
file2=r"E:\TFM\20190222_I\40x\IFBmshift\03after_shift.tif" # before image
file3=r"E:\TFM\20190222_I\40x\IFBmshift\03bf_before_shift.tif" # bf image



sigma=0.49#poison ratio
young=25536# youngsmodulus
pixelsize = 6.25 / 40
window_size = 64
overlapp = 32
pixel_factor = window_size - overlapp
pixelsize = pixelsize * pixel_factor
arrow_treshold=1









bf_image=imread(file3)

print("calculating deformation")
u,v,x,y,mask,mask_std=calculate_deformation(file1,file2,window_size,overlapp)
print("calculating traktion")
tx,ty=ffttc_traction(u,v,young,pixelsize,bf_image=False,filter="mean")  # u and v shifted internally

tx=tx/1000
ty=ty/1000
t_abs=np.sqrt(tx ** 2 + ty ** 2)
u_shift = u - np.mean(u)
v_shift = v - np.mean(v)

## pre calulating energy on all points
print("calculating strain energies")
energy_points = 0.5 * pixelsize * pixelsize * (np.sqrt((tx * (u_shift / pixel_factor) * pixelsize) ** 2 + (
        ty * (v_shift / pixel_factor) * pixelsize) ** 2)) / 1000
bg = np.percentile(energy_points, 50)





if not os.path.isfile(out_file):
    with open((out_file), "w") as f:
        f.write("filename\tstrain_energy\tcontractile force\n")





### gui






custom_cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#DBDC3E","yellow"])
pixx = np.arange(np.shape(tx)[0])
pixy = np.arange(np.shape(tx)[1])
xv, yv = np.meshgrid(pixy, pixx)
pix = np.vstack((xv.flatten(), yv.flatten())).T


# displaying only every 9th arrow with as mean of sourrounding arrows


tx_show=gaussian(tx,sigma=3)
ty_show=gaussian(ty,sigma=3) # blocks of edge size 3
select_x=((xv-1)%3)==0
select_y=((yv-1)%3)==0
selct_size=t_abs>2
select=select_x*select_y*selct_size
tx_show[~select]=0
ty_show[~select]=0


ind = []
energy = 0
contractile_force=0
plt_str1 = 0
plt_str2 = 0

plt.close("all")
fig = plt.figure(figsize=(12,4))
ax1 = plt.axes([0.45, 0.25, 0.3, 0.6]) # ax image bf
ax2 = plt.axes([0.05, 0.25, 0.3, 0.6]) # ax img traction forces
ax_col=plt.axes([0.72, 0.25, 0.05, 0.6]) # ax for colorbar
ax_col.set_axis_off()
ax_text=plt.axes([0.45, 0.1, 0.1, 0.1]) # ax for colorbar
ax_text.set_axis_off()



im = ax1.imshow(t_abs,cmap="rainbow")
ratio=0.2 # ratio of length of biggest arrow to max axis lenght
scale=ratio *np.max(np.shape(tx))/np.max(np.sqrt(tx_show**2+ty_show**2))# automatic sacleing in dependace of the image size


cbar=fig.colorbar(mappable=im,ax=ax_col)
cbar.set_label('traction forces in kPa')
x_small, y_small = get_xy_for_quiver(tx)
ax1.quiver(x_small, y_small, tx_show*scale, ty_show*scale,scale=1, scale_units='xy')


ax2.imshow(bf_image,cmap="gray")






def onselect1(verts):
    global array, pix, ind, energy, plt_str1,plt_str2,ax_text,energy_points,bg,mask,verts_out,p,mask_zoom,ax1,mask_zoom2
    global contractile_force,center,select,tx_show,ty_show,scale
    verts_out = verts

    # simple interpolation of points, to save computation time
    interpol_factors = np.array([np.shape(tx)[1] / np.shape(bf_image)[1], np.shape(tx)[0] / np.shape(bf_image)[0]])
    verts_array=np.array(verts_out)
    verts_interpol=np.zeros(np.shape( verts_array))
    for i in range(np.shape(verts_interpol)[0]):
        verts_interpol[i,:]=verts_array[i,] * interpol_factors  # element wise multiplikatioon for each row
    p = path.Path(verts_interpol)

    #retrieving mask
    ind = p.contains_points(pix, radius=1)
    print(ind)  # boolean mask
    mask = np.reshape(np.array(ind), (np.shape(tx)))

    #showing mask as overlay
    mask_show=np.zeros(np.shape(tx))+np.nan
    mask_show[mask]=1

    ax1.clear()
    ax1.imshow(t_abs,alpha=0.8,cmap="rainbow")
    x_small, y_small = get_xy_for_quiver(tx)
    ax1.quiver(x_small, y_small, tx_show*scale, ty_show*scale,scale=1, scale_units='xy')
    ax1.imshow(mask_show,alpha=0.4,cmap=custom_cmap1)

    #energy caclulation
    energy = np.sum(energy_points[mask]) - bg * np.sum(mask)
    print("strain energy" ,energy)

    # contractile forces calculation
    contractile_force, proj_x, proj_y, center = contractility(tx*1000, ty*1000, pixelsize, mask)  # needs tx in Pa


    print("contractile force ", contractile_force)
    ax1.plot(center[0],center[1],"or",color="red")

    # writing energy text
    if plt_str1 != 0:
        plt_str1.remove()
    plt_str1 = ax_text.text(0, 0, "strain energy = " + str(np.round(energy, 2))+" pJ")

    if plt_str2 != 0:
        plt_str2.remove()
    plt_str2 = ax_text.text(0, 0.5, "contractile force = " + str(np.round(contractile_force*10**6, 2)) + " µN")

    fig.canvas.draw_idle()


lasso1 = LassoSelector(ax2, onselect1)





# closing button
def close(event):
    plt.close("all")


button_ax_c = plt.axes([0.7, 0.1, 0.1, 0.075])
button_ax_c.xaxis.set_ticks_position('none')
button_ax_c.yaxis.set_ticks_position('none')
button_ax_c.set_xticks([])
button_ax_c.set_yticks([])
b_close = Button(button_ax_c, 'close')
b_close.on_clicked(close)

# saving button

def save(event):
    global energy, out_file
    with open((out_file), "a") as f:
        f.write(file1 + "\t" + str(np.round(energy, 2)) +"\t"+ str(np.round(-contractile_force*10**6, 2))+"\n")
    print("saved line " + file1 + "\t" + str(np.round(energy, 2)) +"\t"+ str(np.round(-contractile_force*10**6, 2))+"\n")


button_ax_s = plt.axes([0.8, 0.1, 0.1, 0.075])
button_ax_s.xaxis.set_ticks_position('none')
button_ax_s.yaxis.set_ticks_position('none')
button_ax_s.set_xticks([])
button_ax_s.set_yticks([])
b_save = Button(button_ax_s, 'save')
b_save.on_clicked(save)








