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
from matplotlib.widgets import LassoSelector
from matplotlib import path
import os
import time


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


folder="/media/user/GINA1-BK/traktion_force_microscopy/analysis_with_custom_program/"

tx_filter=np.load(os.path.join(folder,"tx_filter.npy"))/1000
ty_filter=np.load(os.path.join(folder,"ty_filter.npy"))/1000

u=np.load(os.path.join(folder,"u.npy"))   ## dx
v=np.load(os.path.join(folder,"v.npy"))    ##dy
## one could further filter these arrows with a maens filter


u_shift=u-np.mean(u)
v_shift=v-np.mean(v)

# selecting mask for eneergy calculation
fig = plt.figure()
plt.imshow(np.sqrt(tx_filter**2+ty_filter**2))
plt.colorbar()
x_small,y_small=get_xy_for_quiver(tx_filter)
plt.quiver(x_small,y_small,tx_filter,-ty_filter)
ax1=plt.gca()

pixx = np.arange(np.shape(u_shift)[0])
pixy=np.arange(np.shape(u_shift)[1])
xv, yv = np.meshgrid(pixx ,pixy)
pix = np.vstack( (xv.flatten(), yv.flatten()) ).T
ind=[]

def onselect(verts):
    global array, pix, ind
    p = path.Path(verts)
    ind = p.contains_points(pix, radius=1)
    print(ind)   # boolean mask
    fig.canvas.draw_idle()

lasso = LassoSelector(ax1, onselect)
plt.waitforbuttonpress()
plt.waitforbuttonpress()
mask=np.reshape(np.array(ind),np.shape(u_shift))



# energy points on all pixels E_p=0.5*t*deformation , Energy then by integration (sum)
#what exactely is calculated here??
# sum of all absolute values of energies?
energy_points=0.5*pixelsize*pixelsize*(np.sqrt((tx_filter*(u_shift/pixel_factor)*pixelsize)**2+(ty_filter*(v_shift/pixel_factor)*pixelsize)**2))/1000
# devision of u and v py pixel_factor to get pixel shift in terms of current number of pixels

### reconsider the units please


bg=np.percentile(energy_points,50) # calculating background
#unit: pJ



#bg=0.2
#plt.figure()
#plt.hist(energy_points.flatten(),bins=100)### usually very sharp background
#integration

energy= np.sum(energy_points[mask])-bg*np.sum(mask)
# unit is allegidy pj???
print(energy)




















'''
# potting
plt.figure()
plt.imshow(np.sqrt(v**2+u**2))
plt.colorbar()

x_small,y_small=get_xy_for_quiver(u)
plt.quiver(x_small,y_small,u,v)
'''