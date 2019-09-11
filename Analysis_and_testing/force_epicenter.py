import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.insert(0, '/media/user/GINA1-BK/Andreas-Python/tracktion_force_microscopy')
from TFM_functions import *
from TFM_functions import *

def get_xy_for_quiver(u):
    xs=np.zeros(np.shape(u))
    for i in range(np.shape(u)[0]):
        xs[i,:]=np.arange(0,np.shape(u)[1],1)
    ys = np.zeros(np.shape(u))
    for j in range(np.shape(u)[1]):  ## is inverted in other skript...
        ys[:, j] = np.arange(0,np.shape(u)[0], 1)
    return xs, ys



pixelsize_org = 6.25 / 40
window_size = 64
overlapp = 32
pixel_factor = window_size - overlapp
pixelsize_defo =  pixelsize_org * pixel_factor # in µm
pixelsize_defo =1




folder="/home/user/Desktop/WTshift/01/"
#folder="/home/user/Desktop/test_for_negative_numbers/01/"
tx=np.load(os.path.join(folder,"tx.npy"))
ty=np.load(os.path.join(folder,"ty.npy"))
#mask=np.ones(np.shape(ty))>0



tx_filter = np.zeros(np.shape(tx))
tx_filter[mask] = tx[mask]

ty_filter = np.zeros(np.shape(ty))
ty_filter[mask] = ty[mask]

tract_abs=np.sqrt(tx_filter**2+ty_filter**2)

area=(pixelsize_defo*(10**-6))**2 # in meter
fx = tx_filter*area### calculating forces by multiplyin with area # forces then in newton
fy = ty_filter*area



## calculate forces wiht area??
#mask=np.ones(np.shape(tx))


x,y=get_xy_for_quiver(tx)


bx=np.sum(x*(tract_abs**2)+fx*(tx_filter*fx+ty_filter*fy))  ## meaning of ths vector??
by=np.sum(y*(tract_abs**2)+fy*(tx_filter*fx+ty_filter*fy))  #

axx=np.sum(tract_abs**2+fx**2)
axy=np.sum(fx*fy) # redundant
ayy=np.sum(tract_abs**2+fy**2)
#ayx=np.sum(tx*ty)

A=np.array([[axx,axy],[axy,ayy]])
b=np.array([bx,by]).T

# solve equation system:
#center*[bx,by]=[[axx,axy],
#                [axy,ayy]]
# given by A*[[1/bx],
#             [1/by]
center=np.matmul(np.linalg.inv(A),b)


# vector projection to origin

dist_x=center[0]-x
dist_y=center[1]-y
dist_abs=np.sqrt(dist_y**2+dist_x**2)
proj_abs=(fx*(dist_x)+fy*(dist_y))/dist_abs   ## minus because of image axis inversion(?)
contractile_force=np.sum(proj_abs)

#project_vectors
proj_x=proj_abs*dist_x/dist_abs
proj_y=proj_abs*dist_y/dist_abs


#np.sum(np.sqrt(proj_x**2+proj_y**2))
#--> gives slightly diffrent results because "negative values"


#contractile_force, proj_x, proj_y,center=contractility(tx,ty,pixelsize_defo,mask)

print(np.sum(proj_x))
print(np.sum(proj_y))

## illustration of projection
fig = plt.figure()
plt.imshow(np.sqrt((tx / 1000) ** 2 + (ty / 1000) ** 2),alpha=0.3)
plt.plot(center[0],center[1],color="red")
cbar = plt.colorbar()
cbar.set_label("traktion forces in kPa")
x1, y1 = get_xy_for_quiver(tx)
ratio=0.02 # ratio of length of biggest arrow to max axis lenght
scale=ratio *np.max(np.shape( proj_x))/np.max(np.sqrt( proj_x**2+ proj_y**2))# automatic sacleing in dependace of the image size

plt.quiver(x1, y1, proj_x * scale, proj_y  *scale, angles="xy", scale=1, scale_units='xy', width=0.002)
plt.arrow(75, -14, 30 * scale, 0, head_width=0,facecolor="black", edgecolor="black",
          clip_on=False, head_starts_at_zero=False, overhang=2)
plt.text(75, -10, "30 kPa Force", color="black")
plt.plot(center[0],center[1],"or",color="red")
plt.imshow(mask,alpha=0.1)
plt.text(10, 70, "contractile force = "+str(np.round(contractile_force*10**6,2)), color="black")


## illustration of distance arrows
fig = plt.figure()
plt.imshow(proj_abs,alpha=0.3)
plt.plot(center[0],center[1],color="red")
cbar = plt.colorbar()
cbar.set_label("projection")
x1, y1 = get_xy_for_quiver(tx)
ratio=0.02 # ratio of length of biggest arrow to max axis lenght
scale=ratio *np.max(np.shape( dist_x))/np.max(np.sqrt( dist_x**2+ dist_y**2))# automatic sacleing in dependace of the image size

plt.quiver(x1, y1, dist_x * scale, dist_y  *scale, angles="xy", scale=1, scale_units='xy', width=0.002)
plt.arrow(75, -14, 30 * scale, 0, head_width=0,facecolor="black", edgecolor="black",
          clip_on=False, head_starts_at_zero=False, overhang=2)
plt.text(75, -10, "30 kPa Force", color="black")
plt.plot(center[0],center[1],"or",color="red")
plt.imshow(mask,alpha=0.1)
plt.text(10, 70, "contractile force = "+str(np.round(contractile_force*10**6,2)), color="black")










## traktion forces
fig = plt.figure()
plt.imshow(np.sqrt((tx / 1000) ** 2 + (ty / 1000) ** 2))
plt.plot(center[0],center[1],color="red")
cbar = plt.colorbar()
cbar.set_label("traktion forces in kPa")
x1, y1 = get_xy_for_quiver(tx)
ratio=0.2 # ratio of length of biggest arrow to max axis lenght
scale=ratio *np.max(np.shape(tx))/np.max(np.sqrt( tx**2+ ty**2))# automatic sacleing in dependace of the image size

plt.quiver(x1, y1,tx*scale, ty  * scale, scale=1,angles="xy", scale_units='xy', width=0.002)
plt.arrow(75, -14, 30 / scale, 0, head_width=0, facecolor="black", edgecolor="black",
          clip_on=False, head_starts_at_zero=False, overhang=2)
plt.text(75, -10, "30 kPa Force", color="black")
plt.plot(center[0],center[1],"or",color="red")



