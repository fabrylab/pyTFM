import numpy as np
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
import sys
sys.path.insert(0, '/media/user/GINA1-BK/Andreas-Python/tracktion_force_microscopy/')
from TFM_functions import *

from skimage.filters import gaussian
def get_xy_for_quiver(u):
    xs=np.zeros(np.shape(u))
    for i in range(np.shape(u)[0]):
        xs[i,:]=np.arange(0,np.shape(u)[1],1)
    ys = np.zeros(np.shape(u))
    for j in range(np.shape(u)[1]):  ## is inverted in other skript...
        ys[:, j] = np.arange(0,np.shape(u)[0], 1)
    return xs, ys



def force_from_center(Force_magnitude,pos1,pos2):
    new_center = np.array([np.mean([pos1[0], pos2[0]]), np.mean([pos1[1], pos2[1]])])  # new exact center position
    F1 = (pos1 - new_center) * Force_magnitude / np.linalg.norm(pos1 - new_center)  # first force as vector
    F2 = (pos2 - new_center) * Force_magnitude / np.linalg.norm(pos2 - new_center)  # second  force
    return F1,F2
#plt.close("all")

folder="/home/user/Desktop"


## gel parameter
sigma=0.49#poison ratio
young=25536# youngs modulus in Pa??
res=2 ## resolution of applied force grid

pixelsize =5 # in µm typicle value for second image
pixel_number=200 # number of pixel ain x and y direction, only squared matrices
edge_dist=0.1# distance for masking when caclulating contractile force as fraction. edges have inherrent background signals
# force paramteres
# two opsoing forces with distnace , centered at center and turned at angel "orientation"
center=np.array([100.5,100.5]) #in pixel
orientation=np.pi*0.25     ## only 45° angles are useful???   --> other wise problems with forces set between pixels
distance= 30 #in µm
Force_magnitude=2.8 #5000 *(10**-9) ## sum of all applied force


n=0 # number of forces
for i in np.arange(0,1.0,1/res):
    for j in np.arange(0, 1.0, 1/res):
        if not (i-int(i)==0 or j-int(j)==0):
            n+=1
n*=2

def_abs= np.zeros((pixel_number,pixel_number))
def_x=np.zeros((pixel_number,pixel_number))
def_y=np.zeros((pixel_number,pixel_number))



## setting up opsing foces
area=(pixelsize*(10**-6))**2
traction=Force_magnitude / area # corresponding force in Pa
print("input traction = "+str(np.round(traction))+" Pa" )

pos1=center+np.array((np.cos(orientation)*distance/2, np.sin(orientation)*distance/2)) # calculating two points from the center
pos2=center-np.array((np.cos(orientation)*distance/2, np.sin(orientation)*distance/2))



F_pos=[]
F=[]
for i in np.arange(0,1.0,1/res):
    for j in np.arange(0, 1.0, 1/res):
        if not (i-int(i)==0 or j-int(j)==0):
            pos1s = np.array([np.round(pos1[0]) + i, np.round(pos1[1]) + j])
            pos2s = np.array([np.round(pos2[0]) + i, np.round(pos2[1]) + j])
            F1,F2=force_from_center(Force_magnitude/n,pos1s,pos2s)
            F.extend([F1,F2])
            F_pos.extend([pos1s,pos2s])



A=(1+sigma)/(np.pi*young)

## iterating over forces, and summing up deformation
# deformation  in µm
dist_x=np.zeros((pixel_number,pixel_number))
dist_y=np.zeros((pixel_number,pixel_number))


def_x_pos,def_y_pos=get_xy_for_quiver(def_x)# positions for quiver

for f,f_pos in tqdm(zip(F,F_pos),total=len(F)):

    x=(def_x_pos-f_pos[0])*pixelsize*(10**-6)
    y=(def_y_pos-f_pos[1])*pixelsize*(10**-6)
    r=np.sqrt(x**2+y**2)
    K1=((1-sigma)*(r**2)+sigma*(x**2))*(A/(r**3))
    K2= sigma*x*y*(A/(r**3))
    K3=  ((1-sigma)*(r**2)+sigma*(y**2))*(A/(r**3))
    def_x+=f[0]*K1+f[1]*K2
    def_y+=f[0]*K2+f[1]*K3



def_abs= np.sqrt(def_x**2+def_y**2)




# transform deformation to pixels
#
def_x_pix=def_x/(pixelsize*(10**-6))## deformation in terms of pixel
def_y_pix=def_y/(pixelsize*(10**-6))
def_abs_pix=def_abs/(pixelsize*(10**-6))


tx,ty=ffttc_traction(def_x_pix,def_y_pix,pixelsize,pixelsize,young,sigma=0.49,bf_image=False,filter="gaussian")
t_abs=np.sqrt(tx**2+ty**2)
np.save(os.path.join(folder,"u_simu.npy"),def_x)
np.save(os.path.join(folder,"v_simu.npy"),def_y)
np.save(os.path.join(folder,"tx_simu.npy"),tx)
np.save(os.path.join(folder,"ty_simu.npy"),ty)

mask=np.ones(tx.shape)<0 ## contractile force on full image, this is not to robust
mask[int(np.shape(tx)[0]*edge_dist):-int(np.shape(tx)[0]*edge_dist),int(np.shape(tx)[1]*edge_dist):-int(np.shape(tx)[1]*edge_dist)]=True

contractile_force,proj_x,proj_y,center=contractility(tx,ty,pixelsize,mask)

#show_quiver(proj_x,proj_y,scale_ratio=0.2)



# sum of input forces as estimate of contractility
e=np.sum(np.sqrt(np.sum(np.array(F)**2,axis=1))) # simple sum of forces


# exact cotractility with force epicenter method
'''
fx=np.zeros((pixel_number*res,pixel_number*res))# 10 is resolution
fy=np.zeros((pixel_number*res,pixel_number*res))
mask_f=np.ones(np.shape(fx))>0
for f,f_pos in zip(F,F_pos):
    fx[int(f_pos[1] * res), int(f_pos[0] * res)] = f[0]/(area) ## contractility functions needs pressures
    fy[int(f_pos[1] * res), int(f_pos[0] * res)] = f[1]/(area)

e_contract, proj_x1, proj_y1, center_org = contractility(fx, fy, pixelsize, mask_f)
'''













# plotting deformation field

def_x_pos,def_y_pos=get_xy_for_quiver(def_x)# positions for quiver
plt.figure()
plt.imshow(def_abs,cmap="jet")

scale_ratio=0.1
scale = scale_ratio * np.max(np.shape(def_x)) /Force_magnitude  # automatic scaleing in dependance of the image size
for f,f_pos in zip(F,F_pos):
    plt.arrow(f_pos[0], f_pos[1],f[0] *scale, f[1] *scale, head_width=0.8, facecolor="red", width=0.01
              , edgecolor="red", head_starts_at_zero=False, overhang=0, head_length=0.5)

plt.title("simulation of deformation field")
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
cbar=plt.colorbar()
cbar.set_label("displacement in m",fontsize=15)

scale_ratio=0.1
scale = scale_ratio * np.max(np.shape(def_y)) / np.max(def_abs)  # automatic scaleing in dependance of the image size
Q=plt.quiver(def_x_pos,def_y_pos,
             def_x*scale,def_y*scale,headwidth=3,angles='xy', scale_units='xy',
             scale=1,width=0.0025,color="black")

#
#plt.arrow(85, -10, 3*scale, 0, head_width=0, facecolor="red", edgecolor="green",
#          clip_on=False, head_starts_at_zero=False,overhang=0)
#plt.text(86,-2,"displacement arrows\nlength 3 Pixel",color="green")





























# traction forces
plt.figure()
plt.title("prediction of traction forces")
plt.imshow(t_abs/1000,cmap="jet")
plt.plot(center[0],center[1],"or",color="red")

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
cbar=plt.colorbar()
cbar.set_label("displacement in pixel",fontsize=15)

cbar.set_label("Traction in kPa")
x_small,y_small=get_xy_for_quiver(tx)
scale_ratio=0.2
scale = scale_ratio * np.max(np.shape(tx)) / np.max(t_abs/1000)
plt.quiver(x_small,y_small,(tx/1000)*scale ,(ty/1000)*scale , angles='xy', scale_units='xy', color="red",
           scale=1,headwidth = 4, width=0.002)

plt.arrow(85, -10, 1*scale, 0, head_width=0, facecolor="red", edgecolor="red",
          clip_on=False, head_starts_at_zero=False, overhang=2)
plt.text(86,-2,"tractiont arrows\nlength 1 kPa",color="red")






## illustration of projection
plt.figure()
plt.imshow(t_abs/1000,alpha=0.3)
plt.plot(center[0],center[1],"or",color="red")
#plt.plot(center_org[0]/res,center_org[1]/res,"+",color="green")
cbar = plt.colorbar()
cbar.set_label("traktion forces in kPa")
x1, y1 = get_xy_for_quiver(tx)
ratio=0.2 # ratio of length of biggest arrow to max axis lenght
scale=ratio *np.max(np.shape( proj_x))/np.max(np.sqrt( proj_x**2+ proj_y**2))# automatic sacleing in dependace of the image size

plt.quiver(x1, y1, proj_x * scale, proj_y  *scale, scale=1,angles="xy", scale_units='xy', width=0.002)
#plt.arrow(75, -14, 0.1 * scale, 0, head_width=0, facecolor="black", edgecolor="black",
 #         clip_on=False, head_starts_at_zero=False, overhang=2)
#plt.text(75, -10, "0.1 kPa Force", color="black")
#plt.plot(center[0],center[1],"or",color="red")
plt.imshow(mask,alpha=0.1)
plt.text(np.shape(def_abs)[0]*0.1, np.shape(def_abs)[0]*0.7, "contractile force prediction = "+str(np.round(contractile_force*10**9,2))+"nN", color="black")
#ecpect valu of the contractile force:

plt.text(np.shape(def_abs)[0]*0.1, np.shape(def_abs)[0]*0.8, "sum of input forces= "+str(np.round(e*10**9,2))+"nN", color="black")
#plt.text(np.shape(def_abs)[0]*0.1, np.shape(def_abs)[0]*0.9, "contractile force of iput forces= "+str(np.round(e_contract*10**9,2))+"nN", color="black")






