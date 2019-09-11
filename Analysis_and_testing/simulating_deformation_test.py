import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.insert(0, '/home/user/Desktop/Andreas-Python/tracktion_force_microscopy/')
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


plt.close("all")
# single force

sigma=0.49#poison ratio
young=25536# youngsmodulus
pixelsize =1

def_abs= np.zeros((100,100))
def_x=np.zeros((100,100))
def_y=np.zeros((100,100))


# three oposing forces: ## make better image
#F_pos_list=[(54.5,43.5),(33.5,59.5),(57.5,56.5)]
#F_list=[[-10000,10000],[10000,-10000],[-10000,-10000]] ##? unit ## traktion from cells would be e.g  at 60000



# two oposing foreces
F_pos_list=[(60.5,20.5),(25.5,55.5)]
F_list=[[-10000,10000],[10000,-10000]] ##? unit ## traktion from cells would be e.g  at 60000


#F_pos_list=[(20.5,60.5),(40.5,30.5)]
#F_list=[[10000,-10000],[-10000,+10000]] ##? unit ## traktion from cells would be e.g  at 60000
#F_pos_list=[(65.5,60.5),(35.5,30.5)]
#F_list=[[10000,10000],[-10000,-10000]]
# three forces at each side
F_pos_list=[(55.5,50.5),(50.5,50.5),(55.5,45.5),(35.5,30.5),(30.5,30.5),(35.5,25.5)]
F_list=[[10000,10000],[10000,10000],[10000,10000],[-10000,-10000],[-10000,-10000],[-10000,-10000]] ##? unit ## traktion from cells would be e.g  at 60000


F_list=[np.array(x) for x in F_list]

A=(1+sigma)/(np.pi*young)
## iterating over forces, and summing up deformation
dist_test=np.zeros(np.shape(def_abs))
for F,F_pos in zip(F_list,F_pos_list):
    for i in range(np.shape(def_abs)[0]):
        for j in range(np.shape(def_abs)[1]):
            x=i-F_pos[1]  ## switched because arrays...
            y=j-F_pos[0]
            r=np.sqrt(x**2+y**2)
            K=(A/(r**3))*np.array([[(1-sigma)*(r**2)+sigma*(x**2), sigma*x*y],
                       [sigma*x*y, (1-sigma)*(r**2)+sigma*(y**2)]])

            defo=np.matmul(K,np.transpose(F))
            def_abs[i,j] += np.sqrt(defo[0]**2+defo[1]**2)
            def_x[i, j] +=defo[0]
            def_y[i, j] +=defo[1]


def_x_pos,def_y_pos=get_xy_for_quiver(def_x)# positions for quiver
folder="/home/user/Desktop/"
#p.save(os.path.join(folder,"u_simu.npy"),def_x)
#np.save(os.path.join(folder,"v_simu.npy"),def_y)


tx,ty=ffttc_traction(def_x,def_y,young,pixelsize,pixel_factor=1,sigma=0.49,bf_image=False,filter="gaussian")
#np.save(os.path.join(folder,"tx_simu.npy"),tx)
#np.save(os.path.join(folder,"ty_simu.npy"),ty)




sum_traction=np.sum(np.array(F_list))
sum_traction_abs=np.sum(np.abs(np.array(F_list)))
print(sum_traction_abs)
sum_traction_pred=np.sum(tx+ty)
sum_traction_pred_abs=np.sum(np.abs(tx)+np.abs(ty))
print(sum_traction_pred_abs)


# rms caclualtion:

#gausssian blurr
############## spilit in x and y ##### wrong appart from that.
blur=np.zeros(np.shape(def_abs))
for F_pos in F_pos_list:
    blur_pos=[int(x) for x in [np.floor(F_pos[0]),np.ceil(F_pos[0])+1,np.floor(F_pos[1]),np.ceil(F_pos[1])+1]]
    blur[blur_pos[2]:blur_pos[3],blur_pos[0]:blur_pos[1]] = 10000 / 4


blur=gaussian(blur,sigma=2)

scale=0.0005
for F,F_pos in zip(F_list,F_pos_list):
    plt.arrow(F_pos[0],F_pos[1],F[0]*scale,F[1]*scale,head_width=1,facecolor="red"
              ,edgecolor="red",head_starts_at_zero=False,overhang=2)


blur=blur[0:np.shape(tx)[0],0:np.shape(tx)[1]]
rms=np.sum((blur-tx)**2/(np.size(blur)))
print("rms=",rms)




#mean blur:


ind_filter=np.array([x%2 ==0 for x in np.arange(0,np.size(tx),1)])
ind_filter_arr=np.reshape(ind_filter,np.shape(tx))

ind_filter_2=np.array([x%2 ==0 for x in np.arange(0,np.size(def_x),1)])
ind_filter_arr_2=np.reshape(ind_filter_2,np.shape(def_x))




mask=np.ones(np.shape(tx))>0
pixelsize=1
contractile_force,proj_x,proj_y,center=contractility(tx,ty,pixelsize,mask)
# expected value:
#e=np.sqrt(2*10000**2)*2*(pixelsize*10**-6)





# deformation field !!!! in typical pixel size of bf-image!!!!
plt.figure()
plt.imshow(def_abs,cmap="jet")
plt.title("simulation of deformation field")
cbar=plt.colorbar()
cbar.set_label("displacement in pixel")

for F,F_pos in zip(F_list,F_pos_list):
    plt.arrow(F_pos[0],F_pos[1],F[0]/1000,F[1]/1000,head_width=1,facecolor="red"
              ,edgecolor="red",head_starts_at_zero=False,overhang=2)
scale_ratio=0.2
scale = scale_ratio * np.max(np.shape(def_y)) / np.max(def_abs)  # automatic scaleing in dependance of the image size
Q=plt.quiver(def_x_pos,def_y_pos,
             def_x*scale,def_y*scale,headwidth=3,angles='xy', scale_units='xy',
             scale=1,width=0.0025,color="black")


plt.arrow(85, -10, 3*scale, 0, head_width=0, facecolor="red", edgecolor="green",
          clip_on=False, head_starts_at_zero=False, overhang=2)
plt.text(86,-2,"displacement arrows\nlength 3 Pixel",color="green")





# traction forces
plt.figure()
plt.title("prediction of traction forces")
plt.imshow(np.sqrt(tx**2+ty**2)/1000,cmap="jet")
plt.plot(center[0],center[1],"or",color="red")
cbar=plt.colorbar()
cbar.set_label("Traction in kPa")
x_small,y_small=get_xy_for_quiver(tx)
scale_ratio=0.2
scale = scale_ratio * np.max(np.shape(tx)) / np.max(np.sqrt(tx**2+ty**2)/1000)
plt.quiver(x_small,y_small,(tx/1000)*scale ,(ty/1000)*scale , angles='xy', scale_units='xy', color="red",
           scale=1,headwidth = 4, width=0.002)

plt.arrow(85, -10, 1*scale, 0, head_width=0, facecolor="red", edgecolor="red",
          clip_on=False, head_starts_at_zero=False, overhang=2)
plt.text(86,-2,"tractiont arrows\nlength 1 kPa",color="red")






## illustration of projection
plt.figure()
plt.imshow(np.sqrt((tx / 1000) ** 2 + (ty / 1000) ** 2),alpha=0.3)
plt.plot(center[0],center[1],color="red")
cbar = plt.colorbar()
cbar.set_label("traktion forces in kPa")
x1, y1 = get_xy_for_quiver(tx)
ratio=0.05 # ratio of length of biggest arrow to max axis lenght
scale=ratio *np.max(np.shape( proj_x))/np.max(np.sqrt( proj_x**2+ proj_y**2))# automatic sacleing in dependace of the image size

plt.quiver(x1, y1, proj_x * scale, proj_y  *scale, scale=1,angles="xy", scale_units='xy', width=0.002)
plt.arrow(75, -14, 1 * scale, 0, head_width=0, facecolor="black", edgecolor="black",
          clip_on=False, head_starts_at_zero=False, overhang=2)
plt.text(75, -10, "1 kPa Force", color="black")
plt.plot(center[0],center[1],"or",color="red")
plt.imshow(mask,alpha=0.1)
plt.text(10, 70, "contractile force = "+str(np.round(contractile_force*10**9,2))+"nm", color="black")
#ecpect valu of the contractile force:
e=np.sqrt(2*(10000**2))*(((1*(10**-6))**2))*2
plt.text(10, 80, "contractile force = "+str(np.round(e*10**9,2))+"nm", color="black")




#some results:
# median filter doesnt work (many "background/noise" signals)

# intersting result: contratile energy depends on cell diameter because of projection
#--> could also project on line