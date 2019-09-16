import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '/home/user/Desktop/Andreas-Python/tracktion_force_microscopy/')
from TFM_functions import *


#tx=np.array([[10000,-10000],[10000,-10000]])/4
#ty=np.array([[10000,10000],[-10000,-10000]])/4



#tx=np.array([[10000,0],[0,-10000]])/4
#ty=np.array([[10000,0],[0,-10000]])/4




tx=np.array([[0,0,0,0,0],
             [0, 0, 0, 0,0],
             [10000, 0, 0, 0,-10000],
                [0, 0, 0, 0,0],
             [0,0,0,0,0]])/4
ty=np.array([[0,0,0,0,0],
             [0, 0, 0, 0,0],
             [0, 0, 0, 0,0],
            [0, 0, 0, 0,0],
             [0,0,0,0,0]])/4

contractility_excepted=((np.sqrt((10000**2))*2)/4)  *((10**-6)**2) ### atlast this is correct
mask=np.ones(np.shape(ty))>0
pixelsize=1#*(10**6) ## maybe some rounding errors??? or because  plot starts on negative numbers???// check out later
contractile_force, proj_x, proj_y,center=contractility(tx,ty,pixelsize,mask)
t_abs=np.sqrt(tx**2+ty**2)

pixx = np.arange(np.shape(tx)[0])
pixy = np.arange(np.shape(ty)[1])
x, y = np.meshgrid(pixy, pixx)
plt.close("all")
plt.figure()
plt.imshow(t_abs)
plt.quiver(x, y, tx, ty,scale=10000, scale_units='xy',angles='xy')
plt.plot(center[0],center[1],"or",color="red")
plt.text(0,-0.5,"contractility_excepted = "+str(np.round(contractility_excepted,11)))
plt.text(0,-0.75,"contractility_calculated = "+str(np.round(contractile_force,11)))

plt.figure()
plt.imshow(t_abs)
plt.quiver(x, y, proj_x , proj_y  , angles="xy", scale=(10**-4)**2, scale_units='xy')

plt.plot(center[0],center[1],"or",color="red")


plt.text(0,-0.5,"contractility_excepted = "+str(np.round(contractility_excepted,11)))
plt.text(0,-0.75,"contractility_calculated = "+str(np.round(contractile_force,11)))