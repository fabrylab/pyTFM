import openpiv.tools
import openpiv.process
import openpiv.scaling
import openpiv.validation
import openpiv.filters
import matplotlib.pyplot as plt
import numpy as np
import copy
import os
from scipy.ndimage.filters import uniform_filter
from scipy.signal import convolve2d

def std_convoluted(image, N):
    im = np.array(image, dtype=float)
    im2 = im**2
    ones = np.ones(im.shape)

    kernel = np.ones((2*N+1, 2*N+1))
    s = convolve2d(im, kernel, mode="same")
    s2 = convolve2d(im2, kernel, mode="same")
    ns = convolve2d(ones, kernel, mode="same")

    return np.sqrt((s2 - s**2 / ns) / ns)

#folder=r"/media/user/GINA1-BK/traktion_force_microscopy/analysis_with_custom_program"
plt.close("all")
file1="/home/user/Desktop/IFBmshift/03after_shift.tif"
file2="/home/user/Desktop/IFBmshift/03before_shift.tif"
frame_a  = np.array(openpiv.tools.imread( file1 ),dtype="int32")
frame_b  = np.array(openpiv.tools.imread( file2 ),dtype="int32")

plt.close("all")


def get_xy_for_quiver(u):
    xs=np.zeros(np.shape(u))
    for i in range(np.shape(u)[0]):
        xs[i,:]=np.arange(0,np.shape(u)[1],1)
    ys = np.zeros(np.shape(u))
    for j in range(np.shape(u)[1]):  ## is inverted in other skript...
        ys[:, j] = np.arange(np.shape(u)[0],0, -1)
    return xs, ys



window_size=64
overlapp=32
pixel_factor=window_size-overlapp # factor to adijuyt pixelsize




u, v, sig2noise = openpiv.process.extended_search_area_piv( frame_a, frame_b, window_size=window_size, overlap=overlapp,
                                                            dt=1,subpixel_method="gaussian", search_area_size=64, sig2noise_method='peak2peak' )
x, y = openpiv.process.get_coordinates( image_size=frame_a.shape, window_size=window_size, overlap=overlapp)

plt.figure(1)
xs, ys = get_xy_for_quiver(u)



plt.quiver(ys , xs , u, v,scale=100,angles='xy' )

plt.figure(2)
plt.imshow(sig2noise)

u, v, mask = openpiv.validation.sig2noise_val( u, v, sig2noise, threshold = 1.05 )





plt.figure(3)
plt.quiver( x, y, u, v ,scale=100)

#

#filter by local mean in u and v ## add brder conditions
N=7
u_std=np.zeros(np.shape(u))
u_mean=np.zeros(np.shape(u))
for i in range(np.shape(u)[0]):
    for j in range(np.shape(u)[1]):
        u_std[i,j]= np.nanstd(u[int((i-N)*((i-N)>0)):int(i+N),int((j-N)*((j-2)>0)):int(j+N)])
        u_mean[i, j] = np.nanmean(u[int((i-N)*((i-N)>0)):int(i+N),int((j-2)*((j-N)>0)):int(j+N)])



v_std=np.zeros(np.shape(v))
v_mean=np.zeros(np.shape(v))
for i in range(np.shape(u)[0]):
    for j in range(np.shape(v)[1]):
        v_std[i,j]=np.nanstd(v[int((i-N))*int(((i-N)>0)):int(i+N),int((j-N)*((j-N)>0)):int(j+N)])
        v_mean[i,j] = np.nanmean(v[int((i-N)*((i-N)>0)):int(i+N),int((j-N)*((j-N)>0)):int(j+N)])




#u_mean=uniform_filter(u,size=5) ## local mean
#v_mean=uniform_filter(v,size=5) #
#u_std=std_convoluted(u, 5) ## local standrat deviation
#v_std=std_convoluted(v, 5)
std_factor=8
mask_std_u=(u<(u_mean-std_factor*u_std)) + (u>(u_mean+std_factor*u_std))
mask_std_v=(v<(v_mean-std_factor*v_std)) + (v>(v_mean+std_factor*v_std))
mask_std=mask_std_u+mask_std_v
  ### naja
plt.figure()
plt.imshow(mask_std)


### ## validation by replacing extreme outliers by std
# same as in u, v, mask=openpiv.validation.global_std(u, v,std_threshold=2)
def_abs=np.sqrt(u**2+v**2)
m=np.nanmean(def_abs)
std=np.nanstd(def_abs)
plt.figure(4)
plt.hist(def_abs[~np.isnan(def_abs)].flatten(),bins=50)
plt.axvline(m)
plt.axvline(std*15+m)
threshold=std*15+m
mask_out=def_abs>threshold
u[mask_out]=np.nan
v[mask_out]=np.nan


plt.figure(5)
plt.imshow(mask_out)

## same as


plt.figure(6)
plt.quiver( x, y, u, v ,scale=100)



u, v = openpiv.filters.replace_outliers( u, v, method='localmean', max_iter=10, kernel_size=2 )
plt.figure(7)
plt.quiver( x, y, u, v,scale=100)





'''
np.save(os.path.join(folder,"u"),u)
np.save(os.path.join(folder,"v"),-v)   ## invert beacuase of image axis??
np.save(os.path.join(folder,"x"),x)
np.save(os.path.join(folder,"y"),y)
np.save(os.path.join(folder,"mask"),mask)
'''


