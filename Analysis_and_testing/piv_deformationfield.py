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
from TFM_functions import *


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
#/media/user/GINA1-BK/traktion_force_microscopy/IFBmshift/01before_shift.tif
#/media/user/GINA1-BK/traktion_force_microscopy/IFBmshift/01after_shift.tif

#/media/user/GINA1-BK/traktion_force_microscopy/new_data_from_magdalena26_02_19/WTshift/01before_shift.tif
#/media/user/GINA1-BK/traktion_force_microscopy/new_data_from_magdalena26_02_19/WTshift/01after_shift.tif

file1=r"//media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/11bf_before_shift.tif"
file2=r"/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/11after_shift.tif"
frame_a  = np.array(openpiv.tools.imread( file1 ),dtype="int32")
frame_b  = np.array(openpiv.tools.imread( file2 ),dtype="int32")
print(np.shape(frame_a))




def get_xy_for_quiver(u):
    xs=np.zeros(np.shape(u))
    for i in range(np.shape(u)[0]):
        xs[i,:]=np.arange(0,np.shape(u)[1],1)
    ys = np.zeros(np.shape(u))
    for j in range(np.shape(u)[1]):  ## is inverted in other skript...
        ys[:, j] = np.arange(np.shape(u)[0],0, -1)
    return xs, ys




window_size=96   # about 4 times bead as seen in image, 10 actual bead diameter
overlapp=64    # window_size/2   one could generally try higer overlapp
pixel_factor=window_size-overlapp # factor to adjust pixelsize




u, v, sig2noise = openpiv.process.extended_search_area_piv( frame_a, frame_b, window_size=window_size, overlap=overlapp,
                                                            dt=1,subpixel_method="gaussian", search_area_size=window_size, sig2noise_method='peak2peak' )
x, y = openpiv.process.get_coordinates( image_size=frame_a.shape, window_size=window_size, overlap=overlapp)

plt.figure()
xs, ys = get_xy_for_quiver(u)



plt.quiver(xs , ys , u, v,scale=100,angles='xy' )

plt.figure()
plt.imshow(sig2noise)

u, v, mask = openpiv.validation.sig2noise_val( u, v, sig2noise, threshold = 1.05 )



u, v = openpiv.filters.replace_outliers( u, v, method='localmean', max_iter=10, kernel_size=2 )
plt.figure()
plt.quiver( x, y, u, v,scale=100)




'''
np.save(os.path.join(folder,"u"),u)
np.save(os.path.join(folder,"v"),-v)   ## invert beacuase of image axis??
np.save(os.path.join(folder,"x"),x)
np.save(os.path.join(folder,"y"),y)
np.save(os.path.join(folder,"mask"),mask)
'''


