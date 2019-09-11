import openpiv.tools
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import numpy as np
import copy
import os
from scipy.ndimage import zoom
from skimage.filters import rank
from skimage.morphology import cube
from scipy.ndimage.filters import uniform_filter,gaussian_filter

folder="/home/user/Desktop/"


u=np.load(os.path.join(folder,"u_simu.npy"))   ## dx
v=np.load(os.path.join(folder,"v_simu.npy"))    ##dy

#bf_image  = np.array(openpiv.tools.imread( bf_image_file ),dtype="int32")  ## those two need the smae shape(?!)
#frame_a  = np.array(openpiv.tools.imread( file1 ),dtype="int32")

#shift correction by dx-dymean (should generally be zero=






## zoom /interpolate to bright filed image size// confirm if same pixel ratio..
#plt.quiver( x, y, u, v,scale=100)
#fttc, fourrier transform traction force cytometire method input u, v (in same size as b image??)

## simple fourrrier trasnform  2D to get to fourrier space




## bens algortithm:

sigma=0.49 #poission ration
young=25536## youngsmodules
pixelsize=1#6.25/40  ## in µm , here with 40x magnification, depends on the microscope
A=(1+sigma)/(np.pi*young)## constant necessary for K

window_size=64
overlapp=32
pixel_factor=window_size-overlapp
##pixelsize=pixelsize*pixel_factor    ### u is given in terms of pixels of original image


# 0) substracting mean(better name for this step)
u_shift = (u - np.mean(u))  # also try without dis
v_shift = (v - np.mean(v))
## bens algortithm:

# 1)Zeropadding to get squared array with even index number
'''
ax1_length = np.shape(u_shift)[0]  # u and v must have same dimensions
ax2_length = np.shape(u_shift)[1]
max_ind = int(np.max((ax1_length, ax2_length)))
if max_ind % 2 != 0:
    max_ind += 1
u_expand = np.zeros((max_ind, max_ind))
v_expand = np.zeros((max_ind, max_ind))
u_expand[:ax1_length, :ax2_length] = u_shift
v_expand[:ax1_length, :ax2_length] = v_shift
'''


ax1_length = np.shape(u_shift)[0]  # u and v must have same dimensions
ax2_length = np.shape(u_shift)[1]
max_ind = int(np.max((ax1_length, ax2_length)))
if max_ind % 2 != 0:
    max_ind += 1
else:
    max_ind += 2
u_expand = np.zeros((max_ind, max_ind))
v_expand = np.zeros((max_ind, max_ind))
u_expand[max_ind-ax1_length:max_ind, max_ind-ax2_length:max_ind] = u_shift
v_expand[max_ind-ax1_length:max_ind, max_ind-ax2_length:max_ind] = v_shift







## multiply pixel size in here...

 #2) produceing wave vectors   ## why this form
 #form 1:max_ind/2 then -(max_ind/2:1)
kx1=np.array([list(range(0,int(max_ind/2),1)),]*int(max_ind))
kx2=np.array([list(range(-int(max_ind/2),0,1)),]*int(max_ind))
kx=np.append(kx1,kx2,axis=1)
ky=np.transpose(kx)
k = np.sqrt(kx**2 + ky**2)/(pixelsize*max_ind)    # matrix with "relative" distances??#





#np.save("/home/user/Desktop/k_test.npy",k)


#2.) caclulating angle between k and kx with atan 2 function (what is this exactely??)  just if statemments to get
# angel from x1 to x2 in fixed direction... (better explanation)
alpha= np.arctan2(ky,kx)   # angle of ky,kx vector to y axis



alpha[0,0]=np.pi/2 ## why do i need this?-> arctan2(n,0) n-->0 is pi/2 for n positive and -pi/2 for n negative arctan(0,0) has been defined on 0
#np.save("/home/user/Desktop/alpha_test.npy",alpha)
#3) calkulation of K--> Tensor to calculate displacements from Tractions. We calculate inverse of K
#(check if correct inversion by your self)
#K⁻¹=[[kix kid],
#     [kid,kiy]]  ,,, so is "diagonal, kid appears two times

kix = ((k*young) / (2*(1-sigma**2))) *((1 - sigma + sigma * np.cos(alpha)**2))
kiy = ((k*young) / (2*(1-sigma**2))) *((1 - sigma + sigma * np.sin(alpha)**2))
kid = ((k*young) / (2*(1-sigma**2))) *(sigma * np.sin(alpha) * np.cos(alpha))
#np.save("/home/user/Desktop/kid_seg2_test.npy",(sigma * np.sin(alpha) * np.cos(alpha)))
## adding zeros in kid diagonals(?? why do i need this)
kid[:,int(max_ind/2)]=np.zeros(max_ind)
kid[int(max_ind/2),:]=np.zeros(max_ind)


#4) calculate fourrier transform of dsiplacements
u_ft=np.fft.fft2(u_expand*pixelsize*2*np.pi)
v_ft=np.fft.fft2(v_expand*pixelsize*2*np.pi)

#u_ft_2=np.fft.fft2(u_cut*pixelsize*2*np.pi,(82,82))
#v_ft_2=np.fft.fft2(v_cut*pixelsize*2*np.pi,(82,82))


#4.1) calculate tractions in fourrier space T=K⁻¹*U, U=[u,v]  here with individual matrix elelemnts..
tx_ft = kix * u_ft  +  kid * v_ft;
ty_ft = kid * u_ft  +  kiy * v_ft;
##plt.imshow(tx_ft.real)
#4.2) go back to real space
tx=np.fft.ifft2(tx_ft).real
ty=np.fft.ifft2(ty_ft).real


#5.2) cut like in script from ben

## not sure how robust
tx_cut=tx[0:ax1_length,0:ax2_length]
ty_cut=ty[0:ax1_length,0:ax2_length]
#5.3) using filter (mean filter on 3 by 3 form)
tx_filter = gaussian_filter(tx_cut, sigma=1.5)
ty_filter = gaussian_filter(ty_cut, sigma=1.5)


np.save(os.path.join(folder,"tx_filter.npy"),tx_filter)  ## note: x is one shape shorter // fix th
np.save(os.path.join(folder,"ty_filter.npy"),ty_filter) ##



#5.x) zooming out for represetnation purposes exactly to bf_image size
#zoom_factors=[np.shape(bf_image)[0]/np.shape(tx_cut)[0],np.shape(bf_image)[1]/np.shape(tx_cut)[1]]
#tx_zoom_out=zoom(tx_filter,zoom_factors)
#ty_zoom_out=zoom(ty_filter,zoom_factors)





def get_xy_for_quiver(u):
    xs=np.zeros(np.shape(u))
    for i in range(np.shape(u)[0]):
        xs[i,:]=np.arange(0,np.shape(u)[1],1)
    ys = np.zeros(np.shape(u))
    for j in range(np.shape(u)[1]):  ## is inverted in other skript...
        ys[:, j] = np.arange(0,np.shape(u)[0], 1)
    return xs, ys


fig = plt.figure()
plt.imshow(np.sqrt((tx_filter / 1000) ** 2 + (ty_filter / 1000) ** 2), cmap="rainbow")
cbar = plt.colorbar()
cbar.set_label("traktion forces in kPa")
x1, y1 = get_xy_for_quiver(tx_filter)
scale_ratio = 0.2
scale = scale_ratio * np.max(np.shape(tx_filter)) / np.max(
    np.sqrt((tx_filter / 1000) ** 2 + (ty_filter / 1000) ** 2))  # automatic sacleing in dependace of the image size #
plt.quiver(x1, y1, (tx_filter / 1000) * scale, (ty_filter / 1000) * scale, scale=1, scale_units='xy', angles="xy")
plt.arrow(75, -14, 1 * scale, 0, head_width=0, facecolor="black", edgecolor="black",
          clip_on=False, head_starts_at_zero=False, overhang=2)
plt.text(75, -10, "1 kPa Force", color="black")




