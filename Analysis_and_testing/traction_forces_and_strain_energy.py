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

folder="/media/user/GINA1-BK/traktion_force_microscopy/analysis_with_custom_program/"

bf_image_file="/media/user/GINA1-BK/traktion_force_microscopy/new_data_from_magdalena26_02_19/WTshift/09bf_before_shift.tif"
file1="/media/user/GINA1-BK/traktion_force_microscopy/new_data_from_magdalena26_02_19/WTshift/09after_shift.tif"

u=np.load(os.path.join(folder,"u.npy"))   ## dx
v=np.load(os.path.join(folder,"v.npy"))    ##dy
x=np.load(os.path.join(folder,"x.npy"))
y=np.load(os.path.join(folder,"y.npy"))
mask=np.load(os.path.join(folder,"mask.npy"))

bf_image  = np.array(openpiv.tools.imread( bf_image_file ),dtype="int32")  ## those two need the smae shape(?!)
frame_a  = np.array(openpiv.tools.imread( file1 ),dtype="int32")

#shift correction by dx-dymean (should generally be zero=



u_shift=u-np.mean(u)  # also try without dis
v_shift=v-np.mean(v)



## zoom /interpolate to bright filed image size// confirm if same pixel ratio..
#plt.quiver( x, y, u, v,scale=100)
#fttc, fourrier transform traction force cytometire method input u, v (in same size as b image??)

## simple fourrrier trasnform  2D to get to fourrier space




## bens algortithm:

sigma=0.49 #poission ration
young=25536## youngsmodules
pixelsize=6.25/40  ## in µm , here with 40x magnification, depends on the microscope
A=(1+sigma)/(np.pi*young)## constant necessary for K



ax1_length=np.shape(u_shift)[0]   # u and v must have same dimensions
ax2_length=np.shape(u_shift)[1]
#1.0) remove 1 pixel if not  axis length devidable by two) ( is this really necessary)
if ax1_length%2!=0:
    v_cut=v_shift[:-1,:]
    u_cut=u_shift[:-1,:]
else:
    v_cut=v_shift[:,:-1]
    u_cut=u_shift[:,:-1]

#1) zero padding
# add zeros ad array edges to produce square matrix with full contetn
ax1_length=np.shape(u_cut)[0]
ax2_length=np.shape(u_cut)[1]
max_ind=int(np.max((ax1_length,ax2_length)))



u_expand=np.zeros((max_ind,max_ind))
v_expand=np.zeros((max_ind,max_ind))
u_expand[:ax1_length,:ax2_length]=u_cut
v_expand[:ax1_length,:ax2_length]=v_cut
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
alpha= np.arctan2(ky,kx)



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

u_ft_2=np.fft.fft2(u_cut*pixelsize*2*np.pi,(82,82))
v_ft_2=np.fft.fft2(v_cut*pixelsize*2*np.pi,(82,82))


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
tx_filter=uniform_filter(tx_cut,size=3)
ty_filter=uniform_filter(ty_cut,size=3)





#5.x) zooming out for represetnation purposes exactly to bf_image size
zoom_factors=[np.shape(bf_image)[0]/np.shape(tx_cut)[0],np.shape(bf_image)[1]/np.shape(tx_cut)[1]]



tx_zoom_out=zoom(tx_filter,zoom_factors)
ty_zoom_out=zoom(ty_filter,zoom_factors)



# raw out put for comparrison with matlab programm



def get_xy_for_quiver(u):
    xs=np.zeros(np.shape(u))
    for i in range(np.shape(u)[0]):
        xs[i,:]=np.arange(0,np.shape(u)[1],1)
    ys = np.zeros(np.shape(u))
    for j in range(np.shape(u)[1]):  ## is inverted in other skript...
        ys[:, j] = np.arange(0,np.shape(u)[0], 1)
    return xs, ys

plt.close("all")
#plt.figure()
#plt.imshow(np.sqrt(tx_cut**2+ty_cut**2)/1000)
#x_small,y_small=get_xy_for_quiver(tx_cut)
#plt.quiver(x_small,y_small,tx_cut/1000,ty_cut/1000,scale=600)
plt.figure()
plt.imshow(np.sqrt(tx_filter**2+ty_filter**2)/1000)
plt.colorbar()

x_small,y_small=get_xy_for_quiver(tx_cut)
plt.quiver(x_small,y_small,tx_filter/1000,ty_filter/1000,scale=600)



#plt.close("all")

#plt.figure()
#plt.imshow(bf_image,cmap="gray")
fig, ax = plt.subplots()
plt.quiver( x, y, u, v,scale=100)


import matplotlib.font_manager as fm
fontprops = fm.FontProperties(size=15)
scalebar = AnchoredSizeBar(ax.transData,
                           400, '10 pixel displacement', 'upper right',  ##### please fix , got this from looking at np.max(np.sqrt(u**2+v**2))
                           pad=0.1,
                           color='black',
                           frameon=True,
                           size_vertical=30,
                           fontproperties=fontprops)

ax.add_artist(scalebar)

plt.figure()
plt.imshow(np.sqrt(tx_zoom_out**2+ty_zoom_out**2)/1000)
plt.colorbar()
plt.quiver(x[:,:-1],np.flip(y,axis=0)[:,:-1],tx_filter/1000,ty_filter/1000,scale=300)
'''
#plt.figure()
#plt.imshow(np.sqrt(bf_image),alpha=0.5)
#plt.imshow(np.sqrt(tx_zoom_out**2+ty_zoom_out**2)/1000,alpha=0.5)
#plt.colorbar()



'''



# plotting with quiver

#x_quiver=np.array([list(range(0,int(max_ind),1)),]*int(max_ind))
#y_quiver=np.array([list(range(-int(max_ind),0,1)),]*int(max_ind))



### todo:
# good representation
#checkt units and scaling
# check other zooming methods
# also recosider cutting and all that stuff


















'''
not inversed K

### use better method from matlab script
### also use negative values for kx and ky ?????????????
K_matrix=np.zeros((np.shape(u_ft)[0],np.shape(u_ft)[1],2,2))  #### convoultion core fore.... what??
#eigene berechnung von K
## problem with devision by zero
for i in range(np.shape(u_ft)[0]):
    for j in range(np.shape(u_ft)[1]):
        k_2=i ** 2 + j ** 2
        k = np.sqrt(k_2)
        k_3 =  k_2**3
        if j==0 and i==0:
            K_matrix[i, j, :, :]=0
        else:
            K_matrix[i,j,:,:]=(A*2*np.pi/k_3)*np.array([[(1-sigma)*k_2+sigma*(j**2), sigma*i*j],
                                                [sigma*i*j,(1-sigma)*k_2+sigma*(j**2)]])           ## matrix of matirces


#is there a problem with K not squaerd??
'''