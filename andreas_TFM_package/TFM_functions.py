# basic functions for traction force microscopy
import openpiv.tools
import numpy as np
import copy
from scipy.ndimage.filters import uniform_filter,median_filter,gaussian_filter
import openpiv.process
import openpiv.scaling
import openpiv.validation
import openpiv.filters
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.ndimage import zoom
from scipy.ndimage.filters import uniform_filter






def ffttc_traction(u,v,pixelsize1,pixelsize2,young,sigma=0.49,bf_image=False,filter="mean"):
    '''
    fourier transform based calculation of the traction force. U and v must be given  as deformations in pixel. Size of
    these pixels must be the pixelsize (size of a pixel in the deformation field u or v). Note that thePiv deformation
    returns deformation in pixel of the size of pixels in the images of beads before and after.
    If bf_image is provided this script will return a traction field that is zoomed to the size of the brightfield image,
    by interpolation. It is not recommended to use this for any calculations.
    The function can use diffretn filters. Recommended filter is gaussian. Mean filter should yield similar results.

    :param u:deformation field in x direction in pixel of the deformation image
    :param v:deformation field in y direction in pixel of the deformation image
    :param young: youngs modulus in Pa
    :param pixelsize1: pixelsize of the original image, needed because u and v is given as displacment of these pixels
    :param pixelsize2: pixelsize of the deformation image
    :param sigma: posiion ratio of the gel
    :param bf_image: give the brightfield image as an array before cells where removed
    :param filter: str, values: "mean","gaussian","median". Diffrent smoothing methods for the traction field
    :return: tx_filter,ty_filter: traction forces in x and y direction in Pa
    '''

    #0) substracting mean(better name for this step)
    u_shift=(u-np.mean(u))  # shifting to zero mean  (translating to pixelsize of u-image is done later)
    v_shift=(v-np.mean(v))
    ## bens algortithm:

    # 1)Zeropadding to get sqauerd array with even index number
    ax1_length=np.shape(u_shift)[0]   # u and v must have same dimensions
    ax2_length=np.shape(u_shift)[1]
    max_ind=int(np.max((ax1_length,ax2_length)))
    if max_ind%2!=0:
        max_ind+=1

    u_expand=np.zeros((max_ind,max_ind))
    v_expand=np.zeros((max_ind,max_ind))
    u_expand[:ax1_length,:ax2_length]=u_shift
    v_expand[:ax1_length,:ax2_length]=v_shift

    #2) produceing wave vectors   ## why this form
    #form 1:max_ind/2 then -(max_ind/2:1)
    kx1=np.array([list(range(0,int(max_ind/2),1)),]*int(max_ind))
    kx2=np.array([list(range(-int(max_ind/2),0,1)),]*int(max_ind))
    kx=np.append(kx1,kx2,axis=1)*2*np.pi  # fourier transform in this case is defined as
    # F(kx)=1/2pi integral(exp(i*kx*x)dk therefore kx must be expressed as a spatial frequency in distance*2*pi

    ky=np.transpose(kx)
    k = np.sqrt(kx**2 + ky**2)/(pixelsize2*max_ind)
    #np.save("/home/user/Desktop/k_test.npy",k)

    #2.1) calculating angle between k and kx with atan 2 function (what is this exactely??)  just if statemments to get
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
    #u_ft=np.fft.fft2(u_expand*pixelsize1*2*np.pi)
    #v_ft=np.fft.fft2(v_expand*pixelsize1*2*np.pi)   #
    u_ft=np.fft.fft2(u_expand*pixelsize1)
    v_ft=np.fft.fft2(v_expand*pixelsize1)


    #4.1) calculate tractions in fourrier space T=K⁻¹*U, U=[u,v]  here with individual matrix elelemnts..
    tx_ft = kix * u_ft  +  kid * v_ft
    ty_ft = kid * u_ft  +  kiy * v_ft
    ##plt.imshow(tx_ft.real)
    #4.2) go back to real space
    tx=np.fft.ifft2(tx_ft).real
    ty=np.fft.ifft2(ty_ft).real

    #5.2) cut back to oringinal shape
    tx_cut=tx[0:ax1_length,0:ax2_length]
    ty_cut=ty[0:ax1_length,0:ax2_length]

    #5.3) using filter
    if filter=="mean":
        #tx_filter=uniform_filter(tx_cut,size=5)     # this would be non responsive to resolution of the image ### there should be a better way
        #ty_filter=uniform_filter(ty_cut,size=5)
        tx_filter = uniform_filter(tx_cut, size=int(int(np.max((ax1_length,ax2_length)))/16))
        ty_filter = uniform_filter(ty_cut, size=int(int(np.max((ax1_length,ax2_length)))/16))
    if filter == "gaussian":
        tx_filter = gaussian_filter(tx_cut, sigma=int(np.max((ax1_length,ax2_length)))/50)
        ty_filter = gaussian_filter(ty_cut, sigma=int(np.max((ax1_length,ax2_length)))/50)
    if filter == "median":
        tx_filter = median_filter(tx_cut, size=int(int(np.max((ax1_length,ax2_length)))/16))
        ty_filter = median_filter(ty_cut, size=int(int(np.max((ax1_length,ax2_length)))/16))
    if not filter:
        tx_filter = tx_cut
        ty_filter = ty_cut
    #show_quiver(tx_filter,ty_filter)

    #5.x) zooming out for represetnation purposes exactly to bf_image size
    if bf_image:
        zoom_factors=[np.shape(bf_image)[0]/np.shape(tx_cut)[0],np.shape(bf_image)[1]/np.shape(tx_cut)[1]]
        tx_zoom_out=zoom(tx_filter,zoom_factors)
        ty_zoom_out=zoom(ty_filter,zoom_factors)
        return (tx_filter, ty_filter,tx_zoom_out,ty_zoom_out)



    return (tx_filter,ty_filter)


def ffttc_traction_pure_shear(u, v, pixelsize1, pixelsize2, h, young, sigma = 0.49, filter = "mean") -> object:
    '''

     limiting case for h*k==0
    Xavier Trepat, Physical forces during collective cell migration, 2009


    :param u:deformation field in x direction in pixel of the deformation image
    :param v:deformation field in y direction in pixel of the deformation image
    :param young: youngs modulus in Pa
    :param pixelsize1: pixelsize of the original image, needed because u and v is given as displacment of these pixels
    :param pixelsize2: pixelsize of the deformation image
    :param h hight of the membrane the cells lie on, in µm
    :param sigma: poission ratio of the gel
    :param bf_image: give the brightfield image as an array before cells where removed
    :param filter: str, values: "mean","gaussian","median". Diffrent smoothing methods for the traction field.
    :return: tx_filter,ty_filter: traction forces in x and y direction in Pa
    '''



     # 0) substracting mean(better name for this step)
    u_shift = (u - np.mean(u))*pixelsize1  # also try without dis
    v_shift = (v - np.mean(v))*pixelsize1



    ## bens algortithm:

    # 1)Zero padding to get sqauerd array with even index number
    ax1_length = np.shape(u_shift)[0]  # u and v must have same dimensions
    ax2_length = np.shape(u_shift)[1]
    max_ind = int(np.max((ax1_length, ax2_length)))
    if max_ind % 2 != 0:
        max_ind += 1
    u_expand = np.zeros((max_ind, max_ind))
    v_expand = np.zeros((max_ind, max_ind))
    u_expand[max_ind - ax1_length:max_ind, max_ind - ax2_length:max_ind] = u_shift
    v_expand[max_ind - ax1_length:max_ind, max_ind - ax2_length:max_ind] = v_shift   # note: seems to be snummerically slightly diffrent then in normal fnction ??? (dont know why

    # 2) producing wave vectors (FT-space) "
    # form 1:max_ind/2 then -(max_ind/2:1)
    kx1 = np.array([list(range(0, int(max_ind / 2), 1)), ] * int(max_ind))
    kx2 = np.array([list(range(-int(max_ind / 2), 0, 1)), ] * int(max_ind))

    kx = np.append(kx1, kx2, axis=1) / (
                pixelsize2* max_ind)*2*np.pi  # spatial frequencies: 1/wavelength,in 1/µm in fractions of total length

    ky = np.transpose(kx)
    k = np.sqrt(kx ** 2 + ky ** 2)  # matrix with "relative" distances??#

    u_ft = np.fft.fft2(u_expand)
    v_ft = np.fft.fft2(v_expand)



    # 4.1) calculate tractions in fourrier space
    mu=young/(2*(1+sigma))
    tx_ft = mu*u_ft/h
    tx_ft[0, 0] = 0  ## zero frequency would represent force every where(constant), so this is legit??
    ty_ft = mu*v_ft/h
    ty_ft[0, 0] = 0


    # 4.2) go back to real space
    tx = np.fft.ifft2(tx_ft).real  ## maybe devide by 2pi here???
    ty = np.fft.ifft2(ty_ft).real

    # 5.2) cut like in script from ben

    tx_cut = tx[max_ind - ax1_length:max_ind, max_ind - ax2_length:max_ind]
    ty_cut = ty[max_ind - ax1_length:max_ind, max_ind - ax2_length:max_ind]

    #5.3) using filter
    if filter=="mean":
        #tx_filter=uniform_filter(tx_cut,size=5)     # this would be non responsive to resolution of the image ### there should be a better way
        #ty_filter=uniform_filter(ty_cut,size=5)
        tx_filter = uniform_filter(tx_cut, size=int(int(np.max((ax1_length,ax2_length)))/16))
        ty_filter = uniform_filter(ty_cut, size=int(int(np.max((ax1_length,ax2_length)))/16))
    if filter == "gaussian":
        tx_filter = gaussian_filter(tx_cut, sigma=int(np.max((ax1_length,ax2_length)))/50)
        ty_filter = gaussian_filter(ty_cut, sigma=int(np.max((ax1_length,ax2_length)))/50)
    if filter == "median":
        tx_filter = median_filter(tx_cut, size=int(int(np.max((ax1_length,ax2_length)))/16))
        ty_filter = median_filter(ty_cut, size=int(int(np.max((ax1_length,ax2_length)))/16))
    if not filter:
        tx_filter = tx_cut
        ty_filter = ty_cut
    #show_quiver(tx_filter,ty_filter)
    return (tx_filter, ty_filter)


def ffttc_traction_finite_thickness(u, v, pixelsize1, pixelsize2, h, young, sigma = 0.49, filter = "mean") -> object:
    '''
    FTTC with correction for finite substrate thikness according to
    Xavier Trepat, Physical forces during collective cell migration, 2009


    :param u:deformation field in x direction in pixel of the deformation image
    :param v:deformation field in y direction in pixel of the deformation image
    :param young: youngs modulus in Pa
    :param pixelsize1: pixelsize of the original image, needed because u and v is given as displacment of these pixels
    :param pixelsize2: pixelsize of the deformation image
    :param h hight of the membrane the cells lie on, in µm
    :param sigma: poission ratio of the gel
    :param bf_image: give the brightfield image as an array before cells where removed
    :param filter: str, values: "mean","gaussian","median". Diffrent smoothing methods for the traction field.
    :return: tx_filter,ty_filter: traction forces in x and y direction in Pa
    '''

     # 0) substracting mean(better name for this step)
    u_shift = (u - np.mean(u))*pixelsize1  # also try without dis
    v_shift = (v - np.mean(v))*pixelsize1

    ## bens algortithm:
    # 1)Zero padding to get sqauerd array with even index number
    ax1_length = np.shape(u_shift)[0]  # u and v must have same dimensions
    ax2_length = np.shape(u_shift)[1]
    max_ind = int(np.max((ax1_length, ax2_length)))
    if max_ind % 2 != 0:
        max_ind += 1
    u_expand = np.zeros((max_ind, max_ind))
    v_expand = np.zeros((max_ind, max_ind))
    u_expand[max_ind - ax1_length:max_ind, max_ind - ax2_length:max_ind] = u_shift
    v_expand[max_ind - ax1_length:max_ind, max_ind - ax2_length:max_ind] = v_shift   # note: seems to be snummerically slightly diffrent then in normal fnction ??? (dont know why

    # 2) producing wave vectors (FT-space) "
    # form 1:max_ind/2 then -(max_ind/2:1)
    kx1 = np.array([list(range(0, int(max_ind / 2), 1)), ] * int(max_ind),dtype=np.float64)
    kx2 = np.array([list(range(-int(max_ind / 2), 0, 1)), ] * int(max_ind),dtype=np.float64)

    kx = np.append(kx1, kx2, axis=1) / (
                pixelsize2* max_ind)*2*np.pi  # spatial frequencies: 1/wavelength,in 1/µm in fractions of total length

    ky = np.transpose(kx)
    k = np.sqrt(kx ** 2 + ky ** 2)  # matrix with "relative" distances??#
    # np.save("/home/user/Desktop/k_test.npy",k)

    c = np.cosh(k * h)
    s = np.sinh(k * h)
    #gamma = ((3 - 4 * sigma) * (c ** 2) + (1 - 2 * sigma) ** 2 + (k * h * 2 * np.pi) ** 2) / (
                #(3 - 4 * sigma) * s * c + k * h * 2 * np.pi)  ## inf values here because k goes to zero
    gamma = ((3 - 4 * sigma) * (c ** 2) + (1 - 2 * sigma) ** 2 + (k * h) ** 2) / (
                (3 - 4 * sigma) * s * c + k * h)  ## inf values here because k goes to zero
    # 4) calculate fourrier transform of displacements

    u_ft = np.fft.fft2(u_expand)
    v_ft = np.fft.fft2(v_expand)

    '''
    #4.0*) approximation for large h according to this paper
    factor3=young/(2*(1-sigma**2)*k)
    factor3[0,0]=factor3[0,1]
    tx_ft=factor3*(u_ft*((k**2)*(1-sigma)+sigma*(kx**2)) + v_ft*kx*ky*sigma)
    ty_ft=factor3*(v_ft*((k**2)*(1-sigma)+sigma*(ky**2)) + u_ft*kx*ky*sigma)   ## confirmed by comparisson with not corrected programm
    '''

    # 4.1) calculate tractions in fourrier space
    factor1 = (v_ft * kx - u_ft * ky)
    factor2 = (u_ft * kx + v_ft * ky)
    tx_ft = ((-young * ky * c) / (2 * k * s * (1 + sigma))) * factor1 + (
                (young * kx) / (2 * k * (1 - sigma ** 2))) * gamma * factor2
    tx_ft[0, 0] = 0  ## zero frequency would represent force every where(constant), so this is legit??
    ty_ft = ((young * kx * c) / (2 * k * s * (1 + sigma))) * factor1 + (
                (young * ky) / (2 * k * (1 - sigma ** 2))) * gamma * factor2
    ty_ft[0, 0] = 0

    # 4.2) go back to real space
    tx = np.fft.ifft2(tx_ft.astype(np.complex128)).real  ## maybe devide by 2pi here???
    ty = np.fft.ifft2(ty_ft.astype(np.complex128)).real

    # 5.2) cut like in script from ben
    tx_cut = tx[max_ind - ax1_length:max_ind, max_ind - ax2_length:max_ind]
    ty_cut = ty[max_ind - ax1_length:max_ind, max_ind - ax2_length:max_ind]

    #5.3) using filter
    if filter=="mean":
        #tx_filter=uniform_filter(tx_cut,size=5)     # this would be non responsive to resolution of the image ### there should be a better way
        #ty_filter=uniform_filter(ty_cut,size=5)
        tx_filter = uniform_filter(tx_cut, size=int(int(np.max((ax1_length,ax2_length)))/16))
        ty_filter = uniform_filter(ty_cut, size=int(int(np.max((ax1_length,ax2_length)))/16))
    if filter == "gaussian":
        tx_filter = gaussian_filter(tx_cut, sigma=int(np.max((ax1_length,ax2_length)))/50)
        ty_filter = gaussian_filter(ty_cut, sigma=int(np.max((ax1_length,ax2_length)))/50)
    if filter == "median":
        tx_filter = median_filter(tx_cut, size=int(int(np.max((ax1_length,ax2_length)))/16))
        ty_filter = median_filter(ty_cut, size=int(int(np.max((ax1_length,ax2_length)))/16))
    if not filter:
        tx_filter = tx_cut
        ty_filter = ty_cut
    #show_quiver(tx_filter,ty_filter)
    return (tx_filter, ty_filter)

def ffttc_traction_old(u,v,young,pixelsize,pixel_factor,sigma=0.49,bf_image=False,filter="mean"):

    '''
    this function is depricated avoid usage
    :param u:
    :param v:
    :param young:
    :param pixelsize:
    :param pixel_factor:
    :param sigma:
    :param bf_image:
    :param filter:
    :return:
    '''


    #0) substracting mean(better name for this step)
    u_shift=(u-np.mean(u))/pixel_factor   # shifting to zero mean and trasnlating to pixelsize of u-image
    v_shift=(v-np.mean(v))/pixel_factor
    ## bens algortithm:

    # 1)Zeropadding to get sqauerd array with even index number
    ax1_length=np.shape(u_shift)[0]   # u and v must have same dimensions
    ax2_length=np.shape(u_shift)[1]
    max_ind=int(np.max((ax1_length,ax2_length)))
    if max_ind%2!=0:
        max_ind+=1
    u_expand=np.zeros((max_ind,max_ind))
    v_expand=np.zeros((max_ind,max_ind))
    u_expand[:ax1_length,:ax2_length]=u_shift
    v_expand[:ax1_length,:ax2_length]=v_shift

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

    #4.1) calculate tractions in fourrier space T=K⁻¹*U, U=[u,v]  here with individual matrix elelemnts..
    tx_ft = kix * u_ft  +  kid * v_ft;
    ty_ft = kid * u_ft  +  kiy * v_ft;
    ##plt.imshow(tx_ft.real)
    #4.2) go back to real space
    tx=np.fft.ifft2(tx_ft).real
    ty=np.fft.ifft2(ty_ft).real

    #5.2) cut back to oringinal shape
    tx_cut=tx[0:ax1_length,0:ax2_length]
    ty_cut=ty[0:ax1_length,0:ax2_length]

    #5.3) using filter
    if filter=="mean":
        tx_filter=uniform_filter(tx_cut,size=5)
        ty_filter=uniform_filter(ty_cut,size=5)
    if filter == "gaussian":
        tx_filter = gaussian_filter(tx_cut, sigma=1.5)
        ty_filter = gaussian_filter(ty_cut, sigma=1.5)
    if filter == "median":
        tx_filter = median_filter(tx_cut, size=5)
        ty_filter = median_filter(ty_cut, size=5)


    #5.x) zooming out for represetnation purposes exactly to bf_image size
    if bf_image:
        zoom_factors=[np.shape(bf_image)[0]/np.shape(tx_cut)[0],np.shape(bf_image)[1]/np.shape(tx_cut)[1]]
        tx_zoom_out=zoom(tx_filter,zoom_factors)
        ty_zoom_out=zoom(ty_filter,zoom_factors)
        return (tx_filter, ty_filter,tx_zoom_out,ty_zoom_out)

    return(tx_filter,ty_filter)    # because of image axis inversion???



def calculate_deformation(im1,im2,window_size=64,overlapp=32,std_factor=20):

    '''
    Calculation of deformation field using particle image velocimetry (PIV). Recommendations: window_size should be about
    6 time the size of bead. Overlapp should be no less then half of the window_size. Std_factor should be kept as high as
    possibel. Make sure to check for to many exclusions caused by this factor e.g. by looking at the mask_std.
    Side note: returns -v because original v is negative if compared to coordinates of images (y-axis is inverted).


    :param file1: after iamge
    :param file2: before image
    :param window_size: integer, size of interrogation windows for PIV
    :param overlapp: integer, overlapp of interrogation windows for PIV
    :param std_factor: filterng extreme outliers beyond mean (deformation) + std_factor*standard deviation (deforamtion)
    :return:u,v deformation in x and y direction in pixel of the before and after image
            x,y psitions of the deformation fiedl in coordinates of the after and before image
            mask, mask_std  mask of filtered values by signal to noise filtering (piv internal) and filtering for
            extreme outliers
    '''
    #accepting either path to file or image data directly
    if isinstance(im1,str):
        frame_a  = np.array(openpiv.tools.imread( im1 ),dtype="int32")
    elif isinstance(im1,np.ndarray):
        frame_a=im1
    if isinstance(im2, str):
        frame_b  = np.array(openpiv.tools.imread( im2) ,dtype="int32")
    elif isinstance(im2, np.ndarray):
        frame_b = im2




    u, v, sig2noise = openpiv.process.extended_search_area_piv( frame_a, frame_b, window_size=window_size, overlap=overlapp,
                                                                dt=1,subpixel_method="gaussian", search_area_size=window_size, sig2noise_method='peak2peak' )
    x, y = openpiv.process.get_coordinates( image_size=frame_a.shape, window_size=window_size, overlap=overlapp)

    u, v, mask = openpiv.validation.sig2noise_val( u, v, sig2noise, threshold = 1.05 )

    def_abs = np.sqrt(u ** 2 + v ** 2)
    m = np.nanmean(def_abs)
    std = np.nanstd(def_abs)

    threshold = std * std_factor + m
    mask_std = def_abs > threshold
    u[mask_std] = np.nan
    v[mask_std] = np.nan

    u, v = openpiv.filters.replace_outliers( u, v, method='localmean', max_iter=10, kernel_size=2 )
    return(u,-v,x,y,mask,mask_std)     # return -v because of image inverted axis



def get_xy_for_quiver(u):
    '''
    accesoiry function to caclulate grid for plt.quiver. Size of the array will correspond to input u.
    :param u:any array,
    :return:
    '''
    xs=np.zeros(np.shape(u))
    for i in range(np.shape(u)[0]):
        xs[i,:]=np.arange(0,np.shape(u)[1],1)
    ys = np.zeros(np.shape(u))
    for j in range(np.shape(u)[1]):  ## is inverted in other skript...
        ys[:, j] = np.arange(0,np.shape(u)[0], 1)
    return xs, ys


def plotting_deformation(u,v):


    fig=plt.figure()
    def_abs=np.sqrt(u**2+v**2)
    plt.imshow(def_abs)
    cbar=plt.colorbar()
    cbar.set_label("displacement in pixel")
    x1,y1=get_xy_for_quiver(u)
    scale_ratio=0.2
    scale = scale_ratio * np.max(np.shape(u)) / np.max(def_abs)  # automatic scaleing in dependance of the image size
    plt.quiver(x1,y1,u*scale,v*scale,scale=1,headwidth=3, scale_units='xy',angles='xy')

    plt.arrow(75, -14, 10*scale, 0, head_width=0, facecolor="black", edgecolor="black",
              clip_on=False, head_starts_at_zero=False, overhang=2)
    plt.text(75, -10, "10 pixel displacement", color="black")


    plt.close()
    return fig

 #defo image:

def plotting_traction(tx, ty):
    fig = plt.figure()
    plt.imshow(np.sqrt((tx/1000) ** 2 + (ty/1000) ** 2),cmap="rainbow")
    cbar = plt.colorbar()
    cbar.set_label("traktion forces in kPa")
    x1, y1 = get_xy_for_quiver(tx)
    scale_ratio=0.2
    scale = scale_ratio * np.max(np.shape(tx)) / np.max(np.sqrt((tx/1000) ** 2 + (ty/1000) ** 2))  # automatic sacleing in dependace of the image size #
    plt.quiver(x1, y1, (tx/1000) * scale, (ty/1000) * scale, scale=1, scale_units='xy',angles="xy")
    plt.arrow(75, -14, 1 * scale, 0, head_width=0, facecolor="black", edgecolor="black",
              clip_on=False, head_starts_at_zero=False, overhang=2)
    plt.text(75, -10, "1 kPa Force", color="black")
    plt.close()
    return fig


### energy calculation and gui



def contractillity(tx,ty,pixelsize,mask):
    '''
    Calculation of contractile force and force epicenter.Contractile force is the sum of all projection of traction
    forces (in N) towards the force epicenter. The force epicenter is the point that maximizes the contractile force.
    :param tx: traction forces in x direction in Pa
    :param ty: traction forces in y direction in Pa
    :param pixelsize: pixelsize of the traction field
    :param mask: mask of which values to use for calculation
    :return: contractile_force,contractile force in N
             proj_x, projection of traction forces towards the foce epicenter, x component
             proj_y, projection of traction forces towards the foce epicenter, y component
             center, coordinates of the force epicenter
    '''
    mask=mask.astype(bool)
    tx_filter = np.zeros(np.shape(tx))
    tx_filter[mask] = tx[mask]

    ty_filter = np.zeros(np.shape(ty))
    ty_filter[mask] = ty[mask]

    tract_abs = np.sqrt(tx_filter ** 2 + ty_filter ** 2)

    area = (pixelsize * (10 ** -6)) ** 2  # in meter
    fx = tx_filter * area  ### calculating forces by multiplyin with area # forces then in newton
    fy = ty_filter * area

    ## calculate forces wiht area??
    # mask=np.ones(np.shape(tx))

    x, y = get_xy_for_quiver(tx)
    bx = np.sum(x * (tract_abs ** 2) + fx * (tx_filter * fx + ty_filter * fy))  ## meaning of ths vector??
    by = np.sum(y * (tract_abs ** 2) + fy * (tx_filter * fx + ty_filter * fy))  #

    axx = np.sum(tract_abs ** 2 + fx ** 2)
    axy = np.sum(fx * fy)  # redundant
    ayy = np.sum(tract_abs ** 2 + fy ** 2)
    # ayx=np.sum(tx*ty)

    A = np.array([[axx, axy], [axy, ayy]])
    b = np.array([bx, by]).T

    # solve equation system:
    # center*[bx,by]=[[axx,axy],
    #                [axy,ayy]]
    # given by A*[[1/bx],
    #             [1/by]
    center = np.matmul(np.linalg.inv(A), b)

    # vector projection to origin

    dist_x = center[0] - x
    dist_y = center[1] - y
    dist_abs = np.sqrt(dist_y ** 2 + dist_x ** 2)
    proj_abs = (fx * (dist_x) + fy * (dist_y)) / dist_abs
    contractile_force = np.nansum(proj_abs)

    # project_vectors
    proj_x = proj_abs * dist_x / dist_abs
    proj_y = proj_abs * dist_y / dist_abs

    return contractile_force, proj_x, proj_y,center # unit of contractile force is N


def contractile_energy_points(u,v,tx,ty,pixelsize1,pixelsize2):
    pixelsize2*=10**-6 # conversion to m
    pixelsize1*=10**-6
    energy_points = 0.5 * (pixelsize2 ** 2) * (np.sqrt((tx  * u * pixelsize1) ** 2 + (
            ty  * v * pixelsize1) ** 2))

    bg = np.percentile(energy_points, 30)  # value of a background point
    energy_points-=bg
    return energy_points


def plot_deformation_few_arrows(u,v,vmin=0,vmax=False,scale_ratio=0.5):
    '''
    plotting deformation with additional gaussian filter and only showing every 3rd arrow (in x and y direction)
    :param u:
    :param v:
    :param vmin:
    :param vmax:
    :param scale_ratio:
    :return:
    '''

    pixx = np.arange(np.shape(u)[0])
    pixy = np.arange(np.shape(u)[1])
    xv, yv = np.meshgrid(pixy, pixx)
    def_abs = np.sqrt((u ** 2 + v ** 2))
    u_show = gaussian_filter(u, sigma=3)
    v_show = gaussian_filter(v, sigma=3)  # blocks of edge size 3
    select_x = ((xv - 1) % 3) == 0
    select_y = ((yv - 1) % 3) == 0
    selct_size = def_abs > 1
    select = select_x * select_y * selct_size
    u_show[~select] = 0
    v_show[~select] = 0

    if not vmax:
        vmax = np.max(def_abs)
    # defo image:
    fig = plt.figure()
    im=plt.imshow(def_abs, vmin=vmin, vmax=vmax, cmap="rainbow")
    x1, y1 = get_xy_for_quiver(u)

    scale = scale_ratio * np.max(np.shape(u)) / np.max(def_abs)  # automatic sacleing in dependace of the image size

    plt.quiver(x1, y1, u_show * scale, v_show * scale, scale=1, headwidth=3, scale_units='xy', width=0.002)
    plt.arrow(75, -14, 5 * scale, 0, head_width=0, facecolor="black", edgecolor="black",
              clip_on=False, head_starts_at_zero=False, overhang=2)
    plt.text(75, -10, "5 pixel displacement", color="black")
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = plt.colorbar(mappable=im, cax=cax)
    cbar.set_label("displacement in pixel", fontsize=15)


def plot_deformation_all_arrows(u,v,vmin=0,vmax=False,scale_ratio=0.5,bar_length=5):

    pixx = np.arange(np.shape(u)[0])
    pixy = np.arange(np.shape(u)[1])
    xv, yv = np.meshgrid(pixy, pixx)
    def_abs = np.sqrt((u ** 2 + v ** 2))


    if not vmax:
        vmax = np.max(def_abs)
    # defo image:
    fig = plt.figure()
    im=plt.imshow(def_abs, vmin=vmin, vmax=vmax, cmap="rainbow")
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    u_filter=copy.deepcopy(u)
    u_filter[def_abs < vmin]=0
    v_filter=copy.deepcopy(v)
    v_filter[def_abs < vmin] = 0

    x1, y1 = get_xy_for_quiver(u)

    scale = scale_ratio * np.max(np.shape(u)) / np.max(def_abs)  # automatic sacleing in dependace of the image size

    plt.quiver(x1, y1, u_filter * scale, v_filter * scale, scale=1, headwidth=3, scale_units='xy',angles="xy", width=0.002)


    plt.arrow(np.shape(u)[0]*1,0-np.shape(u)[1]*0.1, bar_length * scale, 0, head_width=0, facecolor="black", edgecolor="black",
              clip_on=False, head_starts_at_zero=False, overhang=2)
    plt.text(np.shape(u)[0]*1,0-np.shape(u)[1]*0.05, "%d pixel displacement" % bar_length, fontsize=10, color="black")
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = plt.colorbar(mappable=im, cax=cax)
    cbar.set_label("displacement in pixel", fontsize=15)





def contractile_projection(proj_x,proj_y,tx,ty,mask,center,contractile_force):
    '''
    plotting the porjection of traction froces towards the force epcenter
    :param proj_x:
    :param proj_y:
    :param tx:
    :param ty:
    :param mask:
    :param center:
    :param contractile_force:
    :return:
    '''


    fig = plt.figure()
    custom_cmap1 = LinearSegmentedColormap.from_list("", ["#DBDC3E",
                                                                            "yellow"])
    im=plt.imshow(np.sqrt((tx / 1000) ** 2 + (ty / 1000) ** 2),vmin=0)
    mask_show=np.zeros(np.shape(mask))+np.nan
    mask_show[mask]=1
    plt.imshow(mask_show,alpha=0.5,cmap= custom_cmap1 )
    plt.plot(center[0], center[1],"or", color="red")

    x1, y1 = get_xy_for_quiver(tx)
    ratio = 0.2  # ratio of length of biggest arrow to max axis lenght
    scale = ratio * np.max(np.shape(proj_x)) / np.max(
        np.sqrt(proj_x ** 2 + proj_y ** 2))  # automatic sacleing in dependace of the image size
    plt.quiver(x1, y1, proj_x * scale, proj_y * scale, angles="xy", scale=1, scale_units='xy', width=0.002)

    plt.text(10, 70, "contractile force = " + str(np.round(contractile_force * 10 ** 6, 2)), color="black")

    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = plt.colorbar(mappable=im, cax=cax)
    cbar.set_label('traction forces in kPa')

def gif_mask_overlay(im1,im2,mask):
    im1.setflags(write=1)
    im2.setflags(write=1)

    im1[:,:,0]+=mask.astype("uint8")*100
    im2[:,:,0]+=mask.astype("uint8")*100
    images = [im1,im2]
    return images

