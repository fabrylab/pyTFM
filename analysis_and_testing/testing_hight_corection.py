import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import os
from tqdm import tqdm
import sys
from simulating_deformation import *
from andreas_TFM_package.functions_for_cell_colonie import *
from andreas_TFM_package.utilities import convert_to_int



def contractility_comparison(fx, fy, hs, pixelsize, young, sigma,edge_dist=0.3):
    '''
    produces a plot of contractilities for traction forces prediction assuming finite thikness and infinite thikness.
    Deformations are calculated using the forces fx and fy, assuming finite thikness of hs
    :param fx:
    :param fy:
    :param hs:
    :return:
    '''

    mask = np.ones(np.shape(fx)) < 0  ## contractile force on full image, this is not to robust
    mask[int(fx.shape[0] * edge_dist):-int(fx.shape[0] * edge_dist),
         int(fx.shape[1] * edge_dist):-int(fx.shape[1] * edge_dist)] = True

    tx_h_list = []
    ty_h_list = []
    tx_list = []
    ty_list = []
    contractile_forces_h = []
    contractile_forces = []

    show_quiver(fx, fy, scale_ratio=0.2,headwidth=10)

    for j, h in enumerate(hs):

        def_x_h1, def_y_h1, def_z_h1 = finite_thickness_convolution(fx, fy, pixelsize, h* 10 ** -6, young, sigma=sigma)
        #show_quiver(def_x_h1, def_y_h1)


       # pixx = np.arange(fx.shape[0])
        #pixy = np.arange(fx.shape[1])
       # coords_x, coords_y = np.meshgrid(pixy, pixx)  # matrices of distance
        #def_x_h1, def_y_h1, def_z_h1 = finite_thickness(fx, fy, 10 * 10 ** -6, young, pixelsize, coords_x + 0.001,
                                                       # coords_y + 0.001,
                                                        #sigma=sigma)
        #show_quiver(def_x_h1, def_y_h1)
        #
        def_x_h1_pix = def_x_h1 / (pixelsize)  ## deformation in terms of pixel
        def_y_h1_pix = def_y_h1 / (pixelsize)

        tx_h, ty_h = ffttc_traction_finite_thickness(def_x_h1_pix, def_y_h1_pix, pixelsize1=pixelsize / (10 ** -6),
                                                     pixelsize2=pixelsize / (10 ** -6), young=young, h=h,
                                                     sigma=0.5, filter=None)
        tx, ty = ffttc_traction(def_x_h1_pix, def_y_h1_pix, pixelsize1=pixelsize / (10 ** -6),
                                pixelsize2=pixelsize / (10 ** -6), young=young, sigma=0.5, filter=None)




        # contractile force with projection method
        contractile_force_h, proj_x_h, proj_y_h, center = contractility(tx_h, ty_h, pixelsize / (10 ** -6),
                                                                        mask)  # note sign is arbitrary
        # sum of input forces as estimate of contractility/# sum of the norm of all force vectors
        e = np.sum(np.linalg.norm(np.stack((fx, fy), axis=2), axis=2))
        print(e, contractile_force_h)

        # contractile force with projection method
        contractile_force, proj_x_h, proj_y_h, center = contractility(tx, ty, pixelsize / (10 ** -6),
                                                                      mask)  # note sign is arbitrary
        # sum of input forces as estimate of contractility/# sum of the norm of all force vectors
        e = np.sum(np.linalg.norm(np.stack((fx, fy), axis=2), axis=2))
        print(e, contractile_force)
        tx_h_list.append(tx_h)
        ty_h_list.append(ty_h)
        tx_list.append(tx)
        ty_list.append(ty)
        contractile_forces_h.append(contractile_force_h)
        contractile_forces.append(contractile_force)

        if h == 5: # comparrison of deformations
            show_quiver(def_x_h1, def_y_h1)
            dx,dy=infinite_thickness_convolution(fx, fy, pixelsize, young, sigma=sigma)
            show_quiver(dx,dy)

    # plotting contractilitys
    fig=plt.figure()
    plt.plot(hs, -np.array(contractile_forces), label="non-corrected FFTC")
    plt.plot(hs, -np.array(contractile_forces_h), label="corrected FFTC")
    plt.hlines(e, hs[0], hs[-1], label="sum of input forces")
    plt.xlabel("gel hight in µm")
    plt.ylabel("contractility in N")
    plt.axvline(300, 0,1,linestyle="--",color="green")
    plt.axvline(7, 0, 1,linestyle="--",color="green")
    plt.legend()
    return fig


def contractility_comparison_from_deformation(def_x, def_y, hs, pixelsize1,pixelsize2, young, sigma,mask):
    '''
     compares contractility on real deformation image
    :param fx:
    :param fy:
    :param hs:
    :return:
    '''



    tx_h_list = []
    ty_h_list = []
    tx_list = []
    ty_list = []
    contractile_forces_h = []
    contractile_forces = []

    show_quiver(def_x, def_y, scale_ratio=0.2,headwidth=6,filter=[0,4])

    for j, h in tqdm(enumerate(hs),total=len(hs)):
        tx_h, ty_h = ffttc_traction_finite_thickness(def_x, def_y, pixelsize1=pixelsize1 / (10 ** -6),
                                                     pixelsize2=pixelsize2/ (10 ** -6), young=young, h=h,
                                                     sigma=0.5, filter="gaussian")
        tx, ty = ffttc_traction(def_x, def_y, pixelsize1=pixelsize1 / (10 ** -6),
                                pixelsize2=pixelsize2 / (10 ** -6), young=young, sigma=0.5, filter="gaussian")




        # contractile force with projection method
        contractile_force_h, proj_x_h, proj_y_h, center = contractility(tx_h, ty_h, pixelsize2 / (10 ** -6),
                                                                        mask)  # note sign is arbitrary
        # sum of input forces as estimate of contractility/# sum of the norm of all force vectors

        # contractile force with projection method
        contractile_force, proj_x, proj_y, center = contractility(tx, ty, pixelsize2 / (10 ** -6),
                                                                      mask)  # note sign is arbitrary
        # sum of input forces as estimate of contractility/# sum of the norm of all force vectors

        tx_h_list.append(tx_h)
        ty_h_list.append(ty_h)
        tx_list.append(tx)
        ty_list.append(ty)
        contractile_forces_h.append(contractile_force_h)
        contractile_forces.append(contractile_force)



    # plotting contractilitys
    fig=plt.figure()
    plt.plot(hs, -np.array(contractile_forces)*10**9, label="non-corrected FFTC")
    plt.plot(hs, -np.array(contractile_forces_h)*10**9, label="corrected FFTC")
    plt.xlabel("gel hight in µm")
    plt.ylabel("contractility in Nn")
    plt.legend()
    #plt.gca().invert_yaxis()
    plt.vlines(300,ymin=np.min(-np.array(contractile_forces_h)*10**9),ymax=np.max(-np.array(contractile_forces_h)*10**9)+1000)
    return fig



def showing_multiple_force_predictions(fx, fy, hs, pixelsize, young, sigma,mask):
    '''
    produces a supblot with 20 slots. Shows one kind of prediction in the first two rows and the other kind in the second two rows
    :param mask:
    :return:
    '''
    tx_h_list = []
    ty_h_list = []
    tx_list = []
    ty_list = []
    contractile_forces_h = []
    contractile_forces = []

    for j, i in enumerate(hs):
        h = i * 10 ** -6
        # def_x,def_y= deformation_from_forces_convolution(fx,fy,pixelsize,young,sigma=0.5)
        def_x_h1, def_y_h1, def_z_h1 = finite_thickness_convolution(fx, fy, pixelsize, h, young, sigma=0.5)
        #show_quiver(def_x_h1, def_y_h1)
        # greens_tensor_exact, As_exact = finite_thickenss_convolution_exact_greens_tensor(fx, fy, pixelsize, h, young,
        #                                                                                 sigma=0.5, kernel_size=None)
        # np.save("/media/user/GINA1-BK/data_traktion_force_microscopy/tensors_for_finite_thikness/greens_tensor_exact" + str(i) + ".npy",  greens_tensor_exact)
        # greens_tensor_exact=np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/tensors_for_finite_thikness/greens_tensor_exact" + str(i) + ".npy")
        # def_x_h2,def_y_h2,def_z_h2=finite_thickenss_convolution_only(fx, fy, greens_tensor_exact)
        # show_quiver(def_x_h2, def_y_h2)
        # def_x_pix=def_x/(pixelsize)## deformation in terms of pixel
        # def_y_pix=def_y/(pixelsize)

        def_x_h1_pix = def_x_h1 / (pixelsize)  ## deformation in terms of pixel
        def_y_h1_pix = def_y_h1 / (pixelsize)
        # def_x_h2_pix = def_x_h2 / (pixelsize)  ## deformation in terms of pixel
        # def_y_h2_pix = def_y_h2 / (pixelsize)
        # show_quiver(def_x_h_pix, def_y_h_pix)

        tx_h, ty_h = ffttc_traction_finite_thickness(def_x_h1_pix, def_y_h1_pix, pixelsize1=pixelsize / (10 ** -6),
                                                     pixelsize2=pixelsize / (10 ** -6), young=young, h=h / (10 ** -6),
                                                     sigma=0.5, filter="gaussian")
        tx, ty = ffttc_traction(def_x_h1_pix, def_y_h1_pix, pixelsize1=pixelsize / (10 ** -6),
                                pixelsize2=pixelsize / (10 ** -6), young=young, sigma=0.5, filter="gaussian")
        # tx,ty=ffttc_traction_pure_shear(def_x_h_pix,def_y_h_pix,pixelsize1=pixelsize/(10**-6),pixelsize2=pixelsize/(10**-6),h=h/(10**-6),young=young,sigma=0.5,filter="gaussian")
        # tx, ty = ffttc_traction_finite_thickness(def_x_h2_pix, def_y_h2_pix, pixelsize1=pixelsize / (10 ** -6),
        #                                   pixelsize2=pixelsize / (10 ** -6),h=h/(10**-6), young=young, sigma=0.5, filter="gaussian")

        # contractile force with projection method
        contractile_force_h, proj_x_h, proj_y_h, center = contractility(tx_h, ty_h, pixelsize / (10 ** -6),
                                                                        mask)  # note sign is arbitrary
        # sum of input forces as estimate of contractility/# sum of the norm of all force vectors
        e = np.sum(np.linalg.norm(np.stack((fx, fy), axis=2), axis=2))
        print(e, contractile_force_h)

        # contractile force with projection method
        contractile_force, proj_x_h, proj_y_h, center = contractility(tx, ty, pixelsize / (10 ** -6),
                                                                      mask)  # note sign is arbitrary
        # sum of input forces as estimate of contractility/# sum of the norm of all force vectors
        e = np.sum(np.linalg.norm(np.stack((fx, fy), axis=2), axis=2))
        print(e, contractile_force)
        tx_h_list.append(tx_h)
        ty_h_list.append(ty_h)
        tx_list.append(tx)
        ty_list.append(ty)
        contractile_forces_h.append(contractile_force_h)
        contractile_forces.append(contractile_force)

        # show_quiver(def_x_h,def_y_h)
        # show_quiver(tx,ty,scale_ratio=0.2)

    # show_quiver(def_x,def_y)
    # show_quiver(tx,ty)

    fig, axs = plt.subplots(4, 4)
    fig.set_size_inches(15, 10)
    axs = axs.flatten()
    # plt.title("force reconstruction from deformations in finitely thick mebrane")

    vmax = np.nanmax([[np.sqrt(tx ** 2 + ty ** 2) for tx, ty in zip(tx_h_list, ty_h_list)],
                      [np.sqrt(tx ** 2 + ty ** 2) for tx, ty in zip(tx_list, ty_list)]])
    vmin = np.nanmin([[np.sqrt(tx ** 2 + ty ** 2) for tx, ty in zip(tx_h_list, ty_h_list)],
                      [np.sqrt(tx ** 2 + ty ** 2) for tx, ty in zip(tx_list, ty_list)]])
    ### amke an actually nice plot out of this...

    #
    # plt.figure()
    # plt.plot(hs,-np.array(contractile_forces))
    # plt.plot(hs,-np.array(contractile_forces_h))
    # plt.hlines(e,hs[0],hs[-1])

    factor = 10 ** 11
    for j, (h, tx_h, ty_h, tx, ty, c_h, c) in enumerate(
            zip(hs, tx_h_list, ty_h_list, tx_list, ty_list, contractile_forces_h, contractile_forces)):
        axs[j] = show_quiver(tx_h, tx_h, vmax=vmax, vmin=vmin, scale_ratio=0.3,ax=axs[j])
        axs[j].text(1, 10, "sum of force " + str(np.round(e * factor, 2)), color="white")
        axs[j].text(1, 20, "contractility " + str(np.round(c_h * factor, 2)), color="white")
        axs[j].text(1, 30, "h " + str(np.round(h, 2)), color="white")
        fig,axs[j + 8] = show_quiver(axs[j + 8], tx, tx, vmax=vmax, vmin=vmin, scale_ratio=0.3)
        axs[j + 8].text(1, 10, "sum of force " + str(np.round(e * factor, 2)), color="white")
        axs[j + 8].text(1, 20, "contractility " + str(np.round(c * factor, 2)), color="white")
        axs[j + 8].text(1, 30, "h " + str(np.round(h, 2)), color="white")

    pos = axs[-1].get_position()
    bounds = [pos.bounds[0], pos.bounds[1], pos.width * 0.2, pos.height]
    axs[-1].remove()
    axs[-1] = fig.add_axes(bounds)
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    cb1 = matplotlib.colorbar.ColorbarBase(axs[-1], cmap="viridis",
                                           norm=norm, label="traction forces",
                                           orientation='vertical')
    plt.draw()
    return fig



plt.ioff()
#for i in [0.5,1,2,3,4,5,6,7,8,9,10,12,14,16]:
for i in [5]:
    ## gel parameter
    sigma = 0.5  # poison ratio
    young = 25536  # youngs modulus in Pa??
    dims = (100, 100)
    pixelsize_bf = 0.10225 * (10 ** -6)
    pixelsize = pixelsize_bf * 25  *i
    edge_dist = 0.05
    hs = [5, 6, 8,10, 12,15, 20,30,40,50,60,70,80,90,100,200,300,350]  # in µm

    ##
    ## example value was 1500 Pa on 100X100 image with pixelsize =25*10**-6 m
    f = 1500 * ((pixelsize / i) ** 2)

    fx = np.zeros((100, 100))
    fy = np.zeros((100, 100))
    fx[45, 45] = -f
    fy[45, 45] = -f
    fx[55, 55] = f
    fy[55, 55] = f
    wave_length = 10 * pixelsize * 10 ** 6

    #fig=show_quiver(fx, fy, scale_ratio=0.2, headwidth=10)
    #ax=fig.gca()
    #convert_axis_tick_unit(ax,   pixelsize * 10 ** 6)
    #ax.set_xlabel("coordinates in µm")
    #ax.set_ylabel("coordinates in µm")

    fig=contractility_comparison(fx, fy, hs, pixelsize, young, sigma,edge_dist=0.3)
    plt.savefig("/home/user/Desktop/results/results_29_to03_august/contractility_comparrison_wl%d.png"%int(np.round(wave_length)))



## gel parameter
sigma = 0.5  # poison ratio
young = 25536  # youngs modulus in Pa??
pixelsize_bf = 0.10225 * (10 ** -6)
pixelsize = pixelsize_bf * 11 ## around the correct pixel size
edge_dist = 0.05
hs = [5, 6, 8,10, 12,15, 20,30,40,50,60,70,80,90,100,200,300,350]  # in µm




def_x=np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/out2/11u_200.npy")
def_y=np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/out2/11v_200.npy")

u_ft = np.fft.fft2(def_x * pixelsize * 2 * np.pi)
v_ft = np.fft.fft2(def_y * pixelsize * 2 * np.pi)

#tx, ty = ffttc_traction_finite_thickness(def_x, def_y, pixelsize1=pixelsize_bf / (10 ** -6),
#                                                     pixelsize2=pixelsize/ (10 ** -6), young=young, h=h,
#                                                     sigma=0.5, filter="gaussian")
fig=show_quiver(def_x,def_y,scale_ratio=0.2,filter=[0,3])
ax = fig.gca()
convert_axis_tick_unit(ax, pixelsize * 10 ** 6)
ax.set_xlabel("coordinates in µm")
ax.set_ylabel("coordinates in µm")

mask = np.ones(np.shape(def_x)) < 0
mask[int(def_x.shape[0] * edge_dist):-int(def_x.shape[0] * edge_dist),
        int(def_x.shape[1] * edge_dist):-int(def_x.shape[1] * edge_dist)] = True

fig=contractility_comparison_from_deformation(def_x, def_y, hs, pixelsize1=pixelsize_bf,pixelsize2=pixelsize,
                                              young=young, sigma=sigma,mask=mask)








































## gel parameter
sigma=0.5#poison ratio
young=25536# youngs modulus in Pa??
dims=(100,100)
h=600*10**-6
pixelsize_bf=0.10225*(10**-6)
pixelsize=pixelsize_bf*25
edge_dist=0.3
#
## example value was 1500 Pa on 100X100 image with pixelsize =25*10**-6 m
f=1500*(pixelsize**2)

fx=np.zeros((100,100))
fy=np.zeros((100,100))
fx[45,45]=-f
fy[45,45]=-f
fx[55,55]=f
fy[55,55]=f
##

mask=np.ones(np.shape(fx))<0 ## contractile force on full image, this is not to robust
mask[int(fx.shape[0]*edge_dist):-int(fx.shape[0]*edge_dist),int(fx.shape[1]*edge_dist):-int(fx.shape[1]*edge_dist)]=True
hs=[1,5,30,40,50,60,70] # in µm

showing_multiple_force_predictions(fx, fy, hs, pixelsize, young, sigma,mask)









