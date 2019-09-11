
import matplotlib.animation as animation
from imageio import imread
from skimage.filters import gaussian
import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable
sys.path.insert(0, '/media/user/GINA1-BK/Andreas-Python/tracktion_force_microscopy/TFM_functions.py')
from TFM_functions import *




out_file=r"/home/user/Desktop/01out.txt"    ## output text file
file1=r"/media/user/GINA1-BK/data_traktion_force_microscopy/good_exampel/11after_shift.tif" # after image
file2=r"/media/user/GINA1-BK/data_traktion_force_microscopy/good_exampel/11before_shift.tif" # before image
file3=r"/media/user/GINA1-BK/data_traktion_force_microscopy/good_exampel/11bf_before_shift.tif" # bf image



sigma=0.49 #poison ratio
young=25536 # youngsmodulus in Pa
pixelsize = 6.25 / 40 # µm/pixel pixelsize of the original images
window_size = 64  # windowsize for piv calculation of the deformation field
overlapp = 32 # overlapp for deformation field calculation, should be half of windowsize
pixel_factor = window_size - overlapp
pixelsize_defo = pixelsize * pixel_factor # pixelsize of the deformation field









# writing  new outputfile or apppending to existing file
if not os.path.isfile(out_file):
    with open((out_file), "w") as f:
        f.write("filename\tstrain_energy\tcontractile force\n")



bf_image=imread(file3)
print("calculating deformation")
u,v,x,y,mask,mask_std=calculate_deformation(file1,file2,window_size=window_size,overlapp=overlapp)
print("calculating traktion")
tx,ty=ffttc_traction(u,v,young,pixelsize,pixelsize_defo,bf_image=False,filter="mean")  # u and v shifted internally

np.save("/media/user/GINA1-BK/data_traktion_force_microscopy/temp_data_for_scripts/tx.npy",tx)
np.save("/media/user/GINA1-BK/data_traktion_force_microscopy/temp_data_for_scripts/ty.npy",ty)


t_abs=np.sqrt(tx ** 2 + ty ** 2)
u_shift = u - np.mean(u)    # mean displacement set to zero, used for caclualtion of strain energy
v_shift = v - np.mean(v)

## pre calulating energy on all points
# as the absolute value of deformation*force/force in N
print("calculating strain energies")
energy_points = 0.5 * (pixelsize_defo**2)  * (np.sqrt((tx/1000 * u_shift* pixelsize) ** 2 + (
        ty/1000 * v_shift* pixelsize) ** 2)) / 1000
bg = np.percentile(energy_points, 50) # value of one background point



### gui####

custom_cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#DBDC3E","yellow"])
#x and y coordinates for quiver and region selection
pixx = np.arange(np.shape(tx)[0])
pixy = np.arange(np.shape(tx)[1])
xv, yv = np.meshgrid(pixy, pixx)
pix = np.vstack((xv.flatten(), yv.flatten())).T


# displaying only every 9th arrow with additional filter
# this is purely cosmetic, no influence on calculation
tx_show=gaussian(tx,sigma=3)
ty_show=gaussian(ty,sigma=3) # blocks of edge size 3
select_x=((xv-1)%3)==0 # mask for every third value in x direction
select_y=((yv-1)%3)==0
select_size=t_abs>100  # threshold for displaying traction arrows
select=select_x*select_y*select_size
tx_show[~select]=0# everything except selection set to zero
ty_show[~select]=0



plt.close("all")
fig = plt.figure(figsize=(12,4))
ax1 = plt.axes([0.45, 0.25, 0.3, 0.6]) # ax traction forces
ax2 = plt.axes([0.05, 0.25, 0.3, 0.6]) # ax img  bf

ax_text=plt.axes([0.45, 0.1, 0.1, 0.1])# axes for plotting text
ax_text.set_axis_off()


im = ax1.imshow(t_abs/1000,cmap="rainbow")
ratio=0.2 # ratio of length of biggest arrow to max axis lenght
scale=ratio *np.max(np.shape(tx))/np.max(np.sqrt(tx_show**2+ty_show**2))# automatic scaleing in dependance of the image size
ax1.quiver(xv, yv, tx_show*scale, ty_show*scale,scale=1, scale_units='xy',angles='xy')
ax2.imshow(bf_image,cmap="gray")

# colorbar
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.1)
cbar = plt.colorbar(mappable=im, cax=cax)
cbar.set_label('traction forces in kPa')

# setting  variables for calculation
ind = [] # indices of line sorunding selected region
energy = 0  # strain energy
contractile_force=0 # contractile force
plt_str1 = 0 # str for writing strain energy to plot
plt_str2 = 0 # str for writing contractile energy to plot



def onselect1(verts):
    '''
    main function to handle selection and displaying traction, bf image, selected area. Also caclulation of
     strain energy  and contractile force

    :param verts:
    :return:
    '''

    global array, pix, ind, energy, plt_str1,plt_str2,ax_text,energy_points,bg,mask,verts_out,p,mask_zoom,ax1,mask_zoom2
    global contractile_force,center,select,tx_show,ty_show,scale,lasso,proj_x,proj_y,xv, yv

    verts_out = verts # vertices of curve drawn by lasso selector

    if lasso2.verts: ## slection in brightfiled
    # simple interpolation of points selected in bf image to traction field
        interpol_factors = np.array([np.shape(tx)[1] / np.shape(bf_image)[1], np.shape(tx)[0] / np.shape(bf_image)[0]])
        verts_array=np.array(verts_out)
        verts_interpol=np.zeros(np.shape( verts_array))
        for i in range(np.shape(verts_interpol)[0]):
            verts_interpol[i,:]=verts_array[i,] * interpol_factors  # element wise multiplikatioon for each row
        p = path.Path(verts_interpol) # internal function to build path from vertices
    if lasso1.verts: # selection in fl image
        p = path.Path(verts) # internal function to build path from vertices

    #retrieving mask
    ind = p.contains_points(pix, radius=1) # all points inside of the selected area
    mask = np.reshape(np.array(ind), (np.shape(tx)))

    #showing mask as overlay adn updating the overlay
    mask_show=np.zeros(np.shape(tx))+np.nan
    mask_show[mask]=1
    ax1.clear()
    ax1.imshow(t_abs/1000,alpha=0.8,cmap="rainbow")
    ax1.quiver(xv, yv, tx_show*scale, ty_show*scale,scale=1, scale_units='xy',angles="xy")
    ax1.imshow(mask_show,alpha=0.4,cmap=custom_cmap1)

    #energy caclulation
    # sum of all energy points under the mask substracted with backround
    energy = np.sum(energy_points[mask]) - bg * np.sum(mask)
    print("strain energy" ,energy)

    # contractile forces calculation
    contractile_force, proj_x, proj_y, center = contractility(tx, ty, pixelsize_defo, mask)
    print("contractile force ", contractile_force)
    ax1.plot(center[0],center[1],"or",color="red") # plotting force epicenter

    # writing energy text
    if plt_str1 != 0:
        plt_str1.remove()
    plt_str1 = ax_text.text(0, 0, "strain energy = " + str(np.round(energy, 2))+" pJ")

    # writing contractile force text
    if plt_str2 != 0:
        plt_str2.remove()
    plt_str2 = ax_text.text(0, 0.5, "contractile force = " + str(np.round(contractile_force*10**6, 2)) + " µN")
    fig.canvas.draw_idle()


lasso2 = LassoSelector(ax2, onselect1) # selector on bf image
lasso1 = LassoSelector(ax1, onselect1) # selector on traction field image









# showing deformation
def show_df(event):
    global u,v
    plot_deformation_all_arrows(u, v,vmin=1,scale_ratio=0.2,bar_length=2)


# defining button position and properties
button_ax_df = plt.axes([0.8, 0.06, 0.13, 0.075])
button_ax_df.xaxis.set_ticks_position('none')
button_ax_df.yaxis.set_ticks_position('none')
button_ax_df.set_xticks([])
button_ax_df.set_yticks([])
b_show_df = Button(button_ax_df, 'show deformation')
b_show_df.on_clicked(show_df)





# showing gif of before and after image
ani=0 # animation needs to be stored in external variable to actually show up
def show_gif(event):
    global file1,file2,ani
    fig = plt.figure()
    imgs = [[plt.imshow(plt.imread(file1), animated=True)], [plt.imshow(plt.imread(file2), animated=True)]]

    ani = animation.ArtistAnimation(fig, imgs, interval=180, blit=False,
                           repeat_delay=0)
    plt.show()
    #return ani

# defining button position and properties
button_ax_gif = plt.axes([0.8, 0.16, 0.13, 0.075])
button_ax_gif.xaxis.set_ticks_position('none')
button_ax_gif.yaxis.set_ticks_position('none')
button_ax_gif.set_xticks([])
button_ax_gif.set_yticks([])
b_show_gif = Button(button_ax_gif, 'show gif')
b_show_gif.on_clicked(show_gif)









# saving contractile force and strain energy to text file
# units would be µN and pJ respectively
def save(event):
    global energy, out_file
    with open((out_file), "a") as f:  # appending lines o text file
        f.write(file1 + "\t" + str(np.round(energy, 2)) +"\t"+ str(np.round(-contractile_force*10**6, 2))+"\n")
    print("saved line " + file1 + "\t" + str(np.round(energy, 2)) +"\t"+ str(np.round(-contractile_force*10**6, 2))+"\n")

# defining button position and properties
button_ax_s = plt.axes([0.8, 0.26, 0.13, 0.075])
button_ax_s.xaxis.set_ticks_position('none')
button_ax_s.yaxis.set_ticks_position('none')
button_ax_s.set_xticks([])
button_ax_s.set_yticks([])
b_save = Button(button_ax_s, 'save')
b_save.on_clicked(save)




def show_contractility_projection(event):
    global proj_x, proj_y, tx, ty, center, contractile_force,mask
    contractile_projection(proj_x,proj_y,tx,ty,mask,center,contractile_force)

# defining button position and properties
button_ax_pro = plt.axes([0.8, 0.36, 0.13, 0.075])
button_ax_pro.xaxis.set_ticks_position('none')
button_ax_pro.yaxis.set_ticks_position('none')
button_ax_pro.set_xticks([])
button_ax_pro.set_yticks([])
b_pro = Button(button_ax_pro, 'show_contractility\nprojection')
b_pro.on_clicked(show_contractility_projection)



## some evaluation massages
# message if extremely large deformation was filtered while caclulating deformation field,
# this filtering is only based on standard deviation
# one should be concerned if this is >5
if np.sum(mask_std)>0:
    print("filtered "+str(np.sum(mask_std))+ " values as outliers" )




