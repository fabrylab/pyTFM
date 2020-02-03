'''
plotting and other accesory functions for the Monolayer stress Microscopy

'''
import matplotlib
#matplotlib.use("TkAgg")
#matplotlib.use("Agg") # when running the script to avoid plots popping up
import matplotlib.pyplot as plt
from andreas_TFM_package.TFM_functions import *
from matplotlib import cm  # making list from colormap
from tqdm import tqdm
from matplotlib.ticker import AutoMinorLocator, MultipleLocator,ScalarFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#from matplotlib.figure import Figure
#from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from scipy.optimize import least_squares
from skimage.measure import regionprops

from scipy.sparse.csr import csr_matrix
from scipy.sparse.linalg import spsolve,lsqr
import os
from collections import defaultdict
from scipy.ndimage import binary_erosion


def make_discrete_colorbar():
    '''
    function to make a three pieced colormap, only used for represetatinof loads in labt atlk example.
    :param values:
    :return:
    '''

    cmap = plt.cm.jet  # define the colormap
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # force the first color entry to be grey

    # create the new map
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'Custom cmap', cmaplist, cmap.N)

    # define the bins and normalize

    norm = matplotlib.colors.BoundaryNorm([-1,-0.5,0.5,1], cmap.N)
    return norm


def show_grid(ax):
    ax.grid(True,color="black",alpha=0.3)
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(1))


def hide_ticks(ax, interval):

    for i,label in enumerate(ax.xaxis.get_ticklabels()[:]):
        if (i%interval)==0:
            continue
        else:
            label.set_visible(False)
    for i, label in enumerate(ax.yaxis.get_ticklabels()[:]):
        if (i % interval) == 0:
            continue
        else:
            label.set_visible(False)


def plot_line_stresses(mask,coords,n_stress,shear_stress):
    fig, (ax1, ax2) = plt.subplots(1, 2)
    plt.subplots_adjust(right=0.8)
    ax3 = fig.add_axes([0, 0, 0, 0])
    n_b = list(ax2.get_position().bounds)
    n_b[0] += n_b[2] + 0.02
    n_b[2] = 0.05
    ax3 = fig.add_axes(n_b)

    im1=np.zeros(mask.shape)+np.nan
    im2=np.zeros(mask.shape)+np.nan
    im1[coords[1:-1,0],coords[1:-1,1]]=n_stress
    im2[coords[1:-1, 0], coords[1:-1, 1]]=shear_stress

    min_v = np.nanmin(np.minimum(n_stress, shear_stress))
    max_v = np.nanmax(np.maximum(n_stress, shear_stress))

    ax1.imshow(im1,cmap="jet",vmin=min_v,vmax=max_v,origin="lower")
    ax1.set_title("noraml stress on line")
    ax2.imshow(im2,cmap="jet",vmin=min_v,vmax=max_v,origin="lower")
    ax2.set_title("shear stress on line")



    norm = matplotlib.colors.Normalize(vmin=min_v, vmax=max_v)
    cb1 = matplotlib.colorbar.ColorbarBase(ax3, cmap=matplotlib.cm.get_cmap("jet"),
                                           norm=norm, label="stress",
                                          orientation='vertical')



def plot_continous_boundary_stresses(shape,edge_lines,lines_interpol,min_v,max_v,mask_boundaries=None,plot_t_vecs=False,plot_n_arrows=False,figsize=(10,7),
                                     scale_ratio=0.2,border_arrow_filter=1,cbar_str="line stress in N/µm",vmin=None,vmax=None,
                                    cbar_width="2%",cbar_height="50%",cbar_axes_fraction=0.2,cbar_tick_label_size=20,background_color="white",cbar_borderpad=2.5,linewidth=4,cmap="jet",cbar_style="clickpoints",**kwargs):


    '''
    plotting the line stresses (total transmitted force of cell boundaries), colored by their absolute values
    as continous lines.
    :param shape:
    :param edge_lines:
    :param lines_interpol:
    :param min_v:
    :param max_v:
    :param mask_boundaries:
    :param plot_t_vecs:
    :param plot_n_arrows:
    :param figsize:
    :param scale_ratio:
    :param arrow_filter:
    :param cbar_str:
    :param vmin:  overwrites max_v and min_v if provided
    :param vmax:  overwrites max_v and min_v if provided
    :return:
    '''
    if isinstance(vmin,(float,int)):
        min_v=vmin
    if isinstance(vmax,(float,int)):
        max_v=vmax
    print("plotting cell border stresses")

    fig = plt.figure(figsize=figsize)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    fig.add_axes(ax)
    ax.set_axis_off()

    if not isinstance(mask_boundaries,np.ndarray):
        mask_boundaries=np.zeros(shape)

    # plotting empty array with background as cmap value of zero
    if background_color=="cmap_0":
        cmap_custom = matplotlib.colors.ListedColormap([matplotlib.cm.get_cmap(cmap)(i) for i in range(3)])
    else:
        cmap_custom = matplotlib.colors.ListedColormap([background_color for i in range(3)])
    bounds = [0,1,2]
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap_custom.N)
    im = ax.imshow(np.zeros(shape).astype(int),cmap=cmap_custom,norm=norm)

    # finding global scaling factor for arrows
    all_t_vecs=np.vstack([subdict["t_vecs"] for subdict in lines_interpol.values()])
    if plot_t_vecs:
        scale = scale_for_quiver(all_t_vecs[:, 0], all_t_vecs[:, 1], dims=mask_boundaries.shape,
                                                      scale_ratio=scale_ratio,return_scale=True)

    for line_id,interp in tqdm(lines_interpol.items(),total=len(lines_interpol.values())):
        p_new=interp["points_new"]
        x_new=p_new[:,0]
        y_new=p_new[:,1]
        t_norm=interp["t_norm"]
        t_vecs = interp["t_vecs"]
        n_vecs=interp ["n_vecs"]
        # plotting line segments

        c = matplotlib.cm.get_cmap(cmap)((t_norm - min_v) / (max_v - min_v))  # normalization and creating a color range
        ## see how well that works
        if line_id in edge_lines: # plot lines at the edge
            for i in range(len(x_new)-1):
                plt.plot([x_new[i],x_new[i+1]],[y_new[i],y_new[i+1]],color="gray",linewidth=linewidth)
        else:
            for i in range(len(x_new) - 1):
                plt.plot([x_new[i], x_new[i + 1]], [y_new[i], y_new[i + 1]], color=c[i], linewidth=linewidth)

        # plotting stressvectors
        if plot_t_vecs:
            t_vecs_scale=t_vecs*scale
            for i,(xn,yn,t) in enumerate(zip(x_new,y_new,t_vecs_scale)):
                if i % border_arrow_filter == 0:
                    plt.arrow(xn,  yn,t[0], t[1],head_width=0.5)
        # plotting normal vectors
        if plot_n_arrows:
            for i in range(len(x_new) - 1):
                if i % border_arrow_filter == 0:
                    plt.arrow(x_new[i],  y_new[i],n_vecs[i][0], n_vecs[i][1],head_width=0.5)



    norm = matplotlib.colors.Normalize(vmin=min_v, vmax=max_v)
    if cbar_style == "clickpoints":
        cbaxes = inset_axes(ax, width=cbar_width, height=cbar_height, loc=5, borderpad=cbar_borderpad)
        cbaxes.set_title(cbar_str, color="white")
        cbaxes.tick_params(colors="white",labelsize=cbar_tick_label_size)
        cb1 = matplotlib.colorbar.ColorbarBase(cbaxes, cmap=matplotlib.cm.get_cmap(cmap),
                                               # overrides previours colorbar
                                               norm=norm,
                                               orientation='vertical')
    else:
        cb0=plt.colorbar(im, aspect=20, shrink=0.8,fraction=cbar_axes_fraction) # just exploiting the axis generation by a plt.colorbar
        cb0.outline.set_visible(False)
        cb0.ax.tick_params(labelsize=cbar_tick_label_size)
        cb1 = matplotlib.colorbar.ColorbarBase(cb0.ax, cmap=matplotlib.cm.get_cmap(cmap),
                                               # overrides previours colorbar
                                               norm=norm,
                                               orientation='vertical')
        cb1.ax.tick_params(labelsize=cbar_tick_label_size)
        cb1.ax.set_title(cbar_str, color="black")
    return fig


def check_order(mask,coords):
    plt.figure()
    plt.imshow(mask,origin="lower")
    plt.plot(coords[:,1],coords[:,0],"o")
    for i,(x1,y1) in enumerate(zip(coords[:,0],coords[:,1])):
        plt.text(y1,x1,str(i))
# plotting:
# Pre-processing
def plot_line_stresses(mask,coords,n_stress,shear_stress):
    fig, (ax1, ax2) = plt.subplots(1, 2)
    plt.subplots_adjust(right=0.8)
    ax3 = fig.add_axes([0, 0, 0, 0])
    n_b = list(ax2.get_position().bounds)
    n_b[0] += n_b[2] + 0.02
    n_b[2] = 0.05
    ax3 = fig.add_axes(n_b)

    im1=np.zeros(mask.shape)+np.nan
    im2=np.zeros(mask.shape)+np.nan
    im1[coords[1:-1,0],coords[1:-1,1]]=n_stress
    im2[coords[1:-1, 0], coords[1:-1, 1]]=shear_stress

    min_v = np.nanmin(np.minimum(n_stress, shear_stress))
    max_v = np.nanmax(np.maximum(n_stress, shear_stress))

    ax1.imshow(im1,cmap="jet",vmin=min_v,vmax=max_v,origin="lower")
    ax1.set_title("noraml stress on line")
    ax2.imshow(im2,cmap="jet",vmin=min_v,vmax=max_v,origin="lower")
    ax2.set_title("shear stress on line")



    norm = matplotlib.colors.Normalize(vmin=min_v, vmax=max_v)
    cb1 = matplotlib.colorbar.ColorbarBase(ax3, cmap=matplotlib.cm.get_cmap("jet"),
                                           norm=norm, label="stress",
                                           orientation='vertical')

def check_normal_vectors(mask,coords,n):

    plt.figure()
    plt.imshow(mask,origin="lower",cmap="viridis")
    plt.plot(coords[:,1],coords[:,0],"o")
    for i,(x1,y1) in enumerate(zip(coords[:,1],coords[:,0])):
        plt.text(x1,y1,str(i))
    for vec,p in zip(n,coords):
        plt.arrow(p[1],p[0],vec[0],vec[1],head_width=0.2)




def check_normal_vectors_graph(mask,n,points):
    '''
    n is dictionary point_id:normal vector
    :param mask:
    :param n:
    :param points:
    :return:
    '''
    plt.figure()
    plt.imshow(mask,origin="lower",cmap="viridis")

    for p,vec in n.items():
        plt.arrow(points[p][1],points[p][0],vec[0],vec[1],head_width=0.2)

def check_normal_vectors_array(mask,n_array,origin="lower"):
    '''
    n is a 3d array.axis 0 and 1 represetn xy coordinates and axis 2 is the normal vector
    :param mask:
    :param n:
    :param points:
    :return:
    '''
    plt.figure()
    plt.imshow(mask,origin=origin,cmap="viridis")
    coords=np.where(~np.logical_and(n_array[:,:,0]==0,n_array[:,:,1]==0))

    for px,py in zip(coords[0],coords[1]):
        plt.arrow(py,px,n_array[px,py,0],n_array[px,py,1],head_width=0.2)

def plot_stress_vectors(mask,arrows_array,origin="lower"):
    '''
    n is a 3d array.axis 0 and 1 represetn xy coordinates and axis 2 is the normal vector
    :param mask:
    :param n:
    :param points:
    :return:
    '''
    mask=mask.astype("bool")
    plt.figure()
    abs_im=np.linalg.norm(arrows_array,axis=2)
    abs_im[~mask]=np.nan
    im=plt.imshow(abs_im, origin=origin, cmap="viridis")
    plt.colorbar(im)

    coords=np.where(~np.logical_and(arrows_array[:,:,0]==0,arrows_array[:,:,1]==0))
    for px,py in zip(coords[0],coords[1]):
        plt.arrow(py,px,arrows_array[px,py,0],arrows_array[px,py,1],head_width=0.2)

def plot_nodes(nodes):
    plt.figure(figsize=(7,7))
    plt.xlim(0,np.sqrt(len(nodes)))
    plt.ylim(0, np.sqrt(len(nodes)))
    for n in nodes:
        # plt.scatter(n[1],n[2],s=100,facecolors="none",edgecolors="black")
        plt.text(n[1], n[2], str(int(n[0])), bbox={"boxstyle": "circle", "edgecolor": "black", "facecolor": "none"},
                 color="black")

def plot_grid(nodes,elements,inverted_axis=False,symbol_size=4,arrows=False,image=0): # a connecting lines


    plt.figure(figsize=(7,7))

    if inverted_axis:
        nodes[:,[1,2]]=nodes[:,[2,1]]
        #plt.xlim(np.min(nodes[:, 1]), np.max(nodes[:, 1]))
        #plt.ylim(np.max(nodes[:, 1]), np.min(nodes[:, 2]))

    plt.xlim(np.min(nodes[:, 1])-10, np.max(nodes[:, 1])+10)
    plt.ylim(np.min(nodes[:, 1])-10, np.max(nodes[:, 2])+10)
    if inverted_axis:
        plt.gca().invert_yaxis()

    for n in nodes:
        # plt.scatter(n[1],n[2],s=100,facecolors="none",edgecolors="black")
        plt.text(n[1], n[2], str(int(n[0])),fontsize=symbol_size,alpha=1,ha="center",va="center", bbox={"boxstyle": "circle", "edgecolor": "black", "facecolor": "white","alpha":0.7},
                 color="black")
    coms_x=(nodes[elements[:,3],1]+nodes[elements[:,4],1]+nodes[elements[:,5],1]+nodes[elements[:,6],1])/4
    coms_y=(nodes[elements[:,3],2]+nodes[elements[:,4],2]+nodes[elements[:,5],2]+nodes[elements[:,6],2])/4
    for i,com in enumerate(zip(coms_x,coms_y)):
        plt.text(com[0], com[1], str(i),fontsize=symbol_size,ha="center",va="center", bbox={ "boxstyle": "circle", "edgecolor": "blue", "facecolor": "white","alpha":0.7},
                 color="blue")
    if arrows:
        cols=["black","red","green","blue"]
        for e in elements:
            if e[0]%10==0:
                col=cols[e[0]%4]
                plt.arrow(nodes[e[3]][1], nodes[e[3]][2], nodes[e[4]][1]-nodes[e[3]][1], nodes[e[4]][2]-nodes[e[3]][2], color=col,head_width=0.2)
                plt.arrow(nodes[e[4]][1], nodes[e[4]][2], nodes[e[5]][1]-nodes[e[4]][1], nodes[e[5]][2]-nodes[e[4]][2], color=col,head_width=0.2)
                plt.arrow(nodes[e[5]][1], nodes[e[5]][2], nodes[e[6]][1]-nodes[e[5]][1], nodes[e[6]][2]-nodes[e[5]][2], color=col,head_width=0.2)
                plt.arrow(nodes[e[6]][1], nodes[e[6]][2], nodes[e[3]][1]-nodes[e[6]][1], nodes[e[3]][2]-nodes[e[6]][2], color=col,head_width=0.2)

    else:
        for e in elements:
            plt.plot([nodes[e[3]][1],nodes[e[4]][1]],[nodes[e[3]][2],nodes[e[4]][2]],color="black")
            plt.plot([nodes[e[4]][1], nodes[e[5]][1]], [nodes[e[4]][2], nodes[e[5]][2]],color="black")
            plt.plot([nodes[e[5]][1], nodes[e[6]][1]], [nodes[e[5]][2], nodes[e[6]][2]],color="black")
            plt.plot([nodes[e[6]][1], nodes[e[3]][1]], [nodes[e[6]][2], nodes[e[3]][2]],color="black")






def show_quiver(fx,fy,filter=False,scale_ratio=0.2,headwidth=3,headlength=3,width=0.002,cmap="rainbow"):
    fx=fx.astype("float64")
    fy=fy.astype("float64")
    dims=fx.shape# save dims for use in scaling, otherwise porblems, because filtering will return flatten array

    fig=plt.figure()
    im = plt.imshow(np.sqrt(fx ** 2 + fy ** 2),cmap=cmap)



    if isinstance(filter,list):
        fx,fy,xs,ys=filter_values(fx,fy,abs_filter=filter[0],f_dist=filter[1])
        if scale_ratio:#
            fx, fy=scale_for_quiver(fx,fy, dims, scale_ratio=scale_ratio)
            plt.quiver(xs, ys, fx, fy, scale_units="xy",scale=1, angles="xy",headwidth=headwidth,headlength=headlength,width=width)
        else:
            plt.quiver(xs, ys, fx, fy, scale_units="xy", angles="xy",headwidth=headwidth,headlength=headlength,width=width)
    #elif np.isnan(np.sum(fx)): ## fastest way to check if any nan in this array (or so they say)
     #   ys,xs=np.where(~np.isnan(fx))
    #    plt.quiver(xs, ys, fx[~np.isnan(fx)], fy[~np.isnan(fx)], scale_units="xy", angles="xy")
    else:
        if scale_ratio:  #
            fx, fy = scale_for_quiver(fx, fy, dims, scale_ratio=scale_ratio)
            plt.quiver(fx, fy, scale_units="xy",scale=1, angles="xy",headwidth=headwidth,headlength=headlength,width=width)
        else:
            plt.quiver(fx, fy, scale_units="xy", angles="xy",headwidth=headwidth,headlength=headlength,width=width)




    plt.colorbar(im)
    return fig

def show_quiver_clickpoints(fx,fy,filter=[0,1],scale_ratio=0.2,headwidth=3,headlength=3,width=0.002,figsize=(6.4, 4.8),cbar_str=""
                            ,cmap="rainbow",vmin=None,vmax=None,cbar_axes_fraction=0.2,cbar_tick_label_size=15,background_color="cmap_0",scale=None,cbar_width="2%",cbar_height="50%",cbar_borderpad=2.5,cbar_style="clickpoints",**kwargs):



    fx=fx.astype("float64")
    fy=fy.astype("float64")
    dims=fx.shape# save dims for use in scaling, otherwise porblems, because filtering will return flatten array

    fig=plt.figure(figsize=figsize)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    fig.add_axes(ax)
    ax.set_axis_off()
    map_values=np.sqrt(fx ** 2 + fy ** 2)
    if vmax and not vmin:
        vmin=vmax-1 if np.min(map_values)>vmax else None # other wise the program behaves unexpectedly if the automatically set vmin is smaller the vmax
    im = ax.imshow(map_values,cmap=cmap,vmin=vmin,vmax=vmax)

    # filtering out arrows
    fx_f,fy_f,xs,ys=filter_values(fx,fy,abs_filter=filter[0],f_dist=filter[1])

    # scaling arrows to be reach a certain fraktion of the length of the longer image axis
    if scale_ratio:
        fx_f, fy_f=scale_for_quiver(fx_f,fy_f, dims, scale_ratio=scale_ratio)
        scale=1
    # plotting the arrows
    plt.quiver(xs, ys, fx_f, fy_f, scale_units="xy",scale=1, angles="xy",headwidth=headwidth,headlength=headlength,width=width)
    if cbar_style=="clickpoints":
        cbaxes = inset_axes(ax, width=cbar_width, height=cbar_height, loc=5,borderpad=cbar_borderpad)
        cbaxes.set_title(cbar_str,color="white")
        cbaxes.tick_params(colors="white",labelsize=cbar_tick_label_size)
        plt.colorbar(im, cax=cbaxes)
    else:
        cb=plt.colorbar(im,aspect=20,shrink=0.8,fraction=cbar_axes_fraction)
        cb.ax.tick_params(labelsize=cbar_tick_label_size)

    return fig

def show_edgeline(mask, ax, color="#696969", alpha=0.5 , n=6,plot_inner_line=False):
    '''
    colors a region close to the edge of mask.
    :param values: boolean mask
    :param ax: matplotlib axis object
    :param color: color of the edge coloring
    :param alpha: imshow alpha
    :param n: pixels distance from the ask edge
    :return:
    '''
    colony_area = mask != 0
    edge_area = np.logical_and(colony_area, ~binary_erosion(colony_area, iterations=n))
    edge_show = np.zeros(edge_area.shape)
    edge_show.fill(np.nan)
    edge_show[edge_area] = 1



    cmap_custom1 = matplotlib.colors.ListedColormap([color for i in range(3)])
    bounds = [0, 1, 2]
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap_custom1.N)
    ax.imshow(edge_show, cmap=cmap_custom1, norm=norm, alpha=alpha)

    if plot_inner_line:
        edge_line_inner=np.logical_and(~binary_erosion(colony_area, iterations=n), binary_erosion(colony_area, iterations=n-1))
        edge_inner_show = np.zeros(edge_line_inner.shape)
        edge_inner_show.fill(np.nan)
        edge_inner_show[edge_line_inner] = 1

        cmap_custom2 = matplotlib.colors.ListedColormap(["black" for i in range(3)])
        bounds = [0, 1, 2]
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap_custom2.N)
        ax.imshow(edge_inner_show, cmap=cmap_custom2, norm=norm)



def show_map_clickpoints(values,figsize=(6.4, 4.8),cbar_str=""
                            ,cmap="rainbow",vmin=None,vmax=None,background_color = "cmap_0",cbar_width="2%",cbar_height="50%",
                         cbar_borderpad=2.5,cbar_tick_label_size=15,cbar_axes_fraction=0.2,cbar_style="clickpoints",**kwargs):

    values=values.astype("float64")
    dims=values.shape# save dims for use in scaling, otherwise porblems, because filtering will return flatten array
    values_show = np.zeros(dims)
    values_show .fill(np.nan)
    values_show[values!=0]=values[values!=0]
    fig=plt.figure(figsize=figsize)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    fig.add_axes(ax)
    ax.set_axis_off()

    if vmax and not vmin:
        # other wise the program behaves unexpectedly if the automatically set vmin is smaller the vmax
        vmin = vmax - 1 if np.min(values) > vmax else None

    if background_color == "cmap_0":
        cmap_custom = matplotlib.colors.ListedColormap([matplotlib.cm.get_cmap(cmap)(i) for i in range(3)])
    else:
        cmap_custom = matplotlib.colors.ListedColormap([background_color for i in range(3)])
    bounds = [0, 1, 2]
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap_custom.N)
    im_bg = ax.imshow(np.zeros(dims).astype(int), cmap=cmap_custom, norm=norm)

    im = ax.imshow(values_show,cmap=cmap,vmin=vmin,vmax=vmax)
    if cbar_style == "clickpoints":
        cbaxes = inset_axes(ax, width=cbar_width, height=cbar_height, loc=5, borderpad=cbar_borderpad)
        cbaxes.set_title(cbar_str, color="white")
        cbaxes.tick_params(colors="white",labelsize=cbar_tick_label_size)
        plt.colorbar(im, cax=cbaxes)
    else:
        cb = plt.colorbar(im, aspect=20, shrink=0.8,fraction=cbar_axes_fraction)
        cb.ax.tick_params(labelsize=cbar_tick_label_size)

    # showing the border edge
    #mask = values != 0
    #show_edgeline(mask, ax, color="#696969", alpha=0.7, n=6)

    return fig


def show_quiver_ax(ax,fx,fy,vmin,vmax,filter=False,scale_ratio=0.2):
    '''
    like show quiver but returns manioulates and returns  a pyplot axes object
    :param fx:
    :param fy:
    :param filter:
    :param scale_ratio:
    :return:
    '''
    fx=fx.astype("float64")
    fy=fy.astype("float64")
    dims=fx.shape# save dims for use in scaling, otherwise porblems, because filtering will return flatten array


    im = ax.imshow(np.sqrt(fx ** 2 + fy ** 2),vmin=vmin,vmax=vmax)



    if isinstance(filter,list):
        fx,fy,xs,ys=filter_values(fx,fy,abs_filter=filter[0],f_dist=filter[1])
        if scale_ratio:#
            fx, fy=scale_for_quiver(fx,fy, dims, scale_ratio=scale_ratio)
            ax.quiver(xs, ys, fx, fy, scale_units="xy",scale=1, angles="xy",headwidth=6)
        else:
            ax.quiver(xs, ys, fx, fy, scale_units="xy", angles="xy")
    #elif np.isnan(np.sum(fx)): ## fastest way to check if any nan in this array (or so they say)
     #   ys,xs=np.where(~np.isnan(fx))
    #    plt.quiver(xs, ys, fx[~np.isnan(fx)], fy[~np.isnan(fx)], scale_units="xy", angles="xy")
    else:
        if scale_ratio:  #
            fx, fy = scale_for_quiver(fx, fy, dims, scale_ratio=scale_ratio)
            ax.quiver(fx, fy, scale_units="xy",scale=1, angles="xy")
        else:
            ax.quiver(fx, fy, scale_units="xy", angles="xy")




    #plt.colorbar(im)
    return ax


def show_quiver_mask(fx,fy,mask,filter=False):

    mask=mask.astype(bool)
    fx=fx.astype("float64")
    fy=fy.astype("float64")
    fx_show = np.zeros(fx.shape) +np.nan
    fy_show = np.zeros(fy.shape) +np.nan
    fx_show[mask] = fx[mask]
    fy_show[mask] = fy[mask]

    # coordinates for quiver
    pixx = np.arange(np.shape(fx)[0])
    pixy = np.arange(np.shape(fy)[1])
    xv, yv = np.meshgrid(pixy, pixx)

    fig=plt.figure()
    im = plt.imshow(np.sqrt(fx_show ** 2 + fy_show ** 2))
    if isinstance(filter,list):
        fx,fy,xs,ys=filter_values(fx_show,fy_show,abs_filter=filter[0],f_dist=filter[1])
        plt.quiver(xs, ys, fx, fy, scale_units="xy", angles="xy")
    else:
        plt.quiver(xv[mask],yv[mask],fx_show[mask].flatten(),fy_show[mask].flatten(), scale_units="xy", angles="xy")


    plt.show()
    plt.colorbar(im)
    return fig



def plot_map(ar1,cbar_str="",origin="upper",title="",mask=0,v_range=0,mask_overlay=0):
    '''
    imshow of an array with a color bar
    :return:
    '''

    ar=copy.deepcopy(ar1)
    if isinstance(mask,np.ndarray):
        mask = mask.astype(bool)
        ar[~mask]=np.nan

    fig, ax1 = plt.subplots(1, 1)
    axs = [ax1]
    plt.subplots_adjust(right=0.8)
    ax_cbar = fig.add_axes([0, 0, 0, 0])
    n_b = list(axs[-1].get_position().bounds)
    n_b[0] += n_b[2] + 0.02
    n_b[2] = 0.03
    n_b[1] += 0.2
    n_b[3] -= 0.4

    ax_cbar = fig.add_axes(n_b)
    if isinstance(v_range,tuple):
        min_v=v_range[0]
        max_v=v_range[1]
    else:
        min_v=  np.nanmin(ar)
        max_v=  np.nanmax(ar)



    axs[0].imshow(ar, origin=origin, cmap="jet", vmin=min_v, vmax=max_v)
    axs[0].set_title(title)

    if isinstance(mask_overlay, np.ndarray): # showing overlay mask
        mask_overlay=mask_overlay.astype("float") # cnat set to np.nan otherwise
        mask_overlay[mask_overlay==0]=np.nan
        axs[0].imshow(mask_overlay,alpha=0.75)
    # plotting one color bar
    norm = matplotlib.colors.Normalize(vmin=min_v, vmax=max_v)
    # norm=make_discrete_colorbar()
    cb1 = matplotlib.colorbar.ColorbarBase(ax_cbar, cmap=matplotlib.cm.get_cmap("jet"),
                                           norm=norm, label=cbar_str,
                                           orientation='vertical')
    format=ScalarFormatter()
    format.set_powerlimits((0,0)) # setting scientific notation for color bar ## disable at some point
    cb1 = matplotlib.colorbar.ColorbarBase(ax_cbar, cmap=matplotlib.cm.get_cmap("jet"),
                                           norm=norm, label=cbar_str,format=format,
                                           orientation='vertical')


    return fig




def plot_fields(nodes,fields=[],titles=[],cbar_str=[],grid_lines=False,grid_lines_intervall=10,dims=None,origin="lower",mask=0,mask_overlay=0):
    if not dims:
        l=np.sqrt(len(nodes)).astype(int)
        fs=[np.zeros((l,l)) for i in range (len(fields))] #note []*x produces x copies of the object with same reference
    else:
        fs = [np.zeros(dims) for i in range(len(fields))]  # note []*x produces x copies of the object with same reference

    for i in range(len(fs)):
        fs[i][nodes[:, 2].astype(int), nodes[:, 1].astype(int)]=fields[i]  # inverted because of imshow

        #print(fs[i][nodes[:, 2].astype(int), nodes[:, 1].astype(int)])

    if isinstance(mask,np.ndarray):  # filling non displayed values with nan
        for i,f in enumerate(fs):
            fs[i][~mask]=np.nan

    if isinstance(mask_overlay, np.ndarray):  # preparing overlay mask
        mask_overlay = mask_overlay.astype("float")  # cant set to np.nan otherwise
        mask_overlay[mask_overlay == 0] = np.nan

    # setting plot subplots
    if len(fs)==1:
        fig, ax1= plt.subplots(1, 1)
        axs=[ax1]

    if len(fs) == 2:
        fig, axs = plt.subplots(1, 2)

    if len(fs) == 3:
        fig, axs = plt.subplots(1, 3)

    if len(fs) == 4:
        fig, axs = plt.subplots(2,2)
        axs=axs.flatten()

    plt.subplots_adjust(right=0.8)
    ax_cbar = fig.add_axes([0, 0, 0, 0])
    n_b = list(axs[-1].get_position().bounds)
    n_b[0] += n_b[2] + 0.02
    n_b[2] = 0.03
    n_b[1]+=0.2
    n_b[3] -= 0.4

    ax_cbar = fig.add_axes(n_b)
    plt.grid(grid_lines)

    min_v=  np.nanmin(np.nanmin(np.dstack(fs),axis=2))
    max_v=  np.nanmax(np.nanmax(np.dstack(fs),axis=2))




    # imshowing every array
    for a,f,t in zip(axs,fs,titles):
        a.imshow(f,origin=origin, cmap="jet",vmin=min_v, vmax=max_v)
        a.set_title(t)
        if grid_lines:
            show_grid(a)
            hide_ticks(a, interval=grid_lines_intervall)
        if isinstance(mask_overlay, np.ndarray):
            a.imshow(mask_overlay,alpha=0.75)
    # plotting one color bar
    norm = matplotlib.colors.Normalize(vmin=min_v, vmax=max_v)
   # norm=make_discrete_colorbar()
    cb1 = matplotlib.colorbar.ColorbarBase(ax_cbar, cmap=matplotlib.cm.get_cmap("jet"),
                                           norm=norm, label=cbar_str,
                                          )
    return fig


def correct_torque(fx,fy,mask_area):
    com = regionprops(mask_area.astype(int))[0].centroid  # finding center of mass
    com = (com[1], com[0])  # as x y coorinate

    c_x, c_y = np.meshgrid(range(fx.shape[1]), range(fx.shape[0]))  # arrays with all x and y coordinates
    r = np.zeros((fx.shape[0], fx.shape[1], 2))  # array with all positional vectors
    r[:, :, 0] = c_x  # note maybe its also enough to chose any point as refernece point
    r[:, :, 1] = c_y
    r = r - np.array(com)

    f = np.zeros((fx.shape[0], fx.shape[1], 2),dtype="float64")  # array with all force vectors
    f[:, :, 0] = fx
    f[:, :, 1] = fy
    q = np.zeros((fx.shape[0], fx.shape[1], 2),dtype="float64")  # rotated positional vectors

    def get_torque_angle(p):
        q[:, :, 0] = + np.cos(p) * (f[:, :, 0]) - np.sin(p) * (f[:, :, 1])  # whats the mathematics behind this??
        q[:, :, 1] = + np.sin(p) * (f[:, :, 0]) + np.cos(p) * (f[:, :, 1])
        torque = np.abs(np.nansum(np.cross(r, q, axisa=2, axisb=2))) ## using nna sum to only look at force values in mask
        return torque.astype("float64")

    # plotting torque angle relation ship
    #ps=np.arange(-np.pi/2,np.pi/2,0.01)
    #torques=[get_torque_angle(p)*1000 for p in ps]
    #plt.figure()
    #ticks=np.arange(-np.pi/2,np.pi/2+np.pi/6,np.pi/6)
    #tick_labels=[r"$-\frac{\pi}{2}$",r"$-\frac{\pi}{3}$",r"$-\frac{\pi}{6}$",r"$0$",r"$\frac{\pi}{6}$",r"$\frac{\pi}{3}$",r"$\frac{\pi}{2}$"]
    #plt.xticks(ticks,tick_labels,fontsize=25)
    #plt.yticks(fontsize=15)
    #plt.plot(ps,torques,linewidth=6)
    #plt.gca().spines['bottom'].set_color('black')
    #plt.gca().spines['left'].set_color('black')
    #plt.gca().tick_params(axis='x', colors='black')
    #plt.gca().tick_params(axis='y', colors='black')
    #plt.savefig("/home/user/Desktop/results/thesis/figures/torque_angle.png")



    pstart = 0
    #bounds = ([-np.pi], [np.pi])
    ## just use normal gradient descent??
    eps=np.finfo(float).eps # minimum machine tolerance, for most exact calculation
    p = least_squares(fun=get_torque_angle, x0=pstart, method="lm",
                      max_nfev=100000000, xtol=eps, ftol=eps,gtol=eps, args=())["x"]  # trust region algorithm,

    q[:, :, 0] = + np.cos(p) * (f[:, :, 0]) - np.sin(p) * (f[:, :, 1])  # corrected forces
    q[:, :, 1] = + np.sin(p) * (f[:, :, 0]) + np.cos(p) * (f[:, :, 1])

    return q[:, :, 0],q[:, :, 1], p  # returns corrected forces and rotation angle




def get_torque1(fx,fy,mask_area,return_map=False):
    com = regionprops(mask_area.astype(int))[0].centroid  # finding center of mass
    com = (com[1], com[0])  # as x y coorinate

    c_x, c_y = np.meshgrid(range(fx.shape[1]), range(fx.shape[0]))  # arrays with all x and y coordinates
    r = np.zeros((fx.shape[0], fx.shape[1], 2))  # array with all positional vectors
    r[:, :, 0] = c_x  # note maybe its also enough to chose any point as refernece point
    r[:, :, 1] = c_y
    r = r - np.array(com)

    f = np.zeros((fx.shape[0], fx.shape[1], 2))  # array with all force vectors
    f[:, :, 0] = fx
    f[:, :, 1] = fy
    if return_map:
        return np.cross(r, f, axisa=2, axisb=2)
    else:
        return np.nansum(np.cross(r, f, axisa=2, axisb=2))

def check_unbalanced_forces(fx,fy,mask=None,raise_error=False):
    if not isinstance(mask, np.ndarray):
        mask=np.logical_or(fy!=0,fx!=0)
    torque=get_torque1(fx, fy, mask, return_map=False) # torque of the system
    net_force=np.array([np.sum(fx),np.sum(fy)])
    print("torque = "+ str(torque))
    print("net force = " + str(net_force))
    if raise_error and (torque==0 or np.sum(net_force==0)>0):
        raise Exception("torque or net force is not zero")



def calculate_rotation(a1, a2, mask):
    '''



    hopefully caclualtes the rotaton of a body from a rotaion filed
    :param dx:eithr array with dx values or nodes
    :param dy:
    :return:
    '''
    mask = mask.astype(bool)
    if a1.shape[1] == 5 and a1.shape[0] != a1.shape[1]:  # this would recognize a nodes array

        dx, dy = make_field(a1.astype(int), a2, mask.shape)  # an construct a field from it

    else:
        dx, dy = a1, a2

    r = np.zeros((dx.shape[0], dx.shape[1], 2))  # array with all positional vectors
    d = np.zeros((dx.shape[0], dx.shape[1], 2))


    c_x, c_y = np.meshgrid(range(mask.shape[1]), range(mask.shape[0]))
    r[:, :, 0] = c_x  # note maybe its also enough to chose any point as refernece point
    r[:, :, 1] = c_y
    com=(np.mean(r[:,:,0][mask]),np.mean(r[:,:,1][mask]))
    r = r - np.array(com)
    r[~mask] = 0
    d[:, :, 0][mask] = dx[mask].flatten()  # note maybe its also enough to chose any point as refernece point
    d[:, :, 1][mask] = dy[mask].flatten()
    return np.sum(np.cross(r, d, axisa=2, axisb=2))


def correct_rotation(def_x, def_y, mask):
    '''
    function to apply rigid body translation and rotation to a deformation field, to minimize rotation
    and translation
    :return:
    '''
    mask=mask.astype(bool)
    trans = np.array([np.mean(def_x[mask]), np.mean(def_y[mask])])  # translation
    def_xc1 = def_x - trans[0]  # correction of translation
    def_yc1 = def_y - trans[1]
    ##
    # insert method to calculate angular rotaton??
    ##
    # constructinoon positional vectors and
    r = np.zeros((def_x.shape[0], def_x.shape[1], 2))  # array with all positional vectors
    r_n = np.zeros((def_x.shape[0], def_x.shape[1], 2))
    d = np.zeros((def_x.shape[0], def_x.shape[1], 2))
    d_n = np.zeros((def_x.shape[0], def_x.shape[1], 2))
    c_x, c_y = np.meshgrid(range(def_x.shape[1]), range(def_x.shape[0]))
    r[:, :, 0] = c_x  # note maybe its also enough to chose any point as refernece point
    r[:, :, 1] = c_y
    com = (np.mean(r[:, :, 0][mask]), np.mean(r[:, :, 1][mask])) ## why inverted indices??????????
    r = r - np.array(com)
    r[~mask] = 0
    d[:,:, 0][mask]= def_xc1[mask].flatten()  # note maybe its also enough to chose any point as refernece point
    d[:, :, 1][mask] = def_yc1[mask].flatten()

    # applying rotation
    def rot_displacement(p):
        r_n[:, :, 0] = + np.cos(p) * (r[:, :, 0]) - np.sin(p) * (r[:, :, 1])  # rotation of postional vectors
        r_n[:, :, 1] = + np.sin(p) * (r[:, :, 0]) + np.cos(p) * (r[:, :, 1])
        disp = r - r_n
        return disp  # norma error measure

    # fit to corect rotation
    def displacement_error(p):
        r_n[:, :, 0] = + np.cos(p) * (r[:, :, 0]) - np.sin(p) * (r[:, :, 1])  # rotation of postional vectors
        r_n[:, :, 1] = + np.sin(p) * (r[:, :, 0]) + np.cos(p) * (r[:, :, 1])
        disp = r - r_n
        return np.sum(np.linalg.norm((d[mask] - disp[mask]),axis=1))  # norma error measure

    pstart = -1
    bounds = ([-np.pi], [np.pi])
    ## just use normal gradient descent??
    p = least_squares(fun=displacement_error, x0=pstart, bounds=bounds, method="trf",
                      max_nfev=100000000, xtol=3e-32, ftol=3e-32, gtol=3e-32, args=())["x"]  # trust region algorithm,
    ## note only works if displacement can be reached in "one rotation!!!!
    # get the part of displacement originating form a rotation of p
    d_rot = rot_displacement(p)
    d_n[mask] = d[mask] - d_rot[mask]  # correcting this part of rotation
    return d_n[:, :, 0], d_n[:, :, 1], trans,p


def make_field(nodes, values, dims):
    '''
    function to write e.g loads or deformation data to array
    :param nodes:
    :param values:
    :return:
    '''
    nodes=nodes.astype(int)
    fx = np.zeros(dims)
    fy = np.zeros(dims)
    fx[nodes[:, 2], nodes[:, 1]] = values[:, 0]
    fy[nodes[:, 2], nodes[:, 1]] = values[:, 1]
    return fx, fy


def make_solids_py_values_list(nodes, fx, fy,mask,shape=1):
    '''
    function to create a list of values, eg deformation as needed by solidspy

    :param nodes:
    :param fx:
    :param fy:
    :return:
    '''
    nodes=nodes.astype(int)
    mask=mask.astype(bool)
    if shape==1:
        l = np.zeros((len(nodes)*2))
        l[np.arange(0, len(nodes)*2) % 2 == 0] = fx[mask].flatten()  #ordering in solidspy is x,y..
        l[np.arange(0, len(nodes)*2) % 2 != 0] = fy[mask].flatten()
    else:
        l = np.zeros((len(nodes), 2))
        l[:, 0] = fx[mask][nodes[:, 2], nodes[:, 1]]
        l[:, 1] = fy[mask][nodes[:, 2], nodes[:, 1]]
    return l

def normalizing(img):
    img = img - np.nanmin(img)
    img = img / np.nanmax(img)
    img[img < 0] = 0.0
    img[img > 1] = 1.0
    return img


def get_torque2(nodes,loads):

    nodes=nodes.astype(int)

    k,l=(np.max(nodes[:,1])+1, np.max(nodes[:,2])+1)
    fx = np.zeros((k,l))
    fy = np.zeros((k,l))
    area = np.zeros((k,l),dtype=int)
    fx[nodes[:, 1], nodes[:, 2]] = loads[:, 1]
    fy[nodes[:, 1], nodes[:, 2]] = loads[:, 2]
    area[nodes[:, 1], nodes[:, 2]]=1

    com = regionprops(area)[0].centroid  # finding center of mass
    com = (com[1], com[0])  # as x y coorinate

    c_x, c_y = np.meshgrid(range(fx.shape[1]), range(fx.shape[0]))  # arrays with all x and y coordinates
    r = np.zeros((fx.shape[0], fx.shape[1], 2))  # array with all positional vectors
    r[:, :, 0] = c_x  # note maybe its also enough to chose any point as refernece point
    r[:, :, 1] = c_y
    r = r - np.array(com)

    f = np.zeros((fx.shape[0], fx.shape[1], 2))  # array with all force vectors
    f[:, :, 0] = fx
    f[:, :, 1] = fy
    torque = np.sum(np.cross(r, f, axisa=2, axisb=2))  # note order ju
    return(torque)



def plot_all_sigmas(sigma_max,sigma_min,tau_max,sigma_avg,nodes):


    s_max=np.zeros((int(np.sqrt(len(nodes))),int(np.sqrt(len(nodes)))))
    s_min=copy.deepcopy((s_max))
    t_max=copy.deepcopy((s_max))
    sig_avg=copy.deepcopy((s_max))


    s_max[nodes[:, 2].astype(int), nodes[:, 1].astype(int)] = sigma_max # should be the correct assignement
    s_min[nodes[:, 2].astype(int), nodes[:, 1].astype(int)] = sigma_min
    t_max[nodes[:, 2].astype(int), nodes[:, 1].astype(int)] = tau_max
    sig_avg[nodes[:, 2].astype(int), nodes[:, 1].astype(int)] = sigma_avg




    fig, ((ax1, ax2), (ax3,ax4)) = plt.subplots(2,2)
    plt.subplots_adjust(right=0.8)
    ax5 = fig.add_axes([0, 0, 0, 0])
    n_b = list(ax4.get_position().bounds)
    n_b[0] += n_b[2] + 0.02
    n_b[2] = 0.03
    ax5 = fig.add_axes(n_b)


    min_v=  np.min(np.min(np.dstack([s_max,s_min,t_max,sig_avg]),axis=2))
    max_v=  np.max(np.max(np.dstack([s_max,s_min,t_max,sig_avg]),axis=2))




    ax1.imshow(s_max, origin='lower', cmap="jet",vmin=min_v, vmax=max_v)  ## without any interpolation and such
    ax2.imshow(s_min, origin='lower', cmap="jet",vmin=min_v, vmax=max_v)
    ax3.imshow(t_max,origin='lower', cmap="jet",vmin=min_v, vmax=max_v)#
    ax4.imshow(sig_avg, origin='lower', cmap="jet", vmin=min_v, vmax=max_v)
    ax1.set_title("sigma_max")
    ax2.set_title("sigma_min")
    ax3.set_title("tau_max")
    ax4.set_title("sigma_avg")
    norm = matplotlib.colors.Normalize(vmin=min_v, vmax=max_v)
    cb1 = matplotlib.colorbar.ColorbarBase(ax5, cmap=matplotlib.cm.get_cmap("jet"),
                                           norm=norm, label="stress",
                                           orientation='vertical')


def construct_elements(sqr1,nodes,elements):
    l=np.zeros(len(elements))
    for i,p in enumerate(zip(sqr1[0],sqr1[1])):
        l[i]=np.where((nodes[:,1]==p[0])*(nodes[:,2]==p[1]))[0]# finding corresponding point
    return l






def vizualize_forces_on_areas(tx_resize,ty_resize,areas,mask):
    fig = plt.figure()
    mask_show=np.zeros(mask.shape)+np.nan
    plt.imshow(np.sqrt((tx_resize / 1000) ** 2 + (ty_resize / 1000) ** 2), cmap="rainbow")
    cbar = plt.colorbar()
    cbar.set_label("traktion forces in kPa")
    pixx = np.arange(mask.shape[0])
    pixy = np.arange(mask.shape[1])
    xv, yv = np.meshgrid(pixy, pixx)
    select_x = ((xv - 1) % 50) == 0
    select_y = ((yv - 1) % 50) == 0
    y=select_x[0,:].sum()
    x=select_y[:,0].sum() ## this is correct
    select = select_x * select_y
    tx_show = tx_resize[select]
    ty_show = ty_resize[select]
    x1=np.where(select)[1].reshape((x,y))
    y1=np.where(select)[0].reshape((x,y))
    ## maek beter xs...
    scale_ratio = 0.2
    scale = scale_ratio * np.max(np.shape(tx_resize )) / np.max(
        np.sqrt((tx_show  / 1000) ** 2 + (ty_show  / 1000) ** 2))  # automatic sacleing in dependace of the image size #
    plt.quiver(x1, y1, (tx_show / 1000) * scale, (ty_show / 1000) * scale, scale=1, scale_units='xy', angles="xy",width=0.002)
    for i,(key,value) in enumerate(areas.items()):
        mask_show[value[0]]=i
    plt.imshow(mask_show,alpha=0.8,cmap="magma")
    scale=0.5*10**-4
    for key,value in areas.items():
        plt.arrow(value[4][0],value[4][1],value[5][0]*scale,value[5][1]*scale,head_width=20,color="red")

def find_areas(start_line,lines,i,com_all,invert_direction=False):
    circ_line=[]
    line_ids=[]
    id=i
    line=np.array(start_line)


    # fixing direction:
    v1=line[1]-line[0]
    v2=com_all-line[1] # vector from line end to center of mass
    cross = (np.cross(v2, v1) > 0) * 2 - 1  # gets angle direction
    angle = np.arccos(np.dot(v2, v1) / (np.linalg.norm(v2) * np.linalg.norm(v1))) * cross
    direction_factor=(angle<0)*2-1 # minus one if positive angle towards center, else negative value
    if invert_direction: direction_factor*=-1 # used
    #plt.figure()
    #plt.imshow(mask)
    #plt.arrow(line[0][0], line[0][1], v1[0], v1[1], head_width=20)
    #plt.arrow(line[1][0], line[1][1],  v2[0],  v2[1], head_width=20)

    check=False

    while not check:
        circ_line.append(line)
        line_ids.append(id)

        # logical and operation to find where both coordinates of a point are close
        child_lines1=np.where(np.isclose(line[1], lines)[:, :, 0] * np.isclose(line[1], lines)[:, :, 1] )# lines at the end
        child_lines2=np.unique(child_lines1[0])
        child_lines3=lines[child_lines2,:,:]
        # reoreintating child lines to get uniform direction
        child_lines=np.array([ [l[int(~np.isclose(l[0],line[1]).all())],l[int(np.isclose(l[0],line[1]).all())]] for l in child_lines3 ])
        #finding left most line
        child_vec=child_lines[:,1,:]-child_lines[:,0,:]
        line_vec=line[1]-line[0]
        # filtering
        #angles = np.arcsin(np.abs(child_vec[:,0] * line_vec[1] - child_vec[:,1] * line_vec[0]))
        cross = (np.cross(child_vec, line_vec) > 0) * 2 - 1 # gets angle direction
        angles = np.arccos(np.dot(child_vec, line_vec) / (np.linalg.norm(child_vec,axis=1) * np.linalg.norm(line_vec))) * cross
        #####
        angles[np.isclose(-angles,np.pi)]=np.nan
        angles[np.isclose(angles, np.pi)] = np.nan
        print(angles)
        id = child_lines2[np.nanargmin(angles *  direction_factor)]
        line=child_lines[np.nanargmin(angles *  direction_factor)] ## maybe exclude zero??
        ##reoreintate line if necessary, new point goes first
        check = np.array([np.array([np.isclose(l[0],line[0]),np.isclose(l[1],line[1])]).all()  for l in circ_line]).any()  # chekcing if finisched by comparing sum of points ## doesnt work............
        #plt.arrow(line[0][0], line[0][1], line[1][0] - line[0][0], line[1][1] - line[0][1], head_width=20)
    return circ_line,line_ids
        #plt.text(line[0][0], line[0][1],str(np.round(angles[np.nanargmin(angles)],4)))



def plot_arrows(nodes,x,y,cbar_str=[],scale_ratio=0.05,grid_lines=False,dims=None,origin="lower",title="",mask=0,filter=0,overlay_mask=0):
    '''

    :param nodes:
    :param x:
    :param y:
    :param cbar_str:
    :param scale_ratio:
    :param grid_lines:
    :param dims:
    :param origin:
    :param title:
    :param mask:
    :param filter: list of pararmeters to filter dispplayed values:[minimal abs value, spacing between arrows]
    :return:
    '''

    if not dims:
        l=np.sqrt(len(nodes)).astype(int)
        fi_x=np.zeros((l,l))
        fi_y=np.zeros((l,l))
        dims=(l,l)

    else:
        fi_x = np.zeros(dims) # note []*x produces x copies of the object with same reference
        fi_y = np.zeros(dims)

    fi_x[nodes[:, 2].astype(int), nodes[:, 1].astype(int)]=x  # inverted because of imshow
    fi_y[nodes[:, 2].astype(int), nodes[:, 1].astype(int)]=y
    if isinstance(mask,np.ndarray):
        fi_x[~mask] = np.nan
        fi_y[~mask] = np.nan

    if isinstance(filter,list):
        fi_x_f,fi_y_f,xs,ys=filter_values(fi_x,fi_y,abs_filter=filter[0],f_dist=filter[1])
        fi_x_f,fi_y_f=scale_for_quiver(fi_x_f, fi_y_f, fi_y_f.shape, scale_ratio=scale_ratio)
        min_v = np.nanmin(np.sqrt(fi_x ** 2 + fi_y ** 2))
        max_v = np.nanmax(np.sqrt(fi_x ** 2 + fi_y ** 2))
    else:
        fi_x,fi_y=scale_for_quiver(fi_x,fi_y, fi_x.shape, scale_ratio=scale_ratio)
        min_v = np.nanmin(np.sqrt(fi_x ** 2 + fi_y ** 2))
        max_v = np.nanmax(np.sqrt(fi_x ** 2 + fi_y ** 2))


    fig, ax1 = plt.subplots(1, 1)
    axs = [ax1]
    plt.subplots_adjust(right=0.8)
    ax_cbar = fig.add_axes([0, 0, 0, 0])
    n_b = list(axs[-1].get_position().bounds)
    n_b[0] += n_b[2] + 0.02
    n_b[2] = 0.03
    n_b[1]+=0.2
    n_b[3] -= 0.4

    ax_cbar = fig.add_axes(n_b)
    plt.grid(grid_lines)


    axs[0].imshow(np.sqrt(fi_x ** 2 + fi_y ** 2), origin=origin, cmap="jet")
    #plt.xticks(fontsize=15)
    #plt.yticks(fontsize=15)

    if isinstance(filter, list):
        axs[0].quiver(xs,ys,fi_x_f , fi_y_f , scale_units="xy", scale=1, angles="xy")
    else:
        axs[0].quiver(fi_x, fi_y , scale_units="xy", scale=1, angles="xy")
    axs[0].set_title(title)

    # plotting additional mask if desired
    if isinstance(overlay_mask,np.ndarray):
        ov_mask_show=np.zeros(overlay_mask.shape)+np.nan
        ov_mask_show[overlay_mask]=1
        axs[0].imshow(ov_mask_show,alpha=0.5)

    # plotting one color bar
    norm = matplotlib.colors.Normalize(vmin=min_v, vmax=max_v)
   # norm=make_discrete_colorbar()
    if len(np.unique(np.sqrt(fi_x ** 2 + fi_y ** 2)))>1: # no color bar if only one absolute value
        cb1 = matplotlib.colorbar.ColorbarBase(ax_cbar, cmap=matplotlib.cm.get_cmap("jet"),
                                               norm=norm, label=cbar_str,
                                               orientation='vertical')
    else:
        ax_cbar.remove()



    return fig


def filter_values(ar1,ar2,abs_filter=0,f_dist=3):
    '''
    function to filter out values from an array for better display
    :param ar1:
    :param ar2:
    :param ar:
    :param f_dist: distance betweeen filtered values
    :return:
    '''
    #ar1_show=np.zeros(ar1.shape)+np.nan
    #ar2_show=np.zeros(ar2.shape)+np.nan
    pixx = np.arange(np.shape(ar1)[0])
    pixy = np.arange(np.shape(ar1)[1])
    xv, yv = np.meshgrid(pixy, pixx)

    def_abs = np.sqrt((ar1 ** 2 + ar2 ** 2))
    select_x = ((xv - 1) % f_dist) == 0
    select_y = ((yv - 1) % f_dist) == 0
    select_size = def_abs > abs_filter
    select = select_x * select_y * select_size
    s1 = ar1[select]
    s2 = ar2[select]
    return s1,s2,xv[select], yv[select]


def check_closet_neigbours(points1, points2, assign,mask1=None,mask2=None):
    plt.figure()
    if isinstance(mask1,np.ndarray) and isinstance(mask2,np.ndarray):
        plt.imshow(mask1.astype(int) + mask2.astype(int))

    for p1_ind, p1 in enumerate(points1):
        plt.text(p1[1], p1[0], str(p1_ind))
    for p2_ind, p1_ind in enumerate(assign):
        plt.text(points2[p2_ind][1], points2[p2_ind][0], str(p1_ind))


def get_line(start, end):
    """Bresenham's Line Algorithm
    Produces a list of tuples from start and end
    start end as tupels of points
    """
    # Setup initial conditions
    x1, y1 = start
    x2, y2 = end
    dx = x2 - x1
    dy = y2 - y1

    # Determine how steep the line is
    is_steep = abs(dy) > abs(dx)

    # Rotate line
    if is_steep:
        x1, y1 = y1, x1
        x2, y2 = y2, x2

    # Swap start and end points if necessary and store swap state
    swapped = False
    if x1 > x2:
        x1, x2 = x2, x1
        y1, y2 = y2, y1
        swapped = True

    # Recalculate differentials
    dx = x2 - x1
    dy = y2 - y1

    # Calculate error
    error = int(dx / 2.0)
    ystep = 1 if y1 < y2 else -1

    # Iterate over bounding box generating points between start and end
    y = y1
    points = []
    for x in range(x1, x2 + 1):
        coord = (y, x) if is_steep else (x, y)
        points.append(coord)
        error -= abs(dy)
        if error < 0:
            y += ystep
            error += dx

    # Reverse the list if the coordinates were swapped
    if swapped:
        points.reverse()
    return points

def scale_for_quiver(ar1,ar2,dims,scale_ratio=0.2,return_scale=False):
    scale = scale_ratio * np.max(dims) / np.nanmax(np.sqrt((ar1) ** 2 + (ar2) ** 2))
    if return_scale:
        return scale
    return ar1*scale,ar2*scale


def custom_solver(mat, rhs, mask_area,verbose=False):
    #IBC is "internal boundary condition" contains information about which nodes are fixed and
    # where the unfixed nodes can be found in the rhs vector


    """Solve a static problem [mat]{u_sol} = {rhs}

    Parameters
    ----------
    mat : array
        Array with the system of equations. It can be stored in
        dense or sparse scheme.
    rhs : array
        Array with right-hand-side of the system of equations.

    Returns
    -------
    u_sol : array
        Solution of the system of equations.

    Raises
    ------
    """

    com = regionprops(mask_area.astype(int))[0].centroid  # finding center of mass
    com = (com[1], com[0])  # as x y coorinate

    c_x, c_y = np.meshgrid(range(mask_area.shape[1]), range(mask_area.shape[0]))  # arrays with all x and y coordinates
    r = np.zeros((mask_area.shape[0], mask_area.shape[1], 2))  # array with all positional vectors
    r[:, :, 0] = c_x  # note maybe its also enough to chose any point as refernece point
    r[:, :, 1] = c_y
    r = r - np.array(com)
    len_disp = mat.shape[1]  # length of the disered displacement vector
    zero_disp_x = np.zeros(len_disp)
    zero_disp_y = np.zeros(len_disp)
    zero_torque = np.zeros(len_disp)
    # adding zero displacement condition
    zero_disp_x[::2] = 1
    zero_disp_y[1::2] = 1
    # torque=sum(r1*f2-r2*f1)
    zero_torque[::2] = r[:, :, 1][mask_area.astype(bool)]  # -r2 factor
    zero_torque[1::2] = r[:, :, 0][mask_area.astype(bool)]  # +r1 factor
    add_matrix=np.vstack([zero_disp_x,zero_disp_y,zero_torque])

    # adding zero conditions for force vector and torque
    rhs=np.append(rhs,np.zeros(3))


    if type(mat) is csr_matrix:
        import scipy.sparse
         # convert additional conditions to sparse matrix
        mat=scipy.sparse.vstack([mat,csr_matrix(add_matrix)],format="csr")
        u_sol,error = np.array(lsqr(mat, rhs,atol=10**-10, btol=10**-10,iter_lim=15000,show=verbose))[[0,3]]# sparse least squares solver
    elif type(mat) is np.ndarray:
        # adding to matrix
        mat=np.append(mat,add_matrix,axis=0)

        u_sol,error = np.array(np.linalg.lstsq(mat, rhs))[[0,1]]
    else:
        raise TypeError("Matrix should be numpy array or csr_matrix.")

    return u_sol,error




def plot_deformations(folder):
    files = os.listdir(folder)
    files_dict = defaultdict(dict)
    for f in files:
        if f[2] == "u":
            files_dict[f[:2]]["u"] = f
        if f[2] == "v":
            files_dict[f[:2]]["v"] = f
    for frame, files in files_dict.items():
        u = np.load(os.path.join(folder, files["u"]))
        v = np.load(os.path.join(folder, files["v"]))
        dpi = 200
        fig1 = show_quiver_clickpoints(u, v, filter=[0, int(np.ceil(u.shape[0] / 40))], scale_ratio=0.2,
                                       headwidth=3, headlength=3, width=0.002,
                                       figsize=(2022 / dpi, 2011 / dpi),
                                       cbar_str="deformation\n[pixel]")
        fig1.savefig(os.path.join(folder, frame + "deformation.png"), dpi=200)
