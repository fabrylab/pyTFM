import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
from matplotlib.widgets import RectangleSelector
import openpiv.tools
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import numpy as np
import copy
import os
from scipy.ndimage import zoom
from skimage.filters import rank
from skimage.morphology import cube,label, remove_small_objects
from scipy.ndimage.filters import uniform_filter,median_filter,gaussian_filter
import openpiv.tools
import openpiv.process
import openpiv.scaling
import openpiv.validation
import openpiv.filters
import matplotlib.pyplot as plt
import numpy as np
import copy
import os
import time
from skimage.measure import regionprops
import cv2 as cv
from skimage.morphology import skeletonize,binary_erosion,binary_closing
import clickpoints
from imageio import imread
from skimage.filters import gaussian






from TFM_functions import *


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


sigma=0.49 #poison ratio
young=25536 # youngsmodulus in Pa
pixelsize = 6.25 / 40 # µm/pixel pixelsize of the original images
window_size = 64  # windowsize for piv calculation of the deformation field
overlapp = 32 # overlapp for deformation field calculation, should be half of windowsize
pixel_factor = window_size - overlapp
pixelsize_defo = pixelsize * pixel_factor # pixelsize of the deformation field











file3=r"/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/11bf_before_shift.tif" # bf image
bf_image=imread(file3)

tx=np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/tx.npy")
ty=np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/ty.npy")
mask=np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/mask.npy")
mask_std=np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/mask_std.npy")
u=np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/u.npy")
v=np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/v.npy")



## loading mask and converitng ## ask someone about this
db=clickpoints.DataFile("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/mask_cell_boundary3.cdb")

# not yet used
mask=db.getMask(frame=0).data
mask=binary_closing(mask)
mask=binary_erosion(mask)
mask=binary_erosion(mask)
mask=binary_erosion(mask)
mask=skeletonize(mask)
#plt.figure();plt.imshow(mask)



db.setMarkerType("new_marker",color="#00FF00")
for marker in db.getMarkers(type="marker"):
    print(marker.id)
    id=marker.id
    x=marker.x
    y=marker.y
    image=marker.image
   # db.deleteMarkers(id=id)
    db.setMarker(image=image,id=id+100,x=x,y=y,text=str(id),type="new_marker")




connections=[[2,4],[50,51],[19,20],[2,5],[21,22],[51,52],[48,37],[46,58],[58,44],[46,47],[46,49],[56,2],[51,53],
             [54,57],[57,34],[6,12],[12,11],[12,14],[13,14],[4,3],[23,24],[24,59],[59,26],[59,33],[3,56],[3,8],[4,7],[7,6],[6,5],[5,13],[13,16],
             [16,17], [16,14], [16,17], [17,18], [18,23], [23,22], [18,15],[35,34],[31,30],[39,42], [14,15], [15,19], [19,11], [11,10], [10,9], [9,7], [9,8],
             [8,52], [8,3], [54,53], [53,55], [55,48],
             [49,48], [49,50], [50,51],[58,45], [34,21], [21,33], [32,33], [32,38], [36,38], [36,35], [35,55], [21,20], [20,22],
             [17,25], [25,27], [27,28], [28,43], [42,43],[54,52],[57,10], [42,44], [44,45],[47,48],[36,37],[31,32], [45,47], [43,40], [40,29], [40,39],[30,29],[29,28],[40,39],
             [39,37], [39,41], [41,31], [30,26], [26,27], [24,25], [41,38],[50,60],[60,61],[56,61]]

connections=[sorted(x) for x in connections]

connections_f=np.unique(connections,axis=0)
np.save("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/connections.npy",connections_f)


# making lines:
lines=[]
for points in connections_f:

    markers=db.getMarkers(id=list(points))
    p1=(markers[0].x,markers[0].y)
    p2 = (markers[1].x, markers[1].y)
    lines.append([p1,p2])

db.db.close()

plt.figure()
plt.imshow(mask)
for line in lines:
    plt.plot([line[0][0],line[1][0]],[line[0][1],line[1][1]])

# finding cells from lines:
# strategy: go left at each intersection until your back to your line




def find_areas(start_line,lines):
    circ_line=[]
    lines=np.array(lines)
    line=np.array(start_line)
    check=False

    while not check:
        circ_line.append(line)
        print(line)
        child_lines1=np.where(np.isclose(line[1],lines))# lines at the end
        child_lines2=np.unique(child_lines1[0])
        child_lines3=lines[child_lines2,:,:]
        # reoreintating child lines to get uniform direction
        child_lines=np.array([ [l[np.isclose(l,line[1])],l[~np.isclose(l,line[1])]] for l in child_lines3 ])
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

        line=child_lines[np.nanargmin(angles)] ## maybe exclude zero??
        ##reoreintate line if necessary, new point goes first
        check = np.array([np.isclose(np.sum(l),np.sum(line)) for l in circ_line]).any()  # chekcing if finisched by comparing sum of points

    return circ_line
        #plt.text(line[0][0], line[0][1],str(np.round(angles[np.nanargmin(angles)],4)))
circ_lines=[]
for i in range(len(lines)):
    circ_lines.append(find_areas(lines[i],lines))

plt.figure()
plt.imshow(mask)
for line in circ_lines[4]:
    plt.arrow(line[0][0], line[0][1],  line[1][0]- line[0][0], line[1][1]-line[0][1]  , head_width=20)
    #plt.text(line[0][0], line[0][1], str([np.round(x, 4) for x in np.arctan2(child_vec[:,1],child_vec[:,0])]))
    #plt.text(line[0][0]+10, line[0][1]+10, str([np.round(x, 4) for x in [np.arctan2(line_vec[1],line_vec[0])]]))











tx=tx/1000
ty=ty/1000
u_shift = u - np.mean(u)
v_shift = v - np.mean(v)

## pre calulating energy on all points
print("calculating strain energies")
energy_points = 0.5 * pixelsize * pixelsize * (np.sqrt((tx * (u_shift / pixel_factor) * pixelsize) ** 2 + (
        ty * (v_shift / pixel_factor) * pixelsize) ** 2)) / 1000
bg = np.percentile(energy_points, 50)







### gui
custom_cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#DBDC3E","yellow"])
pixx = np.arange(np.shape(tx)[0])
pixy = np.arange(np.shape(tx)[1])
xv, yv = np.meshgrid(pixy, pixx)
pix = np.vstack((xv.flatten(), yv.flatten())).T
ind = []
energy = 0
contractile_force=0
plt_str1 = 0
plt_str2 = 0

plt.close("all")
fig = plt.figure(figsize=(12,4))
ax1 = plt.axes([0.45, 0.25, 0.3, 0.6]) # ax image bf
ax2 = plt.axes([0.05, 0.25, 0.3, 0.6]) # ax img traction forces
#ax_col=plt.axes([0.65, 0.25, 0.1, 0.6]) # ax for colorbar
#ax_col.set_axis_off()
ax_text=plt.axes([0.45, 0.1, 0.1, 0.1]) # ax for colorbar
ax_text.set_axis_off()



im = ax1.imshow(np.sqrt(tx ** 2 + ty ** 2))
#cbar=fig.colorbar(mappable=im,ax=ax_col)
#cbar.set_label('traction forces in kPa')
x_small, y_small = get_xy_for_quiver(tx)
ax1.quiver(x_small, y_small, tx, ty)


ax2.imshow(bf_image,cmap="gray")






def onselect(verts):
    global array, pix, ind, energy, plt_str1,plt_str2,ax_text,energy_points,bg,mask,verts_out,p,ax1
    global contractile_force,center,im,im_mask
    verts_out = verts

    # simple interpolation of points, to save computation time
    interpol_factors = np.array([np.shape(tx)[1] / np.shape(bf_image)[1], np.shape(tx)[0] / np.shape(bf_image)[0]])
    verts_array=np.array(verts_out)
    verts_interpol=np.zeros(np.shape( verts_array))
    for i in range(np.shape(verts_interpol)[0]):
        verts_interpol[i,:]=verts_array[i,] * interpol_factors  # element wise multiplikatioon for each row
    p = path.Path(verts_interpol)


    #retrieving mask
    ind = p.contains_points(pix, radius=1)
    print(ind)  # boolean mask
    mask = np.reshape(np.array(ind), (np.shape(tx)))

    #showing mask as overlay
    mask_show=np.zeros(np.shape(tx))+np.nan
    mask_show[mask]=1

    ax1.clear()
    im=ax1.imshow(np.sqrt(tx ** 2 + ty ** 2),alpha=0.8)
    x_small, y_small = get_xy_for_quiver(tx)
    ax1.quiver(x_small, y_small, tx, ty)
    im_mask=ax1.imshow(mask_show,alpha=0.4,cmap=custom_cmap1)

    fig.canvas.draw_idle()


lasso = LassoSelector(ax2, onselect,button=1)



p1=[]
p2=[]

def on_press(event):
    global p1,p2,point_draw1,point_draw2,line_draw1,ax1,ax2
    print(event.button,[event.xdata, event.ydata])
    if event.button==3 and event.inaxes==ax1:

        if len(p2) != 0 and len(p1) != 0:
            p1, p2 = [], []
            point_draw1[0].remove()
            point_draw2[0].remove()
            line_draw1.remove()


        if len(p2)==0 and len(p1)!=0:
            p2 = [int(event.xdata), int(event.ydata)]
            point_draw2=ax1.plot(p2[0],p2[1], "or", color="red")
            line_draw1,=ax1.plot([p1[0],p2[0]],[p1[1],p2[1]], color="green")
            fig.canvas.draw_idle()

        if len(p1)==0:
            p1=[int(event.xdata), int(event.ydata)]
            print(p1)
            point_draw1=ax1.plot(p1[0],p1[1],"or",color="red")
            fig.canvas.draw_idle()




cid = fig.canvas.mpl_connect('button_press_event',on_press)


p_set=[]
def on_press(event):
    global p1,p2,point_draw1,point_draw2,line_draw1,ax1,ax2
    print(event.button,[event.xdata, event.ydata])
    if event.button==3 and event.inaxes==ax2:

        if len(p2) != 0 and len(p1) != 0:
            p1, p2 = [], []
            point_draw1[0].remove()
            point_draw2[0].remove()
            line_draw1.remove()


        if len(p2)==0 and len(p1)!=0:
            p2 = [int(event.xdata), int(event.ydata)]
            point_draw2=ax1.plot(p2[0],p2[1], "or", color="red")
            line_draw1,=ax1.plot([p1[0],p2[0]],[p1[1],p2[1]], color="green")
            fig.canvas.draw_idle()

        if len(p1)==0:
            p1=[int(event.xdata), int(event.ydata)]
            print(p1)
            point_draw1=ax1.plot(p1[0],p1[1],"or",color="red")
            fig.canvas.draw_idle()




cid = fig.canvas.mpl_connect('button_press_event',on_press)









def caclulate_normal_force(event):
    global p1,p2,mask,ax1,ax2,tx,ty,x_small,y_small,_custom_cmap1,mask_show,mask_l,im_mask,points
    global f1,f2,b,n1,n2,mid_point

    # split mask
    points=get_line(p1,p2)
    for point in points:
        mask[point[1],point[0]]=False





    mask_l=label(mask,connectivity=1)
    areas=np.array([r.area for r in regionprops(mask_l)])
    v=np.partition(areas, -2)[-2] # finds second largest value
    mask_l=remove_small_objects(mask_l,min_size=v)

    mask_show2 = np.zeros(np.shape(tx)) + np.nan
    mask_show2[mask_l==1] = 1
    mask_show2[mask_l >1] = 2

    im_mask.remove()
    plt.figure()
    plt.imshow(mask_show2)
    plt.figure()
    plt.imshow(mask_l)
    im_mask = ax1.imshow(mask_show2, alpha=0.4, cmap=custom_cmap1)
    fig.canvas.draw_idle()

    mask1=np.zeros(np.shape(mask))
    mask1[mask_l == 1] = 1 # two parts of mask  also dont ignore
    mask2 = np.zeros(np.shape(mask))
    mask2[mask_l > 1] = 1 # two parts of mask  also dont ignore


    ## arrows of total force in one region
    sum_tx1 = np.sum(tx[mask1>0])
    sum_ty1 = np.sum(ty[mask1>0])
    sum_tx2 = np.sum(tx[mask2 > 0])
    sum_ty2 = np.sum(ty[mask2 > 0])
    com1,com2 = [(int(r.centroid[0]),int(r.centroid[1])) for r in regionprops(mask_l)]

    scale = 0.9
    ax1.arrow(com1[1], com1[0], sum_tx1 * scale, sum_ty1 * scale,head_width=5, facecolor="red", edgecolor="red",
              clip_on=False, head_starts_at_zero=False)
    ax1.arrow(com2[1], com2[0], sum_tx2 * scale, sum_ty2 * scale, head_width=5,
              facecolor="red", edgecolor="red",
              clip_on=False, head_starts_at_zero=False)

    print(com1,com2)


    # arows normal to vboundary
    b=np.array(p1)-np.array(p2)#boundaray as vector
    f1=np.array([sum_tx1, sum_ty1]) # forces as vectors
    f2=np.array( [sum_tx2, sum_ty2])

    n1=f1-(b*np.dot(b,f1)/(np.linalg.norm(b)**2)) # forces as normal vectors
    ################## this could make sense, but please confirm...########################
    n2=f2-(b*np.dot(b,f2)/(np.linalg.norm(b)**2))
    ### doesnt make much sense?? negative values





    mid_point=np.mean(np.array([p1,p2]),axis=0)
    scale = 1
    ax1.plot(mid_point[0],mid_point[1],"or",color="red")
    ax1.arrow(mid_point[0],mid_point[1],n1[0]* scale, n1[1]* scale, head_width=2, facecolor="green", edgecolor="green",
              clip_on=False, head_starts_at_zero=False)
    ax1.text(mid_point[0]+n1[0]* scale,mid_point[1]+n1[1]* scale,str(np.round(n1,2)),color="white")
    ax1.arrow(mid_point[0],mid_point[1],n2[0]* scale, n2[1]* scale,  head_width=2,
              facecolor="green", edgecolor="green",
              clip_on=False, head_starts_at_zero=False)
    ax1.text(mid_point[0] + n2[0] * scale, mid_point[1] + n2[1] * scale,str(np.round(n2,2)),color="white")
    ax1.arrow(mid_point[0], mid_point[1], b[0] * scale, b[1] * scale, head_width=2,
              facecolor="blue", edgecolor="green",
              clip_on=False, head_starts_at_zero=False)









button_force = plt.axes([0.7, 0.1, 0.1, 0.075])
button_force.xaxis.set_ticks_position('none')
button_force.yaxis.set_ticks_position('none')
button_force.set_xticks([])
button_force.set_yticks([])
b_force = Button(button_force, 'caclulate normal force')
b_force.on_clicked(caclulate_normal_force)








'''
def toggle_selector(event): ## would alow to toggle a selector
    print('Key pressed.')
    if event.key in ['Q', 'q'] and selector.active:
        print('RectangleSelector deactivated.')
        selector.set_active(False)
    if event.key in ['A', 'a'] and not selector.active:
        print('RectangleSelector activated.')
        selector.set_active(True)


fig.canvas.mpl_connect('key_press_event', toggle_selector)

'''
