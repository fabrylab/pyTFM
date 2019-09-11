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
from scipy.ndimage.measurements import center_of_mass
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
from scipy.spatial import ConvexHull

from matplotlib.path import Path



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

    markers=db.getMarkers(id=list(points),type="marker")
    p1 = (markers[0].x,markers[0].y)
    p2 = (markers[1].x, markers[1].y)
    lines.append([p1,p2])
lines=np.array(lines)
db.db.close()

# finding cells from lines:
# strategy: go left at each intersection until your back to your line


## getting all edge points and removing edge lines as start points for algorithm





com_all=np.array([np.mean(lines[:,:,0]),np.mean(lines[:,:,1])])

#better stretegy: always turn towards center of mass of the whole thing
#edge_points=ConvexHull(np.reshape(lines,(int(len(lines)*2),2))).vertices
#edge_points=np.reshape(lines,(int(len(lines)*2),2))[edge_points]
#filter=np.zeros(np.shape(lines))
#for point in edge_points:
#    filter+=np.array(np.isclose(point,lines),dtype=int)
#filter=np.mean(np.mean(filter,axis=1),axis=1)
#
circ_lines=[]
ids=[]
for i in range(len(lines)):
    #if not filter[i]>0:
    arrea_line,line_id=find_areas(lines[i], lines,i,com_all)
    circ_lines.append(np.array(arrea_line))
    ids.append(line_id)





#filtering unique values


com=[]
for l in circ_lines:
    px=l[:,:,0]# points in x
    py = l[:, :, 1]
    com.append(np.array([np.mean(np.unique(px)),np.mean(np.unique(py))])) # center of mass of each area
com=np.array(com)
com,com_filter=np.unique(com,axis=0,return_index=True)# returns only inices of the uniq values
circ_lines=[circ_lines[i] for i in com_filter]
ids=[ids[i] for i in com_filter]

plt.figure()
plt.imshow(mask)
for line in circ_lines[3]:  ## check out this problematic region
    plt.arrow(line[0][0], line[0][1], line[1][0]- line[0][0],  line[1][1]-line[0][1], head_width=20)
    #plt.text(line[0][0], line[0][1], str([np.round(x, 4) for x in np.arctan2(child_vec[:,1],child_vec[:,0])]))
    #plt.text(line[0][0]+10, line[0][1]+10, str([np.round(x, 4) for x in [np.arctan2(line_vec[1],line_vec[0])]]))

for id,center in zip(com_filter,com):
    plt.text(center[0],center[1],str(id))


#solve with arrays only?

### getting area for each polygone


x, y = np.meshgrid(np.arange(mask.shape[1]), np.arange(mask.shape[0])) # make a canvas with coordinates
x, y = x.flatten(), y.flatten()

def find_area(line,x,y,mask):   ## rename.....
    points = np.vstack((x, y)).T
    p = Path(line[:, 1])  # make a polygon
    grid = p.contains_points(points)
    area_mask = grid.reshape(mask.shape)  #
    return area_mask


def get_point_ids(line_ids):
    point_ids=[]
    for id in line_ids:
        point_ids.append(connections_f[id])
    point_ids=np.array(point_ids)
    point_ids=np.unique(point_ids)
    point_ids.sort()
    return point_ids












areas={}
for i,(line_group,line_ids) in enumerate(zip(circ_lines,ids)):   ## optimze??
    line_ids.sort()
    area_mask=find_area(line_group,x,y,mask)
    com2 = center_of_mass(area_mask)
    com2=np.array([com2[1],com2[0]]) # changeing because of x y switch in image
    areas[i]=[area_mask,line_group,line_ids,get_point_ids(line_ids),com2]   ## also might be a bit to big...






## loading traction force data
# traction forces in these areas
# data from WT2, image 11 in TFM_cell_patches folder
pixelsize = 6.25 / 40 # µm/pixel pixelsize of the original images
window_size = 64  # windowsize for piv calculation of the deformation field
overlapp = 32 # overlapp for deformation field calculation, should be half of windowsize
pixel_factor = window_size - overlapp
pixelsize_defo = pixelsize * pixel_factor # pixelsize of the deformation field and traction force image
tx=np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/temp_data_for_scripts/tx.npy")
ty=np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/temp_data_for_scripts/ty.npy")


#### need to use interpolation to adjust coordinates
tx_resize=cv.resize(tx,dsize=(mask.shape[1],mask.shape[0]),interpolation=cv.INTER_LINEAR)
ty_resize=cv.resize(ty,dsize=(mask.shape[1],mask.shape[0]),interpolation=cv.INTER_LINEAR) ## looks good



complete_area=np.zeros(mask.shape) ## works reasonably fast....
for i in areas.keys():
    complete_area+=areas[i][0]
complete_area=complete_area>0
com_total=center_of_mass(complete_area)
T_total=np.array([np.mean(tx_resize[complete_area]),np.mean(ty_resize[complete_area])]) # correction for total Traction forces

tx_cor=np.zeros(tx_resize.shape)    # cuting out relevant part and normalize all force in cell patch to zero
tx_cor[complete_area]=tx_resize[complete_area]-T_total[0]
ty_cor=np.zeros(ty_resize.shape)
ty_cor[complete_area]=ty_resize[complete_area]-T_total[1]
check=np.sum(ty_cor) # not exactely zero ?? is thsi problematic??
# caclualting force for each area:
for area_id in areas.keys():
    T=np.array([np.sum(tx_cor[areas[area_id ][0]]),np.sum(ty_cor[areas[area_id ][0]])])
    areas[area_id].append(T)

vizualize_forces_on_areas(tx_resize,ty_resize,areas,mask)
'''
plt.figure()
plt.imshow(np.sqrt(tx_resize**2+ty_resize**2)/1000)
plt.colorbar()
complete_area_show=np.zeros(complete_area.shape)+np.nan
complete_area_show[complete_area]=1500

plt.imshow(complete_area_show,alpha=0.4,cmap="gnuplot")
for edge in graph_dijkstra.edges:
    plt.plot([points[edge.start][0],points[edge.end][0]],[points[edge.start][1],points[edge.end][1]],"--",color="grey",alpha=0.2)
k=6
for i in range(len(paths[k]) - 1):
    plt.plot([points[paths[k][i]][0], points[paths[k][i + 1]][0]], [points[paths[k][i]][1], points[paths[k][i + 1]][1]],linewidth=2)
for i in range(len(border_points) - 1):
    plt.plot([points[border_points[i]][0], points[border_points[i + 1]][0]], [points[border_points[i]][1], points[border_points[i + 1]][1]],color="green",alpha=0.6)
com1=np.zeros(2)
for area_ids in paths_dict[k][3]["side1"]:
    com1+=areas[area_ids][4]

com1 /=len(paths_dict[k][3]["side1"])
com2=np.zeros(2)
for area_ids in paths_dict[k][3]["side2"]:
    com2+=areas[area_ids][4]
com2/=len(paths_dict[k][3]["side2"])
F1=-paths_dict[k][4]
F2=-F1
scale=0.00003
plt.arrow(com1[0],com1[1],F1[0]*scale,F1[1]*scale,head_width=50,width=10,color="black")
plt.arrow(com2[0],com2[1],F2[0]*scale,F2[1]*scale,head_width=50,width=10,color="black")


plt.axis("off")
'''



import pickle


with open('/media/user/GINA1-BK/data_traktion_force_microscopy/temp_data_for_scripts/areas.pickle', 'wb') as handle:
    pickle.dump(areas, handle, protocol=pickle.HIGHEST_PROTOCOL)

### make everything to class later...
#class line:
 #   def __init__(self,id,coordinates):
#        self.id=id
#        self.points=coordinates

