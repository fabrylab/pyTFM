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
from tqdm import tqdm
from scipy.ndimage import zoom
from skimage.filters import rank
from skimage.morphology import cube,label, remove_small_objects
from scipy.ndimage.filters import uniform_filter,median_filter,gaussian_filter
import openpiv.tools
from scipy.ndimage.morphology import binary_fill_holes
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
import itertools
from skimage.morphology import skeletonize,binary_erosion,binary_closing
import clickpoints
from imageio import imread
from skimage.filters import gaussian
from TFM_functions import *
from scipy.spatial import ConvexHull
from functions_for_cell_colonie import *
from collections import deque, namedtuple
# we'll use infinity as a default distance to nodes.

def _all_simple_paths_graph(G, source, target, cutoff=None):
    if cutoff < 1:
        return
    visited = [source]
    stack = [iter(G[source])]
    while stack:
        children = stack[-1]
        child = next(children, None)
        if child is None:
            stack.pop()
            visited.pop()
        elif len(visited) < cutoff:
            if child == target:
                yield visited + [target]
            elif child not in visited:
                visited.append(child)
                stack.append(iter(G[child]))
        else: #len(visited) == cutoff:
            if child == target or target in children:
                yield visited + [target]
            stack.pop()
            visited.pop()

def check_duplication(paths):
    ps=[sorted(list(l)) for l in paths]
    u_paths,index,counts=np.unique(ps,return_counts=True,return_index=True)
    n_duplicates=np.sum(counts>1)
    return n_duplicates

connections=[[2,4],[50,51],[19,20],[2,5],[21,22],[51,52],[48,37],[46,58],[58,44],[46,47],[46,49],[56,2],[51,53],
             [54,57],[57,34],[6,12],[12,11],[12,14],[13,14],[4,3],[23,24],[24,59],[59,26],[59,33],[3,56],[3,8],[4,7],[7,6],[6,5],[5,13],[13,16],
             [16,17], [16,14], [16,17], [17,18], [18,23], [23,22], [18,15],[35,34],[31,30],[39,42], [14,15], [15,19], [19,11], [11,10], [10,9], [9,7], [9,8],
             [8,52], [8,3], [54,53], [53,55], [55,48],
             [49,48], [49,50], [50,51],[58,45], [34,21], [21,33], [32,33], [32,38], [36,38], [36,35], [35,55], [21,20], [20,22],
             [17,25], [25,27], [27,28], [28,43], [42,43],[54,52],[57,10], [42,44], [44,45],[47,48],[36,37],[31,32], [45,47], [43,40], [40,29], [40,39],[30,29],[29,28],[40,39],
             [39,37], [39,41], [41,31], [30,26], [26,27], [24,25], [41,38],[50,60],[60,61],[56,61]]
connections=[sorted(x) for x in connections]
connections=np.array(connections)
connections_f=np.unique(connections,axis=0)


# representation as graph: (could be better written like this in the begining)

graph1={i:[] for i in np.unique(connections.flatten())}

for i,j in connections:
    graph1[i].append(j)
    graph1[j].append(i)
for key,value in graph1.items():
    graph1[key]=np.unique(value)

ps=list(_all_simple_paths_graph(graph1_cent, 46, 25, cutoff=73)) ## could also be lsess cutoff??

all_paths=[]
for i,e_p in enumerate(bps_comb):
    print(i)
    all_paths.extend(list(_all_simple_paths_graph(graph1_cent, e_p[0], e_p[1], cutoff=73)))


for i in all_paths[1:20]:
    plot_path(i,mask)

check_duplication(all_paths)

paths_dict={}
Ax=np.zeros(((len(all_paths),len(central_lines_p))))

for i,pas in tqdm(enumerate(all_paths),total=len(all_paths)):   #### needs some debugging ,is force correct???
    ids=get_line_ids_from_path(pas,central_lines_p,return_mask=True) ## this is index on central_lines_  ## pass ints here already!!
    Ax[i,:]=ids
Ax=Ax.astype(int)
#np.save("/media/user/GINA1-BK/data_traktion_force_microscopy/temp_data_for_scripts/Ax.npy",Ax)
Ax=np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/temp_data_for_scripts/Ax.npy")
def f(x):
    return(np.max(np.where(x!=0)[0]))
R=np.linalg.qr(Ax.T)[1]
col_len=np.apply_along_axis(f,0,R)
a=(col_len[1:] - col_len[:-1])>0 # finding all "blocks "
a=np.insert(a,0,True)
a[-1]=True  ### ?????
independent=Ax[a]
np.linalg.matrix_rank(independent)

#i=np.random.choice(list(range(Ax.shape[0])),100000,replace=False)
#ax_part=Ax[i]
#np.linalg.matrix_rank(ax_part)


import numpy as np

def row_echelon(A):
    """ Return Row Echelon Form of matrix A """

    # if matrix A has no columns or rows,
    # it is already in REF, so we return itself
    r, c = A.shape
    if r == 0 or c == 0:
        return A

    # we search for non-zero element in the first column
    for i in range(len(A)):
        if A[i,0] != 0:
            break
    else:
        # if all elements in the first column is zero,
        # we perform REF on matrix from second column
        B = row_echelon(A[:,1:])
        # and then add the first zero-column back
        return np.hstack([A[:,:1], B])

    # if non-zero element happens not in the first row,
    # we switch rows
    if i > 0:
        ith_row = A[i].copy()
        A[i] = A[0]
        A[0] = ith_row

    # we divide first row by first element in it
    A[0] = A[0] / A[0,0]
    # we subtract all subsequent rows with first row (it has 1 now as first element)
    # multiplied by the corresponding element in the first column
    A[1:] -= A[0] * A[1:,0:1]

    # we perform REF on matrix from second row, from second column
    B = row_echelon(A[1:,1:])

    # we add first row and first (zero) column, and return
    return np.vstack([A[:1], np.hstack([A[1:,:1], B]) ])

A = np.array([[4, 7, 3, 8],
              [8, 3, 8, 7],
              [2, 9, 5, 3]], dtype='float')
A= np.array([[0, 0, 1, 1,1],
              [0, 1, 1, 0,1],
              [0, 1, 1, 0,1],
              [0, 0, 1, 1,0]], dtype='int')
row_echelon(A)
## search for leading 1

A_part=Ax[:1000].T
row_echelon(A_part)


a1=np.sum(R==0,axis=0)[1:] # finding all "blocks "
a2=np.sum(R==0,axis=0)[:-1]
a3=(a1-a2)<0
a3=np.insert(a3,0,True)
independent=R[:,a3]
independent=Ax[:100][a3]
##https://math.stackexchange.com/questions/748500/how-to-find-linearly-independent-columns-in-a-matrix



import numpy as np
import pickle

with open("/media/user/GINA1-BK/data_traktion_force_microscopy/temp_data_for_scripts/Ax_rref.npy", "wb") as f:
       pickle.dump(id, f)







# solve for x and y direction independently???
paths_dict={}
Ax=np.zeros(((len(all_paths),len(central_lines_p))))
Fx=np.zeros(len(all_paths))
Fy=np.zeros(len(all_paths))
for i,pas in tqdm(enumerate(all_paths),total=len(all_paths)):   #### needs some debugging ,is force correct???
    assigned=assigne_areas(pas,points,edges_dict,areas)[1]
    F1=traction_from_areas(assigned["side1"],areas)
    ids=get_line_ids_from_path(pas,central_lines_p,return_mask=True) ## this is index on central_lines_p
    paths_dict[i] = [pas, F1, ids, assigned, F1]

    #Fx[i]=F1[0]
    Fx[i] = F1[0]
    Fy[i]=F1[1]

    Ax[i,:]=ids

with open("/media/user/GINA1-BK/data_traktion_force_microscopy/temp_data_for_scripts/Fx", "wb") as f:
    pickle.dump(Fx, f)
with open("/media/user/GINA1-BK/data_traktion_force_microscopy/temp_data_for_scripts/Fy", "wb") as f:
    pickle.dump(Fy, f)
with open("/media/user/GINA1-BK/data_traktion_force_microscopy/temp_data_for_scripts/paths_dict", "wb") as f:
    pickle.dump(paths_dict, f)

x_sol=np.linalg.lstsq(Ax,Fx)


## mal teile ausprobieren