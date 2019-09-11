import numpy as np
import clickpoints
import matplotlib.pyplot as plt

from skimage.morphology import skeletonize,remove_small_holes,remove_small_objects,label,binary_dilation
from scipy.ndimage.morphology import binary_fill_holes
from scipy.spatial import cKDTree
from collections import defaultdict
from functions_for_cell_colonie import *
# reading mask


def mask_to_graph(mask):
    graph = defaultdict(list)

    points = np.array(np.where(mask)).T
    point_tree = cKDTree(points)  # look up table for nearest neigbours ??
    d = np.sqrt(2)  # maximal allowd distance
    for i, p in enumerate(points):
        neighbours = point_tree.query_ball_point(p, d)
        neighbours.remove(i)  # removing the point itself from list of its neighbours
        graph[i].extend(neighbours)
    return graph


def show_points(ps,mask):
    plt.figure()
    plt.imshow(mask)
    plt.plot(ps[:,1],ps[:,0],"or")
def graph_to_mask(mask,graph,points):
    m=np.zeros(mask.shape)
    ps=np.array([y for x in list(graph.values()) for y in x]) #flattening list of point ids
    ps_coord=points[ps] # getting coordinates
    m[ps_coord[:,0],ps_coord[:,1]]=1 #writing points
    return m
def normal_vector(coords): # works confirmed!
    '''
    :param coords: order set of coordinates
    :param intervall: intevall for calculating the normal vector, 1 means point 1 and 3 are used for support for vector on 2
    :return: n list of normal vectors
    '''
    coords = np.append(coords[-1:],coords,axis=0)  # adding values to make looping possible
    coords = np.append(coords,coords[:1],axis=0)
    vec = coords[:-(2)] - coords[2:]


    n = copy.deepcopy(vec)  # vec[:,[1,0]]
    n[:, 1] = -n[:, 1]
    n = n / np.linalg.norm(n, axis=1)[:, np.newaxis]

    n=np.append(n[1:],n[:1],axis=0)    # fixing arrangement of values
    return n

def normal_vector_from_graph(graph,points): ## implement larger window (rather difficult)  ##

    ### problems with orientation of normal vectoirs, although not importatn if we look at e.g. absoulte values of forces...
    '''
    iterating through all values of the graph dict, caclulating connecting vector between directly adjacent points

    :param graph: dictionary with connections
    :param points: list of points, where the key of the graph refferes to the row index
    :return: n: dictionary of normal vectors at each relevant point
    '''
    n={}
    for key,values in graph.items():
        if len(values)==2:
            n_vec=points[values[0]]-points[values[1]] # neighbouring points
            n_vec[1]=-n_vec[1]
            n[key]=n_vec # writing to dictionary
    return n







def normal_vector_interv(coords,interval=1): # intervall thing doesnt work , vector on the wrong point???
    '''
    :param coords: order set of coordinates
    :param intervall: intevall for calculating the normal vector, 1 means point 1 and 3 are used for support for vector on 2
    :return: n list of normal vectors
    '''
    coords = np.append(coords[-interval:],coords,axis=0)  # adding values to make looping possible
    coords = np.append(coords,coords[:interval],axis=0)
    vec = coords[:-(interval+1)] - coords[(interval+1):]


    n = copy.deepcopy(vec)  # vec[:,[1,0]]
    n[:, 1] = -n[:, 1]
    n = n / np.linalg.norm(n, axis=1)[:, np.newaxis]

    n=np.append(n[interval:],n[:interval],axis=0)    # fixing arrangement of values
    return n



def n_shear_stress(coords,stress_tensor,n):
    '''

    :param coords: ordered set of coordinates (in the same order as n)
    :param stress_tensor: 2d array of all stress tensor
    :param n: list of normal vetors at points coords
    :return:
    '''
    stress_vec=[]
    for nv,tens in zip(n,stress_tensor[coords[1:-1,0],coords[1:-1,1]]): # do this fully array based
        stress_vec.append(np.matmul(tens,nv))
    stress_vec=np.array(stress_vec)

    n_stress=[]
    shear_stress=[]
    for nv, st in zip(n,stress_vec):   ## aslo matrix
        n_stress.append(np.dot(nv,st))
        shear_stress.append(np.sqrt(np.dot(st,st)-(np.dot(nv,st)**2)))
    #n_stress=np.matmul(stress_vec,n.T) # only diagonla elements are relevant heres???# equivalent to dot products for al pairs
    n_stress = np.array(n_stress)
    shear_stress= np.array(shear_stress)
    return n_stress,shear_stress



def remove_endpoint(graph,ep):
    print(graph[ep])
    new_p = graph[ep]  # neighours of this point

    check_line_break = np.array([len(graph[p])==2 for p in new_p]).any() and len(new_p)>1 #check if any removal of new_p would cause single and
    # if new_p is not just another loose end, this deals with special cases when the line hit another line
    if check_line_break:
        return
    for p in new_p:
        graph[p].remove(ep)  # remove all further references in the graph
    graph.pop(ep)  # deleting loose endpoint from graph, gives empty list

    # check if removing new point would produce line brek


    if len(new_p)<2:
        remove_endpoint(graph,new_p[0]) # would stop if it encounters point with three neighbours
    else:
        return

### is ok, but sometimes leaves a single point to much aslo could fail under certain conditions...


db=clickpoints.DataFile("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/mask_cell_boundary3.cdb","r")



mask=db.getMask(frame=0).data
# removing small holes
mask=remove_small_holes(mask,100)
mask=remove_small_objects(label(mask),1000)>0 # removing other small bits

#plt.figure()
#plt.imshow(label(mask))
mask=skeletonize(mask)
#plt.figure()
#plt.imshow(label(mask))


# construc a graph with all points using algoriths to identfy nearest neigbours
graph=defaultdict(list)

points=np.array(np.where(mask)).T
point_tree = cKDTree(points)  #look up table for nearest neigbours ??
d=np.sqrt(2) # maximal allowd distance
for i,p in enumerate(points):
    neighbours=point_tree.query_ball_point(p,d)
    neighbours.remove(i) # removing the point itself from list of its neighbours
    graph[i].extend(neighbours)



end_points=np.where(np.array([len(v) for v in graph.values()])<2)[0]
points_2=np.array(list(graph.keys())) # all keys as array
eps=points_2[end_points] # keys of endpoints

for ep in eps:     ## do i even need this???
    print(ep)
    if ep in list(graph.keys()):
        remove_endpoint(graph,ep)
show_points(points[eps],mask) #showing endpoints before removel
plt.figure()
plt.imshow(graph_to_mask(mask,graph,points)) # showing mask after removal of endpoints



mask=graph_to_mask(mask,graph,points) # rewriting mask with deleted points


###########  naja
# simple identification  of all cell borders:
# is abit slow but works good // probalbl better with purely graph theory based approach
borders={} # dictionary to contain set of points reffering to a border
mask2=binary_fill_holes(mask)
mask3=mask2-mask # masking the cells

mask3_l=label(mask3,connectivity=1) #labeling each cell

#points_to_flatt array map:
points_fl=(points[:,0])+mask.shape[0]*(points[:,1])

for l in np.unique(mask3_l)[1:]: # getting each cell border by binary dialtion of one pixel
    m=(mask3_l==l)*1       ## this is mathematically correct?!
    edge=binary_dilation(m)-m
    ps=np.where(edge)
    ps_fl=(ps[0])+mask.shape[0]*(ps[1]) #points_to_flatt array map:
    borders[l]= np.where(np.isin(points_fl, ps_fl))[0] #ids as used in the graph



border_p=borders[l]
graph_part={} # extarcting relevant part of the graph
for p in border_p:
    graph_part[p]=[k for k in graph[p] if k in border_p]  # also removing neighbours if the dont belong there

# np.array([len(graph_part[p]) > 1 for p in graph_part.keys()]).all() check if circular
order=find_path_circular(graph_part, border_p[0], border_p[0], path=[])[:-1] # orderd set of points on circle
points[order]



n=normal_vector(points[order])
check_normal_vectors(mask,points[order],n)



n_shear_stress(points[order],stress_tensor,n)




#ps=points[order] #reordering
#check_order(mask,ps)
# calculating normal vectors for all border points










'''
#showing endpoints
end_points=np.where(np.array([len(v) for v in graph.values()])<2)[0]
points_2=np.array(list(graph.keys()))
ps=points[points_2[end_points]]




'''
