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
inf = float('inf')
Edge = namedtuple('Edge', 'start, end, cost')

class path_error(Exception):
    pass
def make_edge(start, end, cost=1):
  return Edge(start, end, cost)

def make_neighbours_dict(data_list):
    neighbours_dict = defaultdict(list)
    for i, j,dist in data_list:
        neighbours_dict[i].append(j)
        neighbours_dict[j].append(i)
    for key, value in neighbours_dict.items():
        neighbours_dict[key] = list(set(value))
    return neighbours_dict


class Graph:
    def __init__(self, edges):
        # let's check that the data is right
        wrong_edges = [i for i in edges if len(i) not in [2, 3]]
        if wrong_edges:
            raise ValueError('Wrong edges data: {}'.format(wrong_edges))

        self.edges = [make_edge(*edge) for edge in edges]
        self.neighbours_dict=make_neighbours_dict(edges) #already because of @prperty solved
    @property
    def vertices(self):
        return set(
            sum(
                ([edge.start, edge.end] for edge in self.edges), []
            )
        )

    def get_node_pairs(self, n1, n2, both_ends=True):
        if both_ends:
            node_pairs = [[n1, n2], [n2, n1]]
        else:
            node_pairs = [[n1, n2]]
        return node_pairs

    def remove_edge(self, n1, n2, both_ends=True):
        node_pairs = self.get_node_pairs(n1, n2, both_ends)
        edges = self.edges[:]
        for edge in edges:
            if [edge.start, edge.end] in node_pairs:
                self.edges.remove(edge)
        if n1 in  self.neighbours_dict.keys():
            if n2 in self.neighbours_dict[n1]: # removing if existingg, necessary because of biothe ends, might already be removed
                self.neighbours_dict[n1].remove(n2)
            if len( self.neighbours_dict[n1])==0:
                self.neighbours_dict.pop(n1)
        if both_ends:
            if n2 in self.neighbours_dict.keys():
                if n1 in self.neighbours_dict[n2]:
                    self.neighbours_dict[n2].remove(n1)
                if len(self.neighbours_dict[n2])==0:
                    self.neighbours_dict.pop(n2)
            # removing empty keys
        #self.neighbours_dict={ key : values for key,values in self.neighbours_dict.items() if len(values)>0}



    def add_edge(self, n1, n2, cost=1, both_ends=True):
        node_pairs = self.get_node_pairs(n1, n2, both_ends)
        for edge in self.edges:
            if [edge.start, edge.end] in node_pairs:
                return ValueError('Edge {} {} already exists'.format(n1, n2))

        self.edges.append(Edge(start=n1, end=n2, cost=cost))
        self.neighbours_dict[n1].append(n2)  # default dict allows appending por creating in one go
        if both_ends:
            self.edges.append(Edge(start=n2, end=n1, cost=cost))
            self.neighbours_dict[n2].append(n1)
    @property
    def neighbours(self):
        neighbours = {vertex: set() for vertex in self.vertices}
        for edge in self.edges:
            neighbours[edge.start].add((edge.end, edge.cost))

        return neighbours

    def dijkstra(self, source, dest):
        assert source in self.vertices, 'Such source node doesn\'t exist'
        distances = {vertex: inf for vertex in self.vertices}
        previous_vertices = {
            vertex: None for vertex in self.vertices
        }
        distances[source] = 0
        vertices = self.vertices.copy()

        while vertices:
            current_vertex = min(
                vertices, key=lambda vertex: distances[vertex])
            vertices.remove(current_vertex)
            if distances[current_vertex] == inf:
                break
            for neighbour, cost in self.neighbours[current_vertex]:
                alternative_route = distances[current_vertex] + cost
                if alternative_route < distances[neighbour]:
                    distances[neighbour] = alternative_route
                    previous_vertices[neighbour] = current_vertex

        path, current_vertex = deque(), dest
        while previous_vertices[current_vertex] is not None:
            path.appendleft(current_vertex)
            current_vertex = previous_vertices[current_vertex]
        if path:
            path.appendleft(current_vertex)
        return path


    def random_path(self, source, dest, path=[]):
        assert source in self.vertices, 'Such source node doesn\'t exist'
        path = path + [source]
        if source == dest:
            return path
        new_nodes=self.neighbours_dict[source]
        np.random.shuffle(new_nodes)
        for node in new_nodes: # dest not needed
            if node not in path:
                newpath = self.random_path(node, dest, path)
                if newpath: return newpath
        return None




def find_path(graph, start, end, path=[]):
    path = path + [start]
    if start == end:
        return path
    if not start in graph.keys():
        return None
    for node in graph[start]:
        if node not in path:
            newpath = find_path(graph, node, end, path)
            if newpath: return newpath
    return None
def plot_graph(graph,points,mask):
    plt.figure()
    plt.imshow(mask)
    for node1,neighbours in graph.items():
        for node2 in neighbours:
            plt.plot([points[node1][0],points[node2][0]],[points[node1][1],points[node2][1]])

def check_duplication(paths):
    ps=[sorted(list(l)) for l in paths]
    u_paths,index,counts=np.unique(ps,return_counts=True,return_index=True)
    n_duplicates=np.sum(counts>1)
    return n_duplicates

# finding all connecting paths:
def plot_path_with_lines_from_single_id(line_ids,mask):
    plt.figure()
    plt.imshow(mask)
    for l_ids in line_ids:
        edge=lines[l_ids]
        plt.arrow(edge[0][0], edge[0][1], edge[1][0] - edge[0][0], edge[1][1] - edge[0][1], head_width=20)


def plot_graph_from_object(graph,points,mask):
    plt.figure()
    plt.imshow(mask)
    plt.axis("off")
    for edge in graph.edges:
        plt.plot([points[edge.start][0],points[edge.end][0]],[points[edge.start][1],points[edge.end][1]])



from collections import defaultdict
a = defaultdict(lambda: "default", key="some_value")


def find_path_custom(graph,start,end, tried_direc=defaultdict(lambda: []),path=[]):
    path=path+[start] # adding new node
    node=path[-1] # setting current node to explore
    print(node)
    print(path)
    if end in path:
        return path,tried_direc
    for nodes in graph[node]:# exploring

        if nodes not in path and nodes not in tried_direc[node]:
            #print(nodes)
            res=find_path_custom(graph,nodes,end,tried_direc,path)
            if res:
                return res,tried_direc

        if len(path) in tried_direc.keys():
            tried_direc[node].append(node)  #
        else:
            tried_direc[node]=[]



def get_point_ids_from_lines(line,return_single=False):
    points_arr=np.array([np.insert(value,0,key) for key,value in points.items()])
    if return_single:
        p=[]
        for edge in line:
            check=(points_arr[:, (1, 2)] - edge[0]) == 0
            check=check[:,0]*check[:,1]
            p.extend(np.where(check)[0])
            check = (points_arr[:, (1, 2)] - edge[1]) == 0
            check = check[:, 0] * check[:, 1]
            p.extend(np.where(check)[0])
            # removing strickly duplicate
        p=np.array(p)
        p[:-1]-p[1:]
        filter=(p[:-1]-p[1:])!=0
        filter=np.insert(filter,0,True)
        p=p[filter]
        points_selec=points_arr[:,0][p].astype(int) # getting actual point ids
        return points_selec
    else:
        l=[]
        for edge in line:
            check=(points_arr[:, (1, 2)] - edge[0]) == 0
            check=check[:,0]*check[:,1]
            p1=np.where(check)[0]
            check = (points_arr[:, (1, 2)] - edge[1]) == 0
            check = check[:, 0] * check[:, 1]
            p2=np.where(check)[0]
            # removing strickly duplicate
            l.append([points_arr[:,0][int(p1)],points_arr[:,0][int(p2)]])

        l=np.array(l,dtype=int)
        return np.array(l)


def plot_path(path,mask):
    plt.figure()
    plt.imshow(mask)
    for i in range(len(path)-1):
        plt.plot([points[path[i]][0], points[path[i+1]][0]], [points[path[i]][1],points[path[i+1]][1]])
def plot_path_from_line_id(l_ids,mask,lines_p):
    plt.figure()
    plt.imshow(mask)
    for l in l_ids:
        line=lines_p[l]
        plt.plot([points[line[0]][0], points[line[1]][0]], [points[line[0]][1],points[line[1]][1]])


def plot_points(points_id,mask,points):
    plt.figure()
    plt.imshow(mask)
    for i,pi in enumerate(points_id):
        plt.plot(points[pi][0],points[pi][1],"o")
        plt.text(points[pi][0]+4,points[pi][1]+4,str(i))



def plot_path_with_lines(line,mask):
    plt.figure()
    plt.imshow(mask)
    for edge in line:
        plt.arrow(edge[0][0], edge[0][1], edge[1][0] - edge[0][0], edge[1][1] - edge[0][1], head_width=20)


def get_line_ids_from_path(path,lines_p,return_mask=True):
    if return_mask:
        m=np.zeros(len(lines_p))
        for i in range(len(path) - 1):
            point_pair = [path[i], path[i + 1]]
            point_pair.sort()
            m+=(lines_p[:, 0] == point_pair[0]) * (lines_p[:, 1] == point_pair[1])
        m=m>0
        return m


    else:
        ids=[]
        for i in range(len(path)-1):
            point_pair=[path[i],path[i+1]]
            point_pair.sort()

            ids.append(np.where((((lines_p[:,0]==point_pair[0]) * (lines_p[:,1]==point_pair[1])) +
                                ((lines_p[:, 1] == point_pair[0]) * (lines_p[:, 0] == point_pair[1])))>0 )[0])

        return np.array(ids).flatten()



def vizualize_assignement(side_assginement,points,areas,path,line_dict1,edges_dict,mask):
    plt.figure()   ## this seems to work
    plt.imshow(mask)
    for i in range(len(path)-1):
        plt.plot([points[path[i]][0], points[path[i+1]][0]], [points[path[i]][1],points[path[i+1]][1]])
    for i,area_id in enumerate(areas.keys()):
        com=areas[area_id][4]
        plt.text(com[0],com[1],str(i))
        if side_assginement[i]:
            plt.plot(com[0],com[1],"o",color="red")
        else:
            plt.plot(com[0], com[1],"o", color="green")
    for v in line_dict1.values():
        plt.arrow(v[1][0],v[1][1],v[0][0],v[0][1],head_width=20)





def assigne_areas(path,points,edges_dict,areas):   ##### just clean this up

    line_dict1 = {}
    line_vecs = []
    line_points = []

    # line_points_association# to witch lines do the points belong
    for i in range(len(path) - 1):
        ps = edges_dict[path[i], path[i + 1]]
        line_points.extend(ps)
        vec = np.array(edges_dict[path[i], path[i + 1]][-1]) - np.array(edges_dict[path[i], path[i + 1]][0])
        line_vecs.extend([vec] * len(ps))
        line_dict1[path[i], path[i + 1]] = [vec, ps[0], ps[-1]]

    line_points = np.array(line_points)  # some points appera twice but thats ok ?
    line_vecs = np.array(line_vecs)  # some points appera twice but thats ok ?

    mask_edge_points = np.sum(np.abs(line_points[1:] - line_points[:-1]),
                              axis=1) == 0  # finding critical edge point regions
    mask_edge_points = np.insert(mask_edge_points, 0, False)
    for i, edge_pos in enumerate(np.where(mask_edge_points)[0]):
        vec_e = line_points[edge_pos + 1] - line_points[edge_pos - 2]
        line_vecs[edge_pos - 1] = vec_e
        line_vecs[edge_pos] = vec_e
        line_dict1[i] = [vec_e, line_points[edge_pos - 2], line_points[edge_pos + 1]]

    side_assginement = []

    for area_id in areas.keys():
        com = areas[area_id][4]  # getting center of mass of area
        p1 = np.argmin(np.linalg.norm(np.abs(line_points - com), axis=1))  # closeset point on line
        ## little question: why does not sum of absolute values work???
        line_vec = line_vecs[p1]
        vec_com = com - line_points[p1]  # vector to center of mass
        side_assginement.append(np.cross(line_vec, vec_com) > 0)  # orientation via cross product
    side_assginement = np.array(side_assginement)
    assigned_areas ={"side1":np.array(list(areas.keys()))[side_assginement],
                     "side2":np.array(list(areas.keys()))[~side_assginement]}
    return side_assginement,assigned_areas
    #vizualize_assignement(side_assginement, points, areas, path, line_dict1, edges_dict, mask)  # plotting assignement


##### loading mask, refining and setting line segments
db=clickpoints.DataFile("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/mask_cell_boundary3.cdb")
# not yet used
mask=db.getMask(frame=0).data
mask=binary_closing(mask)
mask=binary_erosion(mask)
mask=binary_erosion(mask)
mask=binary_erosion(mask)
mask=skeletonize(mask)
#plt.figure();plt.imshow(mask)

# getting full area:
mask_tract=binary_fill_holes(mask)


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
np.save("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/connections.npy",connections_f)


# list point coordinates as dict:
points={id:np.array([db.getMarkers(id=id)[0].x,db.getMarkers(id=id)[0].y]) for id in np.unique(connections.flatten())}

# making lines:
lines=[]
lines_p=[] # line as list of point ids
for ps in connections_f:

    markers=db.getMarkers(id=list(ps),type="marker")
    p1 = (markers[0].x,markers[0].y)
    p2 = (markers[1].x, markers[1].y)
    lines.append([p1,p2])
    lines_p.append(ps)
lines=np.array(lines)
lines_p=np.array(lines_p)
edges_dict={}
for i,j in lines_p:
    edges_dict[i, j]=get_line(points[i].round().astype(int),points[j].round().astype(int))
    edges_dict[j, i] = get_line(points[j].round().astype(int),points[i].round().astype(int))


db.db.close()



# representation as graph: (could be better written like this in the begining)

graph1={i:[] for i in np.unique(connections.flatten())}

for i,j in connections:
    graph1[i].append(j)
    graph1[j].append(i)
for key,value in graph1.items():
    graph1[key]=np.unique(value)

### finding(the shortest path between two points)
### desired data format : node 1 node 2 length of edge
data_list=[]
for p1,ps in graph1.items():
    for p2 in ps:
        dist=np.linalg.norm(points[p1]-points[p2])
        data_list.append((p1, p2, dist))
        data_list.append((p2, p1, dist))

graph2 = Graph(data_list)

path=graph2.dijkstra(56, 28)



# calculating tractions on two side of the path



### loading areas object from other script:

# area dictionary ith associated lines from other script
import pickle
with open('/media/user/GINA1-BK/data_traktion_force_microscopy/temp_data_for_scripts/areas.pickle', 'rb') as handle:
    areas = pickle.load(handle)

##identifying areas to the left and right of deviding line








## setting up a system of equations:
#write as s=A*e:
#s sum of forces over cut trough cell colony
#A matrix of coefficients, where coefficent simply denotes the presence or absence of an edge
# e: forces over individual edge
# goal: get m (number of edges) independent cuts, then solve:
#e=s*A-1






# finding all necessary equations:
# strategie: find points at outer edge
# iterate through all other lines: find connection betweeen two random end node pairs, including the lines

#1. find edges:
# using same function to identify al areas
start_point=np.unravel_index(np.argmax(lines),lines.shape) # maximum coordinate must lie at the border..
com_all=np.array([np.mean(lines[:,:,0]),np.mean(lines[:,:,1])])
border_line,border_line_id=find_areas(lines[start_point[0]],lines,start_point[0],com_all,invert_direction=True)
border_points=get_point_ids_from_lines(border_line,return_single=True)
border_points2=np.unique(np.array([p for p in border_points if len(graph1[p])>2])) ## borderpoints without points that have no internal connections

border_line_id=get_point_ids_from_lines(border_line,return_single=False)
#plot_path(border_points,mask)
#plot_path_with_lines(border_line,mask)
central_points=list(set(list(points.keys()))-set(border_points))  # "diffrence" of borderpoints and alll points
#2. iterating through all central lines
##




# reducing to central lines
border_lines_id=np.zeros(len(lines_p))>0
for i in border_line_id:
    border_lines_id += (lines_p[:, 0] == i[0]) * (lines_p[:, 1] == i[1])
    border_lines_id += (lines_p[:, 1] == i[0]) * (lines_p[:, 0] == i[1])
border_lines_id=np.where(border_lines_id)[0] ## indices of lines at border
central_lines_id=list(set(list(range(len(lines_p))))-set(border_lines_id)) ## indices of all lines
central_lines_p=lines_p[central_lines_id]
# reducing the graph
graph1_cent=copy.deepcopy(graph1)
for l_id in border_lines_id:
    graph1_cent[lines_p[l_id][0]]=graph1_cent[lines_p[l_id][0]][graph1_cent[lines_p[l_id][0]]!=lines_p[l_id][1]] #removing edges in both directions
    graph1_cent[lines_p[l_id][1]]=graph1_cent[lines_p[l_id][1]][graph1_cent[lines_p[l_id][1]]!=lines_p[l_id][0]] #removing edges in both directions
graph1_cent={k: v for k, v in graph1_cent.items() if len(v)!=0} # removing empty keys= nodes at border with no internal connection
# removing nodes at border with no internal connection



plot_graph(graph1,points,mask)
plot_graph(graph1_cent,points,mask)






def remove_path_form_graph(graph,path): ## write this into the class
    for p in range(len(path)-1):
        graph.remove_edge(path[p],path[p+1],both_ends=True)

def add_path_form_graph(graph, path):
    for p in range(len(path) - 1):
        graph.add_edge(path[p], path[p + 1],cost=np.linalg.norm(points[path[p]] - points[path[p+1]]), both_ends=True)

def remove_node_from_graph(graph_obj,graph_dict,nodes):
    for node in nodes:
        for p in graph_dict[node]:
            graph_obj.remove_edge(p,node,both_ends=True)

def add_node_from_graph(graph_obj,graph_dict,nodes):
    for node in nodes:
        for p in graph_dict[node]:
            graph_obj.add_edge(p,node,cost=np.linalg.norm(points[node] - points[p]),both_ends=True)

def get_optimal_point_pair(l,start,end,points):
    l1=np.array([l[0],l[1]]) # look up for minimum
    l2=np.array([start,end])

    a1=np.array([points[l[0]],points[l[1]]])
    a2=np.array([points[start],points[end]])
    dist=np.linalg.norm(a1[None,:]-a2[:,None],axis=2) # distance between all points in a1 and a2
    # matrix is
    #[[l0-a21,l1-a21],
     #[l0-a22,l1-a22]]
    # with just two combintations
    dist1,dist2=dist[0,0]+dist[1,1],dist[0,1]+dist[1,0]
    if dist1<dist2:
        return l[0],start,l[1],end
    else:
        return l[0], end, l[1], start








graph_as_list=[]
for l in lines_p[central_lines_id]:
    dist=np.linalg.norm(points[l[0]]-points[l[1]])
    graph_as_list.append((l[0],l[1],dist))
    graph_as_list.append((l[1], l[0], dist))
graph_dijkstra=Graph(graph_as_list) ## format needed for shortest path algorithme
graph_di_temp=copy.deepcopy(graph_dijkstra)




def trying_paths(l,start,end, graph_di_temp,points,graph1_cent,reversed_dir=False):

    if reversed_dir:
        l[1], end, l[0], start = get_optimal_point_pair(l, start, end, points)  # finding optimal pair to connect
    else:
        l[0], start, l[1], end = get_optimal_point_pair(l, start, end, points)

    #print("l1",l[1], "end",end, "l0",l[0], "start",start)
    remove_node_from_graph(graph_di_temp, graph1_cent, [l[1]])  # avoid passing over other node
    if start not in graph_di_temp.vertices or l[0] not in graph_di_temp.vertices: #### eventuially just use try and except for thiss.......
        return None

    path1 = graph_di_temp.dijkstra(start, l[0])
    #plot_path(path1,mask)
    #plt.title("path1")
    add_node_from_graph(graph_di_temp, graph1_cent, [l[1]])  # re adding the node, to allow start at this node

    remove_node_from_graph(graph_di_temp, graph1_cent, path1)  # avoid using any edges that are the same as before (might) as well use node here??

    if not path1: # if finding first part is not posiible
        #print("warning, disconnected graph")
        return None
    if end not in graph_di_temp.vertices or l[1] not in graph_di_temp.vertices:  ## happens sometimes if path is blocked
        #print(i, l, start, end)
        #plot_graph_from_object(graph_di_temp, points, mask)
        #plot_path(path1, mask)
        return None


    path2 = graph_di_temp.dijkstra(l[1], end)
    #plot_path(path2, mask)
    #plt.title("path2")
    if not path2: # if finding second part is not possible
        return None
    path1 += path2
    return path1

def trying_paths_random(l,start,end, graph_dijkstra,points,graph1_cent,paths_binary,ids,paths_set,reversed_dir=False):

    for k in range(50):

        print("k = ",k)

        graph_di_temp = copy.deepcopy(graph_dijkstra)
        if reversed_dir:
            l[1], end, l[0], start = get_optimal_point_pair(l, start, end, points)  # finding optimal pair to connect
        else:
            l[0], start, l[1], end = get_optimal_point_pair(l, start, end, points)

        #print("l1",l[1], "end",end, "l0",l[0], "start",start)
        remove_node_from_graph(graph_di_temp, graph1_cent, [l[1]])  # avoid passing over other node
        if start not in graph_di_temp.vertices or l[0] not in graph_di_temp.vertices: #### eventuially just use try and except for thiss.......
            return None

        path1 = graph_di_temp.random_path(start, l[0])
        #plot_path(path1,mask)
        #plt.title("path1")
        add_node_from_graph(graph_di_temp, graph1_cent, [l[1]])  # re adding the node, to allow start at this node

        remove_node_from_graph(graph_di_temp, graph1_cent, path1)  # avoid using any edges that are the same as before (might) as well use node here??

        if not path1: # if finding first part is not posiible
            #print("warning, disconnected graph")
            continue
        if end not in graph_di_temp.vertices or l[1] not in graph_di_temp.vertices:  ## happens sometimes if path is blocked
            #print(i, l, start, end)
            #plot_graph_from_object(graph_di_temp, points, mask)
            #plot_path(path1, mask)
            continue


        path2 = graph_di_temp.random_path(l[1], end)
        #plot_path(path2, mask)
        #plt.title("path2")

        if not path2: # if finding second part is not possible
            continue
        path1 += path2
        tqdm.write(str(check_path(path1, paths_binary, ids, paths_set)))
        if check_path(path1,paths_binary,ids,paths_set):
            return path1
    return None

def find_path_ran(i,l,border_points2,graph_dijkstra,points,mask,graph1_cent,paths_binary,ids,paths_set,bps_comb):
    np.random.shuffle(bps_comb)
    for i,bps in enumerate(bps_comb):
        print("bps_count = ",i)
        start,end=bps
        if l[0] in border_points2:  # selecting either random start poin or "close border point and random end point)
            for k in range(50):
                start = l[0]
                end = np.random.choice(border_points2[border_points2 != l[0]], 1, replace=False)[0]
                path1 = graph_dijkstra.random_path(start, end)
                if check_path(path1,paths_binary,ids,paths_set):
                    return path1


        elif l[1] in border_points2:
            for k in range(50):
                start = l[1]
                end = np.random.choice(border_points2[border_points2 != l[1]], 1, replace=False)[0]
                path1 = graph_dijkstra.random_path(start, end)
                if check_path(path1,paths_binary, ids, paths_set):
                    return path1
        else:
            path1=trying_paths_random(l,start,end, graph_dijkstra,points,graph1_cent,paths_binary,ids,paths_set,reversed_dir=False)
            if path1:
                return path1
            path1 = trying_paths_random(l,start,end, graph_dijkstra,points,graph1_cent,paths_binary,ids,paths_set,reversed_dir=True)
            if path1:
                return path1

    return None

def check_path(path1,paths_binary,ids,paths_set):
    if not ids is None:
        mp_temp = np.array(paths_binary + ids)
    if len(paths_binary)==0:
        return True
    if path1 is None:
        return False
    return ( (len(paths_binary) - np.linalg.matrix_rank(mp_temp)) == 0 and (set(path1) not in paths_set))


def get_bp_pairs(border_points2):
    bps_comb=list(itertools.combinations(border_points2, 2))# all combinations of borderpoints
    for l in bps_comb:
        if l[0]==l[1]:
            bps_comb.pop(l)
    return np.array(bps_comb)

bps_comb=get_bp_pairs(border_points2)
       # plot_points([start,l[0],l[1],end],mask,points)
paths=[]
paths_set=[] # no directional information
paths_binary=[] # representation as 01
ids=None

for i,l in enumerate(lines_p[central_lines_id]): #iterating only through central lines
    print("line_seg = ",i)
    path1=find_path_ran(i,l,border_points2,graph_dijkstra,points,mask,graph1_cent,paths_binary,ids,paths_set,bps_comb)
    #print(b)
    #print(path1)
    if path1:
        ids = get_line_ids_from_path(path1, central_lines_p, return_mask=True)  ## this is index on central_lines_p
        paths.append(path1)

        paths_binary.append(ids)
        paths_set.append(set(path1))
    else:
        print(" no paths found for ",i)
print(len(paths))
print(np.linalg.matrix_rank(np.array(paths_binary)))

#checking paths:
if None in paths:
    raise path_error("path contains nones")
print(check_duplication(paths)," duplicates found")


for i in paths:
    plot_path(i,mask)
# using path to construct equation system


def traction_from_areas(ids,areas):
    T=np.zeros(2)
    for id in ids:
        T+=areas[id][5]
    return T


# solve for x and y direction independently???
paths_dict={}
Ax=np.zeros(((len(paths),len(central_lines_p))))
Fx=np.zeros(len(paths))
Fy=np.zeros(len(paths))
for i,pas in enumerate(paths):   #### needs some debugging ,is force correct???
    assigned=assigne_areas(pas,points,edges_dict,areas)[1]
    F1=traction_from_areas(assigned["side1"],areas)
    ids=get_line_ids_from_path(pas,central_lines_p,return_mask=True) ## this is index on central_lines_p
    paths_dict[i] = [pas, F1, ids, assigned, F1]

    #Fx[i]=F1[0]
    Fx[i] = np.linalg.norm(F1)
    Fy[i]=F1[1]

    Ax[i,:]=ids
Ax=Ax.astype(int)
#np.linalg.solve(Ax,Fx)
###### whats with the side problem (minus or plus signe???)

#e_v=np.isclose(np.linalg.eig(Ax)[0],0)
x_sol=np.linalg.lstsq(Ax,Fx)[0]
test=np.matmul(Ax,x_sol[0])
Fx-test
# system of equations:
#F=A*X



Q,R = np.linalg.qr(Ax) # qr decomposition of A
Qf = np.dot(Q.T,Fx) # computing Q^T*b (project b onto the range of A)
x_qr = np.linalg.solve(R,Qf) # solving R*x = Q^T*b


# problem: 60 42 32 18

## prducing the required cuts



import random
import sympy as symp
import itertools

#A=np.array([[0,1,1],[1,0,1],[1,1,0],[1,1,1]])
#b=[0,1,2,3,4,5]  # indices representing list
#indices=list(itertools.combinations(b,2))

#l=[]
#for i in indices:
#    empty = [0] * 6
#    for j in i:
#        empty[j]=1
#    l.append(empty)



## testing out if it is generally posiible
l=[[0,1,1,0,1],[0,1,1,1,0],[1,1,0,0,0],[0,0,0,1,1],[1,0,1,0,1],[1,0,1,1,0]]
all_l_comb=list(itertools.combinations(l,6))

##
#alternative parametrization

l=[[1,1,0,0,0],[-1,0,-1,-1,0],[0,-1,1,0,1],[0,0,0,1,-1],[-1,0,-1,0,-1],[0,-1,1,1,0]]


all_l_comb=list(itertools.combinations(l,6))

no_inverse=0
for k in all_l_comb:
    A=np.array(k)
    A=symp.Matrix(A)
    #if len(A.rref()[1])<5:
        #no_inverse+=1
#print(no_inverse/len(all_l_comb))
# simple example:

ar1,ar2,ar3,ar4=[1,1,-1,-1]

A=np.array(list(itertools.combinations(l,6))[0])
F=np.array([[ar1],[ar2],[ar3],[ar4],[ar2+ar4],[ar3+ar4]])
F=F.flatten().astype(int)
#F=A*x
for sel in list(itertools.combinations([0,1,2,3,4,5],5)):
    sel=np.array(sel)
    try:
        print(np.linalg.solve(A[sel],F[sel]))
    except:
        pass
##### ask christoph about thisssssssssssss

np.linalg.lstsq(A[:5],F[:5])

x1=np.array([1,0,0,0,1])
x2=np.linalg.lstsq(A[:5],F[:5])[0]
np.matmul(A,x1.T)
np.matmul(A,x2.T)
'''


### testing for area assignement


pa=paths[60]
line_dict1={}
line_vecs=[]
line_points=[]

#line_points_association# to witch lines do the points belong
for i in range(len(pa)-1):
    ps=edges_dict[pa[i],pa[i+1]]
    line_points.extend(ps)
    vec=np.array(edges_dict[pa[i],pa[i+1]][-1])-np.array(edges_dict[pa[i],pa[i+1]][0])
    line_vecs.extend([vec]*len(ps))
    line_dict1[pa[i],pa[i+1]]=[vec,ps[0],ps[-1]]

line_points=np.array(line_points) # some points appera twice but thats ok ?
line_vecs=np.array(line_vecs) # some points appera twice but thats ok ?


mask_edge_points=np.sum(np.abs(line_points[1:]-line_points[:-1]),axis=1)==0 # finding critical edge point regions
mask_edge_points=np.insert(mask_edge_points,0,False)
for i,edge_pos in enumerate(np.where(mask_edge_points)[0]):
    vec_e=line_points[edge_pos+1]-line_points[edge_pos-2]
    line_vecs[edge_pos-1]=vec_e
    line_vecs[edge_pos]=vec_e
    line_dict1[i]=[vec_e,line_points[edge_pos-2],line_points[edge_pos+1]]



side_assginement=[]
for area_id in areas.keys():
    com=areas[area_id][4] # getting center of mass of area
    p1=np.argmin(np.linalg.norm(np.abs(line_points-com),axis=1))  # closeset point on line
    ## little question: why does not sum of absolute values work???
    line_vec= line_vecs[p1]
    vec_com=com-line_points[p1]# vector to center of mass
    side_assginement.append(np.cross(line_vec,vec_com)>0)  # orientation via cross product
side_assginement=np.array(side_assginement)

vizualize_assignement(side_assginement,points,areas,pa,line_dict1,edges_dict,mask) # plotting assignement

'''


'''

plt.figure()   ## this seems to work
mask_show=copy.deepcopy(mask)*np.nan
for n,p in enumerate(line_points):
    mask_show[p[1],p[0]]=n
plt.imshow(mask_show)
plt.colorbar()
plt.plot(com[0],com[1],"o")
plt.plot(p1[0],p1[1],"o")

'''





'''
# splitting along line
  ## please put this into a class
graph_nodes={}
for p1,ps in graph1.items():
    for p2 in ps:
        graph_nodes[(p1,p2)]=get_line(np.round(points[p1]).astype(int),np.round(points[p2]).astype(int))  #dictionary with all points on lines of nodes

#point one one full path
path_points=[]
## adding a bit extended points at the end and the begining
for i in range(len(path1)-1):

    path_points.extend(graph_nodes[(path1[i],path1[i+1])])

## adding a bit extended points at the end and the begining
b_p=points[path1[0]]-points[path1[1]] + points[path1[0]]
e_p=points[path1[-1]] - (points[path1[-2]]-points[path1[-1]])
path_points.extend(get_line(np.round(b_p).astype(int),np.round(points[path1[0]]).astype(int)))
path_points.extend(get_line(np.round(e_p).astype(int),np.round(points[path1[-1]]).astype(int)))
path_points=np.array(path_points)


'''

'''

mask_path=copy.deepcopy(mask_tract)*1


mask_path[path_points[:,1],path_points[:,0]]=0
plt.figure()
plt.imshow(mask_path)
plt.plot(e_p[0],e_p[1],"o")
plt.plot(b_p[0],b_p[1],"o")
## need to fit a curev in ther to extend two the sides?? ---- use fits at the end for everything!!!
# maybe just sratighten out a it     ##### or better solution in general
#### use areas and connection with lines.....

## just simple example to get areas
areas=label(mask_path,connectivity=1)
areas=remove_small_objects(areas,min_size=1000,connectivity=1)
print(np.unique(areas))

coords=[r.coords for r in regionprops(areas)]
coms=[r.centroid for r in regionprops(areas)]
a1=np.zeros(np.shape(mask))
a1[coords[0][:,0],coords[0][:,1]]=1
a1=a1>0
a2=np.zeros(np.shape(mask))
a2[coords[1][:,0],coords[1][:,1]]=1
a2=a2>0


f1=np.array([np.sum(tx_resize[a1]*((pixelsize*10**-6)**2)),np.sum(ty_resize[a1]*((pixelsize*10**-6)**2))])
f2=np.array([np.sum(tx_resize[a2]*((pixelsize*10**-6)**2)),np.sum(ty_resize[a2]*((pixelsize*10**-6)**2))])


plt.figure()
plt.imshow(a1*1+2*a2)
scale=10**9
plt.arrow(coms[0][1],coms[0][0],f1[1]*scale,f1[0]*scale,head_width=20) # check arrow orientation again
plt.arrow(coms[1][1],coms[1][0],f2[1]*scale,f2[0]*scale,head_width=20)

'''