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
from scipy.spatial import ConvexHull

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

connections=[[2,4],[50,56],[50,51],[19,20],[2,5],[21,22],[51,52],[48,37],[46,58],[58,44],[46,47],[46,49],[56,2],[51,53],
             [54,57],[57,34],[6,12],[12,11],[12,14],[13,14],[4,3],[23,24],[24,59],[59,26],[59,33],[3,56],[3,8],[4,7],[7,6],[6,5],[5,13],[13,16],
             [16,17], [16,14], [16,17], [17,18], [18,23], [23,22], [18,15],[35,34],[31,30],[39,42], [14,15], [15,19], [19,11], [11,10], [10,9], [9,7], [9,8],
             [8,52], [8,3], [54,53], [53,55], [55,48],
             [49,48], [49,50], [50,51], [50,56],[58,45], [34,21], [21,33], [32,33], [32,38], [36,38], [36,35], [35,55], [21,20], [20,22],
             [17,25], [25,27], [27,28], [28,43], [42,43],[54,52],[57,10], [42,44], [44,45],[47,48],[36,37],[31,32], [45,47], [43,40], [40,29], [40,39],[30,29],[29,28],[40,39],
             [39,37], [39,41], [41,31], [30,26], [26,27], [24,25], [41,38]]
connections=[sorted(x) for x in connections]
connections=np.array(connections)
connections_f=np.unique(connections,axis=0)
np.save("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/connections.npy",connections_f)


# list point coordinates as dict:
points={id:np.array([db.getMarkers(id=id)[0].x,db.getMarkers(id=id)[0].y]) for id in np.unique(connections.flatten())}

# making lines:
lines=[]
for ps in connections_f:

    markers=db.getMarkers(id=list(ps),type="marker")
    p1 = (markers[0].x,markers[0].y)
    p2 = (markers[1].x, markers[1].y)
    lines.append([p1,p2])
lines=np.array(lines)
db.db.close()


# representation as graph: (could be better written like this in the begining)

graph={i:[] for i in np.unique(connections.flatten())}

for i,j in connections:
    graph[i].append(j)
    graph[j].append(i)
for key,value in graph.items():
    graph[key]=np.unique(value)

# finding all connecting paths:
def find_path(graph, start, end, path=[]):
    path = path + [start]
    if start == end:
        return path
    if not start in graph.keys():
        return None

    new_nodes=graph[start]
    np.random.shuffle(new_nodes)
    for node in new_nodes:  # random arangement for now... could also think about shortest path....
        if node not in path:
            newpath = find_path(graph, node, end, path)
            if newpath: return newpath
    return None



def plot_path(path):
    plt.figure()
    plt.imshow(mask)
    for i in range(len(path)-1):
        plt.plot([points[path[i]][0], points[path[i+1]][0]], [points[path[i]][1],points[path[i+1]][1]])

path=find_path(graph,56,28,path=[])
plot_path(path)




from collections import deque, namedtuple
# we'll use infinity as a default distance to nodes.
inf = float('inf')
Edge = namedtuple('Edge', 'start, end, cost')


def make_edge(start, end, cost=1):
  return Edge(start, end, cost)


class Graph:
    def __init__(self, edges):
        # let's check that the data is right
        wrong_edges = [i for i in edges if len(i) not in [2, 3]]
        if wrong_edges:
            raise ValueError('Wrong edges data: {}'.format(wrong_edges))

        self.edges = [make_edge(*edge) for edge in edges]

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

    def add_edge(self, n1, n2, cost=1, both_ends=True):
        node_pairs = self.get_node_pairs(n1, n2, both_ends)
        for edge in self.edges:
            if [edge.start, edge.end] in node_pairs:
                return ValueError('Edge {} {} already exists'.format(n1, n2))

        self.edges.append(Edge(start=n1, end=n2, cost=cost))
        if both_ends:
            self.edges.append(Edge(start=n2, end=n1, cost=cost))

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



### desired data format : node 1 node 2 length of edge
data_list=[]
for p1,ps in graph.items():
    for p2 in ps:
        dist=np.linalg.norm(points[p1]-points[p2])
        data_list.append((p1, p2, dist))

graph = Graph(data_list)

path=graph.dijkstra(56, 28)

plot_path(path)









#### previous try
from collections import defaultdict
class Graph:
    def __init__(self):
        self.nodes = set()
        self.edges = defaultdict(list) # just sets a defaul value
        self.distances = {}

    def add_node(self, value):
        self.nodes.add(value)

    def add_edge(self, from_node, to_node, distance):
        self.edges[from_node].append(to_node)        # generates new key entry if not already existing (standart for all dictionaries)
        #self.edges[to_node].append(from_node) # only useful in directerd graph
        self.distances[(from_node, to_node)] = distance # distances is more relevant part



def dijsktra(graph, initial):
    visited = {initial: 0}
    path = {}

    nodes = set(graph.nodes)

    while nodes: # runs always??
        min_node = None
        for node in nodes:
            if node in visited:
                if min_node is None:
                    min_node = node
                elif visited[node] < visited[min_node]:
                    min_node = node

        if min_node is None:
            break

        nodes.remove(min_node)
        current_weight = visited[min_node]

        for edge in graph.edges[min_node]:
            weight = current_weight + graph.distance[(min_node, edge)]
            if edge not in visited or weight < visited[edge]:
                visited[edge] = weight
                path[edge] = min_node

    return visited, path


graph_c=Graph() # setting up graph object
for p1,ps in graph.items():
    graph_c.add_node(p1)
    for p2 in ps:
        dist=np.linalg.norm(points[p1]-points[p2])
        graph_c.add_edge(p1,p2,dist)




