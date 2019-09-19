# handling cell boundaries and masks by applying graph theory
import numpy as np
from scipy.spatial import cKDTree
from collections import defaultdict
from collections import deque, namedtuple
import copy
import matplotlib.pyplot as plt

def graph_to_mask(graph,points,dims):
    m=np.zeros(dims)
    ps=np.array([y for x in list(graph.values()) for y in x]) #flattening list of point ids
    ps_coord=points[ps] # getting coordinates
    m[ps_coord[:,0],ps_coord[:,1]]=1 #writing points
    return m


def remove_endpoint(graph,ep):
    '''
    recursive function to remove dead ends in a graph starting from point ep. Ep has one neighbour.
    Function stops if it hits a point with 3 neigbours or the removal of a point would cause the appearence of two more
    loose lines.
    :param graph: graph as a dictionary
    :param ep: start point
    :return:
    '''

    #print(graph[ep])
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


def find_line_segement(graph,start,path=[],left_right=0):
    '''
    recursive function, finds path from a point going from the first or second neigbour until
    it reaches an intersection (point with three neighbours)
    :param graph: graph as a dictionary
    :param start: start point
    :param path: path as list of nodes
    :param left_right: define wich neighbour from start to explore (0 or 1
    :return:
    '''
    if len(graph[start]) >2: # stop if intersection (point with 3 neighbours is reached
        return path   # returns the function before next recursion
    ## other wise there is some otherlapp

    path = path + [start]
    new_ps=np.array(graph[start]) ## next points
    if len(path)==1:
        new_p=new_ps[left_right]  #just choose one point
    else:
        new_p = new_ps[new_ps != path[-2]][0]  # next pint that wasnt the previous point

    # recursive function
    newpath = find_line_segement(graph, new_p, path)  # next step
    if newpath:
        return newpath # return if recursion is completed


def mask_to_graph(mask):
    '''
    converts a binary mask to a  graph (dictionary of neighbours)
    Neighbours are identified by cKDTree method
    :param mask:
    :return:
    '''
    graph = defaultdict(list)
    points = np.array(np.where(mask)).T
    point_tree = cKDTree(points)  # look up table for nearest neigbours ??
    d = np.sqrt(2)  # maximal allowd distance
    for i, p in enumerate(points):
        neighbours = point_tree.query_ball_point(p, d)
        neighbours.remove(i)  # removing the point itself from list of its neighbours
        graph[i].extend(neighbours)
    return graph,points


def identify_line_segments(graph,points): ## could be  abit improved , sometimes a few points at the edges are left out...
    '''
        function to identify all line segents(representing individual cell boundaries. Segments are returnt as a dictionary
        with an id as key and a list of points (referring to the points array) that are included in the line. The
        points are in correct order already
    :param graph:
    :param points:
    :return: dictionary with orderd points  in the line
    '''



    lines_dict={}
    n=0 # counter in while loop
    all_points=list(range(len(points)))
    intersect_ps=[key for key,values in graph.items() if len(values)>2] # fining intersection points

    remaining=set(all_points)-set(intersect_ps) # remaining point ids
    while len(remaining)>0: # stop when all points are assigned to intersect or line segment
        start=next(iter(remaining)) # first point of remaining points
        line_seg=list(reversed(find_line_segement(graph,start=start,left_right=0)))[:-1] + find_line_segement(graph,start=start,left_right=1) # find line segment
        remaining-=set(line_seg)# updating remaining  list
        ## maybe reomve this codition??
        if len(line_seg)>1: # avoid single point lines, which are not usefull (cant really get normal vectors from them....)
            lines_dict[n]=line_seg
            n=n+1

        if n>20000:  # expectation if loop should get suck
            raise Exception("found more than 20000 cell borders; something went wrong")


    # plot to confirm correct lines
    #plt.figure()
    #plt.imshow(graph_to_mask(graph,points,mask_boundaries.shape))
    #for seg_ps in lines_dict.values():
    #    for i, p in enumerate(seg_ps):
    #        plt.plot(points[p, 1], points[p, 0], "o")
    #        plt.text(points[p, 1], points[p, 0], str(i))


    return lines_dict

def find_path(graph, start, end, path=[]):
    '''
    finds a path (not the shortest one) through a graph from start to end node
    :param graph:
    :param start:
    :param end:
    :param path:
    :return:
    '''
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
def find_path_circular(graph, start):
    '''
    function to find the order in a circular graph , of one
    :param graph:
    :param start:
    :return:
    '''
    graph2=copy.deepcopy(graph)
    end = graph2[start][0]  # one of the neigbouring points
    graph2[end].remove(start)  # removing connection between them
    graph2[start].remove(end)
    path = find_path(graph2, start, end)
    return path

def points_to_graph(points):
    graph = defaultdict(list)
    point_tree = cKDTree(points)  # look up table for nearest neigbours ??
    d = np.sqrt(2)  # maximal allowd distance
    for i, p in enumerate(points):
        neighbours = point_tree.query_ball_point(p, d)
        neighbours.remove(i)  # removing the point itself from list of its neighbours
        graph[i].extend(neighbours)
    return graph


###### dijacstra algortihm to fuind shortest path
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


def plot_graph(graph,points,mask,number_nodes=False):
    plt.figure()
    plt.imshow(mask)
    for node1,neighbours in graph.items():
        for node2 in neighbours:
            plt.plot([points[node1][1],points[node2][1]],[points[node1][0],points[node2][0]])
    if number_nodes:
        for node in graph.keys():
               plt.text(points[node][1], points[node][0],str(node))


def find_neighbor_lines(graph, start_ps,other_endpoint,own_points,end_points, visited=[],neighbours=[]):
    '''
    recursive function to find neighbouring line. Explores the graph around the endpoint of a line. Notes the id of
    the line if it hits another line. Doesnt explore any points beyond the endpoints of lines.
    it reaches an intersection (point with three neighbours)
    :param graph: graph as a dictionary
    :param  start_ps: start point as a list with len == 1
    :param own_points: all points in own line
    :param end_points_lines: # list of all endpoints
    :param visited: list of visited nodes. Is filled during recursion
    :param neighbours: id of neighbouring lines

    :return: visited: list of visited nodes
    :return: neighbours
    '''
    visited = visited + start_ps # update visited list  ## start must already be a list

    next_ps =   [graph[s] for s in start_ps]  ## next points
    next_ps= [p for ps in next_ps for p in ps] # flatten
    next_ps=list(np.unique(next_ps)) # removing duplication

    # avoid conecting the two endpoints on own line directly
    # only true in first "iteration layer"
    if other_endpoint in  next_ps and len(visited) == 1:
        next_ps.remove(other_endpoint)

    # remove if point is in visited list or in own line
    for p in copy.deepcopy(next_ps): ##### change in the list while iterating is not a nice idea-->
        if p in visited or p in own_points:
            next_ps.remove(p)



    # extract if point can be found in other line
    for p in copy.deepcopy(next_ps): ##### change in the list while iterating is not a nice idea--> make a copy
        if p in end_points:
            next_ps.remove(p)
            neighbours.append(p)

    # use other points for next iteration layer:
    if len(next_ps)==0: # stop recursion if no more next points are left
        return visited,neighbours

    visited,neighbours = find_neighbor_lines(graph, next_ps, other_endpoint, own_points,end_points, visited=visited,neighbours=neighbours)
    # return when iteration is finished
    if visited:
        return visited,neighbours

def find_exact_line_endpoints(lines_points, points, graph):
    '''
    function to find the exact meeting points of lines.
    First find the next closes points on neighbouring lines by exploring the graph. Then calcualtes a
    new endpoint as center of mass of these neighbouring points. Results are stored in a seperate dictionaryl to be
    used in spline interpoltaion
    :param lines_points: dictionary with line_id: list of all points in correct order
    :param points: array of point coordinates
    :param graph: dictionary with connectivity of points
    :return: lines_endpoints_com: dictionary with the line_id:[new endpoint at start, new_endpoint at end]
    '''


    end_points = [[ps[0], ps[-1]] for ps in lines_points.values()]  # all end points in lines
    end_points = [p for ps in end_points for p in ps]


    # finding all neighbouring edpoints for one endpoint of a line
    lines_endpoints = {}
    for line, l_points in lines_points.items():
        # points on the line with out both endpoints, other wise the algorithm can connect two endpoints on the same line
        l_points_core=l_points[1:-1]
        end1 = l_points[0]
        end2 = l_points[-1]
        v, neighbours1 = find_neighbor_lines(graph, [end1], end2, l_points_core, end_points, visited=[], neighbours=[])
        v, neighbours2 = find_neighbor_lines(graph, [end2], end1, l_points_core, end_points, visited=[], neighbours=[])
        lines_endpoints[line] = (neighbours1+[end1], neighbours2+[end2]) # also adding own endpoints here
        # note adding endpoints after find_neighbour_lines is easiest
    ## calcualte new endpoints:
    # just center of mass of the endpoints
    # write to new dictionary and use this in splines calculation, without any points in between
    lines_endpoints_com = {}
    for line, endpoints in lines_endpoints.items():
        com1 = np.mean(points[np.array(endpoints[0])], axis=0)
        com2 = np.mean(points[np.array(endpoints[1])], axis=0)
        # line from new end points to old endpoints:
        lines_endpoints_com[line] = com1, com2

    return lines_endpoints_com,lines_endpoints



class Graph:
# stolen from
# https://dev.to/mxl/dijkstras-algorithm-in-python-algorithms-for-beginners-dkc


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