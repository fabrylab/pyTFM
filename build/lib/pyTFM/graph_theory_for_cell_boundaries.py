# handling cell boundaries and masks by applying graph theory
import copy
from collections import defaultdict
from collections import deque, namedtuple
from itertools import chain

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree

inf = float('inf')
Edge = namedtuple('Edge', 'start, end, cost')


class path_error(Exception):
    pass


class FindingBorderError(Exception):
    pass


def graph_to_mask(graph, points, dims):
    m = np.zeros(dims)
    ps = np.array([y for x in list(graph.values()) for y in x])  # flattening list of point ids
    ps_coord = points[ps]  # getting coordinates
    m[ps_coord[:, 0], ps_coord[:, 1]] = 1  # writing points
    return m


def remove_endpoints_wrapper(graph, points):
    graph_cp = copy.deepcopy(graph)
    eps_id = find_endpoints(graph_cp)  # keys of endpoints, this fast and efficient
    removed = []
    for ep in eps_id:  # removing the endpoint and its connections from the graph
        remove_endpoints(graph_cp, ep, removed)
    return graph_cp, points, removed


def check_connectivity(graph, ep):
    '''
    checking if removing a node from a graph changes the connectivity of the neighbouring nodes.
    in other words: are the neighbouring nodes still connected to one another even if I remove the original node.
    The neighouring points must be connected via maximally one other node. (more would represent ana actual hole in the
    sceletonized mask)
    This classifies loose ends and points that can be removed.
    :return:
    '''
    # this works for points 1,2 or nneighbours
    # identifying an endpoint by checking if removing this point changes the connectivity of its neighbouring nodes
    # i.e. if i break a connection between two points by removing ep
    l1_ps = graph[ep]  # first layer of points
    # check if this causes a line break
    l2_ps = [[pi for pi in graph[p] if pi != ep] for p in l1_ps]  # second layer of points if original point was removed
    # third layer of points // don't need to go deeper due to properties of skeletonize
    # also adding points form secnod layer
    l3_ps = [np.unique(list(chain.from_iterable([graph[p] for p in l2_p] + [l2_p]))) for l2_p in l2_ps]
    # check if all points in l1_ps are connected even if ep is removed
    connectivity = all([all([p in sub_group for p in l1_ps]) for sub_group in l3_ps])
    # check for connection between points in layer 1--> no connection means
    # that removing ep introduces a line break
    return connectivity


def remove_endpoints(graph, ep, removed=[]):
    '''
    recursive function to remove dead ends in a graph starting from point ep. Ep has one neighbour.
    Function stops if it hits a point with 3 neigbours or the removal of a point would cause the appearence of two more
    loose lines.
    :param graph: graph as a dictionary
    :param ep: start point
    :return:
    '''
    connectivity = check_connectivity(graph, ep)
    if connectivity:
        nps = graph[ep]
        remove_point_from_graph(graph, ep)  # removes the point and all connections form the graph
        removed.append(ep)
    else:
        return
    for p in nps:  # iterating through further points
        remove_endpoints(graph, p, removed)

    return


def remove_point_from_graph(graph, point):
    '''
    removes a point and all connections to this point from a graph
    :param graph:
    :param point:
    :return:
    '''
    nps = graph[point]  # neighbouring points/nddes
    graph.pop(point)  # removing the node of the graph
    for p in nps:  # removing all connections to this node
        graph[p].remove(point)


def find_endpoints(graph):
    '''
    identifies "loose ends":
     goes through all points and checks if removing them would introduce a line break.
    this is just as fast as checking for the number  of neighbours and then checking the distance of these neighbours
    :param graph:
    :return:
    '''
    eps = [ep for ep in graph.keys() if check_connectivity(graph, ep)]
    return np.array(eps)


def find_dead_end_lines(graph, non_dead_end_points, max_id):
    '''
    finds dead end line segments from their start to the point where they hit a none dead end line.
    The point in the none dead edn line is included
    :param graph:
    :param non_dead_end_points:
    :param points:
    :return:
    '''

    eps_id = find_endpoints(graph)  # keys of endpoints, this fast and efficient
    dead_end_lines = {}
    for ep in eps_id:
        lps = find_path(graph, ep, non_dead_end_points, path=[])
        non_dead_end_points.extend(lps)  # adding points in the newly discoverd line segments to termination points
        if len(lps) > 3:  # filtering single points and very small bits
            max_id += 1
            dead_end_lines[max_id] = lps
    return dead_end_lines, max_id


def find_lines_simple(graph):
    # find all endpoints
    graph_cp = copy.deepcopy(graph)
    lines_points = {}
    i = 0
    while len(graph_cp.keys()) > 0:
        # first ednpoint, if no endpoint the first point
        new_endpoint = next((x for x in iter(graph_cp.keys()) if len(graph_cp[x]) == 1), next(iter(graph_cp.keys())))
        line = find_path_to_endpoint(graph_cp, new_endpoint, path=[], first=True)  # only explores one direction
        for p in line:
            remove_point_from_graph(graph_cp, p)
        if len(line) > 2:
            lines_points[i] = line
            i += 1
        if i > 10000:
            raise FindingBorderError("found more than 100000 lines; something went wrong")
            break
    return lines_points


def find_path(graph, start, end, path=[]):

    '''
    recursive function
    finds a path (not necessarily the shortest one) through a graph from start to an end node (not necessarily the closest one).

    :param graph: dict, graph
    :param start: int, start point, must be a key in the graph
    :param end: list, list of endpoints. when any endpoint is reach the path search is stopped
    :param path: list, all nodes visited on the way from start to the first endpoint
    :return:
    '''
    path = path + [start]
    if start in end:
        return path
    if not start in graph.keys():
        return None
    for node in graph[start]:
        if node not in path:
            newpath = find_path(graph, node, end, path)
            if newpath:
                return newpath
    return None  # only partial path


def find_path_to_endpoint(graph, start, path=[], first=False):

    '''
    recursive function
    finds a path to a (not specific) point with only one neighbour

    :param graph: dict, graph
    :param start: int, start point, must be a key in the graph
    :param end: list, list of endpoints. when any endpoint is reach the path search is stopped
    :param path: list, all nodes visited on the way from start to the first endpoint
    :return:
    '''
    path = path + [start]
    if len(graph[start]) < 2 and not first:  # stop if we reached a point with only one neighbour
        return path
    if not start in graph.keys():
        return None
    for node in graph[start]:
        if node not in path:
            newpath = find_path_to_endpoint(graph, node, path)
            if newpath:
                return newpath
    return path  # only partial path


def find_line_segement_recursive(graph, start, path=[], left_right=0):
    '''
    ---> would sometimes cause stack overflow/recursion error
    recursive function
    finds path from a point going from the first or second neigbour until
    it reaches an intersection (point with three neighbours)
    :param graph: graph as a dictionary
    :param start: start point
    :param path: path as list of nodes
    :param left_right: define which neighbour from start to explore (0 or 1
    :return:
    '''
    if len(graph[start]) > 2:  # stop if intersection (point with 3 neighbours is reached
        return path  # returns the function before next recursion
    ## other wise there is some otherlapp

    path = path + [start]
    new_ps = np.array(graph[start])  ## next points
    if len(path) == 1:
        new_p = new_ps[left_right]  # just choose one point
    else:
        new_p = new_ps[new_ps != path[-2]][0]  # next point that wasn't the previous point
    # recursive function
    newpath = find_line_segement_recursive(graph, new_p, path)  # next step

    if newpath:
        return newpath  # return if recursion is completed


def find_line_segement(graph, start, path=None, left_right=0):
    '''
    ---> would sometimes cause stack overflow/recursion error
    recursive function
    finds path from a point going from the first or second neigbour until
    it reaches an intersection (point with three neighbours)
    :param graph: graph as a dictionary
    :param start: start point
    :param path: path as list of nodes
    :param left_right: define which neighbour from start to explore (0 or 1
    :return:
    '''

    # first point
    path = []
    path.append(start)
    new_ps = graph[start]
    new_p = new_ps[left_right]  # just choose one point to the left or right
    # break of if we already hit another intersection
    if len(graph[new_p]) > 2:
        return path
    else:
        path.append(new_p)

    while True:  # stop if intersection (point with 3 neighbours) is reached
        new_ps = graph[new_p]
        new_p = [p for p in new_ps if p not in path]
        new_p = new_p[0]
        # check if the point has more (or less) then 2 neighbours.
        # should be more then 2 to indicate intersection
        if len(graph[new_p]) != 2:
            break
        else:
            path.append(new_p)
    return path






def mask_to_graph(mask, d=np.sqrt(2)):
    '''
    converts a binary mask to a  graph (dictionary of neighbours)
    Neighbours are identified by cKDTree method
    :param mask:
    :param d: maximal allowed distance
    :return:
    '''
    graph = defaultdict(list)
    points = np.array(np.where(mask)).T
    point_tree = cKDTree(points)  # look up table for nearest neigbours ??
    for i, p in enumerate(points):
        neighbours = point_tree.query_ball_point(p, d)
        neighbours.remove(i)  # removing the point itself from list of its neighbours
        graph[i].extend(neighbours)
    return graph, points


def identify_line_segments(graph, points):  #
    '''
        function to identify all line segments (representing individual cell boundaries. Segments are returned as a dictionary
        with an id as key and a list of points (referring to the points array) that are included in the line. The
        points are in correct order already
    :param graph:
    :param points:
    :return: dictionary with orderd points  in the line
    '''

    lines_dict = {}
    n = 0  # counter in while loop
    all_points = list(graph.keys())  # all points in the graph
    intersect_ps = [key for key, values in graph.items() if len(values) > 2]  # finding intersection points
    if len(intersect_ps) == 0:
        raise FindingBorderError("Can't identify internal cell borders.")
    remaining = set(all_points) - set(intersect_ps)  # remaining point ids
    while len(remaining) > 0:  # stop when all points are assigned to intersect or line segment
        start = next(iter(remaining))  # first point of remaining points
        # finding a line segment
        line_seg = list(reversed(find_line_segement(graph, start=start, left_right=0)))[:-1] + find_line_segement(graph,
                                                                                                                  start=start,
                                                                                                                  left_right=1)

        remaining -= set(line_seg)  # updating remaining list

        # avoid single point lines, which are not usefull (cant really get normal vectors from them)
        if len(line_seg) > 1:
            lines_dict[n] = line_seg
            n = n + 1

        if n > 20000:  # expectation if loop should get suck
            raise FindingBorderError("found more than 20000 cell borders; something went wrong")

    # plot to confirm correct lines
    # plt.figure()
    # plt.imshow(graph_to_mask(graph,points,mask_boundaries.shape))
    # for seg_ps in lines_dict.values():
    #    for i, p in enumerate(seg_ps):
    #        plt.plot(points[p, 1], points[p, 0], "o")
    #        plt.text(points[p, 1], points[p, 0], str(i))

    return lines_dict


def find_path_circular(graph, start):
    '''
    function to find the order in a circular graph , of one
    :param graph:
    :param start:
    :return:
    '''
    graph2 = copy.deepcopy(graph)
    end = graph2[start][0]  # one of the neigbouring points
    graph2[end].remove(start)  # removing connection between them
    graph2[start].remove(end)
    path = find_path(graph2, start, [end])
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


def make_edge(start, end, cost=1):
    return Edge(start, end, cost)


def make_neighbours_dict(data_list):
    neighbours_dict = defaultdict(list)
    for i, j, dist in data_list:
        neighbours_dict[i].append(j)
        neighbours_dict[j].append(i)
    for key, value in neighbours_dict.items():
        neighbours_dict[key] = list(set(value))
    return neighbours_dict


def plot_graph(graph, points, mask, number_nodes=False):
    plt.figure()
    plt.imshow(mask)
    for node1, neighbours in graph.items():
        for node2 in neighbours:
            plt.plot([points[node1][1], points[node2][1]], [points[node1][0], points[node2][0]])
    if number_nodes:
        for node in graph.keys():
            plt.text(points[node][1], points[node][0], str(node))


def find_neighbor_lines(graph, start_ps, other_endpoint, own_points, end_points, visited=[], neighbours=[]):
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
    visited = visited + start_ps  # update visited list  ## start must already be a list

    next_ps = [graph[s] for s in start_ps]  ## next points
    next_ps = [p for ps in next_ps for p in ps]  # flatten
    next_ps = list(np.unique(next_ps))  # removing duplication

    # avoid connecting the two endpoints on own line directly
    # only true in first "iteration layer"
    if other_endpoint in next_ps and len(visited) == 1:
        next_ps.remove(other_endpoint)

    # remove if point is in visited list or in own line
    for p in copy.deepcopy(next_ps):  ##### change in the list while iterating is not a nice idea-->
        if p in visited or p in own_points:
            next_ps.remove(p)

    # extract if point can be found in other line
    for p in copy.deepcopy(next_ps):  ##### change in the list while iterating is not a nice idea--> make a copy
        if p in end_points:
            next_ps.remove(p)
            neighbours.append(p)

    # use other points for next iteration layer:
    if len(next_ps) == 0:  # stop recursion if no more next points are left
        return visited, neighbours

    visited, neighbours = find_neighbor_lines(graph, next_ps, other_endpoint, own_points, end_points, visited=visited,
                                              neighbours=neighbours)
    # return when iteration is finished
    if visited:
        return visited, neighbours


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
        l_points_core = l_points[1:-1]
        end1 = l_points[0]
        end2 = l_points[-1]
        v, neighbours1 = find_neighbor_lines(graph, [end1], end2, l_points_core, end_points, visited=[], neighbours=[])
        v, neighbours2 = find_neighbor_lines(graph, [end2], end1, l_points_core, end_points, visited=[], neighbours=[])
        lines_endpoints[line] = (neighbours1 + [end1], neighbours2 + [end2])  # also adding own endpoints here
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

    return lines_endpoints_com, lines_endpoints


class Graph:
    # stolen from
    # https://dev.to/mxl/dijkstras-algorithm-in-python-algorithms-for-beginners-dkc

    def __init__(self, edges):
        # let's check that the data is right
        wrong_edges = [i for i in edges if len(i) not in [2, 3]]
        if wrong_edges:
            raise ValueError('Wrong edges data: {}'.format(wrong_edges))

        self.edges = [make_edge(*edge) for edge in edges]
        self.neighbours_dict = make_neighbours_dict(edges)  # already because of @prperty solved

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
        if n1 in self.neighbours_dict.keys():
            if n2 in self.neighbours_dict[
                n1]:  # removing if existingg, necessary because of biothe ends, might already be removed
                self.neighbours_dict[n1].remove(n2)
            if len(self.neighbours_dict[n1]) == 0:
                self.neighbours_dict.pop(n1)
        if both_ends:
            if n2 in self.neighbours_dict.keys():
                if n1 in self.neighbours_dict[n2]:
                    self.neighbours_dict[n2].remove(n1)
                if len(self.neighbours_dict[n2]) == 0:
                    self.neighbours_dict.pop(n2)
            # removing empty keys
        # self.neighbours_dict={ key : values for key,values in self.neighbours_dict.items() if len(values)>0}

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
        new_nodes = self.neighbours_dict[source]
        np.random.shuffle(new_nodes)
        for node in new_nodes:  # dest not needed
            if node not in path:
                newpath = self.random_path(node, dest, path)
                if newpath: return newpath
        return None
