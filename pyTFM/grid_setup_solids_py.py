from collections import Counter
from contextlib import suppress
from itertools import chain

import solidspy.assemutil as ass
import solidspy.postprocesor as pos
import solidspy.solutil as sol
from pyTFM.stress_functions import *
from pyTFM.utilities_TFM import make_random_discrete_color_range, join_dictionary
from scipy.interpolate import splprep, splev
from scipy.ndimage import binary_closing as binary_clo
from scipy.ndimage.measurements import find_objects
from scipy.ndimage.morphology import binary_fill_holes
from scipy.optimize import least_squares
from scipy.signal import convolve2d
from scipy.sparse.csr import csr_matrix
from scipy.sparse.linalg import lsqr
from skimage.measure import regionprops
from skimage.morphology import skeletonize, remove_small_holes, remove_small_objects, label, binary_dilation


def show_points(ps, mask):
    plt.figure()
    plt.imshow(mask)
    plt.plot(ps[:, 1], ps[:, 0], "or")


def identify_cells(mask_area, mask_boundaries, points):
    '''
    function to identfy cells. Each cell is is a dictionary entry with a list of ids, reffering to
    points.
    :param mask:
    :param area:
    :param mask_boundaries:
    :return:
    '''

    cells = {}  # dictionary containg a list of point idsthat sourround each cell
    cells_area = {}  # dictionary containg a all pixels belonging to that cell as boolean aray
    # points_to_flatt array map:
    # labeling each cell
    m = mask_area.astype(int) - mask_boundaries.astype(int)
    ml = label(m, connectivity=1)
    # creating a list of point coordinates corresponding to a flattend array
    # this will allow easier identification of the id of a point
    points_fl = (points[:, 0]) + mask_area.shape[0] * (points[:, 1])
    sort_ids = np.argsort(points_fl)  # sorting will allow np.searchsorted function;
    points_fl = points_fl[sort_ids]  #

    for i, l in enumerate(
            np.unique(ml)[1:]):  # getting each cell border by binary dilation of one pixel; iterating over each cell
        m_part = (ml == l).astype(bool)  # extracting a cell area
        edge = np.logical_and(binary_dilation(m_part), ~m_part)  # getting the boundary of a cell
        ps = np.array(np.where(edge)).T  # finding coordinates
        ps_fl = (ps[:, 0]) + mask_area.shape[0] * (ps[:, 1])  # convert coordinates to the one of a flat array
        p_ids = sort_ids[np.searchsorted(points_fl,
                                         ps_fl)]  # find indices, where i would need to insert,supposed to be the fastest way
        #  and read index from unsorted list
        cells[i] = p_ids  # save to dictionary
        cells_area[i] = m_part

    ## vizualization
    # creating a(random lsit of hex colors)
    # colors = []
    # for i in range(len(cells.items())):
    #    colors.append('#%06X' % np.random.randint(0, 0xFFFFFF))
    # plt.figure()
    # plt.imshow(mask_boundaries)
    # offset = 0.01
    # for cell, ps in cells.items():
    #    for p in ps:
    #        plt.plot(points[p][1] + offset * cell, points[p][0] + offset * cell, "o", color=colors[cell])

    return cells, cells_area


def spline_interpolation(line, points, k=3, endpoints=None):
    '''
    function takes points from a line and uses spline interpolation to create a smooth representation.

    :param line: list point ids that define a line. Must be incorrect order
    :param points: List of all points. This list is wher the point ids correspond to. yx order
    :param endpoints: additional pair of endpoints in yx order, to include in the spline interpolation
    :return: tck, array defining the spline, with knot position, parameter values and order
            points_new, new position of the points on the line according to spline interploation in x y coorinates !!!
    '''

    x = points[line, 1]  # list x and y coordinate for points on line
    y = points[line, 0]

    if endpoints:  # adding endpoints if provided
        # splprep throws error if two points are exactely identidcal
        # for ecxample for the "dead-end-lines" the endpoint is identical to the first and last point and thus shouldnt
        # be added here
        if not (endpoints[0][1] == x[0] and endpoints[0][0] == y[0]):
            x = np.concatenate([[endpoints[0][1]], x], axis=0)
            y = np.concatenate([[endpoints[0][0]], y], axis=0)
        if not (endpoints[1][1] == x[-1] and endpoints[1][0] == y[-1]):
            x = np.concatenate([x, [endpoints[1][1]]], axis=0)
            y = np.concatenate([y, [endpoints[1][0]]], axis=0)

    k = len(x) - 1 if len(x) <= 3 else 3  # addapt spline order, according to number of points
    tck, u = splprep([x, y], s=10, k=k)  ### parametric spline interpolation
    # fits essentially function: [x,y] =f(t) , t(paramter) is default np.linspace(0,1,len(x))
    # tck is array with x_position of knot, y_position of knot, parameters of the plne, order of the spline
    # ui is paramter given to x,y points, in this case default (values from 0 to 1)
    # s is smoothing factor(?)
    # k is order of the spline, cubic is fine

    if endpoints:  # forcing splines to through endpoints
        tck[1][0][0] = endpoints[0][1]
        tck[1][0][-1] = endpoints[1][1]
        tck[1][1][0] = endpoints[0][0]
        tck[1][1][-1] = endpoints[1][0]
    # new points from spline interpolation (thes points will be used for normal stress vector calculation
    points_new = np.array(splev(u, tck, der=0)).T
    # in xy orientation
    # points_new = np.round(points_new).astype(int)  ## could also use exact values and the interpolate the stress tensor, but thats probably to complicated
    return tck, u, points_new  # points new in xy orientation


def arrange_lines_from_endpoints(cells_lines, lines_endpoints_com):
    '''
    rearranging the order of lines in the line_id list for one one cell. The endpoints for all lines for ne cell are
    extracted and then puzzled together. Each endpoint must occure in two lines.

    :param cells_lines: dictionary cell_id:[line_ids]
    :param lines_endpoints_com: dictionary line_id:[endpoint1,endpoint2] , each endoinpt is an array with x and y coordinates
    :return: cells_lines_new: updated cells_lines dictionary

    '''
    cells_lines_new = {}

    for cell_id, line_ids in cells_lines.items():
        # extracting relevant endpoints
        local_endpoints = {line_id: lines_endpoints_com[line_id] for line_id in line_ids}
        # rearranging endpoints into an suitabel array: axis0: lines axis 1:[endpoint1, endpoint2], axis3: x,y coordinates
        eps = np.array([np.array([value[0], value[1]]) for value in local_endpoints.values()])

        new_line_ids = []  # newly_arranged line_ids
        p_ind1 = 0  # start ids
        p_ind2 = 0
        # iterating through eps, until all endpoints have been visited
        while not np.isnan(eps).all():
            point = copy.deepcopy(eps[p_ind1, p_ind2])  # extracting an endpoint
            eps[p_ind1, p_ind2] = np.nan  # remove the point from array
            # find second occurrence, by taking the norm of the diffrence between the current point and all other points
            # this should be zero
            np_ind1, np_ind2 = np.array(np.where(np.linalg.norm(eps - point, axis=2) == 0)).T[0]
            new_line_ids.append(line_ids[np_ind1])  # note corresponding line_ids
            eps[np_ind1, np_ind2] = np.nan  # remove this point from array
            p_ind1, p_ind2 = np_ind1, np.abs(np_ind2 - 1)  # change to the other end of the line

        cells_lines_new[cell_id] = new_line_ids  # update dictionary

    return cells_lines_new


def find_edge_lines(cells_lines):
    '''
    Finding all lines (cell borders) at the edge of a cell colony. Simply checks if
    a line is associated to only one cell.
    :param cells_lines: dictionary with cell_id:[associated line_ids]
    :return: edge_lines: lsit of line ids at the edge of the cell colony
    '''
    all_lines = np.array(list(chain.from_iterable(cells_lines.values())))  # unpacking all line ids
    counts = Counter(all_lines)  # counting occurences
    edge_lines = [line for line, count in counts.items() if
                  count == 1]  # select if line_id was only associated to one cell
    return edge_lines


def center_of_mass_cells(cells_points, points):
    '''
    calulating the "center of mass" of a cell using only the points at the edge of the cell
    :param cells_points:
    :param points:
    :return:
    '''
    cells_com = {}
    for cell_id, hull_points in cells_points.items():
        cells_com[cell_id] = np.mean(points[hull_points], axis=0)

    return cells_com


def remove_circular_line(allLines_points, lines_endpoints_com, lines_points, lines_endpoints):
    '''
    finds lines that are circular by checking if the first and second endpoint are identical. The lines are
    delted from all input dictionaries
    :param lines_endpoints_com:
    :param lines_points:
    :param lines_endpoints:
    :return:
    '''
    # finding all lines where first and second endpoint is identical
    circular = [l_id for l_id, endpoints in lines_endpoints_com.items() if
                np.linalg.norm(endpoints[0] - endpoints[1]) == 0]
    # print(circular)
    # clearing these lines from the input dictionaries
    for l_id in circular:
        del lines_endpoints_com[l_id]
        del lines_points[l_id]
        del lines_endpoints[l_id]
        del allLines_points[l_id]


def interpolate_cell_area(cells_area, shape):
    cells_area_interpol = {}
    for cell_id, areas in cells_area.items():
        cells_area_interpol[cell_id] = interpolation(areas, shape, min_cell_size=0)
    return cells_area_interpol


def interpolate_points_dict(points_dict, shape_target, shape_orgin):
    points_dict_interp = {}
    for p_id, coords in points_dict.items():
        points_dict_interp[p_id] = (interpolation_single_point(coords[0], shape_target, shape_orgin),
                                    interpolation_single_point(coords[1], shape_target, shape_orgin))
    return points_dict_interp


class Cells_and_Lines:
    # container for cells and lines, assignement of points with them, and assignement with each other
    def __init__(self, mask_boundaries, shape, graph, points):

        # masks, graph and points including dead-end lines // this distinction is mostly due to historic reasons
        self.mask_boundaries_wp = mask_boundaries
        self.inter_shape = shape
        # graph as a dictionary with key = point id, values: ids of neighbouring points
        # any point id is the index in the points array (contains coordinate of these points
        self.graph_wp = graph
        self.points_wp = points
        self.graph, self.points, removed = remove_endpoints_wrapper(self.graph_wp, self.points_wp)
        # masks, graph and points excluding dead-end lines
        self.mask_boundaries = graph_to_mask(self.graph, self.points, mask_boundaries.shape)  # rebuilding the mask

        # interpolate points to the size of the future FEM_grid
        self.points_interpol = interpolation_single_point(self.points, self.inter_shape, self.mask_boundaries.shape)
        # interpolation factors used in the fun
        # self.inerpol_factors=np.array([self.inter_shape[0] /self.mask_boundaries.shape[0], self.inter_shape[1] / self.mask_boundaries.shape[1]])
        # points as dictionary with key=points id, values: points coordinates
        self.points_dict = {i: self.points[i] for i in range(len(self.points))}

        # lines as a dictionary with key=line id, values: ids of containing points (in correct order)
        self.lines_points = identify_line_segments(self.graph, self.points_interpol)

        # cells as a dictionary with key=cell id, values: ids of containing points (not ordered)
        self.cells_points, self.cells_area = identify_cells(self.mask_boundaries,
                                                            binary_fill_holes(self.mask_boundaries), self.points)

        # interpolate the area of individual cells to the size of deformation
        self.cells_area_interpol = interpolate_cell_area(self.cells_area, self.inter_shape)

        # self.points_lines = invert_dictionary(self.lines_points) # point_id:line_id
        self.max_line_id = np.max(list(self.lines_points.keys()))
        self.de_lines_points, self.max_line_id = find_dead_end_lines(self.graph_wp, list(self.graph.keys()),
                                                                     self.max_line_id)
        self.de_endpoints = {key: (self.points[value[0]], self.points[value[-1]]) for key, value in
                             self.de_lines_points.items()}  # using exact endpoints for the dead end lines
        self.allLines_points = join_dictionary(self.lines_points, self.de_lines_points)

        # dictionary with endpoints, needed to completely fill the gaps between all cell_lines
        self.lines_endpoints_com, self.lines_endpoints = find_exact_line_endpoints(self.lines_points, self.points,
                                                                                   self.graph)

        #self.simple_line_plotting(self.allLines_points, subset=np.inf)
        #for p in self.lines_endpoints_com.values():
        #    plt.plot(p[0][1],p[0][0],"o")
        #    plt.plot(p[1][1], p[1][0], "o")

        # removing all lines that are predicted to be circular. Mostly a problem for very short lines
        remove_circular_line(self.allLines_points, self.lines_endpoints_com, self.lines_points, self.lines_endpoints)

        # center of mass of cells, calculated only from the hull points
        self.cells_com = center_of_mass_cells(self.cells_points, self.points)
        # dictionary to associate cells with correct lines, key is cell id, value is line_id
        self.cells_lines = defaultdict(list)
        # dictionary to associate cells with correct lines, key is line id, value is cell

        self.lines_cells = defaultdict(list)
        for l_id, l in self.lines_points.items():
            for c_id, c in self.cells_points.items():
                if l[int(len(l)/2)] in c:
                    self.cells_lines[c_id].append(l_id)
                    self.lines_cells[l_id].append(c_id)
        # using the new endpoints to arrange lines in the correct way
        self.cells_lines = arrange_lines_from_endpoints(self.cells_lines, self.lines_endpoints_com)

        # adding dead end endpoints only now to avoid complications when identifying cells
        self.de_endpoints = {key: (self.points[value[0]], self.points[value[-1]]) for key, value in
                             self.de_lines_points.items()}
        self.lines_endpoints_com = join_dictionary(self.lines_endpoints_com, self.de_endpoints)
        self.lines_endpoints_interpol = interpolate_points_dict(self.lines_endpoints_com, self.inter_shape,
                                                                self.mask_boundaries.shape)

        # list of ids
        self.point_ids = list(self.points_dict.keys())
        self.cell_ids = list(self.cells_points.keys())
        self.line_ids = list(self.allLines_points.keys())
        # list of all lines at the edge of the cell colony
        self.edge_lines = find_edge_lines(self.cells_lines)
        # list of dead end lines
        self.dead_end_lines = list(self.de_lines_points.keys())
        # list of central boundary (non dead-end, non-edge lines)
        self.central_lines = [line_id for line_id in self.allLines_points.keys() if
                              line_id not in self.edge_lines and line_id not in self.dead_end_lines]
        self.n_cells = len(self.cell_ids)
        # list with original line lengths--> later used for interpolation
        self.line_lengths = {key: len(value) for key, value in self.allLines_points.items()}

        # dictionary containing the spline represetation of the points as a parametric function
        # [x,y]=f(u). u is always np.linspace(0,1,"number of points in the line). Use scipy.interpolate.splev
        # to evaluate at other points
        self.lines_splines = defaultdict(list)

        # dictionary with the normal vectors for each line acording to the spline interpolation. Contains list
        # of the normal vectors , position of these points is listed in lines_spline_points
        # interpolation !!! positions are given in xy order !!!
        self.lines_n_vectors = defaultdict(list)

        # dictionary with example points (where a normal vector originates) from spline interpolation
        self.lines_spline_points = defaultdict(list)

        for line_id, line_ps in self.allLines_points.items():
            endpoints = self.lines_endpoints_interpol[line_id]
            k = len(line_ps) + 2 - 1 if len(line_ps) <= 3 else 3  # addapt spline order, according to number of points
            # spline order must be below nomber of points, so choose 2 if lne has one point + 2 endpoints
            tck, u, points_new = spline_interpolation(line_ps, self.points_interpol, k=k,
                                                      endpoints=endpoints)  # spline interpolation
            self.lines_splines[line_id] = tck  # saving the spline object
            # saving a few oints and n vectors for easy representations/ debugging,
            # these points will not be used further
            n_vectors = normal_vectors_from_splines(u, tck)  # calculating normal vectors as a list
            self.lines_n_vectors[line_id] = n_vectors
            self.lines_spline_points[line_id] = points_new

    def cut_to_FEM_grid(self, FEM_mask):
        '''
        removing any points that lie outside the FEM Grid
        :param FEM_mask:
        :return:
        '''
        self.lines_outside = []
        # spline points is already interpolated to shape of self.mask_area
        for l_id, line_points in list(self.lines_spline_points.items()):
            line_points = np.round(line_points).astype(int)
            if np.sum(~FEM_mask[line_points[:, 1], line_points[:, 0]]) > 5:
                self.lines_spline_points.pop(l_id, None)
                self.allLines_points.pop(l_id, None)
                self.line_lengths.pop(l_id, None)
                self.lines_endpoints_com.pop(l_id, None)
                self.de_lines_points.pop(l_id, None)
                self.lines_points.pop(l_id, None)
                with suppress(ValueError): self.edge_lines.remove(l_id)
                with suppress(ValueError): self.edge_lines.remove(l_id)
                with suppress(ValueError): self.edge_lines.remove(l_id)
                self.lines_outside.append(l_id)

    def filter_small_de_line(self, min_length):
        for l_id in copy.deepcopy(self.dead_end_lines):  # does not filter small line segments around cells -
            if self.line_lengths[l_id] < min_length:
                with suppress(AttributeError): self.allLines_points.pop(l_id, None)
                with suppress(AttributeError): self.lines_endpoints_com.pop(l_id, None)
                with suppress(AttributeError): self.lines_endpoints_interpol.pop(l_id, None)
                with suppress(AttributeError): self.lines_splines.pop(l_id, None)
                with suppress(AttributeError): self.lines_n_vectors.pop(l_id, None)
                with suppress(AttributeError):  self.lines_spline_points.pop(l_id, None)
                with suppress(AttributeError):  self.line_lengths.pop(l_id, None)
                with suppress(AttributeError): self.de_endpoints.pop(l_id, None)
                with suppress(AttributeError): self.de_lines_points.pop(l_id, None)

                with suppress(ValueError, AttributeError): self.line_ids.remove(l_id)  # list
                with suppress(ValueError, AttributeError): self.dead_end_lines.remove(l_id)

    def return_n_array(self, fill_nan=True):
        '''
        writes (normal) vectors in a two dimensional array according to their position from spline interpolation
        :return:
        '''
        if fill_nan:
            n_array = np.zeros((self.mask_boundaries.shape[0], self.mask_boundaries.shape[1], 2)) + np.nan
        else:
            n_array = np.zeros((self.mask_boundaries.shape[0], self.mask_boundaries.shape[1], 2))
        for vecs, ps in zip(self.lines_n_vectors.values(), self.lines_spline_points.values()):
            n_array[ps[:, 1], ps[:, 0]] = vecs
        return n_array

    def vizualize_lines_and_cells(self, sample_factor=1, plot_n_vectors=False):
        '''
        plotting the id of lines for each point, after spline interpolation,
        lines at the edge have a different color
        plotting the id of cells
        plotting_normal_vectors, as predcited by intial spline interpolation
        :param cells_points: factor to reduce line_id labels that are plotted, must be <= 1
        :return:
        '''

        offset = 0.005  # ofset for point id text
        fig = plt.figure()
        plt.imshow(self.mask_boundaries, cmap="jet")
        plt.plot([], [], color="green", label="line_ids")  # legend for line ids
        plt.plot([], [], color="red", label="cells_ids")  # legend for cells ids
        # plotting all points with line id as label
        colors = np.array(["C1", "C2", "C3"])  # choosing colors according to line type
        line_classifier = [self.central_lines, self.edge_lines, self.dead_end_lines]
        for l, ps in self.allLines_points.items():
            ps = np.array(ps)
            color = colors[np.array([l in line_ids for line_ids in line_classifier])][0]
            p_indices = np.array(list(range(len(ps))))  # all indicces
            # randomly choosing a few indices, without first and laast index
            ps_select = np.random.choice(p_indices[1:-1], size=int((len(ps) - 2) * sample_factor), replace=False)
            ps_select = np.append(ps_select, p_indices[np.array([0, -1])])  # adding back first and last index

            plt.plot(self.points[ps][:, 1], self.points[ps][:, 0], "o", color=color)
            for p in ps[ps_select]:  # labeling selected points
                plt.text(self.points[p][1] + 1 * offset * l, self.points[p][0] + 1 * offset * l, s=str(l),
                         color="green")
        # plotting cel id at center of mass of cell
        for cell_id, com in self.cells_com.items():
            plt.text(com[1], com[0], str(cell_id), color="red", fontsize=13)

        if plot_n_vectors:
            for points, n_vectors in zip(self.lines_spline_points.values(), self.lines_n_vectors.values()):
                for n_vec, p in zip(n_vectors,
                                    interpolation_single_point(points, self.mask_boundaries.shape, self.inter_shape)):
                    plt.arrow(p[0], p[1], n_vec[0], n_vec[1], head_width=0.15)
        plt.legend()
        return fig

    def vizualize_splines(self, sample_factor=1, subset=np.inf):
        # sample factor: only every nth point is used for plotting
        plt.figure()
        plt.imshow(self.mask_boundaries)
        colors = make_random_discrete_color_range(len(self.lines_spline_points.keys()))
        for i, (l_id, points) in enumerate(self.lines_spline_points.items()):
            if i > subset:
                break
            points = interpolation_single_point(points, self.mask_boundaries.shape, self.inter_shape)
            plt.plot(points[::sample_factor, 0], points[::sample_factor, 1], color=colors[i])

    def simple_line_plotting(self, lines, subset=np.inf):
        plt.figure()
        plt.imshow(self.mask_boundaries)
        for i, ps in enumerate(lines.values()):
            if i > subset:
                break
            ps2 = self.points[ps]
            if len(ps) > 1:
                plt.plot(ps2[:, 1], ps2[:, 0])
            else:
                plt.plot(ps2[:, 1], ps2[:, 0], "o")


class Cells_and_Lines2(Cells_and_Lines):
    # much simplyfied version. Doesnt find cells and doesn't guarnatee to find long line segments. Can't identify
    # lines at the edge of a colony

    def __init__(self, mask_boundaries, shape, graph, points):
        self.mask_boundaries = mask_boundaries
        self.inter_shape = shape
        # graph as a dictionary with key=point id, values: ids of neighbouring points
        # any point id is the index in the points array (contains coordinate of these points
        self.graph = graph
        self.points = points
        # finding line segments
        self.lines_points = find_lines_simple(self.graph)
        self.points_interpol = interpolation_single_point(self.points, self.inter_shape, self.mask_boundaries.shape)

        self.lines_splines = defaultdict(list)
        self.lines_n_vectors = defaultdict(list)
        self.lines_spline_points = defaultdict(list)

        # spline interpolation
        for line_id, line_ps in self.lines_points.items():
            # spline order must be below nomber of points, so choose 2 if lne has one point + 2 endpoints
            tck, u, points_new = spline_interpolation(line_ps, self.points_interpol)  # spline interpolation
            self.lines_splines[line_id] = tck  # saving the spline object
            # saving a few oints and n vectors for easy representations/ debugging,
            # these points will not be used further
            n_vectors = normal_vectors_from_splines(u, tck)  # calculating normal vectors as a list
            self.lines_n_vectors[line_id] = n_vectors
            self.lines_spline_points[line_id] = points_new

        # categories from other lines object
        self.cell_ids = []
        self.line_ids = list(self.lines_points.keys())
        self.edge_lines = []
        self.dead_end_lines = []
        self.central_lines = self.line_ids
        self.line_lengths = {key: len(value) for key, value in self.lines_points.items()}

        # very rough estimate of cell number
        label_mask, self.n_cells = label(~mask_boundaries, connectivity=1, return_num=True)


def prepare_mask_FEM(mask, shape):
    '''
    this function usese skeletonize to transform the cellboundraies to one pixel width. Loose ends are trimmed by
    converison to a graph and then deleting all nodes with only one neighbour.
    :param mask:
    :param min_cell_size: minimal sze of cells in pixles, any hole below that will be filled up. If none, then
        some estimated value is used
    :return: mask of type bool
    '''
    mask = remove_small_objects(mask.astype(bool), 200).astype(bool)  # removing other small bits
    # interpolating the area to the size of future FEM-grd
    mask_int = interpolation(binary_fill_holes(mask), shape)
    mask_int = binary_fill_holes(mask_int).astype(bool)  # sometime first binary fill holes is not enough

    # removing unsuitable pixels for grid later --> a pixel must be connected to at least more then 2 others
    m = convolve2d(mask_int.astype(int), np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]]), mode="same", boundary="fill",
                   fillvalue=0)  # convoultion will produce 1 if only one direct (distance 1) neighbour exists
    p = np.where(np.logical_and(m == 1, mask_int))  # problematic coordinates

    for x, y in zip(p[0], p[1]):
        new_ps = np.array([[x + 1, y], [x, y + 1], [x - 1, y], [x, y - 1]])  # all relevant neigbouring points
        new_p = new_ps[mask_int[new_ps[:, 0], new_ps[:, 1]]][
            0]  # checking where a possible replacement point could be located this can only have one result
        mask_int[x, y] = False  # removing old point
        mask_int[new_p[0], new_p[1]] = True  # adding new point
    return mask_int


def find_borders(mask, shape, raise_error=True, type="colony", min_length=0):
    #### maybe reintroduce small cell filter
    # removing small bits
    mask = remove_small_objects(mask.astype(bool), 1000).astype(bool)
    # generating borders
    mask_boundaries = skeletonize(mask.astype(int))
    # converting mask to graph object
    graph, points = mask_to_graph(mask_boundaries)
    # finding dead ends: cell borders which don't connect to other cell borders at one end:
    # this is use full to get a clean structure of cell borders, which is later used to identifying the number and area of cells
    # applying remove endpoints multiple times to deal with forking dead ends

    try:
        if type == "colony":
            c_l = Cells_and_Lines(mask_boundaries, shape, graph, points)
            c_l.filter_small_de_line(min_length)
        if type == "cell layer":
            c_l = Cells_and_Lines2(mask_boundaries, shape, graph, points)
            c_l.filter_small_de_line(min_length)
    except (RecursionError, FindingBorderError, IndexError, KeyError) as e:
        print("original error: ", e)
        if raise_error:
            raise FindingBorderError
        else:
            return None

    # c_l.cut_to_FEM_grid(mask_int) # this should never be necessary
    # c_l.vizualize_lines_and_cells(sample_factor=0.2,plot_n_vectors=True)
    # c_l.vizualize_splines(sample_factor=4,subset=200000)
    # c_l.simple_line_plotting(c_l.lines_points)
    ## vizualization of spline interpolation with new line enpoints
    return c_l


def interpolation(mask, dims, min_cell_size=100, dtype=bool):
    #
    # some pre clean up of the mask
    mask = remove_small_holes(mask.astype(bool), min_cell_size)
    mask = remove_small_objects(mask.astype(bool), 1000)  # removing other small bits
    # note: remove_small_objects labels automatically if mask is bool
    coords = np.array(np.where(mask)).astype(float)  # coordinates of all points
    interpol_factors = np.array([dims[0] / mask.shape[0], dims[1] / mask.shape[1]])
    coords[0] = coords[0] * interpol_factors[0]  # interpolating x coordinates
    coords[1] = coords[1] * interpol_factors[1]  # interpolating xy coordinates
    coords = np.round(coords).astype(int)

    coords[0, coords[0] >= dims[0]] = dims[0] - 1  # fixing issue when interpolated object is just at the image border
    coords[1, coords[1] >= dims[1]] = dims[1] - 1

    mask_int = np.zeros(dims)
    mask_int[coords[0], coords[1]] = 1
    mask_int = mask_int.astype(int)
    # filling gaps if we interpolate upwards
    if dims[0] * dims[1] >= mask.shape[0] * mask.shape[1]:
        iter = int(np.ceil(np.max([mask.shape[0] / dims[0], mask.shape[0] / dims[0]])) * 5)  # times 5 is safety factor
        mask_int = binary_clo(mask_int, iterations=10)
        print(iter)
    return mask_int.astype(bool)


def interpolation_single_point(point, shape_target, shape_origin):
    # is also works with 2d arrays of shape(n,2)
    interpol_factors = np.array([shape_target[0] / shape_origin[0], shape_target[1] / shape_origin[1]])
    point_interp = point * interpol_factors
    return point_interp


def alligne_objects(mask1, mask2):
    '''
    function to cut the object from mask2 and set it in an array of mask1.shape, so that the center of the bounding
    box of the object in mask1 and the new array are the same
    :param mask1: array deetermining the output shape and center of the binding box f the output object
    :param mask2: array from which an object is cut out
    :return: mas2_alligne: output array
    '''
    rectangle1 = find_objects(mask1, 1)
    lengths1 = [rectangle1[0][0].stop - rectangle1[0][0].start,
                rectangle1[0][1].stop - rectangle1[0][1].start]  # ä side lengths of hte recangle that was cut out
    rect_center1 = [rectangle1[0][0].start + lengths1[0] / 2, rectangle1[0][1].start + lengths1[1] / 2]

    rectangle2 = find_objects(mask2, 1)
    lengths2 = [rectangle2[0][0].stop - rectangle2[0][0].start,
                rectangle2[0][1].stop - rectangle2[0][1].start]  # ä side lengths of hte recangle that was cut out

    # set in around com of fl image:
    mask2_alligne = np.zeros(np.shape(mask1))

    new_cords1 = [int(rect_center1[0] - lengths2[0] / 2),
                  int(rect_center1[0] + lengths2[0] / 2),
                  int(rect_center1[1] - lengths2[1] / 2),
                  int(rect_center1[1] + lengths2[1] / 2)
                  ]
    mask2_alligne[new_cords1[0]: new_cords1[1], new_cords1[2]:new_cords1[3]] = mask2[rectangle2[0][0], rectangle2[0][1]]
    return mask2_alligne.astype(int)


def cut_mask_from_edge(mask, cut_factor, warn_flag=False, fill=True):
    sum_mask1 = np.sum(mask)
    dims = mask.shape
    inds = [int(dims[0] * cut_factor), int(dims[0] - (dims[0] * cut_factor)), int(dims[1] * cut_factor),
            int(dims[1] - (dims[1] * cut_factor))]
    if fill:  # filling to the original shape
        mask_cut = copy.deepcopy(mask)
        mask_cut[:inds[0], :] = 0
        mask_cut[inds[1]:, :] = 0
        mask_cut[:, :inds[2]] = 0
        mask_cut[:, inds[3]:] = 0
    else:  # new array with new shape
        mask_cut = np.zeros((inds[1] - inds[0], inds[3] - inds[2]))
        mask_cut = mask[inds[0]:inds[1], inds[2]:inds[3]]

    sum_mask2 = np.sum(mask_cut)
    warn = "mask was cut close to image edge" if (sum_mask2 < sum_mask1 and warn_flag) else ""
    return mask_cut, warn


def cut_mask_from_edge_wrapper(cut_factor, mask, parameter_dict, cut=True, warn=""):
    if cut:
        mask, warn = cut_mask_from_edge(mask, cut_factor, parameter_dict["TFM_mode"] == "colony")
    return mask, warn


def FEM_simulation(nodes, elements, loads, mats, mask_area, verbose=False, **kwargs):
    DME, IBC, neq = ass.DME(nodes, elements)  # boundary conditions asembly??
    print("Number of elements: {}".format(elements.shape[0]))
    print("Number of equations: {}".format(neq))

    # System assembly
    KG = ass.assembler(elements, mats, nodes, neq, DME, sparse=True)
    RHSG = ass.loadasem(loads, IBC, neq)

    if np.sum(IBC == -1) < 3:  # 1 or zero fixed nodes/ pure neumann-boundary-condition system needs further constraints
        # System solution with custom conditions
        # solver with constraints to zero translation and zero rotation
        UG_sol, rx = custom_solver(KG, RHSG, mask_area, nodes, IBC, verbose=verbose)

    # System solution with default solver
    else:
        UG_sol = sol.static_sol(KG, RHSG)  # automatically detect sparce matrix
        if not (np.allclose(KG.dot(UG_sol) / KG.max(), RHSG / KG.max())):
            print("The system is not in equilibrium!")

    # average shear and normal stress on the colony area
    UC = pos.complete_disp(IBC, nodes, UG_sol)  # uc are x and y displacements
    E_nodes, S_nodes = pos.strain_nodes(nodes, elements, mats, UC)  # stresses and strains
    stress_tensor = calculate_stress_tensor(S_nodes, nodes, dims=mask_area.shape)  # assembling the stress tensor
    return UG_sol, stress_tensor


def grid_setup(mask_area, f_x, f_y, E=1, sigma=0.5, edge_factor=0):
    '''
    setup of nodes, elements, loads and mats(elastic material properties) lists for solids pys finite elements analysis. Every pixel of
    the provided mask is used as a node. Values from f_x,f_y at these pixels are used as loads. Mats is just
    [E, sigma].
    :param mask_area:
    :param f_x:
    :param f_y:
    :param E:
    :param sigma:
    :return:
    '''

    coords = np.array(np.where(mask_area))  # retrieving all coordintates from the  points  in the mask

    # setting up nodes list:[node_id,x_coordinate,y_coordinate,fixation_y,fixation_x]
    nodes = np.zeros((coords.shape[1], 5))
    nodes[:, 0] = np.arange(coords.shape[1])
    nodes[:, 1] = coords[1]  # x coordinate
    nodes[:, 2] = coords[0]  # y coordinate

    # creating an 2D array, with the node id of each pixel. Non assigned pixel is -1.
    ids = np.zeros(mask_area.shape).T - 1
    ids[coords[0], coords[1]] = np.arange(coords.shape[1], dtype=int)  # filling with node ids

    # fix all nodes that are exactely at the edge of the image (minus any regions close to the image edge that are
    # supposed to be ignored)in the movement direction perpendicular to the edge
    ids_cut, w = cut_mask_from_edge(ids, edge_factor, "", fill=False)
    edge_nodes_horizontal = np.hstack([ids_cut[:, 0], ids_cut[:, -1]]).astype(int)  # upper and lower image edge
    edge_nodes_vertical = np.hstack([ids_cut[0, :], ids_cut[-1, :]]).astype(int)  # left and right image edge
    edge_nodes_horizontal = edge_nodes_horizontal[edge_nodes_horizontal >= 0]
    edge_nodes_vertical = edge_nodes_vertical[edge_nodes_vertical >= 0]
    nodes[edge_nodes_vertical, 3] = -1  # fixed in x direction
    nodes[edge_nodes_horizontal, 4] = -1  # fixed in y direction
    nodes = nodes.astype(int)

    # setting up elements list:[ele_id,element type,reference_to_material_properties,node1,node2,node3,node4]
    # nodes must be list in counter clockwise oder for solidspy reasons
    elements = np.zeros((coords.shape[1], 7))

    # list the square(node,node left,node left down, node down) for each node. These are all posiible square shaped
    # elements, with the coorect orientation

    sqr = [(coords[0], coords[1] - 1), (coords[0] - 1, coords[1] - 1), (coords[0] - 1, coords[1]),
           (coords[0], coords[1])]
    # this produce negative indices, when at the edge of the mask
    # filtering these values
    filter = np.sum(np.array([(s[0] < 0) + (s[1] < 0) for s in sqr]),
                    axis=0) > 0  # logical to find any square with negative coordinates
    sqr = [(s[0][~filter], s[1][~filter]) for s in sqr]  # applying filter
    elements = elements[~filter]  # shortening length of elements list according to the same filter

    # enter node ids in elements, needs counter clockwise arangement
    # check by calling pyTFM.functions_for_cell_colonie.plot_grid(nodes,elements,inverted_axis=False,symbol_size=4,arrows=True,image=0)
    elements[:, 6] = ids[sqr[3]]
    elements[:, 5] = ids[sqr[2]]
    elements[:, 4] = ids[sqr[1]]
    elements[:, 3] = ids[sqr[0]]

    # cleaning up elements with nodes outside of cell area
    elements = np.delete(elements, np.where(elements == -1)[0], 0)
    # setting element id and attributing material and element type properties
    elements[:, 0] = np.arange(len(elements))  # id
    elements[:, 1] = 1  # element type/geometry (squares)
    elements[:, 2] = 0  # elastic properties reference
    elements = elements.astype(int)

    # setting up forces
    loads = np.zeros((len(nodes), 3))
    loads[:, 0] = np.arange(len(nodes))
    loads[:, 1] = f_x[coords[0], coords[1]]
    loads[:, 2] = f_y[coords[0], coords[1]]
    # loads=loads[not edge_nodes,:] ## check if this works/ is necessary

    mats = np.array([[E, sigma]])  # material properties: youngsmodulus, poisson ratio
    return nodes, elements, loads, mats


def prepare_forces(tx, ty, ps, mask):
    f_x = tx * ((ps * (10 ** -6)) ** 2)  # point force for each node from tractions
    f_y = ty * ((ps * (10 ** -6)) ** 2)
    f_x[~mask] = np.nan  # setting all values outside of mask area to zero
    f_y[~mask] = np.nan
    f_x_c1 = f_x - np.nanmean(f_x)  # normalizing traction force to sum up to zero (no displacement)
    f_y_c1 = f_y - np.nanmean(f_y)
    f_x_c2, f_y_c2, p = correct_torque(f_x_c1, f_y_c1, mask)
    return f_x_c2, f_y_c2

def correct_forces(f_x,f_y, mask_area):
    f_x[~mask_area] = np.nan  # setting all values outside of mask area to zero
    f_y[~mask_area] = np.nan
    f_x_c1 = f_x - np.nanmean(f_x)  # normalizing traction force to sum up to zero (no displacement)
    f_y_c1 = f_y - np.nanmean(f_y)
    f_x_c2, f_y_c2, p = correct_torque(f_x_c1, f_y_c1, mask_area)
    return f_x_c2, f_y_c2, p

def correct_torque(fx, fy, mask_area):
    com = regionprops(mask_area.astype(int))[0].centroid  # finding center of mass
    com = (com[1], com[0])  # as x y coordinate

    c_x, c_y = np.meshgrid(range(fx.shape[1]), range(fx.shape[0]))  # arrays with all x and y coordinates
    r = np.zeros((fx.shape[0], fx.shape[1], 2))  # array with all positional vectors
    r[:, :, 0] = c_x  # note maybe its also enough to chose any point as refernece point
    r[:, :, 1] = c_y
    r = r - np.array(com)

    f = np.zeros((fx.shape[0], fx.shape[1], 2), dtype="float64")  # array with all force vectors
    f[:, :, 0] = fx
    f[:, :, 1] = fy
    q = np.zeros((fx.shape[0], fx.shape[1], 2), dtype="float64")  # rotated positional vectors

    def get_torque_angle(p):
        q[:, :, 0] = + np.cos(p) * (f[:, :, 0]) - np.sin(p) * (f[:, :, 1])  # whats the mathematics behind this??
        q[:, :, 1] = + np.sin(p) * (f[:, :, 0]) + np.cos(p) * (f[:, :, 1])
        torque = np.abs(
            np.nansum(np.cross(r, q, axisa=2, axisb=2)))  ## using nna sum to only look at force values in mask
        return torque.astype("float64")

    # plotting torque angle relation ship
    # ps=np.arange(-np.pi/2,np.pi/2,0.01)
    # torques=[get_torque_angle(p)*1000 for p in ps]
    # plt.figure()
    # ticks=np.arange(-np.pi/2,np.pi/2+np.pi/6,np.pi/6)
    # tick_labels=[r"$-\frac{\pi}{2}$",r"$-\frac{\pi}{3}$",r"$-\frac{\pi}{6}$",r"$0$",r"$\frac{\pi}{6}$",r"$\frac{\pi}{3}$",r"$\frac{\pi}{2}$"]
    # plt.xticks(ticks,tick_labels,fontsize=25)
    # plt.yticks(fontsize=15)
    # plt.plot(ps,torques,linewidth=6)
    # plt.gca().spines['bottom'].set_color('black')
    # plt.gca().spines['left'].set_color('black')
    # plt.gca().tick_params(axis='x', colors='black')
    # plt.gca().tick_params(axis='y', colors='black')
    # plt.savefig("/home/user/Desktop/results/thesis/figures/torque_angle.png")

    pstart = 0
    # bounds = ([-np.pi], [np.pi])
    ## just use normal gradient descent??
    eps = np.finfo(float).eps  # minimum machine tolerance, for most exact calculation
    # trust region algorithm,
    # there seems to be a bug when using very small tolerances close to the machine precision limit (eps)
    # in rare cases there is an error. see also https://github.com/scipy/scipy/issues/11572
    try:
        p = least_squares(fun=get_torque_angle, x0=pstart, method="lm",
                          max_nfev=100000000, xtol=eps, ftol=eps, gtol=eps, args=())["x"]
    except KeyError:
        eps *= 5
        p = least_squares(fun=get_torque_angle, x0=pstart, method="lm",
                          max_nfev=100000000, xtol=eps, ftol=eps, gtol=eps, args=())["x"]



    q[:, :, 0] = + np.cos(p) * (f[:, :, 0]) - np.sin(p) * (f[:, :, 1])  # corrected forces
    q[:, :, 1] = + np.sin(p) * (f[:, :, 0]) + np.cos(p) * (f[:, :, 1])

    return q[:, :, 0], q[:, :, 1], p  # returns corrected forces and rotation angle


def get_torque1(fx, fy, mask_area, return_map=False):
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


def check_unbalanced_forces(fx, fy, mask=None, raise_error=False):
    if not isinstance(mask, np.ndarray):
        mask = np.logical_or(fy != 0, fx != 0)
    torque = get_torque1(fx, fy, mask, return_map=False)  # torque of the system
    net_force = np.array([np.sum(fx), np.sum(fy)])
    print("torque = " + str(torque))
    print("net force = " + str(net_force))
    if raise_error and (torque == 0 or np.sum(net_force == 0) > 0):
        raise Exception("torque or net force is not zero")


def make_field(nodes, values, dims):
    '''
    function to write e.g. loads or deformation data to array
    :param nodes:
    :param values:
    :return:
    '''
    nodes = nodes.astype(int)
    fx = np.zeros(dims)
    fy = np.zeros(dims)
    fx[nodes[:, 2], nodes[:, 1]] = values[:, 0]
    fy[nodes[:, 2], nodes[:, 1]] = values[:, 1]
    return fx, fy


def make_solids_py_values_list(nodes, fx, fy, mask, shape=1):
    '''
    function to create a list of values, eg deformation as needed by solidspy

    :param nodes:
    :param fx:
    :param fy:
    :return:
    '''
    nodes = nodes.astype(int)
    mask = mask.astype(bool)
    if shape == 1:
        l = np.zeros((len(nodes) * 2))
        l[np.arange(0, len(nodes) * 2) % 2 == 0] = fx[mask].flatten()  # ordering in solidspy is x,y..
        l[np.arange(0, len(nodes) * 2) % 2 != 0] = fy[mask].flatten()
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


def get_torque2(nodes, loads):
    nodes = nodes.astype(int)

    k, l = (np.max(nodes[:, 1]) + 1, np.max(nodes[:, 2]) + 1)
    fx = np.zeros((k, l))
    fy = np.zeros((k, l))
    area = np.zeros((k, l), dtype=int)
    fx[nodes[:, 1], nodes[:, 2]] = loads[:, 1]
    fy[nodes[:, 1], nodes[:, 2]] = loads[:, 2]
    area[nodes[:, 1], nodes[:, 2]] = 1

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
    return (torque)


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
    com = (np.mean(r[:, :, 0][mask]), np.mean(r[:, :, 1][mask]))
    r = r - np.array(com)
    r[~mask] = 0
    d[:, :, 0][mask] = dx[mask].flatten()  # note maybe its also enough to chose any point as refernece point
    d[:, :, 1][mask] = dy[mask].flatten()
    return np.sum(np.cross(r, d, axisa=2, axisb=2))


# applying rotation
def rot_displacement(p, r, r_n):
    r_n[:, :, 0] = + np.cos(p) * (r[:, :, 0]) - np.sin(p) * (r[:, :, 1])  # rotation of postional vectors
    r_n[:, :, 1] = + np.sin(p) * (r[:, :, 0]) + np.cos(p) * (r[:, :, 1])
    disp = r - r_n
    return disp  # norma error measure


def correct_rotation(def_x, def_y, mask):
    '''
    function to apply rigid body translation and rotation to a deformation field, to minimize rotation
    and translation
    :return:
    '''
    mask = mask.astype(bool)
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
    com = (np.mean(r[:, :, 0][mask]), np.mean(r[:, :, 1][mask]))  ## why inverted indices??????????
    r = r - np.array(com)
    r[~mask] = 0
    d[:, :, 0][mask] = def_xc1[mask].flatten()  # note maybe its also enough to chose any point as refernece point
    d[:, :, 1][mask] = def_yc1[mask].flatten()

    # fit to corect rotation
    def displacement_error(p):
        r_n[:, :, 0] = + np.cos(p) * (r[:, :, 0]) - np.sin(p) * (r[:, :, 1])  # rotation of postional vectors
        r_n[:, :, 1] = + np.sin(p) * (r[:, :, 0]) + np.cos(p) * (r[:, :, 1])
        disp = r - r_n
        return np.sum(np.linalg.norm((d[mask] - disp[mask]), axis=1))  # norma error measure

    pstart = -1
    bounds = ([-np.pi], [np.pi])
    ## just use normal gradient descent??
    p = least_squares(fun=displacement_error, x0=pstart, bounds=bounds, method="trf",
                      max_nfev=100000000, xtol=3e-32, ftol=3e-32, gtol=3e-32, args=())["x"]  # trust region algorithm,
    ## note only works if displacement can be reached in "one rotation!!!!
    # get the part of displacement originating form a rotation of p
    d_rot = rot_displacement(p, r, r_n)
    d_n[mask] = d[mask] - d_rot[mask]  # correcting this part of rotation
    return d_n[:, :, 0], d_n[:, :, 1], trans, p


def find_eq_position(nodes, IBC, neq):
    # based on solidspy.assemutil.loadasem

    nloads = IBC.shape[0]
    RHSG = np.zeros((neq, 2))
    x_points = np.zeros((neq)).astype(bool)  # mask showing which point has x deformation
    y_points = np.zeros((neq)).astype(bool)  # mask showing which point has y deformation
    for i in range(nloads):
        il = int(nodes[i, 0])  # index of the node
        ilx = IBC[il, 0]  # indices in RHSG or fixed nodes, if -1
        ily = IBC[il, 1]
        if ilx != -1:
            RHSG[ilx] = nodes[i, [1, 2]]  # x,y position/ not the orientation
            x_points[ilx] = [True]
        if ily != -1:
            RHSG[ily] = nodes[i, [1, 2]]
            y_points[ily] = [True]

    return RHSG.astype(int), x_points, y_points


def custom_solver(mat, rhs, mask_area, nodes, IBC, verbose=False):
    # IBC is "internal boundary condition" contains information about which nodes are fixed and
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

    len_disp = mat.shape[1]  # length of the  displacement vector
    zero_disp_x = np.zeros(len_disp)
    zero_disp_y = np.zeros(len_disp)
    zero_torque = np.zeros(len_disp)

    com = regionprops(mask_area.astype(int))[0].centroid  # finding center of mass
    com = (com[1], com[0])  # as x y coordinate

    c_x, c_y = np.meshgrid(range(mask_area.shape[1]), range(mask_area.shape[0]))  # arrays with all x and y coordinates
    r = np.zeros((mask_area.shape[0], mask_area.shape[1], 2))  # array with all positional vectors
    r[:, :, 0] = c_x  # Note: maybe its also enough to chose any point as reference point
    r[:, :, 1] = c_y
    # solidspy function that is used to construct the loads vector (rhs)
    nodes_xy_ordered, x_points, y_points = find_eq_position(nodes, IBC, len_disp)
    r = r[nodes_xy_ordered[:, 1], nodes_xy_ordered[:, 0], :]  # ordering r in the same order as rhs
    r = r - np.array(com)

    zero_disp_x[x_points] = 1
    zero_disp_y[y_points] = 1

    # torque=sum(r1*f2-r2*f1)   # TDOD: this is actually zero rotation
    zero_torque[x_points] = r[x_points, 1]  # -r2 factor
    zero_torque[y_points] = -r[y_points, 0]  # r1 factor
    add_matrix = np.vstack([zero_disp_x, zero_disp_y, zero_torque])
    # adding zero conditions for force vector and torque
    rhs = np.append(rhs, np.zeros(3))

    if type(mat) is csr_matrix:
        import scipy.sparse
        # convert additional conditions to sparse matrix
        mat = scipy.sparse.vstack([mat, csr_matrix(add_matrix)], format="csr")
        u_sol, error = \
        np.array(lsqr(mat, rhs, atol=10 ** -12, btol=10 ** -12, iter_lim=200000, show=verbose, conlim=10 ** 12))[
            [0, 3]]  # sparse least squares solver
    elif type(mat) is np.ndarray:
        # adding to matrix
        mat = np.append(mat, add_matrix, axis=0)

        u_sol, error = np.array(np.linalg.lstsq(mat, rhs))[[0, 1]]
    else:
        raise TypeError("Matrix should be numpy array or csr_matrix.")

    return u_sol, error


if __name__ == "__main__":
    import clickpoints

    db = clickpoints.DataFile(
        "/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/mask_cell_boundary3.cdb",
        "r")
    mask = db.getMask(frame=0).data

    # loading traction forces
    t_x = np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/tx.npy")
    t_y = np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/ty.npy")

    ## try interpolation on t_x, ty for higher resolution or whatever

    # interpolation

    # from scipy.misc import imresize
    # mask_int2=imresize(mask,t_x.shape,interp='lanczos') # this succs by the way

    # some use full pre clean up
    mask = remove_small_holes(mask, 100)
    mask = remove_small_objects(label(mask), 1000) > 0  # removing other small bits

    mask_int = interpolation(mask, t_x)
    # plt.figure()
    # plt.imshow(mask)
    # plt.figure()
    # plt.imshow(mask_int)

    # plt.figure()
    # plt.imshow(t_x)

    ## mask data is prepared in exactly the same way as for stress analysis along lines
    mask_area, mask_boundaries = prepare_mask(mask_int)
    # plt.figure()
    # plt.imshow(mask_area)
    # plt.figure()
    # plt.imshow(mask_boundaries)

    # setting up files for finite elements
    # make matrix with ids
    nodes, elements, loads, mats = grid_setup(mask_area, t_x, t_y, 1, 0.3)
    # plot_grid(nodes,elements,inverted_axis=True)  # only use with <1000 nodes

#TODO: implement seperate mask to average stresses and stuff optionally even multiple stresses
#TODO: Move away from recursive functions to identfy lines