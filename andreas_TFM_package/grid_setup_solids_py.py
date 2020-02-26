
from andreas_TFM_package.solids_py_stress_functions import *
from andreas_TFM_package.utilities_TFM import make_random_discrete_color_range, invert_dictionary,join_dictionary
from skimage.morphology import skeletonize,remove_small_holes,remove_small_objects,label,binary_dilation,binary_erosion
from scipy.ndimage.morphology import binary_fill_holes
from scipy.ndimage.measurements import find_objects
from scipy.ndimage import binary_closing as binary_clo
from scipy.signal import convolve2d
from scipy.interpolate import splprep, splev
from itertools import chain
from collections import Counter
from contextlib import suppress



def show_points(ps,mask):
    plt.figure()
    plt.imshow(mask)
    plt.plot(ps[:,1],ps[:,0],"or")

def identify_cells(mask_area,mask_boundaries,points):
    '''
    function to identfy cells. Each cell is is a dictionary entry with a list of ids, reffering to
    points.
    :param mask:
    :param area:
    :param mask_boundaries:
    :return:
    '''


    cells={}  # dictionary containg a list of point idsthat sourround each cell
    cells_area={} # dictionary containg a all pixels belonging to that cell as boolean aray
    # points_to_flatt array map:
    # labeling each cell
    m = mask_area.astype(int) - mask_boundaries.astype(int)
    ml = label(m, connectivity=1)
    # creating a list of point coordinates corresponding to a flattend array
    # this will allow easier identification of the id of a point
    points_fl = (points[:, 0]) + mask_area.shape[0] * (points[:, 1])
    sort_ids=np.argsort(points_fl) # sorting will allow np.searchsorted function;
    points_fl=points_fl[sort_ids] #

    for i,l in enumerate(np.unique(ml)[1:]):  # getting each cell border by binary dilation of one pixel; iterating over each cell
        m_part = (ml == l).astype(bool)  #extracting a cell area
        edge = np.logical_and(binary_dilation(m_part), ~m_part) # getting the boundary of a cell
        ps= np.array(np.where(edge)).T # finding coordinates
        ps_fl = (ps[:,0]) + mask_area.shape[0] * (ps[:,1])  # convert coordinates to the one of a flat array
        p_ids=sort_ids[np.searchsorted(points_fl,ps_fl)] # find indices, where i would need to insert,supposed to be the fastest way
        #  and read index from unsorted list
        cells[i]=p_ids # save to dictionary
        cells_area[i]=m_part

    ## vizualization
    # creating a(random lsit of hex colors)
    #colors = []
    #for i in range(len(cells.items())):
    #    colors.append('#%06X' % np.random.randint(0, 0xFFFFFF))
    #plt.figure()
    #plt.imshow(mask_boundaries)
    #offset = 0.01
    #for cell, ps in cells.items():
    #    for p in ps:
    #        plt.plot(points[p][1] + offset * cell, points[p][0] + offset * cell, "o", color=colors[cell])

    return cells,cells_area

def spline_interpolation(line, points,k=3,endpoints=None):
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

    if endpoints: # adding endpoints if provided
        # splprep throws error if two points are exactely identidcal
        # for ecxample for the "dead-end-lines" the endpoint is identical to the first and last point and thus shouldnt
        # be added here
        if not (endpoints[0][1]==x[0] and  endpoints[0][0]==y[0]):
            x = np.concatenate([[endpoints[0][1]], x],axis=0)
            y = np.concatenate([[endpoints[0][0]], y],axis=0)
        if not (endpoints[1][1] == x[-1] and endpoints[1][0] == y[-1]):
            x = np.concatenate([x, [endpoints[1][1]]],axis=0)
            y = np.concatenate([y, [endpoints[1][0]]],axis=0)


    k = len(x) - 1 if len(x) <= 3 else 3  # addapt spline order, according to number of points
    tck, u = splprep([x, y], s=10, k=k)  ### parametric spline interpolation
    # fits essentially function: [x,y] =f(t) , t(paramter) is default np.linspace(0,1,len(x))
    # tck is array with x_position of knot, y_position of knot, parameters of the plne, order of the spline
    # ui is paramter given to x,y points, in this case default (values from 0 to 1)
    # s is smoothing factor(?)
    # k is order of the spline, cubic is fine

    if endpoints: # forcing splines to through endpoints
        tck[1][0][0] = endpoints[0][1]
        tck[1][0][-1] = endpoints[1][1]
        tck[1][1][0] = endpoints[0][0]
        tck[1][1][-1]  =endpoints[1][0]
    # new points from spline interpolation (thes points will be used for normal stress vector calculation
    points_new = np.array(splev(u, tck, der=0)).T
    # in xy orientation
    #points_new = np.round(points_new).astype(int)  ## could also use exact values and the interpolate the stress tensor, but thats probably to complicated
    return tck, u, points_new # points new in xy orientation


def arrange_lines_from_endpoints(cells_lines,lines_endpoints_com):

    '''
    rearranging the order of lines in the line_id list for one one cell. The endpoints for all lines for ne cell are
    extracted and then puzzled together. Each endpoint must occure in two lines.

    :param cells_lines: dictionary cell_id:[line_ids]
    :param lines_endpoints_com: dictionary line_id:[endpoint1,endpoint2] , each endoinpt is an array with x and y coordinates
    :return: cells_lines_new: updated cells_lines dictionary

    '''
    cells_lines_new={}

    for cell_id,line_ids in cells_lines.items():
        # extracting relevant endpoints
        local_endpoints = {line_id:lines_endpoints_com[line_id] for line_id in line_ids}
        # rearranging endpoints into an suitabel array: axis0: lines axis 1:[endpoint1, endpoint2], axis3: x,y coordinates
        eps=np.array([np.array([value[0],value[1]]) for value in local_endpoints.values()])

        new_line_ids=[] #newly_arranged line_ids
        p_ind1=0 # start ids
        p_ind2=0
        # iterating through eps, until all endpoints have been visited
        while not np.isnan(eps).all():
            point = copy.deepcopy(eps[p_ind1, p_ind2]) #extracting an endpoint
            eps[p_ind1, p_ind2] = np.nan # remove the point from array
            # find second occurrence, by taking the norm of the diffrece between the current point and all other points
            # this should be zero
            np_ind1,np_ind2=np.array(np.where(np.linalg.norm(eps - point,axis=2)==0)).T[0]
            new_line_ids.append(line_ids[np_ind1]) # note corresponding line_ids
            eps[np_ind1, np_ind2] = np.nan  # remove this point from array
            p_ind1, p_ind2= np_ind1,np.abs(np_ind2-1)# change to the other end of the line

        cells_lines_new[cell_id]=new_line_ids # update dictionary

    return cells_lines_new


def find_edge_lines(cells_lines):
    '''
    Finding all lines (cell borders) at the edge of a cell colony. Simply checks if
    a line is associated to only one cell.
    :param cells_lines: dictionary with cell_id:[associated line_ids]
    :return: edge_lines: lsit of line ids at the edge of the cell colony
    '''
    all_lines=np.array(list(chain.from_iterable(cells_lines.values()))) # unpacking all line ids
    counts=Counter(all_lines) # counting occurences
    edge_lines=[line for line,count in counts.items() if count==1] # select if line_id was only associated to one cell
    return edge_lines

def center_of_mass_cells(cells_points,points):
    '''
    calulating the "center of mass" of a cell using only the points at the edge of the cell
    :param cells_points:
    :param points:
    :return:
    '''
    cells_com={}
    for cell_id,hull_points in cells_points.items():
        cells_com[cell_id]=np.mean(points[hull_points],axis=0)

    return cells_com

def remove_circular_line(lines_endpoints_com, lines_points, lines_endpoints):
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
    #print(circular)
    # clearing these lines from the input dictionaries
    for l_id in circular:
        del lines_endpoints_com[l_id]
        del lines_points[l_id]
        del lines_endpoints[l_id]



def interpolate_cell_area(cells_area,shape):
    cells_area_interpol={}
    for cell_id, areas in cells_area.items():
        cells_area_interpol[cell_id]=interpolation(areas,shape,min_cell_size=0)
    return cells_area_interpol

def interpolate_points_dict(points_dict,shape_target,shape_orgin):
    points_dict_interp={}
    for p_id,coords in points_dict.items():
        points_dict_interp[p_id]=(interpolation_single_point(coords[0],shape_target,shape_orgin),
                                  interpolation_single_point(coords[1],shape_target,shape_orgin))
    return points_dict_interp






class Cells_and_Lines:
    # container for cells and lines, assignement of points with them, and assignement with each other
    def __init__(self, mask_boundaries, shape, graph,points):

        # masks, graph and points including dead-end lines // this distinction is mostly due to historic reasons
        self.mask_boundaries_wp = mask_boundaries
        self.inter_shape = shape
        # graph as a dictionray with key=ponit id, values: ids of neighbouring points
        # any point id is the index in the points array (contains coordinate of these points
        self.graph_wp = graph
        self.points_wp = points
        self.graph, self.points,removed=remove_endpoints_wrapper(self.graph_wp, self.points_wp)
        # masks, graph and points excluding dead-end lines
        self.mask_boundaries = graph_to_mask(self.graph, self.points, mask_boundaries.shape)  # rebuilding the mask

        # interpolate points to the size of the future FEM_grid
        self.points_interpol=interpolation_single_point(self.points, self.inter_shape,self.mask_boundaries.shape)
        # interpolation factors used in the fun
        #self.inerpol_factors=np.array([self.inter_shape[0] /self.mask_boundaries.shape[0], self.inter_shape[1] / self.mask_boundaries.shape[1]])
        # points as dictionary with key=points id, values: points coordinates
        self.points_dict={i:self.points[i] for i in range(len(self.points))}

        # lines as a dictionary with key=line id, values: ids of containg points (in correct order)
        self.lines_points=identify_line_segments(self.graph)

        # cells as a dictionary with key=cell id, values: ids of containing points (not ordered)
        self.cells_points, self.cells_area=identify_cells(self.mask_boundaries, binary_fill_holes(self.mask_boundaries), self.points)

        # interpolate the area of individual cells to the size of deformation
        self.cells_area_interpol=interpolate_cell_area(self.cells_area,self.inter_shape)

        # self.points_lines = invert_dictionary(self.lines_points) # point_id:line_id
        self.max_line_id = np.max(list(self.lines_points.keys()))
        self.de_lines_points, self.max_line_id = find_dead_end_lines(self.graph_wp, list(self.graph.keys()),
                                                                    self.max_line_id)
        self.de_endpoints={key:(self.points[value[0]],self.points[value[-1]]) for key,value in self.de_lines_points.items()}# using exact endpoints for the dead end lines
        self.allLines_points=join_dictionary(self.lines_points,self.de_lines_points)

        # dictionary with endpoints, needed to completely fill the gaps between all cell_lines
        self.lines_endpoints_com,self.lines_endpoints = find_exact_line_endpoints(self.lines_points, self.points, self.graph)

        # removing all lines that are predicted to be circular. Mostly a problem for very short lines
        remove_circular_line(self.lines_endpoints_com,self.lines_points,self.lines_endpoints)


        # dictionary with both dead end lines and circular lines


        # center of mass of cells, calculated only from the hull points
        self.cells_com=center_of_mass_cells(self.cells_points,self.points)
        # dictionary to associate cells with correct lines, key is cell id, value is line_id
        self.cells_lines= defaultdict(list) ## improve by using endpoints for cell_line_association
        # dictionary to associate cells with correct lines, key is line id, value is cell


        self.lines_cells=defaultdict(list)
        for l_id, l in self.lines_points.items():
            for c_id, c in self.cells_points.items():
                if l[0] in c:
                    self.cells_lines[c_id].append(l_id)
                    self.lines_cells[l_id].append(c_id)

        # using the new endpoints to arrange lines in the correct way
        self.cells_lines = arrange_lines_from_endpoints(self.cells_lines, self.lines_endpoints_com)


        # adding dead end endpoints only now to avoid complications when identifying cells
        self.de_endpoints = {key: (self.points[value[0]], self.points[value[-1]]) for key, value in
                             self.de_lines_points.items()}
        self.lines_endpoints_com=join_dictionary(self.lines_endpoints_com,self.de_endpoints)
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
        self.central_lines = [line_id for line_id in self.allLines_points.keys() if line_id not in self.edge_lines and line_id not in self.dead_end_lines]
        self.n_cells=len(self.cell_ids)
        # list with original line lengths--> later used for interpolation
        self.line_lengths={key:len(value) for key, value in self.allLines_points.items()}

        # dictionary containing the spline represetation of the points as a parametric function
        #[x,y]=f(u). u is always np.linspace(0,1,"number of points in the line). Use scipy.interpolate.splev
        # to evaluate at other points
        self.lines_splines = defaultdict(list)

        # dictionary with the normal vectors for each line acording to the spline interpolation. Contains list
        # of the normal vectors , position of these points is listed in lines_spline_points
        # interpolation !!! positions are given in xy order !!!
        self.lines_n_vectors= defaultdict(list)

        # dictionary with example points (where a normal vector originates) from spline interpolation
        self.lines_spline_points = defaultdict(list)

        for line_id, line_ps in self.allLines_points.items():
            endpoints = self.lines_endpoints_interpol[line_id]
            k = len(line_ps)+2 - 1 if len(line_ps) <= 3 else 3 # addapt spline order, according to number of points
            # spline order must be below nomber of points, so choose 2 if lne has one point + 2 endpoints
            tck, u, points_new=spline_interpolation(line_ps, self.points_interpol,k=k,endpoints=endpoints) #spline interpolation
            self.lines_splines[line_id]=tck # saving the spline object
            # saving a few oints and n vectors for easy representations/ debugging,
            # these points will not be used further
            n_vectors=normal_vectors_from_splines(u, tck) # calculating normal vectors as a list
            self.lines_n_vectors[line_id]=n_vectors
            self.lines_spline_points[line_id] = points_new

    def cut_to_FEM_grid(self, FEM_mask):
        '''
        removing any points that line outside the FEM Grid
        :param FEM_mask:
        :return:
        '''
        self.lines_outside = []
        # spline points is already interpolated to shape of self.mask_area
        for line_id, line_points in list(self.lines_spline_points.items()):
            line_points = np.round(line_points).astype(int)
            if np.sum(~FEM_mask[line_points[:, 1], line_points[:, 0]]) > 5:
                self.lines_spline_points.pop(line_id, None)
                self.allLines_points.pop(line_id, None)
                self.line_lengths.pop(line_id, None)
                self.lines_endpoints_com.pop(line_id, None)
                self.de_lines_points.pop(line_id, None)
                self.lines_points.pop(line_id, None)
                with suppress(ValueError): self.edge_lines.remove(line_id)
                with suppress(ValueError): self.edge_lines.remove(line_id)
                with suppress(ValueError): self.edge_lines.remove(line_id)
                self.lines_outside.append(line_id)
    def return_n_array(self,fill_nan=True):
        '''
        writes (normal) vectors in a two dimensional array according to their position from spline interpolation
        :return:
        '''
        if fill_nan:
            n_array = np.zeros((self.mask_boundaries.shape[0], self.mask_boundaries.shape[1], 2)) + np.nan
        else:
            n_array = np.zeros((self.mask_boundaries.shape[0], self.mask_boundaries.shape[1], 2))
        for vecs, ps in zip(self.lines_n_vectors.values(),self.lines_spline_points.values()):

            n_array[ps[:, 1], ps[:, 0]] = vecs
        return n_array

    def vizualize_lines_and_cells(self,sample_factor=1,plot_n_vectors=False):
        '''
        plotting the id of lines for each point, after spline interpolation,
        lines at the edge have a different color
        plotting the id of cells
        plotting_normal_vectors, as predcited by intial spline interpolation
        :param cells_points: factor to reduce line_id labels that are plotted, must be <= 1
        :return:
        '''

        offset = 0.005 # ofset for point id text
        fig=plt.figure()
        plt.imshow(self.mask_boundaries,cmap="jet")
        plt.plot([], [], color="green",label="line_ids") #legend for line ids
        plt.plot([], [], color="red", label="cells_ids") #legend for cells ids
        # plotting all points with line id as label
        colors=np.array(["C1","C2","C3"]) # choosing colors according to line type
        line_classifier=[self.central_lines,self.edge_lines,self.dead_end_lines]
        for l, ps in self.allLines_points.items():
            ps=np.array(ps)
            color=colors[np.array([l in line_ids for line_ids in line_classifier])][0]
            p_indices=np.array(list(range(len(ps)))) #all indicces
            # randomly choising a few indices, without first and laast index
            ps_select=np.random.choice(p_indices[1:-1],size=int((len(ps)-2)*sample_factor),replace=False)
            ps_select=np.append(ps_select,p_indices[np.array([1,-1])]) # adding back first and last index

            plt.plot(self.points[ps][:,1], self.points[ps][:,0],"o",color=color)
            for p in ps[ps_select]: # labeling selected points
                plt.text(self.points[p][1] + 1 * offset * l, self.points[p][0] + 1 * offset * l, s=str(l),color="green")
        # plotting cel id at center of mass of cell
        for cell_id,com in self.cells_com.items():
            plt.text(com[1],com[0],str(cell_id),color="red",fontsize=13)

        if plot_n_vectors:
            for points,n_vectors in zip(self.lines_spline_points.values(),self.lines_n_vectors.values()):
                for n_vec,p in zip(n_vectors,interpolation_single_point(points, self.mask_boundaries.shape, self.inter_shape)):
                    plt.arrow(p[0], p[1], n_vec[0], n_vec[1], head_width=0.15)
        plt.legend()
        return fig

    def vizualize_splines(self,sample_factor=1,subset=np.inf):
        # sample factor: only every nth point is used for plotting
        plt.figure()
        plt.imshow(self.mask_boundaries)
        colors=make_random_discrete_color_range(len(self.lines_spline_points.keys()))
        for i,(l_id,points) in enumerate(self.lines_spline_points.items()):
            if i>subset:
                break
            points=interpolation_single_point(points,self.mask_boundaries.shape,self.inter_shape)
            plt.plot(points[::sample_factor,0],points[::sample_factor,1],color=colors[i])

    def simple_line_plotting(self,lines,subset=np.inf):
        plt.figure()
        plt.imshow(self.mask_boundaries)
        for i,ps in enumerate(lines.values()):
            if i>subset:
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
        self.lines_points=find_lines_simple(self.graph)
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
        self.central_lines= self.line_ids
        self.line_lengths={key:len(value) for key, value in self.lines_points.items()}

        # very rough estimate of cell number
        label_mask,self.n_cells=label(~mask_boundaries,connectivity=1,return_num=True)


def prepare_mask_FEM(mask,shape):
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
    mask_int = interpolation(binary_fill_holes(mask) , shape)
    mask_int = binary_fill_holes(mask_int).astype(bool) # sometime first binary fill holes is not enough

    # removing unsuitable pixels for grid later --> a pixel must be connected to at least more then 2 others
    m=convolve2d(mask_int.astype(int),np.array([[0,1,0],[1,0,1],[0,1,0]]),mode="same",boundary="fill",  fillvalue=0) # convoultion will produce 1 if only one direct (distance 1) neighbour exists
    p=np.where(np.logical_and(m==1,mask_int)) # problematic coordinates

    for x,y in zip(p[0],p[1]):
        new_ps=np.array([[x+1,y],[x,y+1],[x-1,y],[x,y-1]])# all relevant neigbouring points
        new_p=new_ps[mask_int[new_ps[:,0],new_ps[:,1]]][0]  # checking where a possible replacement point could be located this can only have one result
        mask_int[x,y]=False # removing old point
        mask_int[new_p[0],new_p[1]]=True # adding new point
    return mask_int





def find_borders(mask,shape,raise_error=True,type="colony"):
 #### maybe reintroduce small cell filter
    # removing small bits
    mask = remove_small_objects(mask.astype(bool), 1000).astype(bool)
    # generating borders
    mask_boundaries = skeletonize(mask.astype(int))
    # converting mask to graph object
    graph,points = mask_to_graph(mask_boundaries)
    # finding dead ends: cell borders wich dont connect to other cell borders at one end:
    # this is use full to get a clean structure of cell borders, which is later used to identifying the number and area of cells
    # applying remove endpoints multiple times to deal with forking dead ends

    try :
        if type=="colony":
            c_l=Cells_and_Lines(mask_boundaries,shape, graph,points)
        if type=="cell layer":
            c_l = Cells_and_Lines2(mask_boundaries, shape, graph, points)
    except (RecursionError,FindingBorderError) as e:
        print("original error: ", e)
        if raise_error:
            raise FindingBorderError
        else:
            return None

    #c_l.cut_to_FEM_grid(mask_int) # this should never be necessary
    #c_l.vizualize_lines_and_cells(sample_factor=0.2,plot_n_vectors=True)
    #c_l.vizualize_splines(sample_factor=4,subset=200000)
    #c_l.simple_line_plotting(c_l.lines_points)
    ## vizualization of spline interpolation with new line enpoints
    return c_l


def interpolation(mask, dims, min_cell_size=100,dtype=bool):
    #
    # some pre clean up of the mask
    mask = remove_small_holes(mask.astype(bool), min_cell_size)
    mask = remove_small_objects(mask.astype(bool), 1000) # removing other small bits
    # note: remove_small_objects labels automatically if mask is bool
    coords = np.array(np.where(mask)).astype(float) # coordinates of all points
    interpol_factors = np.array([dims[0] / mask.shape[0], dims[1] / mask.shape[1]])
    coords[0] = coords[0] * interpol_factors[0]  # interpolating x coordinates
    coords[1] = coords[1] * interpol_factors[1]  # interpolating xy coordinates
    coords = np.round(coords).astype(int)

    coords[0, coords[0] >= dims[0]] = dims[0]-1 # fixing issue when interpolated object is just at the image border
    coords[1, coords[1] >= dims[1]] = dims[1]-1

    mask_int = np.zeros(dims)
    mask_int[coords[0], coords[1]] = 1
    mask_int = mask_int.astype(int)
    #filling gaps if we interpolate upwards
    if dims[0]*dims[1]>=mask.shape[0]*mask.shape[1]:
        iter=int(np.ceil(np.max([mask.shape[0]/dims[0],mask.shape[0]/dims[0]]))*5) # times 5 is safety factor
        mask_int=binary_clo(mask_int,iterations=10)
        print(iter)
    return mask_int.astype(bool)

def interpolation_single_point(point,shape_target,shape_origin):
    # is also works with 2d arrays of shape(n,2)
    interpol_factors = np.array([shape_target[0] /shape_origin[0], shape_target[1] / shape_origin[1]])
    point_interp=point*interpol_factors
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

def cut_mask_from_edge(mask,cut_factor,warn_flag=False):
    mask_cut=copy.deepcopy(mask)
    sum_mask1=np.sum(mask_cut)
    dims=mask_cut.shape
    inds=[int(dims[0]*cut_factor),int(dims[0]-(dims[0]*cut_factor)),int(dims[1]*cut_factor),int(dims[1]-(dims[1]*cut_factor))]
    mask_cut[:inds[0], :] = 0
    mask_cut[inds[1]:, :] = 0
    mask_cut[:, :inds[2]] = 0
    mask_cut[:, inds[3]:] = 0
    sum_mask2 = np.sum(mask_cut)
    warn="mask was cut close to image edge" if (sum_mask2<sum_mask1 and warn_flag) else ""
    return mask_cut, warn

def cut_mask_from_edge_wrapper(cut_factor,mask,parameter_dict,cut=True,warn=""):
    warn_flag = parameter_dict["cut_instruction"][parameter_dict["FEM_mode"]]
    if cut:
        mask, warn = cut_mask_from_edge(mask,cut_factor,parameter_dict["TFM_mode"]=="colony")
    return mask, warn

def grid_setup(mask_area, f_x, f_y, E, sigma,edge_factor):
    '''
    setup of nodes, elements, loads and mats(elastic material properties) lists for solids pys finite elments analysis. Every pixel of
    the provided mask is used as a node. Values from f_x,f_y at these pixels are used as loads. Mats is just
    [E, sigma].
    :param mask_area:
    :param f_x:
    :param f_y:
    :param E:
    :param sigma:
    :return:
    '''

    coords = np.array(np.where(mask_area)) # retrieving all coordintates from the  points  in the mask

    # creating an 2D array, with the node id of each pixel. Non assigned pixel is -1.
    ids = np.zeros(mask_area.shape).T - 1
    ids[coords[1], coords[0]] = np.arange(coords.shape[1], dtype=int) # filling with node ids

    # setting up nodes list:[node_id,x_coordinate,y_coordinate,fixation_y,fixation_x]
    nodes = np.zeros((coords.shape[1], 5))
    nodes[:, 0] = np.arange(coords.shape[1])
    nodes[:, 1] = coords[1]  # x coordinate
    nodes[:, 2] = coords[0]  # y coordinate

    # fix all nodes that are exactely at the edge of the image (minus any regions close to the image edge that are
    # supposed to be ignored)in the movement direction perpendicular to the edge
    ids_cut,w=cut_mask_from_edge(ids,edge_factor,"")
    edge_nodes_horizontal = np.hstack([ids_cut[:, 0], ids_cut[:, -1]]).astype(int) # upper and lower image edge
    edge_nodes_vertical = np.hstack([ids_cut[0, :], ids_cut[-1, :]]).astype(int) # left and right image edge
    edge_nodes_horizontal = edge_nodes_horizontal[edge_nodes_horizontal >= 0]
    edge_nodes_vertical = edge_nodes_vertical[edge_nodes_vertical >= 0]

    nodes[edge_nodes_vertical, 3] = -1  # fixed in x direction
    nodes[edge_nodes_horizontal, 4] = -1  # fixed in y direction
    nodes = nodes.astype(int)




    # setting up elements list:[ele_id,element type,reference_to_material_properties,node1,node2,node3,node4]
    # nodes must be list in counter clockwise oder for solidspy reasons
    elements = np.zeros((coords.shape[1], 7))

    #list the square(node,node left,node left down, node down) for each node. These are all posiible square shaped
    # elements, with the coorect orientation

    sqr = [(coords[1], coords[0] - 1),(coords[1] - 1, coords[0] - 1),(coords[1] - 1, coords[0]),(coords[1], coords[0])]
    # this produce negative indices, when at the edge of the mask
    # filtering these values
    filter=np.sum(np.array([(s[0]<0)+(s[1]<0) for s in sqr]),axis=0)>0 # logical to find any square with negative coordinates
    sqr=[(s[0][~filter],s[1][~filter]) for s in sqr] # applying filter
    elements=elements[~filter] # alos shortening length of elements list according to the same filter



    elements[:, 6] = ids[sqr[0][0], sqr[0][1]]  ## zusammenfassen
    elements[:, 5] = ids[sqr[1][0], sqr[1][1]]
    elements[:, 4] = ids[sqr[2][0], sqr[2][1]]
    elements[:, 3] = ids[sqr[3][0], sqr[3][1]]



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
    #loads=loads[not edge_nodes,:] ## check if this works/ is necessary


    # filtering all nodes that are not in elemnte(e.g some singel points at the edge) using sets
    ## hopefully fixed on in prepare mask, so we dnt need this

    #set1=set(np.unique(elements[:,[3,4,5,6]].flatten()))
    #set2=set(nodes[:, 0])
    #dif=set2.difference(set1) # finding complement
    #nodes_id=copy.deepcopy(nodes) #"nan padded nodes, for easy retrieval of ids.. mybe find better solution
    #nodes_id[list(dif)]=np.nan
    #nodes=np.delete(nodes,list(dif),axis=0)
    #loads = np.delete(loads, list(dif), axis=0)


    mats = np.array([[E, sigma]])  # material properties: youngsmodulus, poisson ratio
    return nodes, elements, loads, mats


if __name__ == "__main__":
    import clickpoints
    db=clickpoints.DataFile("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/mask_cell_boundary3.cdb","r")
    mask=db.getMask(frame=0).data

    # loading traction forces
    t_x=np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/tx.npy")
    t_y=np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/test_data/ty.npy")




    ## try interpolation on t_x, ty for higher resolution or whatever



    # interpolation

    #from scipy.misc import imresize
    #mask_int2=imresize(mask,t_x.shape,interp='lanczos') # this succs by the way

    # some use full pre clean up
    mask = remove_small_holes(mask, 100)
    mask = remove_small_objects(label(mask), 1000) > 0  # removing other small bits


    mask_int=interpolation(mask,t_x)
    #plt.figure()
    #plt.imshow(mask)
    #plt.figure()
    #plt.imshow(mask_int)

    #plt.figure()
    #plt.imshow(t_x)

    ## mask data is prepared in exactly the same way as for stress analysis along lines
    mask_area, mask_boundaries=prepare_mask(mask_int)
    #plt.figure()
    #plt.imshow(mask_area)
    #plt.figure()
    #plt.imshow(mask_boundaries)

    # setting up files for finite elements
    # make matrix with ids
    nodes,elements,loads,mats=grid_setup(mask_area,t_x,t_y,1,0.3)
    # plot_grid(nodes,elements,inverted_axis=True)  # only use with <1000 nodes
