
from andreas_TFM_package.solids_py_stress_functions import *
from skimage.morphology import skeletonize,remove_small_holes,remove_small_objects,label,binary_dilation,binary_erosion
from scipy.ndimage.morphology import binary_fill_holes
from scipy.ndimage.measurements import find_objects
from scipy.ndimage import binary_closing as binary_clo
from scipy.signal import convolve2d
from scipy.interpolate import splprep, splev
from itertools import chain
from collections import Counter

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
    m = mask_area - mask_boundaries
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
        x= np.concatenate([[endpoints[0][1]],x,[endpoints[1][1]]])
        y = np.concatenate([[endpoints[0][0]], y, [endpoints[1][0]]])




    tck, u = splprep([x, y], s=10, k=k)  ### parametric spline interpolation
    # fits essentially function: [x,y] =f(t) , t(paramter) is default np.linspace(0,1,len(x))
    # tck is array with x_position of knot, y_position of knot, parameters of the plne, order of the spline
    # ui is paramter given to x,y points, in this case default (values from 0 to 1)
    # s is smootjing factor(?)
    # k is order of the spline, cubic is fine



    if endpoints: # forcing splines to through endpoints
        tck[1][0][0]=endpoints[0][1]
        tck[1][0][-1] =endpoints[1][1]
        tck[1][1][0] =endpoints[0][0]
        tck[1][1][-1] =endpoints[1][0]



    points_new = np.array(splev(u, tck,
                                der=0)).T  # new points from spline interpolation (thes points will be used for normal stress vector calculation
    # in xy orientation
    points_new = np.round(points_new).astype(int)  ## could also use exact values and the interpolate the stress tensor, but thats probably to complicated
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

class Cells_and_Lines:
    # containter for cells and lines, assignement of points with them, and assignement with each other
    def __init__(self, mask_area,mask_boundaries):


        self.mask_area=mask_area
        self.mask_boundaries=mask_boundaries

        # graph as a dictinonray with key=ponit id, values: ids of neighbouring points
        # any point id is the index in the points array (contains coordinate of these points
        self.graph,self.points=mask_to_graph(mask_boundaries)

        # points as dictionary with key=points id, values: points coordinates
        self.points_dict={i:self.points[i] for i in range(len(self.points))}

        # lines as a dictinonray with key=line id, values: ids of containg points (in correct order)
        self.lines_points=identify_line_segments(self.graph,self.points)

        # dictionary with endpoints, needed to completely fill the gaps between all cell_lines
        self.lines_endpoints_com,self.lines_endpoints = find_exact_line_endpoints(self.lines_points, self.points, self.graph)

        # removing all lines that are predicted to be circular. Mostly a problem for very short lines
        remove_circular_line(self.lines_endpoints_com,self.lines_points,self.lines_endpoints)

        # cells as a dictinonray with key=cell id, values: ids of containing points (not ordered)
        self.cells_points, self.cells_area=identify_cells(mask_area, mask_boundaries, self.points)
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

        # list of ids
        self.point_ids = list(self.points_dict.keys())
        self.line_ids=list(self.cells_points.keys())
        self.cell_ids=list(self.lines_points.keys())
        # list of all lines at the edge of the cell colony
        self.edge_lines = find_edge_lines(self.cells_lines)

        # dictionary containing the spline represetation of the points as a parametric function
        #[x,y]=f(u). u is always np.linspace(0,1,"number of points in the line). Use scipy.interpolate.splev
        # to evaluate at other points
        self.lines_splines = defaultdict(list)

        # dictionary with the normal vectors for each line acording to the spline interpolation. Contains list
        # of the normal vectors , position of these points is listed in lines_spline_points
        # interpolation !!! positions are given in xy order !!!
        self.lines_n_vectors= defaultdict(list)

        # dictinary with example points (where a normal vector originates) from spline interpolation
        self.lines_spline_points = defaultdict(list)

        for  line_ids, line_ps in self.lines_points.items():
            endpoints= self.lines_endpoints_com[line_ids]
            k = len(line_ps)+2 - 1 if len(line_ps) <= 3 else 3 # addapt spline order, according to number of points
            # spline order must be below nomber of points, so choose 2 if lne has one point + 2 endpoints
            tck, u, points_new=spline_interpolation(line_ps, self.points,k=k,endpoints=endpoints) #spline interpolation
            self.lines_splines[line_ids]=tck # saving the spline obejct
            # saving a few oints and n vectors for easy representations/ debugging,
            # these points will not be used further
            n_vectors=normal_vectors_from_splines(u, tck) # calculating normal vectors as a list
            self.lines_n_vectors[line_ids]=n_vectors
            self.lines_spline_points[line_ids] = points_new

    # self.filter_single_endpoints()
    #def filter_single_endpoints(self):
        #'''
        #finds lines that are circular. this is a worka round. Remove if interpolation /mask shrinking isue is resolved
        #:param self:
        #:return:
        #'''
        #single_endpoints=[line_id for line_id,(ep1,ep2) in self.lines_endpoints_com.items() if np.linalg.norm(ep1-ep2)<0.2]
        #self.lines_endpoints_com={line_id:value for line_id,value in self.lines_endpoints_com.items() if line_id not in single_endpoints}
        #self.lines_endpoints = {line_id: value for line_id, value in self.lines_endpoints.items() if line_id not in single_endpoints}
        #self.lines_points = {line_id: value for line_id, value in self.lines_points.items() if line_id not in single_endpoints}


        plt.figure()
        plt.imshow(self.mask_boundaries)
        for l_id,points in self.lines_points.items():
            for p in points:
                plt.text(self.points[p][1],self.points[p][0],str(l_id))




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
        colors = []
        for i in range(len(self.lines_points.keys())):
            colors.append('#%06X' % np.random.randint(0, 0xFFFFFF))
        offset = 0.005 # ofset for point id text
        fig=plt.figure()
        plt.imshow(self.mask_boundaries,cmap="jet")

        plt.plot([],[],color="green",label="line_ids") #legend for line ids
        plt.plot([], [], color="red", label="cells_ids") #legend for cells ids

        # plotting all points with line id as label
        for l, ps in self.lines_spline_points.items():
            if l in self.edge_lines:
                color="orange"
            else:
                color = "green"
            p_indices=np.array(list(range(len(ps)))) #all indicces
            # randomly choising a few indices, without first and laast index
            ps_select=np.random.choice(p_indices[1:-1],size=int((len(ps)-2)*sample_factor),replace=False)
            ps_select=np.append(ps_select,p_indices[np.array([1,-1])]) # adding back first and last indeyx

            plt.plot(ps[:, 0], ps[:, 1]) #plotting all points
            for p in ps[ps_select]: # labeling selected points
                plt.text(p[0] + 1 * offset * l, p[1] + 1 * offset * l, s=str(l),color=color)
        # plotting cel id at center of mass of cell
        for cell_id,com in self.cells_com.items():
            plt.text(com[1],com[0],str(cell_id),color="red",fontsize=13)

        if plot_n_vectors:
            for points,n_vectors in zip(self.lines_spline_points.values(),self.lines_n_vectors.values()):
                for n_vec,p in zip(n_vectors,points):
                    plt.arrow(p[0], p[1], n_vec[0], n_vec[1], head_width=0.15)

        plt.legend()
        return fig


    #def line_stresses(self,n_array,stress_tensor):
    #    '''
    #    calculates line stresses and nomr of the line stresses form the normal vectors as an array.
    #    :param n_array:
    #    :param stress_tensor:
    #    :return:
    #    '''
    #
     #   self.ls_array=calculate_stress_vector(n_array, stress_tensor)
    #   self.ls_norm_array=np.linalg.norm(self.ls_array,axis=2)




def prepare_mask(mask,min_cell_size=None):
    '''
    this function usese skeletonize to transform the cellboundraies to one pixel width. Loose ends are trimmed by
    converison to a graph and then deleting all nodes with only one neighbour.
    :param mask:
    :param min_cell_size: minimal sze of cells in pixles, any hole below that will be filled up. If none, then
        some estimated value is used
    :return:
    '''

    if not min_cell_size:
        min_cell_size=mask.shape[0]*mask.shape[1]/1000 ## is this robust????

    min_cell_size=min_cell_size if min_cell_size>2 else 2  # force to be at least 2
    mask = remove_small_holes(mask.astype(bool), min_cell_size)
    mask = remove_small_objects(mask.astype(bool), min_cell_size*10) > 0  # removing other small bits
    mask = skeletonize(mask)

    # converting mask to graph object
    graph,points = mask_to_graph(mask)

    # removing dead ends
    end_points = np.where(np.array([len(v) for v in graph.values()]) < 2)[0]
    points_2 = np.array(list(graph.keys()))  # all keys as array
    eps = points_2[end_points]  # keys of endpoints
    for ep in eps:
        if ep in list(graph.keys()):
            remove_endpoint(graph, ep)


    mask_boundaries = graph_to_mask( graph, points, mask.shape)  # rewriting mask with deleted points

    mask_area = binary_fill_holes(mask_boundaries)


       # removing unsuitable pixels for grid later  # not so very smoooth  /maybe improve
    m=convolve2d(mask_area,np.array([[0,1,0],[1,0,1],[0,1,0]]),mode="same",boundary="fill",  fillvalue=0) # convoultion will produce 1 if only one direct (distance 1) neighbour exists
    p=np.where(np.logical_and(m==1,mask_boundaries)) # problematic coordinates

    for x,y in zip(p[0],p[1]):
        new_ps=np.array([[x+1,y],[x,y+1],[x-1,y],[x,y-1]])# all relevant neigbouring points
        new_p=new_ps[mask_area[new_ps[:,0],new_ps[:,1]]][0]  # checking where a possible replacement point could be located this can only have one result
        mask_boundaries[x,y]=False # removing old point
        mask_boundaries[new_p[0],new_p[1]]=True # adding new point
        mask_area[x,y]=False
        mask_area[new_p[0], new_p[1]] = True

    # updating the garph (coul be done more efficinetly, but this is pretty fast anyway
    graph,points=mask_to_graph(mask_boundaries)

    c_l=Cells_and_Lines(mask_area,mask_boundaries)

    #c_l.vizualize_lines_and_cells(sample_factor=0.5,plot_n_vectors=True)


    ## vizualization of spline interpolation with new line enpoints





    return mask_area, mask_boundaries,c_l


def prepare_mask_spread(mask,dil_factor=1.2):
    '''
   same as prepare_mask, but adds an adddtional area arund the actual mask.
   this function is inteded for tessting only.

    :param mask:
    :return:
    '''

    min_cell_size = mask.shape[0] * mask.shape[1] / 1000  ## is this robust????
    mask = remove_small_holes(mask, min_cell_size)
    mask = remove_small_objects(label(mask), min_cell_size * 10) > 0  # removing other small bits
    mask = skeletonize(mask)


    graph = defaultdict(list)
    points = np.array(np.where(mask)).T
    point_tree = cKDTree(points)  # look up table for nearest neigbours ??
    d = np.sqrt(2)  # maximal allowd distance
    for i, p in enumerate(points):
        neighbours = point_tree.query_ball_point(p, d)
        neighbours.remove(i)  # removing the point itself from list of its neighbours
        graph[i].extend(neighbours)

    end_points = np.where(np.array([len(v) for v in graph.values()]) < 2)[0]
    points_2 = np.array(list(graph.keys()))  # all keys as array
    eps = points_2[end_points]  # keys of endpoints
    for ep in eps:  ## do i even need this???
        if ep in list(graph.keys()):
            remove_endpoint(graph, ep)

    mask_boundaries = graph_to_mask( graph, points,mask.shape)  # rewriting mask with deleted points

    mask_area = binary_fill_holes(mask_boundaries)

    mask_dil = interpolation(mask_area, (int(mask_area.shape[0] * dil_factor), int(mask_area.shape[1] * dil_factor))) # intrerpolate to larger size by a factor
    mask_dil=alligne_objects(mask_area,mask_dil) # cut out and place "ontop" of previous mask

    # removing unsuitable pixels for grid later  for area mask

    m = convolve2d(mask_area, np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]]), mode="same", boundary="fill",
                   fillvalue=0)  # convolution will produce 1 if only one direct (distance 1) neighbour exists
    p = np.where(np.logical_and(m == 1, mask_boundaries))  # problematic coordinates, only at the edge of the amsk

    for x, y in zip(p[0], p[1]):
        new_ps = np.array([[x + 1, y], [x, y + 1], [x - 1, y], [x, y - 1]])  # all relevant neigbouring points
        new_p = new_ps[mask_area[new_ps[:, 0], new_ps[:, 1]]][
            0]  # checking where a possible replacement point could be located this can only have one result
        mask_boundaries[x, y] = False  # removing old point
        mask_boundaries[new_p[0], new_p[1]] = True  # adding new point
        mask_area[x, y] = False
        mask_area[new_p[0], new_p[1]] = True

    graph, points = mask_to_graph(mask_boundaries)
    c_l = Cells_and_Lines(mask_area, mask_boundaries)



    #for border mask
    m = convolve2d(mask_dil, np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]]), mode="same", boundary="fill",
                   fillvalue=0)  # convoultion will produce 1 if only one direct (distance 1) neighbour exists
    p = np.where(np.logical_and(m == 1, mask_dil-binary_erosion(mask_dil)))  # problematic coordinate


    mask_dil=mask_dil.astype(bool)
    for x, y in zip(p[0], p[1]):
        new_ps = np.array([[x + 1, y], [x, y + 1], [x - 1, y], [x, y - 1]])  # all relevant neigbouring points
        new_p = new_ps[mask_dil[new_ps[:, 0], new_ps[:, 1]]][0]  # checking where a possible replacement point could be located this can only have one result
        mask_dil[x, y] = False
        mask_dil[new_p[0], new_p[1]] = True





    return mask_area, mask_boundaries,mask_dil, c_l


def interpolation(mask, dims, min_cell_size=100):
    #
    # some pre clean up of the mask
    mask = remove_small_holes(mask.astype(bool), min_cell_size)
    mask = remove_small_objects(mask.astype(bool), 1000) # removing other small bits
    # note: remove_small_objects labels automatically if mask is bool
    coords = np.array(np.where(mask)).astype(float)
    interpol_factors = np.array([dims[0] / np.shape(mask)[0], dims[1] / np.shape(mask)[1]])
    coords[0] = coords[0] * interpol_factors[0]  # interpolting x coordinates
    coords[1] = coords[1] * interpol_factors[1]  # interpolating xy coordinates
    coords = np.round(coords).astype(int)
    # test if coorect
    mask_int = np.zeros(dims)
    mask_int[coords[0], coords[1]] = 1
    mask_int = mask_int.astype(int)
    if dims[0]*dims[1]>=mask.shape[0]*mask.shape[1]:
        iter=int(np.ceil(np.max([mask.shape[0]/dims[0],mask.shape[0]/dims[0]]))*5) # times 5 is safety factor
        mask_int=binary_clo(mask_int,iterations=10)
        print(iter)
    return mask_int


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


def grid_setup(mask_area, f_x, f_y, E, sigma):
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
    nodes = nodes.astype(int)
    # no constraints for any node here



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