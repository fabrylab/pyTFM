
from __future__ import absolute_import, division, print_function
from solidspy import solids_GUI
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import time
import copy
import os
from datetime import datetime as dt
import solidspy.preprocesor as pre
import solidspy.postprocesor as pos
import solidspy.assemutil as ass
import solidspy.solutil as sol
import matplotlib
from pyTFM.plotting import *

folder="/media/user/GINA1-BK/data_traktion_force_microscopy/fem_example2/"
plt.close("all")


t1=time.time()
## generating input files
# regular grid of size n
n=100
c_x,c_y=np.meshgrid(range(int(np.sqrt(n))),range(int(np.sqrt(n))))


nodes=np.zeros((n,5))
nodes[:,0]=np.arange(n)
nodes[:,1]=c_x.flatten()
nodes[:,2]=c_y.flatten()
nodes=nodes.astype(int)

nodes[47,[3,4]]=-1 # constraint to avoid "all rigid body motion" ## do more research on this...
nodes[48,[3,4]]=-1
# choose point on y middle point, but in smaller half of x middle side (maybe ...)

n_el=int((np.sqrt(n)-1)**2) # number of elements
elements=np.zeros((n_el,7))



elements[:,0]=np.arange(n_el) #id
elements[:,1]=1 #element type/geometry (squares)
elements[:,2]=0 # elastic properties reference

sqr=[(c_x[:-1,:-1],c_y[:-1,:-1]), (c_x[:-1,:-1]+1,c_y[:-1,:-1]), (c_x[:-1,:-1]+1,c_y[:-1,:-1]+1),(c_x[:-1,:-1],c_y[:-1,:-1]+1)] #needs "counter clockwise orientation
sqr=[(x.flatten(),y.flatten()) for x,y in sqr]



ids=np.zeros(c_x.shape)   # array where all coordinates are attributed their id
ids[c_x.flatten(),c_y.flatten()]=np.arange(n,dtype=int)




elements[:,3]=ids[sqr[0][0],sqr[0][1]]  ## zusammenfassen
elements[:,4]=ids[sqr[1][0],sqr[1][1]]
elements[:,5]=ids[sqr[2][0],sqr[2][1]]
elements[:,6]=ids[sqr[3][0],sqr[3][1]]
elements=elements.astype(int)
# writing to files


loads=np.zeros((n,3))
loads[:,0]=np.arange(n)


#for i in range(10,110,10): # forces with non zero moments
#    if i<=50:
#        loads[i-10:i-5,2]=-0.001
#    else:
#        loads[i-5:i,2]=0.001




#loads[20:29,2]=-0.001
#loads[70:79,2]=+0.001

x_inds1=[np.arange(5)+i*10 for i in range(10)]
x_inds2=[np.arange(5,10)+i*10 for i in range(10)]
loads[x_inds1, 1] = 0.001
loads[x_inds2, 1] = -0.001
fx=np.zeros((10,10))
fy=np.zeros((10,10))

fx=np.zeros((10,10))
fy=np.zeros((10,10))
fx[nodes[:, 2].astype(int), nodes[:, 1].astype(int)]=loads[:,1]
fy[nodes[:, 2].astype(int), nodes[:, 1].astype(int)]=loads[:,2]
fig=show_quiver(fx, fy,filter=[0,0], scale_ratio=0.05,headwidth=10,headlength=6)
#loads[0,2]=-0.001
#loads[-1,2]= 0.001
#loads[0,1]=-0.001
#loads[-1,1]= 0.001












np.savetxt(os.path.join(folder,"nodes.txt"),nodes,fmt="%i",delimiter="\t")
np.savetxt(os.path.join(folder,"loads.txt"),loads,fmt="%.3f",delimiter="\t")
np.savetxt(os.path.join(folder,"eles.txt"),elements,fmt="%i",delimiter="\t")
t2=time.time()
print("grid assembly ", t2-t1)

start_time = dt.now()
echo = False
compute_strains=True
plot_contours=False



nodes, mats, elements, loads = pre.readin(folder=folder)
if echo:
    pre.echomod(nodes, mats, elements, loads, folder=folder)
DME, IBC, neq = ass.DME(nodes, elements)  # boundary conditions asembly??
# DME
#IBC list of "equations in x y direction for each node, -1 indicates no equations of movement, as completly fixed --> related to boundary conditions
# neq total number of equations
# DME is like list of independent equations per element(square, so just rearanged IBC for elements,side node zeros are at the end because of fixed length( for other element shape=
#
#
# print("Number of nodes: {}".format(nodes.shape[0]))
print("Number of elements: {}".format(elements.shape[0]))
print("Number of equations: {}".format(neq))



# System assembly
t1=time.time()
KG = ass.assembler(elements, mats, nodes, neq, DME,sparse=False)
RHSG = ass.loadasem(loads, IBC, neq)

u=np.linalg.solve(KG,RHSG)
test=np.matmul(KG,u)-RHSG
#np.allclose(np.dot(KG, u), RHSG)



t2=time.time()
print("sys assembly ",t2-t1)

# System solution

t1=time.time()
UG = sol.static_sol(KG, RHSG)  #   nodes added by me
if not (np.allclose(KG.dot(UG) / KG.max(), RHSG / KG.max())):
    print("The system is not in equilibrium!")

t2=time.time()
print("system solution ", t2-t1)

end_time = dt.now()
#print('Duration for system solution: {}'.format(end_time - start_time))

# Post-processing
start_time = dt.now()
UC = pos.complete_disp(IBC, nodes, UG)  # uc are x and y displacements
#UC2 = pos.complete_disp(IBC, nodes, UG2)  # uc are x and y displacements


E_nodes, S_nodes = None, None
if compute_strains:
    E_nodes, S_nodes = pos.strain_nodes(nodes, elements, mats, UC)
    #E_nodes2,S_nodes2=pos.strain_nodes(nodes, elements, mats, UC2)
if plot_contours:
    pos.fields_plot(elements, nodes, UC, E_nodes=E_nodes, S_nodes=S_nodes) # some matplotlib internal function to make it look smoother
end_time = dt.now()



# calculation of principal stresses
## are signes importnatn here???

# chek this and others
## https://wp.optics.arizona.edu/optomech/wp-content/uploads/sites/53/2016/10/OPTI_222_W21.pdf
sig_x=S_nodes[:,0]
sig_y=S_nodes[:,1]
tau_xy=S_nodes[:,2]

sigma_max=(sig_x+sig_y)/2+np.sqrt(((sig_x-sig_y)/2)**2+tau_xy**2)
sigma_min=(sig_x+sig_y)/2-np.sqrt(((sig_x-sig_y)/2)**2+tau_xy**2)
# maximum shear stress
tau_max=np.sqrt(((sig_x-sig_y)/2)**2+tau_xy**2)
# angle of maximal principal stress
phi_n= np.arctan(2*tau_xy/(sig_x-sig_y))/2
# angel of maximal shear stress
phi_shear= np.arctan(-(sig_x-sig_y)/(2*tau_xy))/2
## side note (phi_n-phi_shear)=pi/4 (also immer 45° unterschied
sigma_avg=(sigma_max+sigma_min)/2# average (maximal?)normal stress  ## thi according to buttler paper

stress_tensor=np.zeros((int(np.sqrt(len(nodes))),int(np.sqrt(len(nodes))),2,2))
stress_tensor[nodes[:, 2].astype(int), nodes[:, 1].astype(int),0,0]=S_nodes[:,0] #sigma_x
stress_tensor[nodes[:, 2].astype(int), nodes[:, 1].astype(int),1,1]=S_nodes[:,1] #sigma_y
stress_tensor[nodes[:, 2].astype(int), nodes[:, 1].astype(int),1,0]=S_nodes[:,2] #sigma_yx
stress_tensor[nodes[:, 2].astype(int), nodes[:, 1].astype(int),0,1]=S_nodes[:,2] #sigma_xy



# plottingloads
plot_fields(nodes,fields=[loads[:,1],loads[:,2]],titles=["x_traction_force","y_traction_force"],cbar_str="force",grid_lines=True)
#plot_fields(nodes,fields=[loads[:,1],loads[:,2]],titles=["",""],cbar_str="force",grid_lines=True)


# plotting deformation
plot_fields(nodes,fields=[UC[:,0],UC[:,1]],titles=["x_deformation","y_deformation"],cbar_str="deformation")
#plot_fields(nodes,fields=[UC[:,0],UC[:,1]],titles=["",""],cbar_str="deformation")


# plotting straints
plot_fields(nodes,fields=[E_nodes[:,0],E_nodes[:,1],E_nodes[:,2]],titles=["x_strain","y_strain","xy_strain"],cbar_str="strain")

# plotting stresses
plot_fields(nodes,fields=[S_nodes[:,0],S_nodes[:,1],S_nodes[:,2]],titles=["x_stress","y_stress","xy_stress"],cbar_str="stress")


#plot_fields(nodes,fields=[sigma_avg,sigma_max,sigma_min,tau_max],titles=["average normal stress","max normal stress","min normal stress","max shear stress"],cbar_str="stress")



'''
#example for boundaries
mask=np.zeros((int(np.sqrt(len(nodes))),int(np.sqrt(len(nodes)))))
#mask[[4,4,4,5,6,4,7,7,7,7,7],[5,6,7,4,3,8,4,5,6,]]=1
#mask[[7,7,7,7,6,5],[3,4,5,6,7,8]]=2
# cicular shape test
from skimage.draw import circle,polygon
from skimage.morphology import binary_dilation

poly=polygon([29,35,35,29,21,15,15,21],[15,21,29,35,35,29,21,15])

mask[poly[0],poly[1]]=1
mask=binary_dilation(mask)*1-mask*1
mask[25,14]=0 # needs start point for line recognition
### fit with splines
# not so nice


# get normalvectors by using values 1 and 3      --> that looks rather good
# first sorted list
coords= np.array(np.where(mask)).T # list aof x ycoordinates of the respective poitns
dist=np.linalg.norm(coords[None,:]-coords[:,None],axis=2)  # only make

# making graph from distance matrix:
dist_bool=dist<=np.sqrt(2) #would be distance of diagonal point+
graph_l=defaultdict(list)
start_end=[]
for i,row in enumerate(dist_bool): # make with arrays only?
    p=list(np.where(row)[0])
    p.remove(i) # neighbouring points
    graph_l[i].extend(p)
    if len(p)==1:
        start_end.append(i)
# connecting the other points

#from functions_for_cell_colonie import find_path
order=find_path(graph_l,start_end[0],start_end[1])
coords=coords[order] #reordering
#check_order(mask,coords)


def calculate_line_stresses(coords,stress_tensor):
    '''
'''
    :param coords: orderd set of coordinates of points in line, start and end must be at 0 and -1 potstion respectively
    :param stress_tensor: stress tensor on the complete field
    :return: n_stress,shear_stress, normal and shear stresses for each point, in the smae order (not for first and last point of the line)
    n: list of normal vectors
    stress_vec stress vectors for each point as calculated from cauchy therorem
    '''
'''

    #normal vectors:
    # vector between point 1 and 3
    vec=coords[:-2]-coords[2:]
    n=copy.deepcopy(vec)#vec[:,[1,0]]
    n[:,1]=-n[:,1]
    n=n/np.linalg.norm(n,axis=1)[:, np.newaxis]
    #calcualting normal and sheer stress along the line
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
    return n,stress_vec,n_stress,shear_stress
n,stress_vec,n_stress,shear_stress=calculate_line_stresses(coords,stress_tensor)
plot_line_stresses(mask,coords,n_stress,shear_stress)


check_normal_vectors(mask,coords,n)
check_normal_vectors(mask,coords,stress_vec)







print('Duration for post processing: {}'.format(end_time - start_time))
print('Analysis terminated successfully!')

#plt.show()




'''
