# function to simulate the surface deformation from point forces applied to a
# half sphere. Also includes a function a finetely thik halfsphere is considered.





from tqdm import tqdm
from scipy.integrate import quad
from scipy.special import jv # first order bessel function
import sys
sys.path.insert(0, '/media/user/GINA1-BK/Andreas-Python/tracktion_force_microscopy/')
from pyTFM.TFM_functions import *
from pyTFM.plotting import *
from scipy.ndimage.morphology import distance_transform_edt
from scipy.signal import convolve2d, fftconvolve
from skimage.filters import gaussian
import threading, queue
import time
try:
    from numba import jit, objmode, cuda
except:
    pass

def threading_wrapper(qu,funct,*args,**kwargs):
    '''

    :param qu: queue.Queue object
    :param funct:
    :param args:
    :param kwargs:
    :return:
    '''

    qu.put(funct(*args,**kwargs))

def execute_as_thread(functions,arg_list,kwargs_list):
    #### Python doesnt work that way??? ######
    threads=[]
    queues=[]
    for i,(funct,args,kwargs) in enumerate(zip(functions,arg_list,kwargs_list)):
        qu = queue.Queue()
        new_thread = threading.Thread(target=threading_wrapper, args=(qu,funct)+args,kwargs=kwargs)
        new_thread.start()
        threads.append(new_thread)
        queues.append(qu)
    for thread in threads: # waiting for all threads to finish
        thread.join()
    return [qu.get() for qu in queues]




def reverse_kernel(k): # not necessary because kernel is symmetric anyway
    '''
    "flips kernel by 180 degrees"
    '''
    k=np.flip(k,axis=0)
    k=np.flip(k,axis=1)
    return k


def get_xy_for_quiver(u):
    xs=np.zeros(np.shape(u))
    for i in range(np.shape(u)[0]):
        xs[i,:]=np.arange(0,np.shape(u)[1],1)
    ys = np.zeros(np.shape(u))
    for j in range(np.shape(u)[1]):  ## is inverted in other skript...
        ys[:, j] = np.arange(0,np.shape(u)[0], 1)
    return xs, ys



def force_from_center(Force_magnitude,pos1,pos2):
    new_center = np.array([np.mean([pos1[0], pos2[0]]), np.mean([pos1[1], pos2[1]])])  # new exact center position
    F1 = (pos1 - new_center) * Force_magnitude / np.linalg.norm(pos1 - new_center)  # first force
    F2 = (pos2 - new_center) * Force_magnitude / np.linalg.norm(pos2 - new_center)  # first force
    return F1,F2
#plt.close("all")


def return_force_coordinate_with_offset(fx):
    '''
    deformation exactly at force origin is infinetly large. (point force is not real physical appearance)
    therefore it is necassary to let the force act on points between the actual positions
    :param fx:
    :return:
    '''
    pixx = np.arange(fx.shape[0])
    pixy = np.arange(fx.shape[1])
    dist_x, dist_y = np.meshgrid(pixy, pixx)
    dist_x=dist_x.astype("float64")
    dist_y=dist_y.astype("float64")
    dist_x+=0.5
    dist_y+=0.5
    return dist_x,dist_y



def infinite_thickness_convolution_kernel(fx,fy,pixelsize,young,sigma,kernel_size=None,force_shift=0.5):

    if not kernel_size:
        kernel_size = fx.shape

    pixx = np.arange(kernel_size[1])
    pixy = np.arange(kernel_size[0])
    dist_x, dist_y = np.meshgrid(pixy, pixx)  # matrices of distance
    dist_x = dist_x.astype(float)
    dist_y = dist_y.astype(float)
    # makeing sure coorinates are shifted and distance is never exactly zero ... (maybe just exlude zero????)
    if np.mean(dist_x[0, :]).is_integer() or np.mean(dist_y[:, 0]).is_integer():
        dist_x -= np.mean(dist_x[0, :])   # centering
        dist_y -= np.mean(dist_y[:, 0])
    else:
        dist_x -= np.mean(dist_x[0, :])  # centering
        dist_y -= np.mean(dist_y[:, 0])

    dist_x = (dist_x+force_shift) * pixelsize
    dist_y = (dist_y+force_shift) * pixelsize
    r=np.sqrt(dist_x**2+dist_y**2)

    # convolution kernels
    A = (1 + sigma) / (np.pi * young)  # a constant factor for boussinesq kernel
    conv_k1=((1 - sigma) * (r ** 2) + sigma * (dist_x ** 2)) * (
                    A / (r ** 3))  # components of the boussinesq kernel (? nameing?)
    conv_k2= sigma * dist_x * dist_y * (A / (r ** 3))
    conv_k3=((1 - sigma) * (r ** 2) + sigma * (dist_y ** 2)) * (A / (r ** 3))

    return ([conv_k1,conv_k2,conv_k3],[A])



def finite_thickness_convolution_greens_tensor(fx, fy,pixelsize, h, young, sigma=0.5,kernel_size=None,force_shift=0):
    '''
    calculating the deformation fieldgiven a tracktion force field using the boussinesq solution for an
     half sphere with finite thickness
     # employs the approximate greens tensor from

     Cell Force Microscopy on Elastic Layers of Finite Thickness
     Rudolf Merkel, Norbert Kirchgeßner, Claudia M. Cesa, and Bernd Hoffman
     Biophysical Journal  Volume 93  November 2007

     This approximation is only spposed to be used for sigma=0.5

    :param fx: forces in x direction
    :param fy: forces iny direction
    :param fx_coords: coordinates of these forces shifted by 0.5 pixel. shift is necessayr because we assume point load
    :param fy_coords:
    :param young: youngs modulus of material
    :param sigma: poisson ratio of material
    :return:
    '''

    # distance for the convolution kernel

    # greens tensor with:
    # [K1,K2,K3]
    # [K4,K5,K6]
    # [K7,K8,K9]
    # returns displacement also in z direction (makes sense because of "pulling" at the buttom of the layer)

    E = young

    if not kernel_size:
        kernel_size=fx.shape

    pixx = np.arange(kernel_size[1])
    pixy = np.arange(kernel_size[0])
    dist_x, dist_y = np.meshgrid(pixy, pixx)  # matrices of distance
    dist_x = dist_x.astype(float)
    dist_y = dist_y.astype(float)
    # makeing sure coorinates are shifted and distance is never exactly zero ... (maybe just exlude zero????)
    if np.mean(dist_x[0,:]).is_integer() or np.mean(dist_y[:, 0]).is_integer():
        dist_x -= np.mean(dist_x[0, :])-0.5  # centering
        dist_y -= np.mean(dist_y[:, 0])-0.5
    else:
        dist_x -= np.mean(dist_x[0,:]) # centering
        dist_y -= np.mean(dist_y[:, 0])

    dist_x = (dist_x+force_shift)* pixelsize
    dist_y = (dist_y+force_shift)* pixelsize

    r=np.sqrt(dist_x**2+dist_y**2)


    ### elements of the greens tensor
    # Usefull parts // this is just an approximation
    # /also only for sigma=0.49
    s = r / h
    mu = E / (2 * (1 + sigma))
    A1 = ((2 - sigma) / (4 * np.pi * mu * h * s)) * (0.12 * np.exp(-0.43 * s) + 0.88 * np.exp(-0.83 * s))
    A2 = -(sigma / (4 * np.pi * mu * h * s)) * (1 + 1.22 * s + 1.31 * (s ** 2.23)) * np.exp(-1.25 * s)
    #A3_b = ((1-2*sigma)/(4*np.pi*mu*h*s))
    A3 = -0.063 * ((np.exp(-0.44 * s) - np.exp(-2.79 * s)) ** 2)#*A3_b
    A4 = ((2 * (1 - sigma)) / (4 * np.pi * mu * h * s)) * (1 + 0.46 * s - 2.50 * (s ** 2.13)) * np.exp(-2.18 * s)

    # components of the greens tensor
    K1 = A1 - ((dist_x ** 2 - dist_y ** 2) / (r ** 2)) * A2
    K2 = -(2 * dist_x * dist_y / (r ** 2)) * A2
    K3 = -(dist_x / r) * A3
    K4 = -(2 * dist_x * dist_y / (r ** 2)) * A2  # is K2
    K5 = A1 + ((dist_x ** 2 - dist_y ** 2) / (r ** 2)) * A2
    K6 = -(dist_y / r) * A3
    K7 = (dist_x / r) * A3
    K8 = (dist_y / r) * A3
    K9 = A4

    return ([K1,K2,K3,K4,K5,K6,K7,K8,K9],[A1,A2,A3,A4])




def finite_thickness_convolution_greens_tensor_sparse(dist_x, dist_y, pixelsize, h, young, sigma=0.5, force_shift=0.5):
    '''
     Greenstensor component for a given list of distances in x and y dircetion
     # employs the approximate greens tensor from

     Cell Force Microscopy on Elastic Layers of Finite Thickness
     Rudolf Merkel, Norbert Kirchgeßner, Claudia M. Cesa, and Bernd Hoffman
     Biophysical Journal  Volume 93  November 2007

     This approximation is only spposed to be used for sigma=0.5

    :param fx_coords: coordinates of these forces shifted by 0.5 pixel. shift is necessayr because we assume point load
    :param fy_coords:
    :param young: youngs modulus of material
    :param sigma: poisson ratio of material
    :return:
    '''

    # distance for the convolution kernel

    # greens tensor with:
    # [K1,K2,K3]
    # [K4,K5,K6]
    # [K7,K8,K9]
    # returns displacement also in z direction (makes sense because of "pulling" at the buttom of the layer)

    E = young
    # makeing sure coordinates are shifted and distance is never exactly zero ... (maybe just exlude zero????)
    dist_x[dist_x % 1 == 0] = dist_x[dist_x%1 == 0] + force_shift
    dist_y[dist_y % 1 == 0] = dist_y[dist_y % 1 == 0] + force_shift

    dist_x = (dist_x)* pixelsize
    dist_y = (dist_y)* pixelsize

    r=np.sqrt(dist_x**2+dist_y**2)


    ### elements of the greens tensor
    # Usefull parts // this is just an approximation
    # /also only for sigma=0.49
    s = r / h
    mu = E / (2 * (1 + sigma))
    A1 = ((2 - sigma) / (4 * np.pi * mu * h * s)) * (0.12 * np.exp(-0.43 * s) + 0.88 * np.exp(-0.83 * s))
    A2 = -(sigma / (4 * np.pi * mu * h * s)) * (1 + 1.22 * s + 1.31 * (s ** 2.23)) * np.exp(-1.25 * s)
    #A3_b = ((1-2*sigma)/(4*np.pi*mu*h*s))
    #A3 = -0.063 * ((np.exp(-0.44 * s) - np.exp(-2.79 * s)) ** 2)#*A3_b
    #A4 = ((2 * (1 - sigma)) / (4 * np.pi * mu * h * s)) * (1 + 0.46 * s - 2.50 * (s ** 2.13)) * np.exp(-2.18 * s)

    # components of the greens tensor
    K1 = A1 - ((dist_x ** 2 - dist_y ** 2) / (r ** 2)) * A2
    K2 = -(2 * dist_x * dist_y / (r ** 2)) * A2
    #K3 = -(dist_x / r) * A3
    K4 = -(2 * dist_x * dist_y / (r ** 2)) * A2  # is K2
    K5 = A1 + ((dist_x ** 2 - dist_y ** 2) / (r ** 2)) * A2
   # K6 = -(dist_y / r) * A3
   # K7 = (dist_x / r) * A3
   # K8 = (dist_y / r) * A3
   # K9 = A4
    return [K1,K2,K4,K5]
 #   return ([K1,K2,K3,K4,K5,K6,K7,K8,K9],[A1,A2,A3,A4])

def finite_thickenss_convolution_exact_greens_tensor(fx, fy,pixelsize, h, young, sigma=0.5,kernel_size=None):
    '''
    calculating the deformation fieldgiven a tracktion force field using the boussinesq solution for an
     half sphere with finite thickness
     # employs the approximate greens tensor from

     Cell Force Microscopy on Elastic Layers of Finite Thickness
     Rudolf Merkel, Norbert Kirchgeßner, Claudia M. Cesa, and Bernd Hoffman
     Biophysical Journal  Volume 93  November 2007

     This is the "exact" solution

    :param fx: forces in x direction
    :param fy: forces iny direction
    :param fx_coords: coordinates of these forces shifted by 0.5 pixel. shift is necessayr because we assume point load
    :param fy_coords:
    :param young: youngs modulus of material
    :param sigma: poisson ratio of material
    :return:
    '''

    # distance for the convolution kernel

    # greens tensor with:
    # [K1,K2,K3]
    # [K4,K5,K6]
    # [K7,K8,K9]
    # returns displacement also in z direction (makes sense because of "pulling" at the buttom of the layer)


    if not kernel_size:
        kernel_size = fx.shape

    pixx = np.arange(kernel_size[1])
    pixy = np.arange(kernel_size[0])
    dist_x, dist_y = np.meshgrid(pixy, pixx)  # matrices of distance
    dist_x = dist_x.astype(float)
    dist_y = dist_y.astype(float)
    # makeing sure coorinates are shifted and distance is never exactly zero ... (maybe just exlude zero????)
    if np.mean(dist_x[0, :]).is_integer() or np.mean(dist_y[:, 0]).is_integer():
        dist_x -= np.mean(dist_x[0, :]) - 0.5  # centering
        dist_y -= np.mean(dist_y[:, 0]) - 0.5
    else:
        dist_x -= np.mean(dist_x[0, :])  # centering
        dist_y -= np.mean(dist_y[:, 0])

    dist_x = dist_x * pixelsize
    dist_y = dist_y * pixelsize

    r = np.sqrt(dist_x ** 2 + dist_y ** 2)


    ### elements of the greens tensor
    mu = young / (2 * (1 + sigma))



    A1 = np.zeros(r.shape)
    A2 = np.zeros(r.shape)
    A3 = np.zeros(r.shape)
    A4=np.zeros(r.shape)
    ## filling the arrays for the green tensor at every r position (could be abit faster i guess with applying symetries...
    # aslo leave out A4


    # functions to be integrated when calcualting the greens tensor
    def integrand_A1(t,s):
        N = (3 - 4 * sigma) * np.exp(-4 * t) + (-24 * sigma + 10 + 4 * t ** 2 + 16 * sigma ** 2) * np.exp(-2 * t) + (
                    3 - 4 * sigma)
        phi1 = (-np.exp(-2 * t) / (N * (1 + np.exp(-2 * t)))) * (
                    2 * (sigma - 1) * t ** 2 + 2 * (1 - sigma) * t + 8 * sigma ** 3 - 20 * sigma ** 2 + 21 * sigma - 8 +
                    (2 * (sigma - 3) * t ** 2 + 2 * (
                                1 - sigma) * t + 8 * sigma ** 3 - 40 * sigma ** 2 + 48 * sigma - 18) * np.exp(-2 * t)
                    + (-4 * sigma ** 2 + 11 * sigma - 6) * np.exp(-4 * t))
        return (phi1 * jv(0, s * t))  # jv(0,t) is zero order besselfunction

    def integrand_A2(t,s):
        N = (3 - 4 * sigma) * np.exp(-4*t)+(-24*sigma+10+4*t**2+16*sigma**2)*np.exp(-2*t)+(3-4*sigma)

        phi2=(-np.exp(-2*t)/(N*(1+np.exp(-2*t))))*(2*(sigma-1)*t**2+2*(1-sigma)*t+8*sigma**3-20*sigma**2+13*sigma-2+
                                           (2*(sigma+1)*t**2+2*(1-sigma)*t+8*(sigma-1)*sigma**2+2)*np.exp(-2*t)
                                           +sigma*(3-4*sigma)*np.exp(-4*t))
        return(phi2*jv(2,s*t)) # jv(2,t) is second order besselfunction

    def integrand_A3(t,s):
        N = (3 - 4 * sigma) * np.exp(-4 * t) + (-24 * sigma + 10 + 4 * t ** 2 + 16 * sigma ** 2) * np.exp(-2 * t) + (
                    3 - 4 * sigma)

        phi3 = ((4 * (1 - sigma) * np.exp(-2 * t)) / N) * (t ** 2 + 2 * (2 * sigma - 1) * (sigma - 1))

        return (phi3 * jv(1, s * t))  # jv(1,t) is first order besselfunction

    def integrand_A4(t,s):
        N = (3 - 4 * sigma) * np.exp(-4 * t) + (-24 * sigma + 10 + 4 * t ** 2 + 16 * sigma ** 2) * np.exp(-2 * t) + (
                    3 - 4 * sigma)

        phi4 = ((2 * (1 - sigma) * np.exp(-2 * t)) / N) * (
                    2 * t * (t + 1) + 8 * sigma ** 2 - 12 * sigma + 5 + (3 - 4 * sigma) * np.exp(-2 * t))

        return (phi4 * jv(0, s * t))  # jv(0,t) is second order besselfunction



    for j,r_s in tqdm(enumerate(r.flatten()),total=r.shape[0]*r.shape[1]):
        s = r_s / h

        # "derivatory part to boussinesq solution"
        A1A=(-1/(2*np.pi*mu*h))*quad(integrand_A1,0,np.inf,args=s)[0]
        A2A=(-1/(2*np.pi*mu*h))*quad(integrand_A2,0,np.inf,args=s)[0]
        A3A=(-1/(2*np.pi*mu*h))*quad(integrand_A3,0,np.inf,args=s)[0]
        A4A=(-1/(2*np.pi*mu*h))*quad(integrand_A4,0,np.inf,args=s)[0]

        # part of the boussinesq solution
        A1B=(2-sigma)/(4*np.pi*mu*h*s)
        A2B=-sigma/(4*np.pi*mu*h*s)
        A3B=(1-2*sigma)/(4*np.pi*mu*h*s)
        A4B=(2*(1-sigma))/(4*np.pi*mu*h*s)
        index=np.unravel_index(j,r.shape)
        A1[index]=A1A+A1B
        A2[index]=A2A+A2B
        A3[index]=A3A+A3B
        A4[index]=A4A+A4B

    # componenets of the greesn tensor
    K1 = A1 - ((dist_x ** 2 - dist_y ** 2) / (r ** 2)) * A2
    K2 = -(2 * dist_x * dist_y / (r ** 2)) * A2
    K3 = -(dist_x / r) * A3
    K4 = -(2 * dist_x * dist_y / (r ** 2)) * A2  # is K2
    K5 = A1 + ((dist_x ** 2 - dist_y ** 2) / (r ** 2)) * A2
    K6 = -(dist_y / r) * A3
    K7 = (dist_x / r) * A3
    K8 = (dist_y / r) * A3
    K9 = A4
    return ([K1,K2,K3,K4,K5,K6,K7,K8,K9],[A1,A2,A3,A4])



def finite_thickenss_convolution_exact_greens_tensor_sparse(dist_x, dist_y, pixelsize, h, young, sigma=0.5, force_shift=0.5):
    '''
    calculating the deformation fieldgiven a tracktion force field using the boussinesq solution for an
     half sphere with finite thickness
     # employs the approximate greens tensor from

     Cell Force Microscopy on Elastic Layers of Finite Thickness
     Rudolf Merkel, Norbert Kirchgeßner, Claudia M. Cesa, and Bernd Hoffman
     Biophysical Journal  Volume 93  November 2007

     This is the "exact" solution

    :param fx: forces in x direction
    :param fy: forces iny direction
    :param fx_coords: coordinates of these forces shifted by 0.5 pixel. shift is necessayr because we assume point load
    :param fy_coords:
    :param young: youngs modulus of material
    :param sigma: poisson ratio of material
    :return:
    '''

    # distance for the convolution kernel

    # greens tensor with:
    # [K1,K2,K3]
    # [K4,K5,K6]
    # [K7,K8,K9]
    # returns displacement also in z direction (makes sense because of "pulling" at the buttom of the layer)

    E = young
    # makeing sure coordinates are shifted and distance is never exactly zero ... (maybe just exlude zero????)
    dist_x[dist_x % 1 == 0] = dist_x[dist_x % 1 == 0] + force_shift
    dist_y[dist_y % 1 == 0] = dist_y[dist_y % 1 == 0] + force_shift

    dist_x = (dist_x) * pixelsize
    dist_y = (dist_y) * pixelsize

    r = np.sqrt(dist_x ** 2 + dist_y ** 2)

    ### elements of the greens tensor
    mu = young / (2 * (1 + sigma))

    A1 = np.zeros(r.shape)
    A2 = np.zeros(r.shape)
   # A3 = np.zeros(r.shape)
    #A4 = np.zeros(r.shape)
    ## filling the arrays for the green tensor at every r position (could be abit faster i guess with applying symetries...
    # aslo leave out A4


    # functions to be integrated when calcualting the greens tensor
    def integrand_A1(t,s):
        N = (3 - 4 * sigma) * np.exp(-4 * t) + (-24 * sigma + 10 + 4 * t ** 2 + 16 * sigma ** 2) * np.exp(-2 * t) + (
                    3 - 4 * sigma)
        phi1 = (-np.exp(-2 * t) / (N * (1 + np.exp(-2 * t)))) * (
                    2 * (sigma - 1) * t ** 2 + 2 * (1 - sigma) * t + 8 * sigma ** 3 - 20 * sigma ** 2 + 21 * sigma - 8 +
                    (2 * (sigma - 3) * t ** 2 + 2 * (
                                1 - sigma) * t + 8 * sigma ** 3 - 40 * sigma ** 2 + 48 * sigma - 18) * np.exp(-2 * t)
                    + (-4 * sigma ** 2 + 11 * sigma - 6) * np.exp(-4 * t))
        return (phi1 * jv(0, s * t))  # jv(0,t) is zero order besselfunction

    def integrand_A2(t,s):
        N = (3 - 4 * sigma) * np.exp(-4*t)+(-24*sigma+10+4*t**2+16*sigma**2)*np.exp(-2*t)+(3-4*sigma)

        phi2=(-np.exp(-2*t)/(N*(1+np.exp(-2*t))))*(2*(sigma-1)*t**2+2*(1-sigma)*t+8*sigma**3-20*sigma**2+13*sigma-2+
                                           (2*(sigma+1)*t**2+2*(1-sigma)*t+8*(sigma-1)*sigma**2+2)*np.exp(-2*t)
                                           +sigma*(3-4*sigma)*np.exp(-4*t))
        return(phi2*jv(2,s*t)) # jv(2,t) is second order besselfunction

    def integrand_A3(t,s):
        N = (3 - 4 * sigma) * np.exp(-4 * t) + (-24 * sigma + 10 + 4 * t ** 2 + 16 * sigma ** 2) * np.exp(-2 * t) + (
                    3 - 4 * sigma)

        phi3 = ((4 * (1 - sigma) * np.exp(-2 * t)) / N) * (t ** 2 + 2 * (2 * sigma - 1) * (sigma - 1))

        return (phi3 * jv(1, s * t))  # jv(1,t) is first order besselfunction

    def integrand_A4(t,s):
        N = (3 - 4 * sigma) * np.exp(-4 * t) + (-24 * sigma + 10 + 4 * t ** 2 + 16 * sigma ** 2) * np.exp(-2 * t) + (
                    3 - 4 * sigma)

        phi4 = ((2 * (1 - sigma) * np.exp(-2 * t)) / N) * (
                    2 * t * (t + 1) + 8 * sigma ** 2 - 12 * sigma + 5 + (3 - 4 * sigma) * np.exp(-2 * t))

        return (phi4 * jv(0, s * t))  # jv(0,t) is second order besselfunction



    for j,r_s in tqdm(enumerate(r.flatten()), total=len(r.flatten())):
        s = r_s / h

        # "derivatory part to boussinesq solution"
        A1A=(-1/(2*np.pi*mu*h))*quad(integrand_A1,0,np.inf,args=s)[0]
        A2A=(-1/(2*np.pi*mu*h))*quad(integrand_A2,0,np.inf,args=s)[0]
        #A3A=(-1/(2*np.pi*mu*h))*quad(integrand_A3,0,np.inf,args=s)[0]
        #A4A=(-1/(2*np.pi*mu*h))*quad(integrand_A4,0,np.inf,args=s)[0]

        # part of the boussinesq solution
        A1B=(2-sigma)/(4*np.pi*mu*h*s)
        A2B=-sigma/(4*np.pi*mu*h*s)
        #A3B=(1-2*sigma)/(4*np.pi*mu*h*s)
        #A4B=(2*(1-sigma))/(4*np.pi*mu*h*s)
        index=np.unravel_index(j,r.shape)
        A1[index]=A1A+A1B
        A2[index]=A2A+A2B
        #A3[index]=A3A+A3B
        #A4[index]=A4A+A4B

    # componenets of the greesn tensor
    K1 = A1 - ((dist_x ** 2 - dist_y ** 2) / (r ** 2)) * A2
    K2 = -(2 * dist_x * dist_y / (r ** 2)) * A2
    #K3 = -(dist_x / r) * A3
    K4 = -(2 * dist_x * dist_y / (r ** 2)) * A2  # is K2
    K5 = A1 + ((dist_x ** 2 - dist_y ** 2) / (r ** 2)) * A2
    #K6 = -(dist_y / r) * A3
    #K7 = (dist_x / r) * A3
    #K8 = (dist_y / r) * A3
    #K9 = A4
    return [K1,K2,K4,K5]

def infinite_thickness_convolution(fx,fy,pixelsize,young,sigma,kernel_size=None,force_shift=0):  ## there should be a much better solution to this...
    '''
    calculating the deformation fieldgiven a tracktion force field using the boussiinesq solution for an
    infinte half sphere of an lineary isotropic elastic material
    :param fx: forces in x direction
    :param fy: forces iny direction
    :param fx_coords: coordinates of these forces shifted by 0.5 pixel. shift is necessary because we assume point load
    :param fy_coords:
    :param young: youngs modulus of material
    :param sigma: poisson ratio of material
    :return:
    '''
    fx=fx.astype(np.float128)
    fy=fy.astype(np.float128)
    # distance for the convolution kernel
    ([conv_k1, conv_k2, conv_k3], [A])=infinite_thickness_convolution_kernel(fx,fy,pixelsize,young,sigma,kernel_size=None,force_shift=force_shift)
    #deformation by convolution with the kernels
    def_x = convolve2d(fx, conv_k1, mode="same", boundary="fill", fillvalue=0) + convolve2d(fy, conv_k2, mode="same", boundary="fill", fillvalue=0)
    def_y = convolve2d(fx, conv_k2, mode="same", boundary="fill", fillvalue=0) + convolve2d(fy, conv_k3, mode="same", boundary="fill", fillvalue=0)

    return def_x,def_y




##### this not working jet #####
#@jit(["float32(float32,float32,float32,float32,float32,float32,float32)"], target="gpu")
def numba_loop(fx, fy, fx_coord, fy_coord,xs,ys, pixelsize, force_shift,out_put_shape,h,E,sigma,u_out,v_out):

    for i,(f_x, f_y, x, y) in enumerate(zip(fx, fy, fx_coord, fy_coord)):
        dist_x = xs - x
        dist_y = ys - y

        # greens tensor, As are central elements of the tensor
        # makeing sure coordinates are shifted and distance is never exactly zero ... (maybe just exlude zero????)
        for i in range(dist_x.shape[0]):
            for j in range(dist_x.shape[1]):
                if dist_x[i,j] == 0:
                    dist_x[i,j] += force_shift

        for i in range(dist_y.shape[0]):
            for j in range(dist_y.shape[1]):
                if dist_y[i,j] == 0:
                    dist_y[i,j] += force_shift

        dist_x = dist_x * pixelsize
        dist_y = dist_y * pixelsize

        r = np.sqrt(dist_x ** 2 + dist_y ** 2)

        ### elements of the greens tensor
        # Usefull parts // this is just an approximation
        # /also only for sigma=0.49
        s = r / h
        mu = E / (2 * (1 + sigma))
        A1 = ((2 - sigma) / (4 * np.pi * mu * h * s)) * (0.12 * np.exp(-0.43 * s) + 0.88 * np.exp(-0.83 * s))
        A2 = -(sigma / (4 * np.pi * mu * h * s)) * (1 + 1.22 * s + 1.31 * (s ** 2.23)) * np.exp(-1.25 * s)
        # A3_b = ((1-2*sigma)/(4*np.pi*mu*h*s))
        # A3 = -0.063 * ((np.exp(-0.44 * s) - np.exp(-2.79 * s)) ** 2)#*A3_b
        # A4 = ((2 * (1 - sigma)) / (4 * np.pi * mu * h * s)) * (1 + 0.46 * s - 2.50 * (s ** 2.13)) * np.exp(-2.18 * s)

        # components of the greens tensor
        K1 = A1 - ((dist_x ** 2 - dist_y ** 2) / (r ** 2)) * A2
        K2 = -(2 * dist_x * dist_y / (r ** 2)) * A2
        # K3 = -(dist_x / r) * A3
        K4 = -(2 * dist_x * dist_y / (r ** 2)) * A2  # is K2
        K5 = A1 + ((dist_x ** 2 - dist_y ** 2) / (r ** 2)) * A2

        def_u = f_x * K1 + f_y * K2
        def_v = f_x * K4 + f_y * K5
        u_out += def_u
        v_out += def_v
        # show_quiver(def_u,def_v)
        with objmode():
            if i % 100 == 0:
                print("###", i)

    return u_out, v_out




def finite_thickness_sparse_numba(fx, fy, fx_coord, fy_coord, pixelsize, h, young, out_put_shape, sigma=0.5, force_shift=0.5):
    '''
    convolution with the greens tensor for finite thikness calculation
    :param fx: forces in x direction
    :param fy: forces iny direction
    :param fx_coords: coordinates of these forces shifted by 0.5 pixel. shift is necessayr because we assume point load
    :param fy_coords:
    :param young: youngs modulus of material
    :param sigma: poisson ratio of material
    :return:
    '''
    u_out = np.zeros(out_put_shape)
    v_out = np.zeros(out_put_shape)
    ys, xs = np.indices(out_put_shape)

    fx = fx.astype(np.float)
    fy = fy.astype(np.float)
    E = young

    u_out,v_out = numba_loop(fx, fy, fx_coord, fy_coord, xs.astype(float), ys.astype(float), pixelsize, force_shift, out_put_shape, h, E, sigma,u_out,v_out)
    return u_out, v_out



def finite_thickness_sparse(fx, fy, fx_coord, fy_coord, pixelsize, h, young, out_put_shape, sigma=0.5, force_shift=0.5, use_exact=False):
    '''
    convolution with the greens tensor for finite thikness calculation
    :param fx: forces in x direction
    :param fy: forces iny direction
    :param fx_coords: coordinates of these forces shifted by 0.5 pixel. shift is necessayr because we assume point load
    :param fy_coords:
    :param young: youngs modulus of material
    :param sigma: poisson ratio of material
    :return:
    '''
    u_out = np.zeros(out_put_shape)
    v_out = np.zeros(out_put_shape)
    ys, xs = np.indices(out_put_shape)

    fx = fx.astype(np.float128)
    fy = fy.astype(np.float128)

    for f_x, f_y, x, y in tqdm(zip(fx, fy, fx_coord, fy_coord), total=len(fx)):
        dist_x = xs - x
        dist_y = ys - y


        # unfortunately the exact variant takes waaaaaaaay to long
        if use_exact:
            [K1, K2, K4, K5] = finite_thickenss_convolution_exact_greens_tensor_sparse(dist_x, dist_y, pixelsize, h, young,
                                                                                 sigma=sigma, force_shift=force_shift)
        else:
            [K1, K2, K4, K5] = finite_thickness_convolution_greens_tensor_sparse(dist_x, dist_y, pixelsize, h, young,
                                                                                 sigma=sigma, force_shift=force_shift)

        def_u = f_x * K1 + f_y * K2
        def_v = f_x * K4 + f_y * K5
        u_out += def_u
        v_out += def_v
        #show_quiver(def_u,def_v)

    return u_out, v_out



def finite_thickness_convolution(fx, fy,pixelsize, h, young, sigma=0.5, kernel_size=None, force_shift=0, method="fftconvolve"):
    '''
    convolution with the greens tensor for finite thikness calculation
    :param fx: forces in x direction
    :param fy: forces iny direction
    :param fx_coords: coordinates of these forces shifted by 0.5 pixel. shift is necessary because we assume point load
    :param fy_coords:
    :param young: youngs modulus of material
    :param sigma: poisson ratio of material
    :return:
    '''
    # returns displacement also in z direction (makes sense because of "pulling" at the buttom of the layer)
    fx = fx.astype(np.float128)
    fy = fy.astype(np.float128)

    # greens tensor, As are central elements of the tensor
    ([K1, K2, K3, K4, K5, K6, K7, K8, K9], [A1, A2, A3, A4]) = finite_thickness_convolution_greens_tensor(fx,
                                                            fy,pixelsize, h, young, sigma=sigma,kernel_size=kernel_size,force_shift=force_shift)
    #deformation by convolution with the kernels
    #v1,v2,v3,v4 = execute_as_thread([convolve2d]*4,[(fx, K1),(fy, K2),(fx, K4),(fx, K4),fy, K5],[{"mode":"same", "boundary":"fill", "fillvalue":0}]*4)

    #def_x=v1+v2
    #def_y=v3+v2
    if method=="convolve2d":
        def_x = convolve2d(fx, K1, mode="same", boundary="fill", fillvalue=0) + convolve2d(fy, K2, mode="same",
                                                                                           boundary="fill", fillvalue=0)
        def_y = convolve2d(fx, K4, mode="same", boundary="fill", fillvalue=0) + convolve2d(fy, K5, mode="same",
                                                                                           boundary="fill", fillvalue=0)
    if method=="fftconvolve": # this is always nearly faster and better??
        def_x = fftconvolve(fx, K1, mode="same") + fftconvolve(fy, K2, mode="same")
        def_y = fftconvolve(fx, K4, mode="same") + fftconvolve(fy, K5, mode="same")

    #def_z = convolv
    #def_z = convolve2d(fx, K7, mode="same", boundary="fill", fillvalue=0) + convolve2d(fy, K8, mode="same", boundary="fill", fillvalue=0)

    return def_x,def_y,#def_z

def finite_thickness_convolution_exact(fx, fy,pixelsize, h, young, sigma=0.5,kernel_size=None):
    '''
       calculating the deformation fieldgiven a tracktion force field using the boussinesq solution for an
        half sphere with finite thickness
        # employs the approximate greens tensor from

        Cell Force Microscopy on Elastic Layers of Finite Thickness
        Rudolf Merkel, Norbert Kirchgeßner, Claudia M. Cesa, and Bernd Hoffman
        Biophysical Journal  Volume 93  November 2007

        This is the "exact solution

       :param fx: forces in x direction
       :param fy: forces iny direction
       :param fx_coords: coordinates of these forces shifted by 0.5 pixel. shift is necessayr because we assume point load
       :param fy_coords:
       :param young: youngs modulus of material
       :param sigma: poisson ratio of material
       :return:
       '''

    # distance for the convolution kernel

    # greens tensor with:
    # [K1,K2,K3]
    # [K4,K5,K6]
    # [K7,K8,K9]
    # returns displacement also in z direction (makes sense because of "pulling" at the buttom of the layer)


    fx = fx.astype(np.float128)
    fy = fy.astype(np.float128)

    ([K1, K2, K3, K4, K5, K6, K7, K8, K9], [A1, A2, A3, A4]) = finite_thickenss_convolution_exact_greens_tensor(fx, fy,pixelsize, h, young, sigma=0.5,kernel_size=None)

    #deformation by convolution with the kernels
    def_x = convolve2d(fx, K1, mode="same", boundary="fill", fillvalue=0) + convolve2d(fy, K2, mode="same", boundary="fill", fillvalue=0)
    def_y = convolve2d(fx, K4, mode="same", boundary="fill", fillvalue=0) + convolve2d(fy, K5, mode="same", boundary="fill", fillvalue=0)
    #def_z = convolve2d(fx, K7, mode="same", boundary="fill", fillvalue=0) + convolve2d(fy, K8, mode="same", boundary="fill", fillvalue=0)


    return def_x,def_y#,def_z



def finite_thickenss_convolution_only(fx,fy,greens_tensor):
    '''
    uses greens_tensor only
    :param greens_tensor: list of tensor components
    :return:
    '''
    K1, K2, K3, K4, K5, K6, K7, K8, K9=greens_tensor
    #deformation by convolution with the kernels
    def_x = convolve2d(fx, K1, mode="same", boundary="fill", fillvalue=0) + convolve2d(fy, K2, mode="same", boundary="fill", fillvalue=0)
    def_y = convolve2d(fx, K4, mode="same", boundary="fill", fillvalue=0) + convolve2d(fy, K5, mode="same", boundary="fill", fillvalue=0)
    def_z = convolve2d(fx, K7, mode="same", boundary="fill", fillvalue=0) + convolve2d(fy, K8, mode="same", boundary="fill", fillvalue=0)
    return def_x,def_y,def_z

def deformation_in_fourrier_space():
    pass


def deformation_by_upsampling(fx, fy, factor, pixelsize=1, sigma=0.5, young=1, h=100, kernel_size=(30,30), method="convolve2d", return_upsampled=False):
   # org_size=fx.shape[0] # only squared shapes
    fxn, fyn = np.zeros((fx.shape[0] * factor, fx.shape[1] * factor)), np.zeros((fx.shape[0]  * factor, fx.shape[1]  * factor))
    # represent g forces by stretches
    fx_pos = np.array(np.where(fx != 0)).T
    fx_stretch = np.array([fx_pos[:, 0] - 0.5, fx_pos[:, 0] + 0.5, fx_pos[:, 1] - 0.5, fx_pos[:, 1] + 0.5]).T
    fy_pos = np.array(np.where(fy != 0)).T
    fy_stretch = np.array([fy_pos[:, 0] - 0.5, fy_pos[:, 0] + 0.5, fy_pos[:, 1] - 0.5, fy_pos[:, 1] + 0.5]).T
    # positions of forces in expanded array
    fxnp = np.array([np.meshgrid(np.linspace(sx1 * factor, sx2 * factor - 1, factor),
                                 np.linspace(sy1 * factor, sy2 * factor - 1, factor)) for sx1, sx2, sy1, sy2 in
                     fx_stretch]).astype(int)
    fynp = np.array([np.meshgrid(np.linspace(sx1 * factor, sx2 * factor - 1, factor),
                                 np.linspace(sy1 * factor, sy2 * factor - 1, factor)) for sx1, sx2, sy1, sy2 in
                     fy_stretch]).astype(int)
    for p, stretch in zip(fx_pos, fxnp):
        fxn[stretch[0], stretch[1]] = fx[p[0], p[1]] / (factor ** 2)
    for p, stretch in zip(fy_pos, fynp):
        fyn[stretch[0], stretch[1]] = fy[p[0], p[1]] / (factor ** 2)

    # calcualte deformation
    u, v = finite_thickness_convolution(fxn, fyn, pixelsize / factor, h, young, sigma=sigma,
                                        kernel_size=kernel_size, method=method)  # somwehat of an approximation
    all_pos = np.array(np.meshgrid(np.arange(fx.shape[0]), np.arange(fx.shape[1] ),indexing="ij")).astype(int) * factor
    ub, vb = u[all_pos[0], all_pos[1]], v[all_pos[0], all_pos[1]]
    if return_upsampled:
        return ub/pixelsize, vb/pixelsize, u/pixelsize, v/pixelsize, fxn, fyn
    else:
        return ub/pixelsize, vb/pixelsize



def deformation_by_upsampling_sparse(fx, fy, factor, pixelsize=1, sigma=0.5, young=1, h=100, use_numba=True,use_exact=False):

    # represent forces by stretches
    fyx_pos = np.array(np.where(np.logical_or(fx != 0, fy != 0))).T
    f_stretch = np.array([fyx_pos[:, 0] - 0.5, fyx_pos[:, 0] + 0.5, fyx_pos[:, 1] - 0.5, fyx_pos[:, 1] + 0.5]).T

    # positions of forces in expanded array
    f_np = np.array([np.meshgrid(np.linspace(sx1, sx2, factor),
                                 np.linspace(sy1, sy2, factor)) for sy1, sy2, sx1, sx2 in f_stretch])

    # value of the new forces
    fx_new = np.array([np.ones(f_np[i][0].shape) * fx[f_pos[0],f_pos[1]]/(factor**2) for i,f_pos in enumerate(fyx_pos)])
    fy_new = np.array([np.ones(f_np[i][0].shape) * fy[f_pos[0],f_pos[1]]/(factor**2) for i,f_pos in enumerate(fyx_pos)])

    # flattening to lsit
    fx_new = fx_new.flatten()
    fy_new = fy_new.flatten()
    x_pos = f_np[:,0,:,:].flatten()
    y_pos = f_np[:,1,:,:].flatten()

    # calculaltion of deformation field
    t1 = time.time()
    if use_numba:
        t1=time.time()
        u, v = finite_thickness_sparse_numba(fx_new, fy_new, x_pos, y_pos, pixelsize, h, young, out_put_shape=fx.shape,
                                       sigma=sigma, force_shift=0.5 / factor)
    else:
        u, v = finite_thickness_sparse(fx_new, fy_new, x_pos, y_pos, pixelsize, h, young, out_put_shape=fx.shape,
                                       sigma=sigma, force_shift=0.5 / factor, use_exact=use_exact)
    t2 = time.time()
    print("#### calculation time",t2-t1)
    # return deformation in pixel
    return u/pixelsize, v/pixelsize




# deformation by fourrier transform
def fourrier_deformation_finite_thickness(tx, ty, pixelsize, h, young, sigma=0.49):
    '''

    :param tx: x tractions [N/m**2]
    :param ty: y tractions [N/m**2]
    :param pixelsize: size of a pixel in m
    :param h: substrate height in m
    :param young: Youngs Modulous in N/m**2
    :param sigma: Poisson's ration (no unit)
    :return: u, v in pixel --> because this is the typical input for most other functions
    '''

    ## bens algortithm:
    # 1)Zero padding to get sqauerd array with even index number
    ax1_length, ax2_length = tx.shape  # u and v must have same dimensions
    max_ind = np.max(tx.shape )
    if max_ind % 2 != 0:
        max_ind += 1
    tx_expand = np.zeros((max_ind, max_ind))
    ty_expand = np.zeros((max_ind, max_ind))
    tx_expand[max_ind - ax1_length:max_ind, max_ind - ax2_length:max_ind] = tx
    ty_expand[max_ind - ax1_length:max_ind, max_ind - ax2_length:max_ind] = ty

    # 2) producing wave vectors (FT-space)
    kx1 = np.array([list(range(0, int(max_ind / 2), 1)), ] * int(max_ind))
    kx2 = np.array([list(range(-int(max_ind / 2), 0, 1)), ] * int(max_ind))
    # fourier transform in this case is defined as
    # F(kx)=1/2pi integral(exp(i*kx*x)dk therefore kx must be expressed as a spatial frequency in distance*2*pi
    kx = np.append(kx1, kx2, axis=1) * 2 * np.pi / (pixelsize * max_ind)
    ky = np.transpose(kx)
    k = np.sqrt(kx ** 2 + ky ** 2)

    # some constants
    r = k * h
    c = np.cosh(r)
    s = np.sinh(r)
    s_c = np.tanh(r) # this formulation is robust (returns correct limit of 1) when k * h gets large
   # gamma = ((3 - 4 * sigma) * (c ** 2) + (1 - 2 * sigma) ** 2 + (r) ** 2) / (
   #         (3 - 4 * sigma) * s * c + r)  ## inf values here because k goes to zero

    # TODO use this formulation in FTTC algorithm
    gamma = ((3 - 4 * sigma) + (((1 - 2 * sigma) ** 2) / (c ** 2)) + ((r ** 2) / (c ** 2)))/ (
            (3 - 4 * sigma) * s_c + r / (c ** 2))

    A = young /(2 * s_c * (1 + sigma))
    B = young * gamma/ (2 * (1 - sigma**2))

    F1 = (ky**2 * A + kx**2 * B )
    F2 = (B-A) * ky * kx
    F3 = (kx**2 * A + ky**2 * B )

    F11 = (-1)*k*F3 / (F2**2 - F1*F3)
    F12 = k*F2/(F2**2 - F1*F3)
    # F12 and F21 are identical
    F22 = (-1)*k*F1 / (F2**2 - F1*F3)

    # 4) calculate fourier transform of displacements
    tx_ft = np.fft.fft2(tx_expand)
    ty_ft = np.fft.fft2(ty_expand)

    u_ft = F11 * tx_ft + F12 *  ty_ft
    v_ft = F12 * tx_ft + F22 *  ty_ft

    # zero frequency would represents deformations everywhere?
    u_ft[0, 0] = 0
    v_ft[0, 0] = 0

    # 4.2) go back to real space
    u = np.fft.ifft2(u_ft.astype(np.complex128)).real
    v = np.fft.ifft2(v_ft.astype(np.complex128)).real

    # 5.2) cut like in script from ben
    u_cut = u[max_ind - ax1_length:max_ind, max_ind - ax2_length:max_ind]
    v_cut = v[max_ind - ax1_length:max_ind, max_ind - ax2_length:max_ind]

   # u_b2, v_b2 = deformation_by_upsampling_sparse(tx * pixelsize**2, ty * pixelsize**2, factor=10, pixelsize=pixelsize, sigma=0.49,
  #                                                young=young,
   #                                               h=h, use_numba=False, use_exact=False)

   # show_quiver(u_cut, v_cut)
   # show_quiver(u_b2, v_b2)

    return u_cut , v_cut


def fourrier_deformation(fx, fy, pixelsize, young, sigma=0.49, h=100, use_finite_height=True):
    # return deformation in pixels

    #conversion of force to tractions
    tx = fx/(pixelsize**2)
    ty = fy/(pixelsize**2)
    if use_finite_height:
        u, v = fourrier_deformation_finite_thickness(tx, ty, pixelsize, h, young, sigma=sigma)
    else:
        u, v = fourrier_deformation_infinite_hight(tx, ty, pixelsize, young, sigma=sigma)
    #show_quiver(u, v)
    return u/pixelsize, v/pixelsize


def fourrier_deformation_infinite_hight(tx, ty, pixelsize, young, sigma=0.49):

    '''
    :param tx: x tractions [N/m**2]
    :param ty: y tractions [N/m**2]
    :param pixelsize: size of a pixel in m
    :param young: Youngs Modulous in N/m**2
    :param sigma: Poisson's ration (no unit)
    :return: u, v in pixel --> because this is the typical input for most other functions
    '''
    #show_quiver(fx, fy)
    # 1)Zero padding to get sqauerd array with even index number
    #

    ax1_length, ax2_length = tx.shape  # u and v must have same dimensions
    max_ind = np.max(tx.shape)
   # max_ind +=2
    if max_ind % 2 != 0:
        max_ind += 1

    tx_expand = np.zeros((max_ind, max_ind))
    ty_expand = np.zeros((max_ind, max_ind))
    tx_expand[max_ind - ax1_length:max_ind, max_ind - ax2_length:max_ind] = tx
    ty_expand[max_ind - ax1_length:max_ind, max_ind - ax2_length:max_ind] = ty

    #fx_expand = np.pad(fx,((int((max_ind-ax1_length)/2),int((max_ind-ax1_length)/2)), (int((max_ind-ax2_length)/2),int((max_ind-ax2_length)/2))))
    #fy_expand = np.pad(fy, ((int((max_ind - ax1_length) / 2), int((max_ind - ax1_length) / 2)),
     #                       (int((max_ind - ax2_length) / 2), int((max_ind - ax2_length) / 2))))

    # 2) produceing wave vectors   ## why this form
    # form 1:max_ind/2 then -(max_ind/2:1)
    kx1 = np.array([list(range(0, int(max_ind / 2), 1)), ] * int(max_ind))
    kx2 = np.array([list(range(-int(max_ind / 2), 0, 1)), ] * int(max_ind))
    # fourier transform in this case is defined as

    kx = np.append(kx1, kx2, axis=1) * 2 * np.pi / (pixelsize * max_ind)
    # F(kx)=1/2pi integral(exp(i*kx*x)dk therefore kx must be expressed as a spatial frequency in distance*2*pi

    ky = np.transpose(kx)
    k = np.sqrt(kx ** 2 + ky ** 2)

    A = (2 + 2 * sigma) / (young * k ** 3)

    f11 = (1 - sigma) * k**2 + sigma * ky**2
    f12 = -sigma * kx * ky  # note the sign error in butler et al 2002 here
    f22 = (1 - sigma) * k**2 + sigma * kx**2

    # not sure if I need this
    f12[:, int(max_ind / 2)] = np.zeros(max_ind)
    f12[int(max_ind / 2), :] = np.zeros(max_ind)

    # 4) calculate fourier transform of displacements
    tx_ft = np.fft.fft2(tx_expand)
    ty_ft = np.fft.fft2(ty_expand)

    u_ft = A * (f11 * tx_ft + f12 * ty_ft)
    v_ft = A * (f12 * tx_ft + f22 * ty_ft)

    # zero frequency would represnt deformations everywhere?
    u_ft[0, 0] = 0
    v_ft[0, 0] = 0

    # 4.2) go back to real space
    u = np.fft.ifft2(u_ft.astype(np.complex128)).real
    v = np.fft.ifft2(v_ft.astype(np.complex128)).real

    # 5.2) cut like in script from ben
    u_cut = u[int((max_ind-ax1_length)/2):max_ind-int((max_ind-ax1_length)/2), int((max_ind-ax2_length)/2):max_ind-int((max_ind-ax2_length)/2)]
    v_cut = v[int((max_ind-ax1_length)/2):max_ind-int((max_ind-ax1_length)/2), int((max_ind-ax2_length)/2):max_ind-int((max_ind-ax2_length)/2)]

    # 0) substracting mean(better name for this step)
    u_shift = (u_cut - np.mean(u_cut))  # shifting to zero mean  (translating to pixelsize of u-image is done later)
    v_shift = (v_cut - np.mean(v_cut))

    return u_cut, v_cut

