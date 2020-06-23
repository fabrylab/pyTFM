from scipy.ndimage import binary_erosion
import sys
sys.path.insert(0,'/home/user/Software/pyTFM/analysis_and_testing/')
from simulating_deformation import *
from playing_with_strains import *
from plotting_evaluation import *
from general_evaluation import *
from pyTFM.TFM_functions import strain_energy_points
import time
import numpy as np
from skimage.filters import gaussian


def deformation_by_upsampling1D(fx, factor, sigma=0.5,E=1):
    ##### this not a good evaluation becuase basic equations are based on 2dimensional boussinesq solution?? #####
    # bousssinesq solution or deformation along one axis, if forces are only found along this axis,
    # this is manly for investigating the "upsampling method"
   # org_size=fx.shape[0] # only squared shapes
    fxn = np.zeros((fx.shape[0] * factor))
    # represent g forces by stretches
    fx_pos = np.array(np.where(fx != 0))[0]
    fx_stretch = np.array([fx_pos - 0.5, fx_pos + 0.5]).T
    fxnp = np.array([np.linspace(sx1 * factor, sx2 * factor - 1, factor) for sx1, sx2 in fx_stretch]).astype(int)
    for p, stretch in zip(fx_pos, fxnp):
        fxn[stretch] = fx[p] / (factor ** 2)
    mid=fxn.shape[0]/2
    if mid==int(mid):
       mid+=0.5
    A=(1+sigma)/np.pi*E
    kernel = np.abs(A/(np.arange(0,fxn.shape[0])-mid))# distances
    u=np.convolve(fxn, kernel, mode="same")

    all_pos = np.array(np.arange(fx.shape[0])).astype(int) * factor
    ub = u[all_pos]
    return ub


def slicing_plot(u,f_b,f_bf,f_nf,f_f):
    line = np.array([np.arange(30, 50), np.arange(30, 50)])
    fxbl = f_b[line[0], line[1]]
    fxbfl =  f_bf[line[0], line[1]]
    ul = u[line[0], line[1]]
    fxfl_nf = f_nf[line[0], line[1]]
    fxfl = f_f[line[0], line[1]]
    plt.figure()
    plt.plot(fxbl, "o", label="f input")
    plt.plot(fxbfl, "o", label="f input filtered")
    plt.plot(ul, "o", label="def")
    plt.plot(fxfl_nf, "o", label="f output un filtered")
    plt.plot(fxfl, "o", label="f output filtered")
    plt.legend()

def compare_strain_energies(u, v, fx_b, fy_b, fx_bf, fy_bf, fx_nf, fy_nf, fx, fy):
    mask = binary_erosion(np.ones((100, 100)), iterations=10).astype(bool)
    energy_points_f = strain_energy_points(u, v, fx_b, fy_b, 1, 1)  # contractile energy at any point
    strain_b = np.sum(energy_points_f[mask])
    energy_points_f = strain_energy_points(u, v, fx_bf, fy_bf, 1, 1)  # contractile energy at any point
    strain_bf = np.sum(energy_points_f[mask])
    energy_points_f = strain_energy_points(u, v, fx_nf, fy_nf, 1, 1)  # contractile energy at any point
    strain_fnf = np.sum(energy_points_f[mask])
    energy_points_f = strain_energy_points(u, v, fx, fy, 1, 1)  # contractile energy at any point
    strain_f = np.sum(energy_points_f[mask])
    print("strain energy\n", "input: ", strain_b, "\ninput filtered: ", strain_bf, "\n output unfiltered: ", strain_fnf, "\noutput filtered: ", strain_f)
    return strain_b,strain_bf, strain_fnf, strain_f

def compare_contracillities(fx_b, fy_b, fx_bf, fy_bf, fx_nf, fy_nf, fx, fy):
    mask = binary_erosion(np.ones((100, 100)), iterations=10).astype(bool)
    contractile_force_b, *others = contractillity(fx_b, fy_b, 1, mask)
    contractile_force_bf, *others = contractillity(fx_bf, fy_bf, 1, mask)
    contractile_force_fnf, *others = contractillity(fx_nf, fy_nf, 1, mask)
    contractile_force_f, *others = contractillity(fx, fy, 1, mask)

    return contractile_force_b, contractile_force_bf, contractile_force_fnf, contractile_force_f

pixelsize=1
h=1000
young=1
mask=np.ones((100,100))


fx_b, fy_b = np.zeros((100,100)), np.zeros((100,100))
l1 = np.array([np.arange(35,45),np.arange(35,45)])
l2 = np.array([np.arange(55,65),np.arange(55,65)])
#l1 = np.array([40,40])
#l2 = np.array([60,60])
fx_b[l1[0],l1[1]]=-1
fy_b[l1[0],l1[1]]=-1
fx_b[l2[0],l2[1]]=1
fy_b[l2[0],l2[1]]=1
#fx_b = gaussian(fx_b,sigma=4)
#fy_b = gaussian(fy_b,sigma=4)

fx_bf = gaussian(fx_b,sigma=1)
fy_bf = gaussian(fy_b,sigma=1)
#show_quiver(fx_b, fy_b)

u_b, v_b = finite_thickness_convolution(fx_b, fy_b, pixelsize, h, young, sigma=0.5,
                                        kernel_size=(20,20))

fx_f_nf, fy_f_nf = traction_wrapper(u_b, v_b, pixelsize, h, young, mask=mask,sigma=0.5,
                              filter=None, fs=None)  # assuming pixelsize == 1
#show_quiver(fx_f_nf, fy_f_nf)
fx_f, fy_f = traction_wrapper(u_b, v_b, pixelsize, h, young, mask=mask,sigma=0.5,
                              filter="gaussian", fs=1)  # assuming pixelsize == 1

print("old method")
strain_b,strain_bf, strain_fnf, strain_f = compare_strain_energies(u_b, v_b, fx_b,fy_b,fx_bf,fy_bf,fx_f_nf,fy_f_nf,fx_f,fy_f)
contractile_force_b,contractile_force_bf, contractile_force_fnf, contractile_force_f = compare_contracillities(fx_b, fy_b,fx_bf,fy_bf, fx_f_nf,
                                                                                             fy_f_nf, fx_f, fy_f)

s_list = []
c_list = []
for i in range(1,15,5):# upsampling method for deformation field
    factor = i
    print(i)
    u_b1, v_b1 = deformation_by_upsampling(fx_b, fy_b, factor, kernel_size=(10*factor, 10*factor))

    # select only specific deformations

    # downsample

    #
    fx_f_nf1, fy_f_nf1 = traction_wrapper(u_b1, v_b1, pixelsize, h, young, mask=mask, sigma=0.5,
                                          filter=None, fs=None)  # assuming pixelsize == 1
    # show_quiver(fx_f_nf1, fy_f_nf1)
    fx_f1, fy_f1 = traction_wrapper(u_b1, v_b1, pixelsize, h, young, mask=mask, sigma=0.5,
                                    filter="gaussian", fs=1)  # assuming pixelsize == 1

    # extracting values to 1d along a single line
    # slicing_plot(u_b,fx_f_nf,fx_f)
    # slicing_plot(u_b1,fx_f_nf1,fx_f1)

    print("new method")
    strain_b1, strain_bf1, strain_fnf1, strain_f1 = compare_strain_energies(u_b1, v_b1, fx_b, fy_b,fx_bf,fy_bf, fx_f_nf1, fy_f_nf1, fx_f1, fy_f1)
    #print(strain_b1, strain_fnf1)2
    s_list.append([strain_b1, strain_fnf1, strain_f1])
    contractile_force_b1,contractile_force_bf1, contractile_force_fnf1, contractile_force_f1 = compare_contracillities(fx_b, fy_b,fx_bf,fy_bf, fx_f_nf1, fy_f_nf1, fx_f1, fy_f1)
    c_list.append([contractile_force_b1, contractile_force_fnf1, contractile_force_f1])

s_list=np.array(s_list)
c_list=np.array(c_list)
plt.figure()
plt.plot(s_list[:,0])
plt.plot(s_list[:,1])
plt.plot(s_list[:,2])

plt.figure()
plt.plot(c_list[:,0])
plt.plot(c_list[:,1])
plt.plot(c_list[:,2])
plt.figure()
plt.bar(list(range(8)),[strain_b,strain_bf, strain_fnf, strain_f, strain_b1,strain_bf1, strain_fnf1, strain_f1])



plt.figure()
plt.bar(list(range(8)),[contractile_force_b, contractile_force_bf, contractile_force_fnf, contractile_force_f, contractile_force_b1, contractile_force_bf1, contractile_force_fnf1, contractile_force_f1 ])

slicing_plot(u_b, fx_b, fx_bf, fx_f_nf, fx_f)
slicing_plot(u_b1, fx_b,fx_bf, fx_f_nf1, fx_f1)

ep1 = strain_energy_points(u_b1, v_b1, fx_f_nf1, fy_f_nf1, 1, 1)
ep2 = strain_energy_points(u_b1, v_b1, fx_f1, fy_f1, 1, 1)
show_quiver(u_b1,v_b1)
show_quiver(fx_f_nf1,fy_f_nf1)
show_quiver(fx_f1,fy_f1)

plt.figure();plt.imshow(ep1);plt.colorbar()
plt.figure();plt.imshow(ep2);plt.colorbar()

### illustrating the method
plt.close("all")
fx = np.random.uniform(-1,1,(10,10))
fy = np.random.uniform(-1,1,(10,10))

fx = np.zeros((10,10))
fy = np.zeros((10,10))
fx[4,4] = 1
fy[4,4] = 1
us, vs, u, v, fxn, fyn  = deformation_by_upsampling(fx, fy, factor=5, pixelsize=1, sigma=0.5, young=1, h=100, kernel_size=(30,30), method="convolve2d", return_upsampled=True)


f, ax = show_quiver(fx,fy, cmap="coolwarm",scale_ratio=0.1, vmax=1.25, width= 0.006,headwidth=7, headaxislength=5, headlength=7)
ax.set_axis_on()
ax.set_yticks(np.arange(-0.5,fx.shape[0],1), minor="true")
ax.set_xticks(np.arange(-0.5,fx.shape[1],1), minor="true")
ax.grid(which="minor", linewidth=1, color="black")

f, ax = show_quiver(us,vs, cmap="coolwarm",scale_ratio=0.05, width= 0.006,headwidth=7, headaxislength=5, headlength=7)
ax.set_axis_on()
ax.set_yticks(np.arange(-0.5,us.shape[0],1), minor="true")
ax.set_xticks(np.arange(-0.5,us.shape[1],1), minor="true")
ax.grid(which="minor", linewidth=1, color="black")

f, ax = show_quiver(fxn, fyn, cmap="coolwarm",scale_ratio=0.02, vmax=1.25, width= 0.003,headwidth=7, headaxislength=5, headlength=7)
ax.set_axis_on()
ax.set_yticks(np.arange(0, u.shape[0],1))
ax.set_xticks(np.arange(0, u.shape[1],1))
ax.set_yticks(np.arange(-0.5,fxn.shape[0],1), minor="true")
ax.set_xticks(np.arange(-0.5,fxn.shape[1],1), minor="true")
ax.grid(which="major", linewidth=1, color="black")

f, ax = show_quiver(u,v, cmap="coolwarm", scale_ratio=0.03, width= 0.003,headwidth=7, headaxislength=5, headlength=7)
ax.set_axis_on()
ax.set_yticks(np.arange(-0.5,u.shape[0],1), minor="true")
ax.set_xticks(np.arange(-0.5,u.shape[1],1), minor="true")
ax.grid(which="minor", linewidth=1, color="black")