from scipy.ndimage import binary_erosion
import sys
sys.path.insert(0,'/home/user/Software/pyTFM/analysis_and_testing/')
from simulating_deformation import *
from playing_with_strains import *
from plotting_evaluation import *
from general_evaluation import *
import time
import numpy as np
from skimage.filters import gaussian

def deformation_by_upsampling(fx, fy, factor, kernel_size=(30,30)):
    org_size=fx.shape[0] # only squared shapes
    fxn, fyn = np.zeros((org_size * factor, org_size * factor)), np.zeros((org_size * factor, org_size * factor))
    # represent g forces by stretches
    fx_pos = np.array(np.where(fx_b != 0)).T
    fx_stretch = np.array([fx_pos[:, 0] - 0.5, fx_pos[:, 0] + 0.5, fx_pos[:, 1] - 0.5, fx_pos[:, 1] + 0.5]).T
    fy_pos = np.array(np.where(fy_b != 0)).T
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
    u, v = finite_thickness_convolution(fxn, fyn, pixelsize / factor, h, young, sigma=0.5,
                                        kernel_size=kernel_size)  # somwhat of an approximation
    all_pos = np.array(np.meshgrid(np.arange(org_size), np.arange(org_size))).astype(int) * factor
    ub, vb = u[all_pos[0], all_pos[1]], v[all_pos[0], all_pos[1]]
    return ub, vb



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
h=1
young=1
mask=np.ones((100,100))


fx_b, fy_b = np.zeros((100,100)), np.zeros((100,100))
l1 = np.array([np.arange(35,45),np.arange(35,45)])
l2 = np.array([np.arange(55,65),np.arange(55,65)])
l1 = np.array([40,40])
l2 = np.array([60,60])
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








# effekt of gausfilter on functions:
f = np.zeros(50)
pos = np.arange(50)
f_pos = 25
f[f_pos] = 1
dist = np.abs(np.arange(50) - 25)
d = 1/dist # boussinesq solution for deormation in x direction if force also goes purely in x direction//point force
# uniform distributed load on one pixel:
from scipy.integrate import dblquad
def int_funct(y,x,sigma=0.5):
    return((1-sigma)/np.sqrt(x**2,y**2))+(sigma*x**2/np.sqrt(x**2+y**2)**3)


d1 = np.log(np.abs((dist+0.5)/(dist-0.5)))
d1[dist==0]= 2* np.log(0.5)*1
plt.figure()
plt.plot(pos,d)
plt.plot(pos,d1)
plt.plot(pos,f,"o")
## show that gaussfiltering alone is problematic:
'''
fi=np.zeros(30)
fi[15]=1


# testing smaller kernel size for convolution to get faster calculation

ts=[]
ubs=[]
vbs=[]
for i in range(2,110,2):
    t1=time.time()
    ub, vb = finite_thickness_convolution(fx_b, fy_b, pixelsize, h, young, sigma=0.5,
                                 kernel_size=(i,i))  # somwhat of an approximation
    t2 = time.time()
    ts.append(t2-t1)
    ubs.append(ub)
    vbs.append(vb)
corsu=[np.corrcoef(u.flatten(),ubs[-1].flatten())[0,1] for u in ubs]
corsv=[np.corrcoef(v.flatten(),vbs[-1].flatten())[0,1] for v in vbs]

show_quiver(ubs[0], vbs[0])
show_quiver(ubs[-1], vbs[-1])

'''