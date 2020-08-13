from scipy.ndimage import binary_erosion
import sys
sys.path.insert(0,'/home/user/Software/pyTFM/analysis_and_testing/')
from simulating_deformation import *
from playing_with_strains import *
from plotting_evaluation import *
from general_evaluation3 import *
import time
import numpy as np
from skimage.filters import gaussian





pixelsize = 1
h = 1000
young = 1
mask=np.ones((100,100))
factor = 5



f = np.zeros(50)
pos = np.arange(50)
f_pos = 25
f[f_pos] = 1

dist_vecs = np.arange(50) - f_pos
dist = np.abs(dist_vecs)
d=1/dist
d[dist==0]=sorted(d.flatten())[-2]*1.2


fx, fy = np.zeros((100,100)), np.zeros((100,100))
#l1 = np.array([np.arange(35,45),np.arange(35,45)])
#l2 = np.array([np.arange(55,65),np.arange(55,65)])
l1 = np.array([40,40])
l2 = np.array([60,60])
fx[l1[0],l1[1]]=-1
fy[l1[0],l1[1]]=-1
fx[l2[0],l2[1]]=1
fy[l2[0],l2[1]]=1
#fx_b = gaussian(fx_b,sigma=4)
#fy_b = gaussian(fy_b,sigma=4)



u, v = deformation_by_upsampling(fx, fy, factor=10, kernel_size=(20*factor, 20*factor))


#fx_filtered = gaussian(fx,sigma=1)
#fy_filtered = gaussian(fy,sigma=1)
#show_quiver(fx_b, fy_b)


f_filtered = gaussian(f, sigma=1)

slice = np.array([np.arange(55, 65),np.arange(55, 65)])
u_slice=u[slice[0],slice[1]]
f_slice=fx[slice[0],slice[1]]
f_slice_filtered=gaussian(f_slice,1)

plt.figure();
plt.bar(slice[0], u_slice, color="C0", width=1, label="deformtaions*")
plt.bar(slice[0], f_slice, color="C1" , width=1, label="force")
plt.bar(slice[0], f_slice*u_slice, color="C2", width=0.5, label="def * force")
plt.ylim((0,1.5))
plt.legend(loc="upper right")
plt.savefig("/home/user/Desktop/backup_from_harddrive/data_traction_force_microscopy/why_strain_energy_is_influenced_by_gaussian1.png")
plt.figure();
plt.bar(slice[0], u_slice, color="C0", width=1, label="deformtaions*")
plt.bar(slice[0], f_slice_filtered, color="C1", width=1, label="gauss-filterd force")
plt.bar(slice[0], f_slice_filtered*u_slice, color="C2", width=0.5, label="def * force")
plt.ylim((0,1.5))
plt.legend(loc="upper right")
plt.savefig("/home/user/Desktop/backup_from_harddrive/data_traction_force_microscopy/why_strain_energy_is_influenced_by_gaussian2.png")





# effekt of gausfilter on functions:

f = np.zeros(50)
pos = np.arange(50)
f_pos = 25
f[f_pos] = 1
dist_vecs = np.arange(50) - 25
dist = np.abs(dist_vecs)
d=1/dist
d[dist==0]=sorted(d.flatten())[-2]
#d1[dist==0]= 2* np.log(0.5)*1
plt.figure()
plt.plot(pos,d)
print(np.sum(f*d))
for i in range(8):
    plt.plot(pos, gaussian(f, sigma=i) * d)
    print(np.sum(gaussian(f, sigma=i) * d))
plt.plot(pos,f,"o")

from scipy.signal import medfilt
from scipy.ndimage import uniform_filter
plt.plot(pos, medfilt(f,kernel_size=3))
print(np.sum(medfilt(f,kernel_size=3) * d))

plt.plot(pos, uniform_filter(f, size=10))
print(np.sum(uniform_filter(f, size=10) * d))

ps=np.linspace(-5,5,100)
x=np.abs(np.linspace(-5,5,100))
y1=1/x
y1[np.argmin(x)]=sorted(y1.flatten())[-2]
y2=-np.abs(x)+5
y3=y1*y2
plt.figure();
plt.plot(ps,y1)
plt.plot(ps,y2)
plt.plot(ps,y3)

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