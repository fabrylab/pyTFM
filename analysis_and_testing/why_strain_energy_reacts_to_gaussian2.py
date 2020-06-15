from scipy.ndimage import binary_erosion
import sys
sys.path.insert(0,'/home/user/Software/pyTFM/analysis_and_testing/')
from simulating_deformation import *
from playing_with_strains import *
from plotting_evaluation import *
from general_evaluation import *
from pyTFM.TFM_functions import strain_energy_points
import time
import matplotlib.pyplot as plt
import numpy as np
from skimage.filters import gaussian


pixelsize=1
h=1000
young=1
mask=np.ones((100,100))
mask = binary_erosion(np.ones((100, 100)), iterations=10).astype(bool)

fx_b, fy_b = np.zeros((100,100)), np.zeros((100,100))
#l1 = np.array([np.arange(35,45),np.arange(35,45)])
#l2 = np.array([np.arange(55,65),np.arange(55,65)])
l1 = np.array([40,40])
l2 = np.array([60,60])
fx_b[l1[0],l1[1]]=-1
fy_b[l1[0],l1[1]]=-1
fx_b[l2[0],l2[1]]=1
fy_b[l2[0],l2[1]]=1

u_b, v_b = deformation_by_upsampling(fx_b, fy_b, factor=5, kernel_size=(50, 50))


fx_f_nf, fy_f_nf = traction_wrapper(u_b, v_b, pixelsize, h, young, mask=mask,sigma=0.5,
                              filter=None, fs=None)  # assuming pixelsize == 1
#show_quiver(fx_f_nf, fy_f_nf)
fx_f, fy_f = traction_wrapper(u_b, v_b, pixelsize, h, young, mask=mask,sigma=0.5,
                              filter="gaussian", fs=1)  # assuming pixelsize == 1


energy_points = strain_energy_points(u_b,v_b, fx_b, fy_b, 1, 1)  # contractile energy at any point
cont1 = np.sum(energy_points[mask])
energy_points = strain_energy_points(u_b, v_b, fx_f_nf, fy_f_nf, 1, 1)  # contractile energy at any point
cont2 = np.sum(energy_points[mask])
energy_points = strain_energy_points(u_b, v_b, fx_f, fy_f, 1, 1)  # contractile energy at any point
cont3 = np.sum(energy_points[mask])

f,a = show_quiver(fx_b,fy_b,cmap="coolwarm")
f.savefig("/home/user/Desktop/backup_from_harddrive/data_traction_force_microscopy/why_strain_energy_is_influenced_by_gaussian3.png")
f,a=show_quiver(fx_f_nf, fy_f_nf,cmap="coolwarm")
f.savefig("/home/user/Desktop/backup_from_harddrive/data_traction_force_microscopy/why_strain_energy_is_influenced_by_gaussian4.png")
f,a=show_quiver(fx_f, fy_f,cmap="coolwarm")
f.savefig("/home/user/Desktop/backup_from_harddrive/data_traction_force_microscopy/why_strain_energy_is_influenced_by_gaussian5.png")
f=plt.figure(figsize=(2,5))
plt.bar([1,2,3],[cont1,cont2,cont3])
plt.xticks([1,2,3],["input", "output unfiltered", "output filtered"])
plt.xticks(rotation=70)
plt.tight_layout()
f.savefig("/home/user/Desktop/backup_from_harddrive/data_traction_force_microscopy/why_strain_energy_is_influenced_by_gaussian6.png")