
from pyTFM.data_analysis import *
import numpy as np
import os



# path to the out.txt text file
file1="/home/user/Desktop/plate1/relaxed/out8.txt"
parameter_dict1,res_dict1 = read_output_file(file1) # reading the file and splitting into parameters and results
n_frames1, values_dict1, frame_list1 =prepare_values(res_dict1,[])


file2="/home/user/Desktop/plate2/beads_yuko_forced/png/out3.txt"
parameter_dict2, res_dict2 = read_output_file(file2)# reading the fie and splitting into parameters and results
n_frames2, values_dict2, frame_list2 = prepare_values(res_dict2,[])

file3="/home/user/Desktop/plate3/relaxed/out1.txt"
parameter_dict3,res_dict3 = read_output_file(file3)# reading the fie and splitting into parameters and results
n_frames3, values_dict3, frame_list3 = prepare_values(res_dict3,[])





norm1 = values_dict1["strain energycolony"] /values_dict1["area_force_measurement"]
norm2 = values_dict2["strain energycolony"] /values_dict2["area_force_measurement"]
norm3 = values_dict3["strain energycolony"] /values_dict3["area_force_measurement"]
scipy_ttest_ind(norm1, norm2)

strain_energy_density1 =  np.mean(norm1)
strain_energy_density2 =  np.mean(norm2)
strain_energy_density3 =  np.mean(norm3)



er1 = np.std(norm1)/(len(norm1)**0.5)
er2 = np.std(norm2)/(len(norm2)**0.5)
er3 = np.std(norm3)/(len(norm3)**0.5)

plt.figure()
plt.bar([0,1,2],[strain_energy_density1, strain_energy_density2, strain_energy_density3], yerr=np.array([er1,er2,er3]) )
plt.xticks([0,1,2],["2000 kPa", "64000 kPa", "16000 kPa"])
plt.ylabel("strain energy density [J/µm]")
plt.xlabel("substrate stiffness")
plt.ylim((0,np.max([strain_energy_density1 + er1,strain_energy_density2+ er2,strain_energy_density3+ er3, ])*1.2))




plate1 = 2000
plate2 = 64000
plate3 = 16000