
from pyTFM.data_analysis import *
import numpy as np
import os


def introduce_cell_number(res_dict, cell_number, max_area):
    res_dict_cp = res_dict.copy()
    for frame, value in res_dict_cp.items():
        cn = cell_number[frame]
        cn_local = [[obj, (area / max_area) * cn] for (obj, area) in value["area_force_measurement"]]
        res_dict[frame]["cell number"] = cn_local


# path to the out.txt text file
file1="/home/user/Desktop/plate1/relaxed/out8.txt"
parameter_dict1,res_dict1 = read_output_file(file1) # reading the file and splitting into parameters and results



file2="/home/user/Desktop/plate2/beads_yuko_forced/png/out3.txt"
parameter_dict2, res_dict2 = read_output_file(file2)# reading the fie and splitting into parameters and results

file3="/home/user/Desktop/plate3/relaxed/out1.txt"
parameter_dict3,res_dict3 = read_output_file(file3)# reading the fie and splitting into parameters and results




max_area = 342000

rename_dict = {"x-1_y-1":"0","x-1_y0":"1","x0_y-1":"2","x0_y0":"3","x0_y1":"4","x1_y0":"5","x1_y1":"6","x-1_y1":"7","x1_y-1":"8"}
ids = ["x0_y0", "x0_y1","x0_y-1", "x1_y0","x-1_y0","x1_y1","x1_y-1","x-1_y1","x-1_y-1"]


cell_number1 = [51, 49, 44, 42, 48, 52, 40, 57, 51]
cell_number2 = [35, 38, 32, 32, 42, 24, 32, 42, 31]
cell_number3 = [31 ,31, 25 ,36, 34, 36, 26, 31, 27]

cell_number1 = {rename_dict[id]:cn for id, cn in zip(ids, cell_number1) }
del cell_number1["0"]
cell_number2 = {rename_dict[id]:cn for id, cn in zip(ids, cell_number2) }
cell_number3 = {rename_dict[id]:cn for id, cn in zip(ids, cell_number3) }

introduce_cell_number(res_dict1, cell_number1, max_area)
introduce_cell_number(res_dict2, cell_number2, max_area)
introduce_cell_number(res_dict3, cell_number3, max_area)

n_frames1, values_dict1, frame_list1 = prepare_values(res_dict1,[])
n_frames2, values_dict2, frame_list2 = prepare_values(res_dict2,[])
n_frames3, values_dict3, frame_list3 = prepare_values(res_dict3,[])



norm1 = values_dict1["strain energycolony"] /values_dict1["cell number"]
norm2 = values_dict2["strain energycolony"] /values_dict2["cell number"]
norm3 = values_dict3["strain energycolony"] /values_dict3["cell number"]
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
plt.ylabel("strain energy [J/cell]")
plt.xlabel("substrate stiffness")
plt.ylim((0,np.max([strain_energy_density1 + er1,strain_energy_density2+ er2,strain_energy_density3+ er3, ])*1.2))




plate1 = 2000
plate2 = 64000
plate3 = 16000