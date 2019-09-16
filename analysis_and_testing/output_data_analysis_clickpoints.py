# example anlaysis of the results

from andreas_TFM_package.TFM_functions_for_clickpoints import units
from andreas_TFM_package.utilities_TFM import *
from andreas_TFM_package.data_analysis import *
from collections import defaultdict
import numpy as np
import os

folder_plots="/media/user/GINA1-BK/data_traktion_force_microscopy/WT_vs_KO_images_10_09_2019/wt_vs_ko_images_Analyzed/results_100_95/"
createFolder(folder_plots)
exclude=[]
# first file
file1="/media/user/GINA1-BK/data_traktion_force_microscopy/WT_vs_KO_images_10_09_2019/wt_vs_ko_images_Analyzed/WTshift/out_100.txt"
parameter_dict1,res_dict1=read_output_file(file1)
n_frames1,values_dict1,mean_values1,std_values1,frame_list1=prepare_values(res_dict1,exclude)
print("WT",np.sum(values_dict1["contractile energy"]))
## n is not always the same!!
# second file
exclude=["01","10"]
file2="/media/user/GINA1-BK/data_traktion_force_microscopy/WT_vs_KO_images_10_09_2019/wt_vs_ko_images_Analyzed/KOshift/out_100.txt"
parameter_dict2,res_dict2=read_output_file(file2)
n_frames2,values_dict2,mean_values2,std_values2,frame_list2=prepare_values(res_dict2,exclude)
print("KO",np.sum(values_dict2["contractile energy"]))

#plt.figure()
#plt.plot(0,0)
#plt.plot(values_dict1["area"],values_dict1["conractility"],"o",color="C1",label="WT")
#plt.plot(values_dict2["area"],values_dict2["conractility"],"o",color="C2",label="KO")
#plt.legend()


#plt.figure()
#plt.plot(0,0)
#plt.plot(values_dict1["area"],values_dict1["contractile energy"],"o",color="C1",label="WT")
#plt.plot(values_dict2["area"],values_dict2["contractile energy"],"o",color="C2",label="KO")
#plt.legend()

# correlation of contractility and contractile energy







all_types=[name for name in values_dict2.keys() if not name.startswith("std ")]
# performing test
t_test_dict=t_test(values_dict1,values_dict2,all_types)

all_types=list(values_dict1.keys())

lables=["WT","KO"]
fig=plot_contractility_correlation(values_dict1,values_dict2,lables,frame_list1,frame_list2)
fig.savefig(os.path.join(folder_plots,"coordinated_contractility_vs_contractile_energy.png"))



types=['avarage normal stress','avarage shear stress']
ylabels=[ty+"\n"+units[ty] for ty in types]
fig=box_plots(values_dict1,values_dict2,lables,t_test_dict,types=types,ylabels=ylabels)
fig.savefig(os.path.join(folder_plots,"stress_measures_on_the_cell_area.png"))

types=['avarage line stress','avarage cell force',
       'avarage cell pressure','avarage cell shear']
ylabels=[ty+"\n"+units[ty] for ty in types]
fig=box_plots(values_dict1,values_dict2,lables,t_test_dict,types=types,ylabels=ylabels)
fig.savefig(os.path.join(folder_plots,"stress_measures_at_cell_borders.png"))

types=['contractility per area','contractile energy per area']
ylabels=[ty+"\n"+units[ty] for ty in types]
fig=box_plots(values_dict1,values_dict2,lables,t_test_dict,types=types,ylabels=ylabels)
fig.savefig(os.path.join(folder_plots,"contractility_contractile_energy.png"))


types=['area']
ylabels=[ty+"\n"+units[ty] for ty in types]
fig=box_plots(values_dict1,values_dict2,lables,t_test_dict,types=types,ylabels=ylabels)
fig.savefig(os.path.join(folder_plots,"cell area.png"))


