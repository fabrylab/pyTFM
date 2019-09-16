from andreas_TFM_package.TFM_functions_for_clickpoints import units
from andreas_TFM_package.utilities_TFM import *
from andreas_TFM_package.data_analysis import *
from collections import defaultdict
import numpy as np


# first file
file="/media/user/GINA1-BK/data_traktion_force_microscopy/WT_vs_KO_images_10_09_2019/wt_vs_ko_images/WTshift/out.txt"
parameter_dict1,res_dict1=read_output_file(file)
n_frames1,values_dict1,mean_values1,std_values1=prepare_values(res_dict1)
## n is not always the same!!
# second file
file="/media/user/GINA1-BK/data_traktion_force_microscopy/WT_vs_KO_images_10_09_2019/wt_vs_ko_images/KOshift/out.txt"
parameter_dict2,res_dict2=read_output_file(file)
n_frames2,values_dict2,mean_values2,std_values2=prepare_values(res_dict2)

# performing test
t_test_dict=t_test(values_dict1,values_dict2)

all_types=list(values_dict1.keys())
lables=["WT","KO"]
# plotting box plots
types=['avarage normal stress']
box_plots(values_dict1,values_dict2,t_test_dict,lables,types)
types=['avarage shear stress']
box_plots(values_dict1,values_dict2,t_test_dict,lables,types)
types=['avarage line stress','avarage cell force','avarage cell pressure','avarage cell shear']
box_plots(values_dict1,values_dict2,t_test_dict,lables,types)
types=['conractility per cell','contractile energy per cell']
box_plots(values_dict1,values_dict2,t_test_dict,lables,types)


types=['avarage normal stress','avarage shear stress','avarage line stress','avarage cell force',
       'avarage cell pressure','avarage cell shear','conractility per cell','contractile energy per cell']
fig=box_plots(values_dict1,values_dict2,t_test_dict,lables,types=types,normalize=True)
fig.savefig("/media/user/GINA1-BK/data_traktion_force_microscopy/WT_vs_KO_images_10_09_2019/wt_vs_ko_images/eval.png")
