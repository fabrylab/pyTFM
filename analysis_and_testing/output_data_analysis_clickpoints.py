#  analysis of the results

from andreas_TFM_package.TFM_functions_for_clickpoints import units
from andreas_TFM_package.utilities_TFM import *
from andreas_TFM_package.data_analysis import *
from collections import defaultdict
import numpy as np
import os


add_per_area_cell(units,add="cell",dict_type="units") # adding per area
add_per_area_cell(units,add="area",dict_type="units") # adding per cell to units list

# setting the output folder for plots
folder_plots="/media/user/GINA1-BK/data_traktion_force_microscopy/WT_vs_KO_images_10_09_2019/wt_vs_ko_images_Analyzed/results_test/"
createFolder(folder_plots) # creating the folder if it doesn't already exist

# first file
exclude=[] # list of frames to be excluded
# path to the out.txt text file
file1="/media/user/GINA1-BK/data_traktion_force_microscopy/WT_vs_KO_images_10_09_2019/wt_vs_ko_images_Analyzed/WTshift/out.txt"
parameter_dict1,res_dict1=read_output_file(file1) # reading the fie and splitting into parameters and results
n_frames1,values_dict1, frame_list1=prepare_values(res_dict1,exclude) # pooling all frames
add_per_area_cell(values_dict1,add="area",dict_type="values") # calculating measures per area

# second file
exclude=["01","10"]  # list of frames to be excluded
# path to the out.txt text file
file2="/media/user/GINA1-BK/data_traktion_force_microscopy/WT_vs_KO_images_10_09_2019/wt_vs_ko_images_Analyzed/KOshift/out.txt"
parameter_dict2,res_dict2=read_output_file(file2)# reading the fie and splitting into parameters and results
n_frames2,values_dict2, frame_list2=prepare_values(res_dict2,exclude) # pooling all frames
add_per_area_cell(values_dict2,add="area",dict_type="values")# calculating measures per area

# list of all names of measures  found in the outputfile
all_types=[name for name in values_dict2.keys() if not name.startswith("std ")]
# performing statistical test (t-test)
t_test_dict=t_test(values_dict1,values_dict2,all_types)
lables=["WT","KO"] # label for the first and second text file; make sure the order is correct


# plotting contractillity vs contractile energy
fig=plot_contractillity_correlation(values_dict1,values_dict2,lables)
#fig=plot_contractility_correlation(values_dict1,values_dict2,lables,frame_list1,frame_list2)   # labeling the frame of each point
fig.savefig(os.path.join(folder_plots,"coordinated_contractility_vs_contractile_energy.png")) # save to output folder

# plotting other measures
# choosing which measures should be displyed
types=['avarage normal stress','avarage shear stress']
# generating labels for the y axis. This uses units stored in a dictionary imported from some where else
ylabels=[ty+"\n"+units[ty] for ty in types] # list of lables on the y axis, must be list of strings with length =len(types)
# plotting box plots, with statistical information. Meaning of the stars:
# ***-> p<0.001 **-> 0.01>p>0.001 *-> 0.05>p>0.01 ns-> p>0.05
fig=box_plots(values_dict1,values_dict2,lables,t_test_dict=t_test_dict,types=types,ylabels=ylabels)
fig.savefig(os.path.join(folder_plots,"stress_measures_on_the_cell_area.png")) # saving to outputfolder



types=['avarage line stress','avarage cell force',
       'avarage cell pressure','avarage cell shear']
ylabels=[ty+"\n"+units[ty] for ty in types]
fig=box_plots(values_dict1,values_dict2,lables,t_test_dict=t_test_dict,types=types,ylabels=ylabels)
fig.savefig(os.path.join(folder_plots,"stress_measures_at_cell_borders.png"))



types=['contractility per area','contractile energy per area']
ylabels=[ty+"\n"+units[ty] for ty in types]
fig=box_plots(values_dict1,values_dict2,lables,t_test_dict=t_test_dict,types=types,ylabels=ylabels)
fig.savefig(os.path.join(folder_plots,"contractility_contractile_energy.png"))


types=['area']
ylabels=[ty+"\n"+units[ty] for ty in types]
fig=box_plots(values_dict1,values_dict2,lables,t_test_dict=t_test_dict,types=types,ylabels=ylabels)
fig.savefig(os.path.join(folder_plots,"cell area.png"))





