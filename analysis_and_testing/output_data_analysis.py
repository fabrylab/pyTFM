# analysis of results from traction force microscopy and mono layer stress microscopy.
# The analysis has been performed on two data bases with 10 or 12 images of cell colonies. One database contains wild
# type cells, the other contains plectin knockout cells.
# First we set an output folder, then we read in the output files for both databases. Next we normalize all measured
# quantities with the area of the cell colony. Then we perform a two sided t-test to check if there is a significant
# difference between any quantity. Lastly we plot a selection of the quantities in box plots.


# importing tools for reading the output file, analysis and plotting
# also imports dictionary with units for quantities
from andreas_TFM_package.data_analysis import *
import numpy as np
import os


# setting the output folder for plots. All plots are saved to this folder.
folder_plots = "/media/user/GINA1-BK/data_traction_force_microscopy/WT_vs_KO_images/lates_analysis"
createFolder(folder_plots) # creating the folder if it doesn't already exist


## reading in results

# reading the first output file
# list of frames to be excluded. values from these frames are not read in. We don't exclude anything for the wildtype,
# but two frames in the ko, due to issues with imaging the beads.
exclude1=[]
# path to the out.txt text file
file1="/media/user/GINA1-BK/data_traction_force_microscopy/WT_vs_KO_images/WTshift/out2.txt"
parameter_dict1,res_dict1=read_output_file(file1) # reading the file and splitting into parameters and results
# pooling all frames: values_dict has "name of the quantity": "list of values for each frame".
# this also returns the number of frames (n_frames) and a list of the label of frame (frame_list). The frame labels are
# ultimately derived from the number at the beginning of the image file names in your database.
# n_frames is the same for all quantities
n_frames1,values_dict1, frame_list1=prepare_values(res_dict1,exclude1)

# second file
exclude2=["01","10"]  # list of frames to be excluded, thes
# path to the out.txt text file
file2="/media/user/GINA1-BK/data_traction_force_microscopy/WT_vs_KO_images/KOshift/out3.txt"
parameter_dict2,res_dict2=read_output_file(file2)# reading the fie and splitting into parameters and results
n_frames2,values_dict2, frame_list2=prepare_values(res_dict2,exclude2) # pooling all frames


## normalizing the quantities by the area of the cell colony

# units is a dictionary with quantity name: unit of the quantity
# its imported from andreas_TFM_package.parameters_and_strings
# here we add an " per area" for every existing entry in unit and update the unit with /m2
units = add_to_units(units, add_name=" per area", add_unit="/m2",exclude=["area", "cell"]) # adding per cell to units list
# you could do the same with per cells..:
# units2 =a dd_to_units(units,add_name=" per cell", add_unit="",exclude=["area", "cell"]) # adding per number of cells

# now we normalize all quantities. We exclude any quantity with a name that contains strings in exclude. Also std_
#(standart deviations) are not normalized.
# all suitable values are devided by values_dict[norm] and get a new name by adding add_name.
values_dict1=normalize_values(values_dict1,norm="area of colony",add_name=" per area",exclude=["area", "cells"]) # calculating measures per area
values_dict2=normalize_values(values_dict2,norm="area of colony",add_name=" per area",exclude=["area", "cells"])# calculating measures per area



# adding coefficeint of variation
folder1=os.path.split(file1)[0]
folder2=os.path.split(file2)[0]
add_mean_normal_stress_cv(values_dict1,exclude1,folder1,n=6)
add_mean_normal_stress_cv(values_dict2,exclude2,folder2,n=6)


## performing statistical analysis

# getting a list of all values that we want to analyze. There is nothing wrong with analyzing using every quantity.
all_types=[name for name in values_dict2.keys() if not name.startswith("std ") and name in values_dict1.keys()]
# performing a two sided t-test comparing values in values_dict1 with values in values_dict2 by a two sided
# independent t-test
t_test_dict=t_test(values_dict1,values_dict2,all_types)




## plotting

# here we produce a few boxplots and compare a set of quantities in each plot.
# additionally we plot contractillity vs the contractile energy.

# label for the first and second text file; make sure the order is correct
lables=["WT","KO"]

# plotting contractillity vs contractile energy.

types=["contractile energy on colony","contractillity on colony"]
#compare_two_values(values_dict1, values_dict2,types, lables, xlabel,ylabel,frame_list1=[], frame_list2=[]
fig=compare_two_values(values_dict1, values_dict2, types, lables,xlabel="contractile energy [J]",ylabel="contractillity [N]")
# You could also add the frame as a lable to each point, if you want to identify them:
#fig=compare_two_values(values_dict1, values_dict2, types, lables,xlabel="contractile energy [J]",
#                       ylabel="contractillity [N]",frame_list1=frame_list1,frame_list2=frame_list2)
# fig=plot_contractillity_correlation(values_dict1,values_dict2,lables,frame_list1,frame_list2)
# saving to output folder
fig.savefig(os.path.join(folder_plots,"coordinated_contractillity_vs_contractile_energy.png"))

# boxplots for the other measures
# choosing which measures should be displayed in this plot
types=['average normal stress colony','average shear stress colony']
# generating labels for the y axis. This uses units stored in a dictionary imported from
# andreas_TFM_package.parameters_and_strings. ylabels must be a list of strings with length=len(types).
# Of cause you can also set labels manually e.g. ylabels=["label1","label2",...].
ylabels=[ty+"\n"+units[ty] for ty in types]

# plotting box plots, with statistical information. Meaning of the stars:
# ***-> p<0.001 **-> 0.01>p>0.001 *-> 0.05>p>0.01 ns-> p>0.05
fig=box_plots(values_dict1,values_dict2,lables,t_test_dict=t_test_dict,types=types,ylabels=ylabels)
# saving to output folder
fig.savefig(os.path.join(folder_plots,"stress_measures_on_the_cell_area.png"))


# same procedure for some other quantities
# measures of cell cell interactions
types=['avarage line stress','avarage cell force',
       'avarage cell pressure','avarage cell shear']
ylabels=[ty+"\n"+units[ty] for ty in types]
fig=box_plots(values_dict1,values_dict2,lables,t_test_dict=t_test_dict,types=types,ylabels=ylabels)
fig.savefig(os.path.join(folder_plots,"stress_measures_at_cell_borders.png"))

# contractillity and contractile energy
types=['contractillity on colony per area','contractile energy on colony per area']
ylabels=[ty+"\n"+units[ty] for ty in types]
fig=box_plots(values_dict1,values_dict2,lables,t_test_dict=t_test_dict,types=types,ylabels=ylabels,plot_legend=False)
fig.savefig(os.path.join(folder_plots,"contractility_contractile_energy.png"))

# only the colony area
types=['area of colony']
ylabels=[ty+"\n"+units[ty] for ty in types]
fig=box_plots(values_dict1,values_dict2,lables,t_test_dict=t_test_dict,types=types,ylabels=ylabels)
fig.savefig(os.path.join(folder_plots,"cell area.png"))


# coefficeint of variation
types=["coefficient of variation mean normal stress"]
ylabels=[ty+"\n"+units[ty] for ty in types]
fig=box_plots(values_dict1,values_dict2,lables,t_test_dict=t_test_dict,types=types,ylabels=ylabels,plot_legend=False)
fig.savefig(os.path.join(folder_plots,"cv.png"))

# plotting the coefficeint of variation of mean normal stress: this is a measure for stress fluctions on the
# whole cell colony area, including both areas inside and between cells








#######
types=['avarage line stress per area']
ylabels=[ty+"\n"+units[ty] for ty in types]
fig=box_plots(values_dict1,values_dict2,lables,t_test_dict=t_test_dict,types=types,ylabels=ylabels)
fig.savefig(os.path.join(folder_plots,"cell area.png"))



types=["contractile energy on colony","avarage line stress per area"]
#compare_two_values(values_dict1, values_dict2,types, lables, xlabel,ylabel,frame_list1=[], frame_list2=[]
fig=compare_two_values(values_dict1, values_dict2, types, lables,xlabel="contractile energy [J]",
                       ylabel="line stress",frame_list1=frame_list1,frame_list2=frame_list2)



types=["contractile energy on colony","avarage cell force per area"]
#compare_two_values(values_dict1, values_dict2,types, lables, xlabel,ylabel,frame_list1=[], frame_list2=[]
fig=compare_two_values(values_dict1, values_dict2, types, lables,xlabel="contractile energy [J]",
                       ylabel="cell froce",frame_list1=frame_list1,frame_list2=frame_list2)


types=["contractile energy on colony","average normal stress colony per area"]
#compare_two_values(values_dict1, values_dict2,types, lables, xlabel,ylabel,frame_list1=[], frame_list2=[]
fig=compare_two_values(values_dict1, values_dict2, types, lables,xlabel="contractile energy [J]",
                       ylabel="avg norma stress",frame_list1=frame_list1,frame_list2=frame_list2)


# normilizing wiht contractile energy
values_dict1["average normal stress per contractile energy"]=values_dict1["average normal stress colony"]/values_dict1['contractile energy on colony']
values_dict2["average normal stress per contractile energy"]=values_dict2["average normal stress colony"]/values_dict2['contractile energy on colony']
values_dict1["contractillity per contractile energy"]=values_dict1["contractillity on colony"]/values_dict1['contractile energy on colony']
values_dict2["contractillity per contractile energy"]=values_dict2["contractillity on colony"]/values_dict2['contractile energy on colony']

values_dict1["avarage line stress per contractile energy"] = values_dict1["avarage line stress"]/values_dict1['contractile energy on colony']
values_dict2["avarage line stress per contractile energy"] = values_dict1["avarage line stress"]/values_dict2['contractile energy on colony']


all_types=[name for name in values_dict2.keys() if not name.startswith("std ")]
t_test_dict=t_test(values_dict1,values_dict2,all_types)
types=["average normal stress per contractile energy","contractillity per contractile energy","avarage line stress per contractile energy"]

ylabels=[ty+"\n"+units[ty] for ty in types]
fig=box_plots(values_dict1,values_dict2,lables,t_test_dict=t_test_dict,types=types,ylabels=ylabels,low_ylim=0)
fig.savefig(os.path.join(folder_plots,"some_normailzations.png"))