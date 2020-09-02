from pyTFM.data_analysis import *

# reading the Wildtype data set. Use your own output text file here
file_WT = r"/home/andy/Software/example_data_for_pyTFM/clickpoints_tutorial/WT/out.txt"
# reading the parameters and the results, sorted for frames and object ids
parameter_dict_WT, res_dict_WT = read_output_file(file_WT)
# pooling all frames and objects together.
n_frames_WT, values_dict_WT, frame_list_WT = prepare_values(res_dict_WT)
# reading the KO data set. Use your own output text file here
file_KO = r"/home/andy/Software/example_data_for_pyTFM/clickpoints_tutorial/KO/out.txt"
parameter_dict_KO, res_dict_KO = read_output_file(file_KO)
n_frames_KO, values_dict_KO, frame_list_KO = prepare_values(res_dict_KO)

# normalizing the strain energy
values_dict_WT["strain energy per area"] = values_dict_WT["strain energy"]/values_dict_WT["area Cell Area"]
values_dict_KO["strain energy per area"] = values_dict_KO["strain energy"]/values_dict_WT["area Cell Area"]
# normalizing the contractility
values_dict_WT["contractility per area"] = values_dict_WT["contractility"]/values_dict_WT["area Cell Area"]
values_dict_KO["contractility per area"] = values_dict_KO["contractility"]/values_dict_WT["area Cell Area"]

# t-test for all value pairs
t_test_dict = t_test(values_dict_WT,values_dict_KO)

lables = ["WT", "KO"] # designations for the two dictionaries that are provided to the box_plots functions
types = ["contractility per area", "strain energy per area"] # name of the measures that are plotted
ylabels = ["contractility per colony area [N/m²]", "strain energy per colony area [J/m²]"] # custom axes labels
# producing a two box plots comparing the strain energy and the contractility in WT and KO
fig_force = box_plots(values_dict_WT, values_dict_KO, lables, t_test_dict=t_test_dict, types=types,
           low_ylim=0, ylabels=ylabels, plot_legend=True)


lables = ["WT", "KO"] # designations for the two dictionaries that are provided to the box_plots functions
types = ["mean normal stress Cell Area", "average magnitude line tension"] # name of the measures that are plotted
ylabels = ["mean normal stress [N/m]", "line tension [N/m]"] #
fig_stress = box_plots(values_dict_WT, values_dict_KO, lables, t_test_dict=t_test_dict, types=types,
          low_ylim=0, ylabels=ylabels, plot_legend=True)


lables = ["WT", "KO"] # designations for the two dictionaries that are provided to the box_plots functions
# name of the measures that are plotted. Must be length 2 for this case.
types = ["contractility per area", "strain energy per area"]
# plotting value of types[0] vs value of types[1]
fig_force2 = compare_two_values(values_dict_WT, values_dict_KO, types, lables,
         xlabel="contractility per colony area [N/m²]", ylabel="strain energy per colony area [J/m²]")

# define and output folder for your figures
folder_plots = r"/home/andy/Software/example_data_for_pyTFM/clickpoints_tutorial/plots/"
# create the folder, if it doesn't already exist
createFolder(folder_plots)
# saving the three figures that were created beforehand
fig_force.savefig(os.path.join(folder_plots, "forces1.png")) # boxplot comparing measures for force generation
fig_stress.savefig(os.path.join(folder_plots, "fig_stress.png")) # boxplot comparing normal stress and line tension
fig_force2.savefig(os.path.join(folder_plots, "forces2.png")) # plot of strain energy vs contractility
