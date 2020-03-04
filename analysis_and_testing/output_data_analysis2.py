from pyTFM.data_analysis import *


folder_plots = "/home/user/Software/pyTFM/example_analysis/plots/"
createFolder(folder_plots) # creating the folder if it doesn't already exist

file_WT="/home/user/Software/pyTFM/example_analysis/WT_analyzed/out.txt"
parameter_dict_WT,res_dict_WT=read_output_file(file_WT)
n_frames_WT,values_dict_WT, frame_list_WT=prepare_values(res_dict_WT)

file_KO="/home/user/Software/pyTFM/example_analysis/KO_analyzed/out.txt"
parameter_dict_KO,res_dict_KO=read_output_file(file_KO)# reading the fie and splitting into parameters and results
n_frames_KO,values_dict_KO, frame_list_KO=prepare_values(res_dict_KO) # pooling all frames


values_dict_WT["strain energy per area"]=values_dict_WT["strain energy on colony"]/values_dict_WT["area of colony"]
values_dict_KO["strain energy per area"]=values_dict_KO["strain energy on colony"]/values_dict_WT["area of colony"]

values_dict_WT["contractillity per area"]=values_dict_WT["contractillity on colony"]/values_dict_WT["area of colony"]
values_dict_KO["contractillity per area"]=values_dict_KO["contractillity on colony"]/values_dict_WT["area of colony"]

# label for the first and second text file; make sure the order is correct
lables=["WT","KO"]
t_test_dict=t_test(values_dict_WT,values_dict_KO)
types=["contractillity per area","strain energy per area"]
ylabels=["contractillity per colony area [N/m²]","strain energy per colony area [J/m²]"]
fig_force=box_plots(values_dict_WT,values_dict_KO,lables,t_test_dict=t_test_dict,types=types,low_ylim=0,ylabels=ylabels,plot_legend=True)
fig_force.savefig(os.path.join(folder_plots,"forces1.png"))
fig_force2=compare_two_values(values_dict_WT, values_dict_KO, types, lables,xlabel="contractillity per colony area [N/m²]",ylabel="strain energy per colony area [J/m²]")
fig_force2.savefig(os.path.join(folder_plots,"forces2.png"))

types=["mean normal stress on colony","average magnitude line tension"]
ylabels=["mean normal stress [N/m]","line tension [N/m]"]
fig_stress=box_plots(values_dict_WT,values_dict_KO,lables,t_test_dict=t_test_dict,types=types,low_ylim=0,ylabels=ylabels,plot_legend=True)
fig_stress.savefig(os.path.join(folder_plots,"fig_stress.png"))

