# functions fro reading in a typical outputfile and performing some statstical analysis as wel as
#data representation

from andreas_TFM_package.TFM_functions_for_clickpoints import units
from andreas_TFM_package.utilities_TFM import *
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.stats import ttest_ind as scipy_ttest_ind

def read_output_file(file_path):

    parameter_dict={}
    res_dict=defaultdict(dict)
    with open(file_path, "r") as f:
        for line in f.readlines():
            elements=line.strip().split("\t")
            if len(elements)>1: # skip if no values can be seperated
                if is_int(elements[0]): # adding to output results if the line starts with the frame
                    res_dict[elements[0]][elements[1]]=float(elements[2])
                else: # else add to parameters dict
                    parameter_dict[elements[0]]=try_float_convert(elements[1])
    return parameter_dict,res_dict
# joing all values from frames

def prepare_values(res_dict,exclude):
    values_dict=defaultdict(list)
    frames=[]
    for frame,subdict in res_dict.items():
        if frame not in exclude:
            frames.append(frame)
            for name ,value in subdict.items():
                values_dict[name].append(value)
    n_frames=len(frames)
    mean_values={name:np.mean(values) for name, values in values_dict.items()}
    std_values={name:np.std(values) for name, values in values_dict.items()}
    return n_frames,values_dict,mean_values,std_values,frames

def t_test(values_dict1,values_dict2,types):
    t_test_dict={}
    for name in values_dict1.keys():
        if name in types:
            t_test_dict[name]= scipy_ttest_ind(values_dict1[name],values_dict2[name])
    return t_test_dict

def set_segnificance_starts(pvalue):
    if pvalue > 0.05:
        stars="ns"
    if pvalue <= 0.05 and pvalue > 0.01:
        stars = "*"
    if pvalue <= 0.01 and pvalue > 0.001:
        stars = "**"
    if pvalue <= 0.001:
        stars = "***"
    return stars
def split_name(name):
    if len(name)>10:
        name_part=name[10:]
        name_part=name_part.replace(" ","\n",1)
        name_return=name[:10]+name_part
    else:
        name_return=copy.deepcopy(name)
    return name_return


def normalize__value_pair(v1,v2): # somewhat weird normalization
    max_abs=np.max(np.abs(np.array(v1+v2)))
    return np.array(v1)/max_abs,np.array(v2)/max_abs

def box_plots(values_dict1,values_dict2,lables,t_test_dict=None,ylabels=[],types=[]):

    if len(types)==0: # use all types
        names=list(values_dict1.keys())
    else:
        names=types
    fig,axs=plt.subplots(1,len(names))
    fig.set_figwidth(len(names)*3)
    axs=make_iterable(axs)
    m_pos=[]
    bps=[]
    for i,(name,ax,label) in enumerate(zip(names,axs,ylabels)):
        bp1=ax.boxplot(values_dict1[name], positions=[0.5], patch_artist=True)
        bp2=ax.boxplot(values_dict2[name], positions=[0.8], patch_artist=True)
        bps.append([bp1,bp2])
        ax.set_ylabel(label)
        m_pos=np.mean([0.5,0.8])
        ax.set_xticks([m_pos])
        name_refined=split_name(name)
        ax.set_xticklabels([name_refined],rotation=60,horizontalalignment="center",multialignment="center")
        ax.ticklabel_format(axis='y', style="sci", scilimits=(-3, 3))

        if isinstance(t_test_dict,dict):
            maxv=np.max(values_dict1[name]+values_dict2[name])
            v_range=maxv-np.min(values_dict1[name]+values_dict2[name])
            ax.plot([0.5,0.5,0.8,0.8],[maxv+v_range*0.05,maxv+v_range*0.1,maxv+v_range*0.1,maxv+v_range*0.05],color="black")  # plotting conecting line
            stars=set_segnificance_starts(t_test_dict[name].pvalue)
            ax.text(m_pos,maxv+v_range*0.15,stars,verticalalignment='center',horizontalalignment='center')

    for bp1,bp2 in bps:
        bp1["boxes"][0].set_facecolor("C1")
        bp2["boxes"][0].set_facecolor("C2")
    #ax=plt.gca()
    #

    plt.plot(0,0,color="C1",label=lables[0])
    plt.plot(0, 0, color="C2", label=lables[1])
    plt.legend(loc="upper right")
    plt.tight_layout()
    #ax.relim()
    #ax.autoscale_view()
    #plt.show()
    return fig


def plot_contractility_correlation(values_dict1, values_dict2, lables, frame_list1=[], frame_list2=[]):
    fig = plt.figure()
    plt.plot(0, 0)  # fixes lower x and y limits
    plt.xlabel("contractility [N]")
    plt.ylabel("contractile energy [J]")
    plt.plot(values_dict1["contractility"], values_dict1["contractile energy"], "o", color="C1", label=lables[0])
    plt.plot(values_dict2["contractility"], values_dict2["contractile energy"], "o", color="C2", label=lables[1])
    if isinstance(frame_list1, list) and isinstance(frame_list2, list):
        for f, v1, v2 in zip(frame_list1, values_dict1["contractility"], values_dict1["contractile energy"]):
            plt.text(v1, v2, f, color="C1")
        for f, v1, v2 in zip(frame_list2, values_dict2["contractility"], values_dict2["contractile energy"]):
            plt.text(v1, v2, f, color="C2")
    plt.ticklabel_format(axis='both', style="sci", scilimits=(-3, 3))
    plt.legend()
    return fig


def full_standard_analysis(res_file1,res_file2,label1,label2,out_folder):


    createFolder(out_folder)
    exclude = []
    # first file
    file1 = res_file1
    parameter_dict1, res_dict1 = read_output_file(file1)
    n_frames1, values_dict1, mean_values1, std_values1, frame_list1 = prepare_values(res_dict1, exclude)
    #print("WT", np.sum(values_dict1["contractile energy"]))
    ## n is not always the same!!
    # second file
    exclude = ["01", "10"]
    file2 = res_file2
    parameter_dict2, res_dict2 = read_output_file(file2)

    n_frames2, values_dict2, mean_values2, std_values2, frame_list2 = prepare_values(res_dict2, exclude)
    #print("KO", np.sum(values_dict2["contractile energy"]))

    all_types = [name for name in values_dict2.keys() if not name.startswith("std ")]
    # performing test
    t_test_dict = t_test(values_dict1, values_dict2, all_types)

    all_types = list(values_dict1.keys())

    lables = [label1,label2]
    plt.ioff()
    fig = plot_contractility_correlation(values_dict1, values_dict2, lables, frame_list1, frame_list2)
    fig.savefig(os.path.join(out_folder, "coordinated_contractility_vs_contractile_energy.png"))
    plt.close(fig)
    types = ['avarage normal stress', 'avarage shear stress']
    ylabels = [ty + "\n" + units[ty] for ty in types]
    fig = box_plots(values_dict1, values_dict2, lables, t_test_dict, types=types, ylabels=ylabels)
    fig.savefig(os.path.join(out_folder, "stress_measures_on_the_cell_area.png"))
    plt.close(fig)
    types = ['avarage line stress', 'avarage cell force',
             'avarage cell pressure', 'avarage cell shear']
    ylabels = [ty + "\n" + units[ty] for ty in types]
    fig = box_plots(values_dict1, values_dict2, lables, t_test_dict, types=types, ylabels=ylabels)
    fig.savefig(os.path.join(out_folder, "stress_measures_at_cell_borders.png"))
    plt.close(fig)
    types = ['contractility per area', 'contractile energy per area']
    ylabels = [ty + "\n" + units[ty] for ty in types]
    fig = box_plots(values_dict1, values_dict2, lables, t_test_dict, types=types, ylabels=ylabels)
    fig.savefig(os.path.join(out_folder, "contractility_contractile_energy.png"))
    plt.close(fig)
    types = ['area']
    ylabels = [ty + "\n" + units[ty] for ty in types]
    fig = box_plots(values_dict1, values_dict2, lables, t_test_dict, types=types, ylabels=ylabels)
    fig.savefig(os.path.join(out_folder, "cell area.png"))
    plt.close(fig)
