# functions fro reading in a typical outputfile and performing some statstical analysis as weel as
#data representation

from andreas_TFM_package.TFM_functions_for_clickpoints import units
from andreas_TFM_package.utilities_TFM import *
from collections import defaultdict
import numpy as np
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
                else: # else add to marmaters dict
                    parameter_dict[elements[0]]=try_float_convert(elements[1])
    return parameter_dict,res_dict
# joing all values from frames

def prepare_values(res_dict):
    values_dict=defaultdict(list)
    frames=[]
    for frame,subdict in res_dict.items():
        frames.append(frame)
        for name ,value in subdict.items():
            values_dict[name].append(value)
    n_frames=len(frames)
    mean_values={name:np.mean(values) for name, values in values_dict.items()}
    std_values={name:np.std(values) for name, values in values_dict.items()}
    return n_frames,values_dict,mean_values,std_values

def t_test(values_dict1,values_dict2):
    t_test_dict={}
    for name in values_dict1.keys():
        t_test_dict[name]= scipy_ttest_ind(values_dict1[name],values_dict2[name])
    return t_test_dict

def normalize__value_pair(v1,v2): # somewhat weird normalization
    max_abs=np.max(np.abs(np.array(v1+v2)))
    return np.array(v1)/max_abs,np.array(v2)/max_abs

def box_plots(values_dict1,values_dict2,t_test_dict,lables,ylabel="",types=[],normalize=False):

    ylabel_use = "relative quantity for each pair" if normalize else ylabel
    fig=plt.figure()
    plt.ylabel(ylabel_use)
    if len(types)==0: # use all types
        names=list(values_dict1.keys())
    else:
        names=types
    m_pos=[]
    bps=[]
    ax=plt.gca()
    for i,name in enumerate(names):
        v1=values_dict1[name]
        v2=values_dict2[name]
        if normalize:
            v1,v2=normalize__value_pair(v1,v2)

        bp1=ax.boxplot(v1, positions=[i+1], patch_artist=True)
        bp2=ax.boxplot(v2, positions=[i + 1.3], patch_artist=True)
        bps.append([bp1,bp2])
        m_pos.append(np.mean([i+1, i+1.3]))
    for bp1,bp2 in bps:
        bp1["boxes"][0].set_facecolor("C1")
        bp2["boxes"][0].set_facecolor("C2")
    #ax=plt.gca()
    #ax.set_xticklabels(names)
    #ax.set_xticks(m_pos)
    plt.xticks(m_pos, names, fontsize=8, rotation=60, horizontalalignment="center", multialignment="right")
    plt.plot(0,0,color="C1",label=lables[0])
    plt.plot(0, 0, color="C2", label=lables[1])
    plt.legend()
    plt.tight_layout()
    #ax.relim()
    #ax.autoscale_view()
    plt.show()
    return fig


