# functions fro reading in a typical outputfile and performing some statstical analysis as wel as
# data representation

from andreas_TFM_package.TFM_functions_for_clickpoints import units
from andreas_TFM_package.utilities_TFM import *
from collections import defaultdict
import numpy as np
# import matplotlib
# matplotlib.use("Agg")  # uncomment if you want to avoid plots showing up while you are working
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind as scipy_ttest_ind


def read_output_file(file_path):
    """
    reading and out put file. The file will be split in a dictionary containing parameters and a dictionary containing
    results. The text file must be tab separated ("\t") were each line is either a parameter line or a results line.
    A parameter line has a string in the first column and a string or number in the second. Column 1 and 2 are assigned
    as key and value pair to "parameter_dict".
    A results line starts with an "integer" or can be converted to one (typically "01","10"..). This number
    represents the frame from which these values originate. The second column is a string with the name of the measured
    quantity. The third column is a float or int with the measured value. Forth and fifth columns can contain strings
    with the unit of the measured quantity and a warning. Both might also be empty and are not added to "res_dict".
    The results are stored in res_dict as a nested dictionary key1:frame(from column 1), key2: name of the measured
    quantity(from column 2), values (from column 3).

    :param file_path: path to the output text file
    :return: parameter_dict,res_dict dictionaries containing the parameters and the analysis results
    """

    parameter_dict = {}  # dictionary to store the parameters
    res_dict = defaultdict(dict)  # nested dictionary to store the results of the analysis
    with open(file_path, "r") as f:
        for line in f.readlines():
            elements = line.strip().split("\t")  # removing tailing "\n" and splitting at "\t"
            if len(elements) > 1:  # skip if no values can be separated
                if is_int(elements[0]):  # adding to output results if the line starts with the frame
                    res_dict[elements[0]][elements[1]] = float(elements[2])
                else:  # else add to parameters dict
                    parameter_dict[elements[0]] = try_float_convert(elements[1])
    return parameter_dict, res_dict


def prepare_values(res_dict, exclude=[]):
    """
    Pooling the values of one quantity for different frames and dropping frames that you want to ignore for
    whatever reason. If you want to certain exclude frames provided a list of them by e.g. exclude=["01","09"]. Exclude
    must be strings.
    :param res_dict: Nested dictionary with key1:frame, key2: name of the measured
    quantity, values. Typically produced from
    ~andreas_TFM_package.data_analysis.read_output_file
    :param exclude:
    :return: values_dict: dictionary with key: name of the measured quantity, values: list of measured values
            n_frames: number of frames that were pooled
            frames: list of frames that were pooled. This can be used to reconstruct from which frame an individual
            value comes from.
    """

    values_dict = defaultdict(list)  # out put dictionary
    frames = []  # list of frames that were pooled
    for frame, subdict in res_dict.items():  # iterating trough the frames
        if frame not in exclude:  # ignoring any frame you want ot exclude
            frames.append(frame)  # save the frame
            for name, value in subdict.items():  # iterating through the measured quantities
                values_dict[name].append(value)  # appending the value for every frame to values_dict
    for key, value in values_dict.items():  # conversion to array. This is necessary to simplify subsequent calculations.
        values_dict[key] = np.array(value)
    n_frames = len(frames)  # number of frames

    # mean_values={name:np.mean(values) for name, values in values_dict.items()} # return the mean of each quantity
    # std_values={name:np.std(values) for name, values in values_dict.items()}  # return the standard deviation  of each quantity
    return n_frames, values_dict, frames


def add_per_area_cell(values_dict, add="cell", dict_type="values"):
    """
    Either normalizing by the area of the cell colony or the number of cells, or generating new units.

    if dict_type =="values" and ad=="area":
    Adds for each existing quantity (except "area" and "n_cells") another quantity with " per area" added to the
    original name. Then calculates the values of this new quantity by dividing the original values by the area
    of the cell colony.
    if dict_type =="values" and ad=="cell":
    Adds for each existing quantity (except "area" and "n_cells") another quantity with " per cell" added to the
    original name. Then calculates the values of this new quantity by dividing the original values by the number of
    cells in the cell colony. Caution the number of cells in the colony can be a little lower then what you have
    drawn, especially if you use low resolution deformation fields.
    if dict_type =="units" and ad=="area"
    Adds an entry to the units dictionary with " per area" added to the original name and "/m2" added to the original
    value.
    if dict_type =="units" and ad=="cell"
    Adds an entry to the units dictionary with " per cell" added to the original name but keeps the value as it was.

    :param values_dict:
    :param add:
    :param dict_type:
    :return:
    """

    assert add in ["cell", "area"]  # checking for correct input
    assert dict_type in ["values", "units"]

    add_strings = [" per cell", " per area"]
    if add == "cell":
        add_str = add_strings[0]
        add_key = "n_cells"
        add_unit = ""
    if add == "area":
        add_str = add_strings[1]
        add_key = "area"
        add_unit = "/m2"

    dict_cp = copy.deepcopy(values_dict)
    if dict_type == "values":
        for key, value in dict_cp.items():
            if not key.startswith("std ") and add_strings[0] not in key and add_strings[
                1] not in key and key != "area" and key != "n_cells":  # ignore entries with standart deviation
                values_dict[key + add_str] = value / values_dict[add_key]
    if dict_type == "units":
        for key, value in dict_cp.items():
            if not key.startswith("std ") and add_strings[0] not in key and add_strings[
                1] not in key and key != "area" and key != "n_cells":  # ignore entries with standart deviation
                values_dict[key + add_str] = value + add_unit


def t_test(values_dict1, values_dict2, types):
    """
    Performing a two-sided t-test with two set of values. The values are provided in the two dictionaries values_dict1 and
    values_dict2. Both must have the same keys. Results are stored in a new dictionary with key: name of the
    compared quantity (simply the key from values_dict1 and 2), value: t-test object generated
    from scpy.stats.ttest_ind. For further uses this object needs the attribute "pvalue".
    see also https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html
    Only the quantities listed in "types" are compared. If you want to compare just one quantity you still have to
    set types as list: types=["quantity 1"].

    :param values_dict1: first set of values. Dictionary with key: name of quantity, values: list of measured values
    :param values_dict2: second set of values. Dictionary with key: name of quantity, values: list of measured values
    :param types: list of quantities ot compare
    :return:  t_test_dict, dictionary with key: name of a quantity, value: t-test object
    """

    t_test_dict = {}  # output dictionary
    for name in types:  # iterating through the names of all quantities that you want to compare
        t_test_dict[name] = scipy_ttest_ind(values_dict1[name], values_dict2[name])  # performing a two-sided t-test
    return t_test_dict


def set_significance_stars(pvalue):
    """
    converting a pvalue to a representation in significance stars.
    The current conversion is # ***-> p<0.001 **-> 0.01>p>0.001 *-> 0.05>p>0.01 ns-> p>0.05. But feel free to change
    that.

    :param pvalue: any float
    :return: stars: string representing the significance level
    """

    if pvalue > 0.05:
        stars = "ns"  # not significant
    if pvalue <= 0.05 and pvalue > 0.01:
        stars = "*"
    if pvalue <= 0.01 and pvalue > 0.001:
        stars = "**"
    if pvalue <= 0.001:
        stars = "***"
    return stars


def split_name(name):
    """
    automatically rearranges a string in two lines if it is to long (>10 characters). The
    new line is inserted at the first blank space (" ") after the 10th chrackter. This is used while plotting
    to make lables look nicer.
    :param name:
    :return:
    """

    if len(name) > 10:  # check the length of the string
        name_part = name[10:]
        name_part = name_part.replace(" ", "\n",
                                      1)  # replace first blank space in the part that is too long with a new line
        name_return = name[:10] + name_part  # rejoining the string
    else:
        name_return = copy.deepcopy(name)  # do nothing if the string is not to long
    return name_return


def box_plots(values_dict1, values_dict2, lables, t_test_dict=None, ylabels=[], types=[]):
    """
    Comparing a list of quantities from two experiments by plotting boxplots and optionally adding statistical
    significance stars. The list of quantities to display is given in "types". The results of the experiments
    are provided with the two dictionaries values_dict1 and values_dict2. If a dictionary t_test_dict containing
    p-values (the p-values must be the attribute "pvalue" of the values of this dictionary),
    the statistical significance of differences between the results of the two experiments are shown with stars.
    Stars are generated from the p-value according to ~andreas_TFM_package.data_analysis.set_significance_stars .

    :param values_dict1: first dictionary with "name of the measured value":array of measured values.
    :param values_dict2: second dictionary with "name of the measured value":array of measured values.
    :param lables: labels describing values_dict1 and values_dict2.
    :param t_test_dict: dictionary containing p values. "name of the measured value":statistical test object
    with attribute "pvalue".
    :param ylabels: list of labels for the y axis. Must be same length as types.
    :param types: list of names of measured values. Must be keys in both values_dict1 and 2
    :return: fig: figure object


    :usage:

    from andreas_TFM_package.data_analysis import *
    ## inputs:
    # two dictionaries containing the measured values.
    #
    # both must have the same keys. Values must be array, but can be of different lenght
    values_dict1 = {'avarage normal stress': np.array([1, 2, 3, 4]),
                    'avarage shear stress': np.array([1, 2, 3, 4])}
    values_dict2 = {'avarage normal stress': np.array([7, 10, 4, 5]),
                    'avarage shear stress': np.array([7, 10, 4, 5])}
    # list of types to be plotted. These must be keys in the dictionaries above
    types = ['avarage normal stress', 'avarage shear stress']
    # list of lables that describe the content of values_dict1 and values_dict2 in this order
    lables = ["KO", "WT"]
    # list of lables on the y axis. Must be same length as  types
    ylabels = ['avarage normal stress in N/m', 'avarage shear stress in N/m']

    ## plotting box plots
    fig = box_plots(values_dict1, values_dict2, lables, types=types, ylabels=ylabels)


    ## adding stastical information
    from scipy.stats import ttest_ind as scipy_ttest_ind
    # creating a dictionary with test results. All keys in the list types must be present.
    # Values must have "pvalue" as attribute
    t_test_dict = {}
    t_test_dict['avarage normal stress'] = scipy_ttest_ind(values_dict1['avarage normal stress'],
                                                           values_dict2['avarage normal stress'])
    t_test_dict['avarage shear stress'] = scipy_ttest_ind(values_dict1['avarage shear stress'],
                                                          values_dict2['avarage shear stress'])

    # plotting box plots, with statistical information. Meaning of the stars:
    # ***-> p<0.001 **-> 0.01>p>0.001 *-> 0.05>p>0.01 ns-> p>0.05
    fig = box_plots(values_dict1, values_dict2, lables, t_test_dict=t_test_dict, types=types, ylabels=ylabels)
    """

    fig, axs = plt.subplots(1, len(types))  # creating subplot with one row and one column per quantity in types
    fig.set_figwidth(len(types) * 3)  # roughly choosing the correct width for the figure
    axs = make_iterable(axs)  # returns [axs] if axs is a single object and not a list. This allows iteration.

    bps = []  # list of boxplot objects. We need this to adjust the color of the boxes.
    # generating boxplots for each quantity in types
    for i, (name, ax, label) in enumerate(zip(types, axs, ylabels)):
        bp1 = ax.boxplot(values_dict1[name], positions=[0.5], patch_artist=True)  # boxlpot of first experiment
        bp2 = ax.boxplot(values_dict2[name], positions=[0.8], patch_artist=True)  # boxlpot of second experiment
        bps.append([bp1, bp2])  # saving boxplot
        ax.set_ylabel(label)  # setting label of y axis from ylabels list
        m_pos = np.mean([0.5, 0.8])  # x position between two boxes
        ax.set_xticks([m_pos])  # add a tick at this position
        name_refined = split_name(name)  # arrange in two lines if the name is to long
        # label this tick with the name of the measure and
        # central alignment and 60 degree rotation of this label
        ax.set_xticklabels([name_refined], rotation=60, horizontalalignment="center", multialignment="center")
        # show numbers on y axis in scientific notation if they are outside of 10**3 to 10**-3
        ax.ticklabel_format(axis='y', style="sci", scilimits=(-3, 3))

        # display statistical significance if t_test_dict was provided
        if isinstance(t_test_dict, dict):
            # plotting a line connecting both boxplots
            # finding suitable y values for the edges of this lines
            maxv = np.max(np.hstack([values_dict1[name], values_dict2[name]]))  # maximum value from both experiments
            minv = np.min(np.hstack([values_dict1[name], values_dict2[name]]))  # minimum value from both experiments
            v_range = maxv - minv  # range of values from both experiments
            # drawing the lines
            ax.plot([0.5, 0.5, 0.8, 0.8],
                    [maxv + v_range * 0.05, maxv + v_range * 0.1, maxv + v_range * 0.1, maxv + v_range * 0.05],
                    color="black")

            # plotting stars to symbolize statistical significance
            # generating the correct number of stars from the p value
            stars = set_significance_stars(t_test_dict[name].pvalue)
            # drawing the stars centrally slightly above the line
            ax.text(m_pos, maxv + v_range * 0.15, stars, verticalalignment='center', horizontalalignment='center')

    # labeling the two different experiments
    # changing the color of the boxes in the boxplots according their experiment
    for bp1, bp2 in bps:
        bp1["boxes"][0].set_facecolor("C1")
        bp2["boxes"][0].set_facecolor("C2")
    # adding a legend  with labels from the "lables" list
    plt.plot(0, 0, color="C1", label=lables[0])  # empty line segments as "anchor" of the legend
    plt.plot(0, 0, color="C2", label=lables[1])
    plt.legend(loc="upper right")  # fixing legend to upper right corner
    plt.tight_layout()  # improving the general plot layout
    return fig


def plot_contractillity_correlation(values_dict1, values_dict2, lables, frame_list1=[], frame_list2=[]):
    """
    plotting contractillity vs contractile energy. This plot is supposed to show differences in the correlation of
    contractillity and contractile energy. Contractillity is the sum of the projections of the traction forces to their
    force epicenter. Thus it is high when contraction originating from a single center. Contractile energy is the sum
    of tracktion forces multiplied with the deformations. Therefore it's independent of orientation. A bunch of individual
    cells contracting has a low contractillity compared to its contractile energy. A colony contracting
    "as a single cell" has a high contractillity compared to its contractile energy.
    If you provide lists of frames for each value with frame_list1 adn frame_list2, each point will by labled by its
    corresponding frame
    :param values_dict1: first dictionary with key: name of a quantity, value: measured values for this quantity,
    must contain the keys "contractility" and "contractile energy"
    :param values_dict2: second dictionary with key: name of a quantity, value: measured values for this quantity,
    must contain the keys "contractility" and "contractile energy"
    :param lables: list of lables describing values_dict1, values_dict2
    :param frame_list1: optional, list of frames corresponding to the values in values_dict1. Must be list of strings.
    :param frame_list2: optional, list of frames corresponding to the values in values_dict2. Must be list of strings.
    :return: fig, figure object
    """

    fig = plt.figure()
    plt.plot(0, 0)  # fixes lower x and y limits at 0 and 0
    plt.xlabel("contractility [N]")  # xlabel
    plt.ylabel("contractile energy [J]")  # label
    # plotting contractillity vs contractile energy in first experiment
    plt.plot(values_dict1["contractility"], values_dict1["contractile energy"], "o", color="C1", label=lables[0])
    # plotting contractillity vs contractile energy in second experiment
    plt.plot(values_dict2["contractility"], values_dict2["contractile energy"], "o", color="C2", label=lables[1])
    # optionally labeling the data points with their corresponding string
    if isinstance(frame_list1, list) and isinstance(frame_list2, list):
        for f, v1, v2 in zip(frame_list1, values_dict1["contractility"], values_dict1["contractile energy"]):
            plt.text(v1, v2, f, color="C1")
        for f, v1, v2 in zip(frame_list2, values_dict2["contractility"], values_dict2["contractile energy"]):
            plt.text(v1, v2, f, color="C2")
    # show numbers on y axis in scientific notation if they are outside of 10**3 to 10**-3
    plt.ticklabel_format(axis='both', style="sci", scilimits=(-3, 3))
    plt.legend()  # adding a legend
    return fig


def full_standard_analysis(res_file1, res_file2, label1, label2, out_folder):
    """
    current anaylyis of an experiment. The results from two samples in res_file1 and res_file2 (text files
    produced from the TFM clickpoints addon) are compared for several quantities. Each quantity is compared with two box
    plots and a t-test. Figures are saved to the output folder. The folder is generated if it is not already
    existing. You need to label each input text file with a string (label1 and label2). E.g. label1="KO", label2="WT".
    :param res_file1: path to first output text file
    :param res_file2: path to second output text file
    :param label1: label of the first output textfile
    :param label2: label of the second output textfile
    :param out_folder: output folder for plots
    :return:
    """

    ### note_ import units from andreas_TM_package.TFM-functions_for_clickpoints
    add_per_area_cell(units, add="cell", dict_type="units")  # adding per area
    add_per_area_cell(units, add="area", dict_type="units")  # adding per cell to units list

    # setting the output folder for plots
    createFolder(out_folder)  # creating the folder if it doesn't already exist

    # first file
    exclude = []  # list of frames to be excluded
    # path to the out.txt text file

    parameter_dict1, res_dict1 = read_output_file(
        res_file1)  # reading the fie and splitting into parameters and results
    n_frames1, values_dict1, frame_list1 = prepare_values(res_dict1, exclude)  # pooling all frames
    add_per_area_cell(values_dict1, add="area", dict_type="values")  # calculating measures per area

    # second file
    exclude = ["01", "10"]  # list of frames to be excluded
    # path to the out.txt text file

    parameter_dict2, res_dict2 = read_output_file(
        res_file2)  # reading the fie and splitting into parameters and results
    n_frames2, values_dict2, frame_list2 = prepare_values(res_dict2, exclude)  # pooling all frames
    add_per_area_cell(values_dict2, add="area", dict_type="values")  # calculating measures per area

    # list of all names of measures  found in the outputfile
    all_types = [name for name in values_dict2.keys() if not name.startswith("std ")]
    # performing statistical test (t-test)
    t_test_dict = t_test(values_dict1, values_dict2, all_types)
    lables = [label1, label2]  # label for the first and second text file; make sure the order is correct

    # plotting contractillity vs contractile energy
    fig = plot_contractillity_correlation(values_dict1, values_dict2, lables)
    # fig=plot_contractility_correlation(values_dict1,values_dict2,lables,frame_list1,frame_list2)   # labeling the frame of each point
    fig.savefig(
        os.path.join(out_folder, "coordinated_contractility_vs_contractile_energy.png"))  # save to output folder

    # plotting other measures
    # choosing which measures should be displyed
    types = ['avarage normal stress', 'avarage shear stress']
    # generating labels for the y axis. This uses units stored in a dictionary imported from some where else
    ylabels = [ty + "\n" + units[ty] for ty in
               types]  # list of lables on the y axis, must be list of strings with length =len(types)
    # plotting box plots, with statistical information. Meaning of the stars:
    # ***-> p<0.001 **-> 0.01>p>0.001 *-> 0.05>p>0.01 ns-> p>0.05
    fig = box_plots(values_dict1, values_dict2, lables, t_test_dict=t_test_dict, types=types, ylabels=ylabels)
    fig.savefig(os.path.join(out_folder, "stress_measures_on_the_cell_area.png"))  # saving to outputfolder

    types = ['avarage line stress', 'avarage cell force',
             'avarage cell pressure', 'avarage cell shear']
    ylabels = [ty + "\n" + units[ty] for ty in types]
    fig = box_plots(values_dict1, values_dict2, lables, t_test_dict=t_test_dict, types=types, ylabels=ylabels)
    fig.savefig(os.path.join(out_folder, "stress_measures_at_cell_borders.png"))

    types = ['contractility per area', 'contractile energy per area']
    ylabels = [ty + "\n" + units[ty] for ty in types]
    fig = box_plots(values_dict1, values_dict2, lables, t_test_dict=t_test_dict, types=types, ylabels=ylabels)
    fig.savefig(os.path.join(out_folder, "contractility_contractile_energy.png"))

    types = ['area']
    ylabels = [ty + "\n" + units[ty] for ty in types]
    fig = box_plots(values_dict1, values_dict2, lables, t_test_dict=t_test_dict, types=types, ylabels=ylabels)
    fig.savefig(os.path.join(out_folder, "cell area.png"))


