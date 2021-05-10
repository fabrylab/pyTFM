# functions fro reading in a typical outputfile and performing some statstical analysis as wel as
# data representation

import warnings
from collections import defaultdict
from functools import partial

# import matplotlib
# matplotlib.use("Agg")  # uncomment if you want to avoid plots showing up while you are working
import matplotlib.pyplot as plt
import numpy as np
from pyTFM.utilities_TFM import *
from scipy.ndimage import binary_erosion
from scipy.stats import ttest_ind as scipy_ttest_ind

plt.rcParams["axes.edgecolor"] = "#696969"
plt.rcParams["axes.labelcolor"] = "#696969"
plt.rcParams["xtick.color"] = "#696969"
plt.rcParams["ytick.color"] = "#696969"


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
    res_dict = defaultdict(lambda: defaultdict(list))  # nested dictionary to store the results of the analysis
    with open(file_path, "r") as f:
        for line in f.readlines():
            elements = line.strip().split("\t")  # removing tailing "\n" and splitting at "\t"
            if len(elements) > 1:  # skip if no values can be separated
                if is_int(elements[0]):  # adding to output results if the line starts with the frame
                    res_dict[elements[0]][elements[2]].append([elements[1], try_float_convert(elements[3])])  # no warnings
                elif len(elements) > 1:  # else add to parameters dict
                    parameter_dict[elements[0]] = try_float_convert(elements[1])
    return parameter_dict, res_dict


def prepare_values(res_dict, exclude=[]):
    """
    Pooling the values of one quantity for different frames and dropping frames that you want to ignore for
    whatever reason. If you want to certain exclude frames provided a list of them by e.g. exclude=["01","09"]. Exclude
    must be strings.
    :param res_dict: Nested dictionary with key1:frame, key2: name of the measured
    quantity, values. Typically produced from
    ~pyTFM.data_analysis.read_output_file
    :param exclude:
    :return: values_dict: dictionary with key: name of the measured quantity, values: list of measured values
            n_frames: number of frames that were pooled
            frames: list of frames that were pooled. This can be used to reconstruct from which frame an individual
            value comes from.
    """
    ## TODO: implement support for mutlipe objects in a singel image


    values_dict = defaultdict(list)  # out put dictionary
    frames = []  # list of frames that were pooled
    for frame, subdict in res_dict.items():  # iterating trough the frames
        if frame not in exclude:  # ignoring any frame you want ot exclude
            frames.append(frame)  # save the frame
            for name, value in subdict.items():  # iterating through the measured quantities
                values_dict[name].extend([v[1] for v in value])  # appending the value for every frame to values_dict
                # cells are just thrown in together
    for key, value in values_dict.items():  # conversion to array. This is necessary to simplify subsequent calculations.
        values_dict[key] = np.array(value)
    n_frames = len(frames)  # number of frames

    # mean_values={name:np.mean(values) for name, values in values_dict.items()} # return the mean of each quantity
    # std_values={name:np.std(values) for name, values in values_dict.items()}  # return the standard deviation  of each quantity
    return n_frames, values_dict, frames


def normalize_values(values_dict, norm, add_name, exclude=["area", "cells"]):
    """
    Normalize all quantities in "values_dict" by the quantity "norm". Norm must be a key in "values_dict".
    The normalized quantities are added to the old dictionary with a new key as "old key" +"add_name".
    Quantities whos name contains any string in exclude are not normalized. Also quantities that start with
    "std " are ignored.
    :param values_dict:
    :param add:
    :param dict_type:
    :return:
    """

    new_values = defaultdict(partial(np.ndarray, 0))
    for key, value in values_dict.items():
        new_values[key] = value  # copying to new dict
        if not key.startswith("std ") and not any([e in key for e in exclude]) and not key == norm:
            new_values[key + add_name] = value / values_dict[norm]  # adding normalization to new dict

    return new_values


def filter_nans(values_dict1, values_dict2, name):
    '''
    filters nans and warns if nans encountered. Also checks if either of the values is empty.
    :return:
    '''
    v1 = values_dict1[name]
    v2 = values_dict2[name]
    # warning if and stoping if array is empty
    if v1.size == 0 or v2.size == 0:
        warnings.warn("No values in %s. Check output file or calculations." % name)
        return v1, v2, True
    # filtering out nans
    # warning if nans are encountered. This should usually not be the case
    if np.isnan(np.concatenate([v1, v2])).any():
        warnings.warn("nans encounterd in %s. There could be something wrong with the calcualtion." % name)
    # filtering out nans
    v1 = v1[~np.isnan(v1)]
    v2 = v2[~np.isnan(v2)]
    return v1, v2, False


def add_to_units(units_dict, add_name, add_unit, exclude=["area", "cells"]):
    """
    Adds add_name to any key and add_unit to their value except for keys that contain the strings in exclude.
    Also keys with "std " are ignored.

    :param values_dict:
    :param add:
    :param dict_type:
    :return:
    """

    new_units = defaultdict(str)
    for key, value in units_dict.items():
        new_units[key] = value  # copying to new dict
        if not key.startswith("std ") and not any([e in key for e in exclude]):
            new_units[key + add_name] = value + add_unit  # adding new unit to dictionary
    return new_units


def t_test(values_dict1, values_dict2, types=None):
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
    if not isinstance(types, list):
        types = list(set(values_dict1.keys()).intersection(values_dict2.keys()))
    t_test_dict = {}  # output dictionary
    for name in types:  # iterating through the names of all quantities that you want to compare
        if check_types(values_dict1, values_dict2, name):
            try:
                v1, v2, empty = filter_nans(values_dict1, values_dict2, name)
                t_test_dict[name] = scipy_ttest_ind(v1, v2)  # performing a two-sided t-test
            except TypeError:
                pass
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

    if len(name) > 20:  # check the length of the string
        name_part = name[20:]
        name_part = name_part.replace(" ", "\n",
                                      1)  # replace first blank space in the part that is too long with a new line
        name_return = name[:20] + name_part  # rejoining the string
    else:
        name_return = copy.deepcopy(name)  # do nothing if the string is not to long
    return name_return


def box_plots(values_dict1, values_dict2, lables, t_test_dict=None, ylabels=[], types=[], low_ylim=None,
              plot_legend=True):
    """
    Comparing a list of quantities from two experiments by plotting boxplots and optionally adding statistical
    significance stars. The list of quantities to display is given in "types". The results of the experiments
    are provided with the two dictionaries values_dict1 and values_dict2. If a dictionary t_test_dict containing
    p-values (the p-values must be the attribute "pvalue" of the values of this dictionary),
    the statistical significance of differences between the results of the two experiments are shown with stars.
    Stars are generated from the p-value according to ~pyTFM.data_analysis.set_significance_stars .

    :param values_dict1: first dictionary with "name of the measured value":array of measured values.
    :param values_dict2: second dictionary with "name of the measured value":array of measured values.
    :param lables: labels describing values_dict1 and values_dict2.
    :param t_test_dict: dictionary containing p values. "name of the measured value":statistical test object
    with attribute "pvalue".
    :param ylabels: list of labels for the y axis. Must be same length as types.
    :param types: list of names of measured values. Must be keys in both values_dict1 and 2
    :return: fig: figure object


    :usage:

    from pyTFM.data_analysis import *
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
        v1, v2, empty = filter_nans(values_dict1, values_dict2, name)
        maxv = np.max(np.hstack([v1, v2]))  # maximum value from both experiments
        minv = np.min(np.hstack([v1, v2]))  # minimum value from both experiments
        v_range = maxv - minv  # range of values from both experiments
        if empty:
            continue
        bp1 = ax.boxplot(v1, positions=[0.5], patch_artist=True)  # boxlpot of first experiment
        bp2 = ax.boxplot(v2, positions=[0.8], patch_artist=True)  # boxlpot of second experiment
        bps.append([bp1, bp2])  # saving boxplot
        ax.set_ylabel(label)  # setting label of y axis from ylabels list
        m_pos = np.mean([0.5, 0.8])  # x position between two boxes
        ax.set_xticks([m_pos])  # add a tick at this position
        name_refined = split_name(name)  # arrange in two lines if the name is to long
        # label this tick with the name of the measure and
        # central alignment and 60 degree rotation of this label
        ax.set_xticklabels("")
        # ax.set_xticklabels([name_refined], rotation=60, horizontalalignment="center", multialignment="center")
        # show numbers on y axis in scientific notation if they are outside of 10**3 to 10**-3
        ax.ticklabel_format(axis='y', style="sci", scilimits=(-3, 3))
        ax.tick_params(axis="y", labelsize=20)

        # display statistical significance if t_test_dict was provided
        if isinstance(t_test_dict, dict):
            # plotting a line connecting both boxplots
            # finding suitable y values for the edges of this lines
            # drawing the lines
            ax.plot([0.5, 0.5, 0.8, 0.8],
                    [maxv + v_range * 0.05, maxv + v_range * 0.1, maxv + v_range * 0.1, maxv + v_range * 0.05],
                    color="black")

            # plotting stars to symbolize statistical significance
            # generating the correct number of stars from the p value
            stars = set_significance_stars(t_test_dict[name].pvalue)
            # drawing the stars centrally slightly above the line
            ax.text(m_pos, maxv + v_range * 0.15, stars, verticalalignment='center', horizontalalignment='center')
        low_ylim1 = minv if low_ylim == None else low_ylim
        ax.set_ylim(bottom=low_ylim1)
    # labeling the two different experiments
    # changing the color of the boxes in the boxplots according their experiment
    for bp1, bp2 in bps:
        bp1["boxes"][0].set_facecolor("C1")
        bp2["boxes"][0].set_facecolor("C2")
    # adding a legend  with labels from the "lables" list
    plt.plot(0, 0, color="C1", label=lables[0])  # empty line segments as "anchor" of the legend
    plt.plot(0, 0, color="C2", label=lables[1])
    if len(types) == 1:
        plt.xlim((0, 2))  # work around to get a nicer legend position if only onw type is plotted
    if plot_legend:
        plt.legend(loc="upper right")  # fixing legend to upper right corner
    plt.tight_layout()  # improving the general plot layout
    return fig


def check_types(values_dict1, values_dict2, types):
    types = make_iterable(types)
    # checks if any value in the dicts is empty listed in types and throws an appropriate warning
    checks = [len(values_dict1[t]) == 0 for t in types] + [len(values_dict2[t]) == 0 for t in types]
    warn_messages = [t + " empty in dictionary 1" for t in types] + [t + " empty in dictionary 2" for t in types]
    for c, warns in zip(checks, warn_messages):
        if c:
            warnings.warn(warns)
    if any(checks):
        return False
    else:
        return True


def compare_two_values(values_dict1, values_dict2, types, lables, xlabel, ylabel, frame_list1=[], frame_list2=[]):
    """

    Mostly used for:

    ###
    update documentation
    ###

    plotting contractillity vs contractile energy. This plot is supposed to show differences in the correlation of
    contractillity and contractile energy. Contractillity is the sum of the projections of the traction forces to their
    force epicenter. Thus it is high when contraction originating from a single center. Contractile energy is the sum
    of tracktion forces multiplied with the deformations. Therefore it's independent of orientation. A bunch of individual
    cells contracting has a low contractillity compared to its contractile energy. A colony contracting
    "as a single cell" has a high contractillity compared to its contractile energy.
    If you provide lists of frames for each value with frame_list1 adn frame_list2, each point will by labled by its
    corresponding frame
    :param values_dict1: first dictionary with key: name of a quantity, value: measured values for this quantity,
    must contain the keys in types
    :param values_dict2: second dictionary with key: name of a quantity, value: measured values for this quantity,
    must contain the keys in types
    :param types: list with the name of the two quantities that you want to compare. Must be length 2
    :param lables: list of lables describing values_dict1, values_dict2
    :param xlabel: label for x axis
    :param ylabel: label for y axis
    :param frame_list1: optional, list of frames corresponding to the values in values_dict1. Must be list of strings.
    :param frame_list2: optional, list of frames corresponding to the values in values_dict2. Must be list of strings.
    :return: fig, figure object
    """

    fig = plt.figure()
    plt.plot(0, 0)  # fixes lower x and y limits at 0 and 0
    plt.xlabel(xlabel)  # xlabel
    plt.ylabel(ylabel)  # ylabel
    # plotting contractillity vs contractile energy in first experiment
    if not check_types(values_dict1, values_dict2, types):
        return

    plt.plot(values_dict1[types[0]], values_dict1[types[1]], "o", color="C1", label=lables[0], markersize=10,
             markeredgewidth=1, markeredgecolor="black")
    # plotting contractillity vs contractile energy in second experiment
    plt.plot(values_dict2[types[0]], values_dict2[types[1]], "o", color="C2", label=lables[1], markersize=10,
             markeredgewidth=1, markeredgecolor="black")
    # optionally labeling the data points with their corresponding string
    if isinstance(frame_list1, list) and isinstance(frame_list2, list):
        for f, v1, v2 in zip(frame_list1, values_dict1[types[0]], values_dict1[types[1]]):
            plt.text(v1, v2, f, color="C1")
        for f, v1, v2 in zip(frame_list2, values_dict2[types[0]], values_dict2[types[1]]):
            plt.text(v1, v2, f, color="C2")
    # show numbers on y axis in scientific notation if they are outside of 10**3 to 10**-3
    plt.gca().ticklabel_format(axis='both', style="sci", scilimits=(-3, 3))
    plt.gca().tick_params(axis="both", labelsize=20)
    plt.legend()  # adding a legend
    plt.tight_layout()
    return fig


def edge_erosion(stress_tensor, n=7):
    '''
    returns stress tensor with nans replacing n pixels at the cell edge
    :param stress_tensor:
    :param n: ignoring n pixels form the edge of the cell area
    :return:
    '''
    cell_area = stress_tensor[:, :, 0, 0] != 0
    cell_area = binary_erosion(cell_area, iterations=n)
    stress_tensor[~cell_area] = np.nan
    return stress_tensor


def cv_folder(folder, exclude, n=6):
    '''
    reads stress tensors from a folder and calculates the coefficient of variation of the mean normal stress. frames listed
    in exclude will be excluded.
    :param folder: folder with the stress tensors
    :param n: ignore n pixels form the colony edge for the calculation
    :param exclude:  list of frames. Must be the exact string at the begin of the stress tensor filename. E.g ["05","12"]

    :return
    '''
    files = [os.path.join(folder, f) for f in os.listdir(folder) if "stress_tensor.npy" in f]
    cvs = []
    for f in files:
        if not any([os.path.split(f)[1].startswith(e) for e in exclude]):  # checking if frame is in exclude list
            stress_tensor = np.load(f)
            stress_tensor = edge_erosion(stress_tensor, n=n)  # erosion of the edges
            mean_normal_stress = (stress_tensor[:, :, 0, 0] + stress_tensor[:, :, 1, 1]) / 2  # mean normal stress
            cvs.append(
                np.nanstd(mean_normal_stress) / np.nanmean(mean_normal_stress))  # standard deviation of a single colony
    return np.array(cvs)


def add_mean_normal_stress_cv(values_dict1, exclude, folder, n=6):
    '''
    adds the coefficeint of variation to the results dictionary
    :param values_dict1: dictionary, where to which the values are added
    :param folder: folder of the output. This Folder needs to contain the stress tensor files
    :param n: ignore n pixels form the colony edge for the calculation
    :param exclude:  list of frames. Must be the exact string at the begin of the stress tensor filename. E.g ["05","12"]
    :return:
    '''
    values_dict1["coefficient of variation mean normal stress"] = cv_folder(folder, exclude, n=n)


def full_standard_analysis(res_file1, res_file2, label1, label2, out_folder, units):
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

    # setting the output folder for plots. All plots are saved to this folder.
    createFolder(out_folder)  # creating the folder if it doesn't already exist

    ## reading in results

    # reading the first output file
    # list of frames to be excluded. values from these frames are not read in. We don't exclude anything for the wildtype,
    # but two frames in the ko, due to issues with imaging the beads.
    exclude = []
    # path to the out.txt text file

    parameter_dict1, res_dict1 = read_output_file(
        res_file1)  # reading the file and splitting into parameters and results
    # pooling all frames: values_dict has "name of the quantity": "list of values for each frame".
    # this also returns the number of frames (n_frames) and a list of the label of frame (frame_list). The frame labels are
    # ultimately derived from the number at the beginning of the image file names in your database.
    # n_frames is the same for all quantities
    n_frames1, values_dict1, frame_list1 = prepare_values(res_dict1, exclude)

    # second file
    exclude = ["01", "10"]  # list of frames to be excluded, thes
    # path to the out.txt text file

    parameter_dict2, res_dict2 = read_output_file(
        res_file2)  # reading the fie and splitting into parameters and results
    n_frames2, values_dict2, frame_list2 = prepare_values(res_dict2, exclude)  # pooling all frames

    ## normalizing the quantities by the area of the cell colony

    # units is a dictionary with quantity name: unit of the quantity
    # its imported from pyTFM.parameters_and_strings
    # here we add an " per area" for every existing entry in unit and update the unit with /m2
    units = add_to_units(units, add_name=" per area", add_unit="/m2",
                         exclude=["area", "cell"])  # adding per cell to units list
    # you could do the same with per cells..:
    # units2 =a dd_to_units(units,add_name=" per cell", add_unit="",exclude=["area", "cell"]) # adding per number of cells

    # now we normalize all quantities. We exclude any quantity with a name that contains strings in exclude. Also std_
    # (standart deviations) are not normalized.
    # all suitable values are devided by values_dict[norm] and get a new name by adding add_name.
    values_dict1 = normalize_values(values_dict1, norm="area of colony", add_name=" per area",
                                    exclude=["area", "cells"])  # calculating measures per area
    values_dict2 = normalize_values(values_dict2, norm="area of colony", add_name=" per area",
                                    exclude=["area", "cells"])  # calculating measures per area

    ## performing statistical analysis

    # getting a list of all values that we want to analyze. There is nothing wrong with analyzing using every quantity.
    all_types = [name for name in values_dict2.keys() if not name.startswith("std ")]
    # performing a two sided t-test comparing values in values_dict1 with values in values_dict2 by a two sided
    # independent t-test
    t_test_dict = t_test(values_dict1, values_dict2, all_types)

    ## plotting

    # here we produce a few boxplots and compare a set of quantities in each plot.
    # additionally we plot contractillity vs the contractile energy.

    # label for the first and second text file; make sure the order is correct
    lables = [label1, label2]

    # plotting contractillity vs contractile energy.

    types = ["contractile energy on cell colony", "contractillity on cell colony"]
    # compare_two_values(values_dict1, values_dict2,types, lables, xlabel,ylabel,frame_list1=[], frame_list2=[]
    fig = compare_two_values(values_dict1, values_dict2, types, lables, xlabel="contractile energy [J]",
                             ylabel="contractillity [N]")
    # You could also add the frame as a lable to each point, if you want to identify them:
    fig = compare_two_values(values_dict1, values_dict2, types, lables, xlabel="contractile energy [J]",
                             ylabel="contractillity [N]", frame_list1=frame_list1, frame_list2=frame_list2)
    # fig=plot_contractillity_correlation(values_dict1,values_dict2,lables,frame_list1,frame_list2)
    # saving to output folder
    fig.savefig(os.path.join(out_folder, "coordinated_contractillity_vs_contractile_energy.png"))

    # boxplots for the other measures
    # choosing which measures should be displayed in this plot
    types = ['average normal stress colony', 'average shear stress colony']
    # generating labels for the y axis. This uses units stored in a dictionary imported from
    # pyTFM.parameters_and_strings. ylabels must be a list of strings with length=len(types).
    # Of cause you can also set labels manually e.g. ylabels=["label1","label2",...].
    ylabels = [ty + "\n" + units[ty] for ty in types]

    # plotting box plots, with statistical information. Meaning of the stars:
    # ***-> p<0.001 **-> 0.01>p>0.001 *-> 0.05>p>0.01 ns-> p>0.05
    fig = box_plots(values_dict1, values_dict2, lables, t_test_dict=t_test_dict, types=types, ylabels=ylabels)
    # saving to output folder
    fig.savefig(os.path.join(out_folder, "stress_measures_on_the_cell_area.png"))

    # same procedure for some other quantities
    # measures of cell cell interactions
    types = ['avarage line stress', 'avarage cell force',
             'avarage cell pressure', 'avarage cell shear']
    ylabels = [ty + "\n" + units[ty] for ty in types]
    fig = box_plots(values_dict1, values_dict2, lables, t_test_dict=t_test_dict, types=types, ylabels=ylabels)
    fig.savefig(os.path.join(out_folder, "stress_measures_at_cell_borders.png"))

    # contractillity and contractile energy
    types = ['contractillity on cell colony per area', 'contractile energy on cell colony per area']
    ylabels = [ty + "\n" + units[ty] for ty in types]
    fig = box_plots(values_dict1, values_dict2, lables, t_test_dict=t_test_dict, types=types, ylabels=ylabels)
    fig.savefig(os.path.join(out_folder, "contractility_contractile_energy.png"))

    # only the colony area
    types = ['area']
    ylabels = [ty + "\n" + units[ty] for ty in types]
    fig = box_plots(values_dict1, values_dict2, lables, t_test_dict=t_test_dict, types=types, ylabels=ylabels)
    fig.savefig(os.path.join(out_folder, "cell area.png"))
