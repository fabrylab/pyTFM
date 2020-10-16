# general usefull function
import copy
import os
import warnings
from collections import defaultdict
import natsort
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage.filters import gaussian_filter


class suppress_warnings():
    def __init__(self, warning_type):
        self.warning_type = warning_type

    def __enter__(self):
        warnings.filterwarnings("ignore", category=self.warning_type)

    def __exit__(self, type, value, traceback):
        warnings.filterwarnings("default", category=self.warning_type)


def make_iterable(value):
    if not hasattr(value, '__iter__') or isinstance(value, str):
        return [value]
    else:
        return value


def convert_none_str(x):
    if isinstance(x, str):
        if x == "None":
            return None
    return x


def make_iterable_args(value):
    # in order to unpack array as one value we need need [array]

    if not hasattr(value, '__iter__') or isinstance(value, str) or isinstance(value, np.ndarray):
        return [value]
    else:
        return value


def convert_axis_tick_unit(ax, factor):
    plt.draw()
    plt.pause(0.001)  # needs to wait for drawing first
    xticks = [convert_to_int(l.get_text()) for l in list(ax.get_xticklabels())]
    ax.set_xticklabels([str(np.round(l * factor)) for l in xticks])
    yticks = [convert_to_int(l.get_text()) for l in list(ax.get_yticklabels())]
    ax.set_yticklabels([str(np.round(l * factor)) for l in yticks])


def convert_str_none(string):
    # converts "None" or "none" to None
    return string if not (string == "None" or string == "none") else None


def make_rank_list(values):

    '''
    produce a list containing the corresponding rank of input values. Ties
    get the same rank. This is used for obtaining the correct sort indices
    (frames in the cdb database) for the list of frames.
    Sorting is performed by the natsort package, wich should recognize signs and scientific notation
    '''

    unique_values = set(values)
    unique_values = natsort.natsorted(unique_values,alg=natsort.REAL)
    unique_values_dict = {value: rank for rank, value in enumerate(unique_values)}  # value:rank of the frame
    rank_list = [unique_values_dict[value] for value in values]
    return rank_list, unique_values_dict


def invert_dictionary(d):
    d_inv = defaultdict(list)
    for key, values in d.items():
        for v in make_iterable(values):
            d_inv[v].append(key)
    return d_inv


def round_flexible(n, digits=2):
    '''
    returns a number rounded to 2 positions after its firs segnificant position
    7.1242*10**-7 --> 7.12*10**-7
    7.1242*10**9 --> 7.12*10**9  and so on
    :param n: float
    :return:
    '''
    if not (isinstance(n, float) or isinstance(n, int)) or n == 0 or np.isnan(n) or np.isinf(n):
        return n
    else:
        rounding_decimal = -int(np.floor(np.log10(np.abs(n)))) + digits
    return np.round(n, rounding_decimal)


def round_flexible_str(n, digits=2, sci_limit=3):
    '''
    returns a number rounded to 2 positions after its firs segnificant position
    7.1242*10**-7 --> 7.12*10**-7
    7.1242*10**9 --> 7.12*10**9  and so on
    :param n: float
    :return:
    '''
    if not (isinstance(n, float) or isinstance(n, int)) or n == 0 or np.isnan(n) or np.isinf(n):
        return n
    else:
        first_pos = -int(np.floor(np.log10(np.abs(n))))
        rounding_decimal = first_pos + digits
        n_round = np.round(n, rounding_decimal)
        if np.abs(rounding_decimal) > sci_limit:
            s = "%" + ".%de" % digits
            string = s % n_round
        elif rounding_decimal >= 0:
            s = "%" + ".%df" % rounding_decimal
            string = s % n_round
        else:
            string = str(n_round)
    return string


def split_path_with_os(folder):
    if not os.path.split(folder)[1] == "":
        parts = [os.path.split(folder)[1]]
        remaining = [os.path.split(folder)[0]]
    else:
        remaining1 = os.path.split(folder)[0]
        parts = [os.path.split(remaining1)[1]]
        remaining = [os.path.split(remaining1)[0]]

    while True:
        path_part = os.path.split(remaining[-1])[1]
        if path_part == "":
            break
        parts.append(path_part)
        remaining.append(os.path.split(remaining[-1])[0])

    return parts


def gaussian_with_nans(arr1, sigma="auto"):
    '''
    applies a gaussian to an array with nans, so that nan values are ignored.
    :param arr1:
    :param sigma: either a string or a float. Sigma as a string: "auto" or "auto-n",n will be applied as a factor
    to determine sigma based on the largest array-axis. Sigma as a float: will directly be used as a sigma for the
    gaussian filter// doesn't depend on axis lengths.
    :return:
    '''
    if isinstance(sigma, str):
        if "auto" in sigma:
            try:
                sigfactor = float(sigma.split("-")[1])
            except IndexError as e:
                sigfactor = 1
                print(e)
            sigma = np.max(arr1.shape) / (500 * sigfactor)  # appropriate sigma

    # not sure why this works....
    arr_zeros = copy.deepcopy(arr1)
    arr_ones = copy.deepcopy(arr1)
    arr_zeros[np.isnan(arr1)] = 0  # array where all nans are replaced by zeros
    arr_ones[~np.isnan(arr1)] = 1
    arr_ones[np.isnan(arr1)] = 0  # array where all nans are replaced by zeros and all other values are replaced by ones

    filter_zeros = gaussian_filter(arr_zeros, sigma=sigma)  # gaussian filter applied to both arrays
    filter_ones = gaussian_filter(arr_ones, sigma=sigma)
    filter_final = filter_zeros / filter_ones  # devision cancles somehow the effect of nan positions
    filter_final[np.isnan(arr1)] = np.nan  # refilling original nans
    return filter_final


def make_display_mask(mask):
    '''
    converts a boolean mask to a mask with 1 and np.nans
    :param mask: np.ndarray with dtype bool or int or suitable float
    :return:
    '''
    mask_show = copy.deepcopy(mask).astype(bool)
    mask_show[~mask_show] = np.nan

    return mask


def find_prefix(n):
    '''
    finds an apropriate prefix (nano, giga, milli..) and returns the rounded number
    before the prefex and the prefix as a string
    :param n: float
    :return:
    '''
    dicct = {-12: "p", -9: "n", -6: "µ", -3: "m", 0: "", 3: "k", 6: "M", 9: "G"}
    exponent = (int(np.floor(np.log10(np.abs(n)) / 3)) * 3)
    # limit to maximal coverd symbols
    exponent = -12 if exponent < -12 else exponent
    exponent = 9 if exponent > 9 else exponent

    n_new = n / 10 ** exponent
    symbol = dicct[exponent]
    return np.round(n_new, 2), symbol


def convert_to_int(a):
    '''
    converts a string to integer. this function is necessary to account for - signe at the beggining in weird formatting
    :param a: string representing an integer
    :return:
    '''

    try:
        n = int(a[0])  # checks if minus signe is present
    except ValueError:
        n = -int(a[1:])
        return n
    return n


def try_int_strip(string):
    try:
        return int(string)
    except ValueError:
        return string.strip("'")


def squeeze_list(l):
    if len(l) == 1:
        if isinstance(l[0], list):
            l = l[0]
    return l


def createFolder(directory):
    '''
    function to create directories if they dont already exist
    '''
    try:
        if not os.path.exists(directory):
            nd = os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

    return directory


def exclude_by_key(d, ex_list):
    ex_dict = {key: values for key, values in d.items() if
               key not in ex_list}
    return ex_dict


def join_dictionary(d1, d2, update_keys=False):
    if update_keys:
        d3 = update_keys(d1, d2)
    else:
        d3 = d2
    return {**d1, **d3}
    # note:z = {**x, **y} and "update" are nices tricks here


def update_keys(d1, d2):
    # recounts keys of d2 by starting with last key of d1. keys must all be integers
    d3 = copy.deepcopy(d2)
    max_key = np.max(list(d1.keys()))
    for key, value in d2.items():
        max_key += 1
        d3[max_key] = value
    return d3


def produce_index_array(u):
    u_indices = np.zeros(np.shape(u))
    for i, j in zip(range(np.shape(u)[0]), range(np.shape(u)[1])):
        u_indices[i, j] = i + j
    return u_indices


def ndargmin(array):
    return np.unravel_index(np.nanargmin(array), array.shape)


def make_random_discrete_color_range(size):
    colors = []
    for i in range(size):
        colors.append('#%06X' % np.random.randint(0, 0xFFFFFF))
    return colors


# function to convert none object from re.search to empty string

def get_group(s=None, group_number=1):
    '''
    This function is used to return empty strings ("") froma re.search object if the object
    is empty. Also it reads out a list of group objects at ones and can read all existing 
    groups. Note that group 0 is per default the entire input string, thus not interesting 
    for us. Also can retrieve all groups
    parameters:
    s- a re.search object
    group_number- integer, specifies which group of the re.search should be returned
    returns:
    string of the group of the re.search object
    '''

    if s is None:  # if no match was found
        return ''

    if isinstance(group_number, str):  # all groups
        if group_number == "all":
            return s.groups()

    if isinstance(group_number, list):  # multiple group ids
        groups = []
        for gn in group_number:
            try:
                groups.append(s.group(gn))
            except IndexError:  # append empty string if group not found
                groups.append('')
        return groups

    try:  # the particular group was not found
        group = s.group(group_number)
    except IndexError:
        return ''

    return s.group(group_number)  # single group


def find_non_nan_region(padded_track):
    '''
    parameters:
    padded_track: Nan-padded Track e.g. from getTracksNanpadded from a clickpointsdata base.
    returns: mask for the non Nan region of the track

    '''
    bound1, bound2 = np.where(~np.isnan(padded_track))[0][0], np.where(~np.isnan(padded_track))[0][
        -1]  # find last and first non nan value
    non_nan_padded_track = np.zeros((len(padded_track),))
    non_nan_padded_track[bound1:bound2] = 1
    non_nan_padded_track = non_nan_padded_track > 0
    return non_nan_padded_track


def is_int(s):
    '''
    checks if string can be converted to int
    :param s:
    :return:
    '''
    try:
        int(s)
        return True
    except ValueError:
        return False


def except_error(func, error, print_error=True, return_v=False, **kwargs):  # take functino and qkwarks
    '''
    wraper to handle errors and return false if the exception is encountered
    :param func:
    :param error:
    :param kwargs:
    :param return_values:
    :return:
    '''

    try:
        values = func(**kwargs)
    except error as e:
        if print_error:
            print(e)
        return return_v
    return values


def try_float_convert(s):
    '''
    checks if string can be converted to int
    :param s:
    :return:
    '''
    try:
        return float(s)
    except ValueError:
        return s


def unpack_list(li):
    if not (isinstance(li, list) or isinstance(li, np.ndarray)):
        return li, ""
    elif len(li) == 1:
        return li[0], ""
    else:
        return li[0], li[1]


def flattten_nested_dict(dict1):
    k_v_list = []
    for k1, v1 in dict1.items():
        if isinstance(v1, (dict, defaultdict)):
            for k2, v2 in v1.items():
                if isinstance(v2, (dict, defaultdict)):
                    for k3, v3 in v1.items():
                        if isinstance(v3, (dict, defaultdict)):
                            for k4, v4 in v3.items():
                                k_v_list.append([k1, k2, k3, k4, v4])
                        else:
                            k_v_list.append([k1, k2, k3, v3])
                else:
                    k_v_list.append([k1, k2, v2])
        else:
            k_v_list.append([k1, v1])

    return k_v_list
