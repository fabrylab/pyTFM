# general usefull function
import matplotlib as plt
import numpy as np
import os
import copy

def make_iterable(value):
    if not hasattr(value, '__iter__') or isinstance(value,str):
        return [value]
    else:
        return value

def make_paramters_dict_tfm(**kwargs):
    '''
    adding default values and "pixelsize2"
    :param kwargs:
    :return:
    '''
    paramter_dict=copy.deepcopy(kwargs)
    if "pixelsize" in paramter_dict.keys():
        paramter_dict["pixelsize_beads_image"]=paramter_dict["pixelsize"]
        del paramter_dict["pixelsize"]



    if not "pixelsize_def_image" in kwargs.keys():
        paramter_dict["pixelsize_def_image"]=paramter_dict["pixelsize_beads_image"]*(paramter_dict["window_size"]-paramter_dict["overlapp"])

    if not "std_factor" in paramter_dict.keys():
        paramter_dict["std_factor"]=15

    return paramter_dict

def convert_axis_tick_unit(ax,factor):
    plt.draw()
    plt.pause(0.001)  # needs to wait for drawing first
    xticks = [convert_to_int(l.get_text()) for l in list(ax.get_xticklabels())]
    ax.set_xticklabels([str(np.round(l * factor)) for l in xticks])
    yticks = [convert_to_int(l.get_text()) for l in list(ax.get_yticklabels())]
    ax.set_yticklabels([str(np.round(l * factor)) for l in yticks])


    
def make_rank_list(values,dtype=int):
    '''
    produce a list containing the corresponding rank of input values. Ties
    get the same rank. This is used for optaining the correct sort indices 
    (frames in the cdb database) for the list of frames.
    '''
    
    values_conv=np.array(values,dtype=dtype) #allows only int or string that can be converted
    unique_values= np.unique(values_conv)# already sorted
    unique_values_dict={value:rank for rank,value in enumerate(unique_values)} # value:rank of the frame
    rank_list = [unique_values_dict[value] for value in values_conv]
    return rank_list

def round_flexible(n):
    '''
    returns a number rounded to 2 positions after its firs segnificant position
    7.1242*10**-7 --> 7.12*10**-7
    7.1242*10**9 --> 7.12*10**9  and so on
    :param n: float
    :return:
    '''

    rounding_decimal=-int(np.floor(np.log10(np.abs(n)))) + 2
    return np.round(n,rounding_decimal)
def find_prefix(n):
    '''
    finds an apropriate prefix (nano, giga, milli..) and returns the rounded number
    before the prefex and the prefix as a string
    :param n: float
    :return:
    '''
    dicct = {-12:"p",-9:"n",-6:"µ",-3:"m",0:"", 3:"k", 6: "M", 9: "G"}
    exponent=(int(np.floor(np.log10(np.abs(n)) / 3)) * 3)
    # limit to maximal coverd symbols
    exponent= -12 if exponent < -12 else exponent
    exponent = 9 if exponent >9 else exponent

    n_new=n / 10 ** exponent
    symbol=dicct[exponent]
    return np.round(n_new,2),symbol

def convert_to_int(a):

    '''
    converts a string to integer. this function is necessary to account for - signe at the beggining in weird formatting
    :param a: string representing an integer
    :return:
    '''

    try:
        n=int(a[0])  # checks if minus signe is present
    except ValueError:
        n = -int(a[1:])
        return n
    return n




def createFolder(directory):
    '''
    function to create directories if they dont already exist
    '''
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)




def produce_index_array(u):
    u_indices=np.zeros(np.shape(u))
    for i,j in zip(range(np.shape(u)[0]),range(np.shape(u)[1])):
        u_indices[i,j]=i+j
    return u_indices

def ndargmin(array):
    return np.unravel_index(np.nanargmin(array), array.shape)



def createFolder(directory):
    '''
    function to create directories, if they dont already exist
    '''
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)



#function to convert none object from re.search to empty string

def get_group(s=None,group_number=1):

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

    if isinstance(group_number,str): # all groups
        if group_number=="all":
            return s.groups()
   
    if isinstance(group_number,list): # multiple group ids
        groups=[]
        for gn in group_number:
            try:
                groups.append(s.group(gn))
            except IndexError: # append empty string if group not found
                groups.append('') 
        return groups

    try: # the particular group was not found
        group=s.group(group_number)
    except IndexError:
        return ''

    return s.group(group_number) # single group



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
    if not (isinstance(li,list) or isinstance(li,np.ndarray)):
        return li,""
    elif len(li)==1:
        return li[0],""
    else:
        return li[0],li[1]





