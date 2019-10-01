


## generating gifs of shift corrected "before" and "after" images.
# the programm will search trough all sub folders in the given folder,
#output is saved in the "folder" as single gifs

# beware: sometimes files on usb drive are saved as hidden. check  with ls commadn

# data structure:
#....
import os
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt
from skimage.feature import register_translation
from scipy.ndimage.interpolation import shift
import os
import re
from correcting_frame_shift_gifs import *
from correcting_frame_shift import *

import imageio


parser = ArgumentParser(description='script for shifting two flourescent images of beads\n'
                                    'this script will shift the after image and\n'
                                    'cut the after and the brightfield before image.\n'
                                    'you can also produce gifs.\n'
                                    'output is saved in the -dest folder.\n '
                        "example usage:\n"
                        "python '/media/user/GINA1-BK/Andreas-Python/tracktion_force_microscopy/correcting_frame_shift_single_pictures_shell.py'\n"
                                    " -dest /home/user/Desktop/20190216/SNU\n"
                                    " -mode 'gifs'")


parser.add_argument("-dest",  dest="folder",default=os.getcwd(),help=("folder, in which the script will search images.\n"
                                                                        "Put the path in quotation marks"))
parser.add_argument("-mode",  dest="mode",default="images",help=("chose weather gifs or single images should be produced.\n"
                                                                 "Possible values are 'images' and 'gifs'.\n"
                                                                 "Use quotation marks"))
parser.add_argument("-selector",  dest="selector",default="",help=("add a negative selector.\n"
                                                                        "The programm will not search in a folder containing this string.\n"
                                                                   "Use quotation marks"))

args = parser.parse_args()
folder=args.folder
selector=args.selector
mode=args.mode


if mode=="images":
    cut_images(folder, selector)
if mode=="gifs":
    cut_images_gifs(folder, selector)

