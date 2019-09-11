
from __future__ import division, print_function
import sys
sys.path.insert(0,"/media/user/GINA1-BK/Andreas-Python/tracktion_force_microscopy")
from TFM_functions_for_clickpoints import *  # must be on top because of some matplotlib backend issues
#from utilities import  get_group,createFolder
import numpy as np

import os
import json
from qtpy import QtCore, QtGui, QtWidgets
from clickpoints.includes.QtShortCuts import AddQLineEdit, AddQComboBox, QInputNumber
from clickpoints.includes import QtShortCuts
import qtawesome as qta
import re
import clickpoints

#folder="/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/out_clickpoints"
#db=clickpoints.DataFile(os.path.join(folder,"database3.cdb"),"r") # genarating new clickpoints file

# calculating deformation with piv


# define default dicts, enables .access on dicts
class dotdict(dict):
    def __getattr__(self, attr):
        if attr.startswith('__'):
            raise AttributeError
        return self.get(attr, None)

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


class Addon(clickpoints.Addon):
    finished = QtCore.Signal()

    def __init__(self, *args, **kwargs):
        clickpoints.Addon.__init__(self, *args, **kwargs)

        """ get or set options """
        self.layer_dict = dict((i, l.name) for i, l in enumerate(list(self.db.getLayers())))
        """ get or set options """
        # layers = ([l.name for l in self.db.getLayers()])

        #self.addOption(key="segmentation_th", display_name="sort_layers", default=False, value_type="bool", tooltip="")

        """ GUI Widgets"""
        # set the title and layout
        self.setWindowTitle("TFM")
        self.setWindowIcon(qta.icon("fa.compress"))
        self.setMinimumWidth(400)
        self.setMinimumHeight(200)
        self.layout = QtWidgets.QVBoxLayout(self)

        # buttons
        self.button_def = QtWidgets.QPushButton("deformation")
        self.button_tra = QtWidgets.QPushButton("tracktion forces")
        self.button_FEM = QtWidgets.QPushButton("FEM analysis")
        self.button_def.clicked.connect(self.def_button)
        self.button_tra.clicked.connect(self.button_tra)
        self.button_FEM.clicked.connect(self.button_FEM)
        self.layout.addWidget(self.button_def)
        self.layout.addWidget(self.button_tra)
        self.layout.addWidget(self.button_FEM)

        self.box=QtWidgets.QComboBox("box")
        self.box.addItems(["Marker Count", "Marker Positions", "Track Positions"])
        self.layout.addWidget(self.box)
        self.combo_style = AddQComboBox(self.layout, "Mode", values=["Marker Count", "Marker Positions", "Track Positions"])


    def def_button(self):
       pass

    def button_tra(self):
        pass

    def button_FEM(self):
        pass



    def processAll(self, save_properties=True):
        """
        Iterates over all available images in the current CP project,
        perform segmentation, extract region properties.
        """

        q_images = self.db.getImages()

        results = []
        for q_img in q_images:
            self.q_img = q_img
            print("Batch processing Image Nr %d" % q_img.sort_index)
            self.updateSegmentation(qimg=q_img)



    """ DEFAULT ADDON FUNCTIONS """
    def buttonPressedEvent(self):
        # show the addon window when the button in ClickPoints is pressed
        self.show()


'''

sigma=0.49#poison ratio
young=25536# youngsmodulus
pixelsize = 4.09 / 40 #µm/pixel
window_size=100   # in pixel , keep at 10 µm
overlapp=60   # window_size/2
std_factor=15
pixel_factor = window_size - overlapp
h=300

parameter_dict=make_paramters_dict_tfm( sigma=sigma, young=young, pixelsize1=pixelsize, window_size=window_size, overlapp=overlapp,
                         std_factor=std_factor, h=h, pixel_factor = window_size - overlapp)
# calculating the deformation field and adding to data base


apply_to_all_frames(db,parameter_dict,analysis_function=deformation) # calculation of deformations on all frames
apply_to_all_frames(db,parameter_dict,analysis_function=tracktion_force) # calculation of deformations on all frames


sigma=0.49#poison ratio
pixelsize = 4.09 / 40 #µm/pixel
parameter_dict = {"sigma":0.49,"pixelsize":4.09 / 40}
folder="/media/user/GINA1-BK/data_traktion_force_microscopy/TFM_cell_patches/WT2shift/out_clickpoints"

apply_to_all_frames(db,parameter_dict,analysis_function=FEM_analysis) # calculation of deformations on all frames
'''