
from __future__ import division, print_function
from andreas_TFM_package.TFM_functions_for_clickpoints import *  # must be on top because of some matplotlib backend issues
#from TFM_functions_for_clickpoints import * local import

#from utilities import  get_group,createFolder
import os
from qtpy import QtCore, QtGui, QtWidgets
from clickpoints.includes.QtShortCuts import AddQLineEdit, AddQComboBox, QInputNumber
from clickpoints.includes import QtShortCuts
import qtawesome as qta
import clickpoints
import asyncio
import threading

import warnings
warnings.simplefilter(action='ignore', category=RuntimeWarning)




def add_parameter_from_list(labels, dict_keys, default_parameters, layout,grid_line_start,con_func):
    params = {}
    param_labels = {}
    for i, (label,dict_key) in enumerate(zip(labels,dict_keys)):
        p = QtWidgets.QLineEdit()  # box to enter number
        p.setFixedWidth(120)
        p.setText(str(default_parameters[dict_key]))
        p.textChanged.connect(con_func)
        p_l = QtWidgets.QLabel()  # text
        p_l.setText(label)
        layout.addWidget(p, grid_line_start + i, 1)  # adding to layout
        layout.addWidget(p_l, grid_line_start + i, 0)
        params[dict_key] = p
        param_labels[dict_key] = p_l
    last_line = grid_line_start + i
    return params, param_labels, last_line

def read_all_paramters(parameter_widgets,parameter_dict):
    for p_name,p_widget in parameter_widgets.items():
        if isinstance(p_widget,QtWidgets.QLineEdit):
            parameter_dict[p_name]=float(p_widget.text())
        if isinstance(p_widget, QtWidgets.QComboBox):
            parameter_dict[p_name] = p_widget.currentText()
    return parameter_dict



class Addon(clickpoints.Addon):

    def __init__(self, *args, **kwargs):
        clickpoints.Addon.__init__(self, *args, **kwargs)

        self.outfile_path=os.path.join(os.path.split(self.db._database_filename)[0],"out.txt") # path for output text file
        self.frame_number=self.db.getImageCount()
        self.db_info, self.all_frames = get_db_info_for_analysis(self.db) # information about the path, image dimensions
        self.res_dict=defaultdict(dict) # dictionary that catches all results
        # conncetion between frames in the data base and numbering of images in their filename, list of all frames


        """ GUI Widgets"""
        # set the title and layout
        self.setWindowTitle("TFM")
        self.setWindowIcon(qta.icon("fa.compress"))
        self.setMinimumWidth(400)
        self.setMinimumHeight(200)


        self.layout = QtWidgets.QGridLayout(self)

        # button to start calculation
        self.button_start = QtWidgets.QPushButton("start")
        self.button_start.clicked.connect(self.start_thread)
        self.layout.addWidget(self.button_start,0,0)

         # check_boxes
        self.check_box_def =QtWidgets.QCheckBox("deformation")
        self.check_box_tra = QtWidgets.QCheckBox("tracktion forces")
        self.check_box_FEM = QtWidgets.QCheckBox("FEM analysis")
        self.check_box_contract = QtWidgets.QCheckBox("contractility_measures")
        self.layout.addWidget(self.check_box_def, 0, 1)
        self.layout.addWidget(self.check_box_tra, 1, 1)
        self.layout.addWidget(self.check_box_FEM, 2, 1)
        self.layout.addWidget(self.check_box_contract, 3, 1)


        #

        # # choosing single or all frames
        self.analysis_mode=QtWidgets.QComboBox()
        self.analysis_mode.addItems(["current frame", "all frames"])
        self.layout.addWidget(self.analysis_mode, 4, 1)
        self.analysis_mode_descript = QtWidgets.QLabel()
        self.analysis_mode_descript.setText("analysis mode")
        self.layout.addWidget(self.analysis_mode_descript, 4, 0)



        #### parameters
        self.parameter_dict = default_parameters
        self.parameter_list=["youngs modulus [Pa]","possion ratio","pixel size [µm]","piv overlapp [µm]","piv window size [µm]","gel hight [µm]"]
        self.param_dict_keys=["young","sigma","pixelsize","overlapp","window_size","h"]
        self.parameter_widgets,self.parameter_lables,last_line=add_parameter_from_list(self.parameter_list,
                                                            self.param_dict_keys,self.parameter_dict,self.layout
                                                                                       ,5,self.parameters_changed)
        # drop down for choosing wether to use height correction
        self.use_h_correction = QtWidgets.QComboBox()
        self.use_h_correction.addItems(["finite_thikness", "infinite_thikness"])
        self.use_h_correction.currentTextChanged.connect(self.parameters_changed) # update self.paramter_dict everytime smethong changed
        self.use_h_correction_descr = QtWidgets.QLabel()
        self.use_h_correction_descr.setText("enable height correction")
        self.layout.addWidget(self.use_h_correction,  last_line+1, 1)
        self.layout.addWidget(self.use_h_correction_descr,   last_line + 1, 0)
        self.parameter_widgets["TFM_mode"]=self.use_h_correction # adding to parameters dict
        self.parameter_lables["TFM_mode"] =  self.use_h_correction_descr  # adding to parameters dict
        # adding to layout
        self.setLayout(self.layout)
        self.parameters_changed() # initialize parameters dict

    # reading paramteres and updating the dictionary
    def parameters_changed(self):
        self.parameter_dict = read_all_paramters(self.parameter_widgets,self.parameter_dict)

       
        #write_output_file(self.parameter_dict, None, "parameters", self.outfile_path) # writes paramters to output file..

    # decorator functions to handle diffrent outputs and writing to text file
    def calculate_general_properties(self,frames):
        apply_to_frames(self.db, self.parameter_dict, analysis_function=general_properties, res_dict=self.res_dict,
                        frames=frames, db_info=self.db_info)  # calculation of colony area, number of cells in colony
    def calculate_deformation(self,frames):
            apply_to_frames(self.db, self.parameter_dict,analysis_function=deformation,res_dict=self.res_dict,
                                frames=frames,db_info=self.db_info)  # calculation of deformations

    def calculate_tracktion(self,frames):
            apply_to_frames(self.db, self.parameter_dict,analysis_function=tracktion_force,res_dict=self.res_dict,
                                frames=frames,db_info=self.db_info)  # calculation of tracktion forces

    def calculate_FEM_analysis(self,frames):
            apply_to_frames(self.db, self.parameter_dict, analysis_function=FEM_analysis,res_dict=self.res_dict,
                            frames=frames,db_info=self.db_info)  # calculation of various stress measures

    def calculate_contractile_measures(self,frames):
            apply_to_frames(self.db, self.parameter_dict,analysis_function=get_contractility_contractile_energy,
                 res_dict=self.res_dict,frames=frames,db_info=self.db_info)  # calculation of contractility and contractile energy




    def start(self): # perform all checked calculations
        cdb_frame = self.cp.getCurrentFrame()
        self.frame = get_frame_from_annotation(self.db, cdb_frame) # current frame


        print("parameters:\n",self.parameter_dict)
        self.mode=self.analysis_mode.currentText() # only current frame or all frames
        if self.mode == "current frame": # only current frame
            frames=self.frame
            print("analyzing current frame = ", frames)
        if self.mode == "all frames": # all frames
            frames=self.all_frames
            print("analyzing frames = ", frames)
            write_output_file(self.parameter_dict, "parameters", self.outfile_path)

        self.calculate_general_properties(frames)
        if self.check_box_def.isChecked():     
            self.calculate_deformation(frames)
        if self.check_box_tra.isChecked():
            self.calculate_tracktion(frames)
        if self.check_box_FEM.isChecked():
            self.calculate_FEM_analysis(frames)
        if self.check_box_contract.isChecked():
            self.calculate_contractile_measures(frames)


        write_output_file(self.res_dict, "results", self.outfile_path)  # writing to output file
        # reloding all images that where changed
        #if self.mode=="current frame":
        #        self.cp.reloadImage(frame_index=cdb_frame)
        #else:
        #        for frame_index in range(self.frame_number):
        #                self.cp.reloadImage(frame_index=frame_index )
        
        print("calculation complete")


    def start_thread(self): ## run in a sepearte thread to keep clickpoints gui responsive (ask richie about this)
        x = threading.Thread(target=self.start)
        x.start()
        #x.join()
    def buttonPressedEvent(self):
        # show the addon window when the button in ClickPoints is pressed
        self.show()
    async def run(self): # disable running in other thread other wise "figsize" is somehow overloaded
        pass




