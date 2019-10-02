
from __future__ import division, print_function
from andreas_TFM_package.TFM_functions_for_clickpoints import *  # must be on top because of some matplotlib backend issues
from andreas_TFM_package.parameters_and_strings import tooltips,default_parameters
#from TFM_functions_for_clickpoints import * local import

#from utilities import  get_group,createFolder
import os
from qtpy import QtCore, QtGui, QtWidgets
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
        p.setToolTip(tooltips[dict_key])
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
        self.parameter_dict = default_parameters # loading default parameters
        # guessing the curent analysis mode
        self.parameter_dict["FEM_mode"],undetermined=guess_TFM_mode(self.db_info,self.parameter_dict)
        self.colony_type_old_id =int(self.parameter_dict["FEM_mode"] == "cell layer")

        # connection between frames in the data base and numbering of images in their filename, list of all frames


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
        self.button_start.setToolTip(tooltips["button_start"])
        self.layout.addWidget(self.button_start,0,0)


        # choosing type of cell system

        self.colony_type = QtWidgets.QComboBox()
        self.colony_type.addItems(["colony", "cell layer"])
        self.colony_type.setToolTip(tooltips["colony_type"])
        self.colony_type.setCurrentIndex(self.colony_type_old_id)
        self.colony_type.activated.connect(self.switch_colony_type_mode)  # activated when user changes the the selected text




        self.sub_layout=QtWidgets.QHBoxLayout()
        self.sub_layout.addWidget(self.colony_type)
        self.sub_layout.addStretch()
        self.layout.addLayout(self.sub_layout, 2, 0)

        # filling areas for cell patches
        self.fill_patches_button = QtWidgets.QPushButton("fill cell area")
        self.fill_patches_button.setToolTip(tooltips["fill_patches_button"])
        self.fill_patches_button.clicked.connect(self.fill_patches)

        self.sub_layout2 = QtWidgets.QHBoxLayout()
        self.sub_layout2.addWidget(self.fill_patches_button)
        self.sub_layout2.addStretch()
        self.layout.addLayout(self.sub_layout2, 3, 0)
        self.fill_patches_button.setVisible(self.colony_type_old_id) # ide button for now

        if not undetermined:
            self.switch_colony_type_mode(first=True) # initializing masks according to the first value of "FEM_mode"


         # check_boxes
        self.check_box_def =QtWidgets.QCheckBox("deformation")
        self.check_box_tra = QtWidgets.QCheckBox("traction forces")
        self.check_box_FEM = QtWidgets.QCheckBox("FEM analysis")
        self.check_box_contract = QtWidgets.QCheckBox("contractillity_measures")

        self.check_box_def.setToolTip(tooltips["check_box_def"])
        self.check_box_tra.setToolTip(tooltips["check_box_tra"])
        self.check_box_FEM.setToolTip(tooltips["check_box_FEM"])
        self.check_box_contract.setToolTip(tooltips["check_box_contract"])

        self.layout.addWidget(self.check_box_def, 0, 1)
        self.layout.addWidget(self.check_box_tra, 1, 1)
        self.layout.addWidget(self.check_box_FEM, 2, 1)
        self.layout.addWidget(self.check_box_contract, 3, 1)


        # # choosing single or all frames
        self.analysis_mode=QtWidgets.QComboBox()
        self.analysis_mode.addItems(["current frame", "all frames"])
        self.analysis_mode.setToolTip(tooltips["apply to"])
        self.layout.addWidget(self.analysis_mode, 4, 1)
        self.analysis_mode_descript = QtWidgets.QLabel()
        self.analysis_mode_descript.setText("apply to")
        self.layout.addWidget(self.analysis_mode_descript, 4, 0)



        #### parameters
        self.parameter_labels=["youngs modulus [Pa]","possion ratio","pixel size [µm]","piv overlapp [µm]","piv window size [µm]","gel hight [µm]"]
        self.param_dict_keys=["young","sigma","pixelsize","overlapp","window_size","h"]
        self.parameter_widgets,self.parameter_lables,last_line=add_parameter_from_list(self.parameter_labels,
                                                            self.param_dict_keys,self.parameter_dict,self.layout
                                                                                       ,5,self.parameters_changed)
        # drop down for choosing wether to use height correction
        self.use_h_correction = QtWidgets.QComboBox()
        self.use_h_correction.addItems(["finite_thickness", "infinite_thickness"])
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

    def calculate_traction(self,frames):
            apply_to_frames(self.db, self.parameter_dict,analysis_function=traction_force,res_dict=self.res_dict,
                                frames=frames,db_info=self.db_info)  # calculation of traction forces

    def calculate_FEM_analysis(self,frames):
            apply_to_frames(self.db, self.parameter_dict, analysis_function=FEM_full_analysis,res_dict=self.res_dict,
                            frames=frames,db_info=self.db_info)  # calculation of various stress measures

    def calculate_contractile_measures(self,frames):
            apply_to_frames(self.db, self.parameter_dict,analysis_function=get_contractillity_contractile_energy,
                 res_dict=self.res_dict,frames=frames,db_info=self.db_info)  # calculation of contractility and contractile energy

    # switching the type of cell colony
    def switch_colony_type_mode(self,first=False):
        if not first:
            choice = QtWidgets.QMessageBox.question(self, 'continue',
                                                "This will delete all previous mask. Do you want to coninue?",
                                                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        else:
            choice=QtWidgets.QMessageBox.Yes

        if choice == QtWidgets.QMessageBox.Yes:
            # new mask types
            self.parameter_dict["FEM_mode"] =self.colony_type.currentText()
            self.colony_type_old_id=self.colony_type.currentIndex()
            if self.colony_type.currentText()=="cell layer":
                setup_masks(self.db,self.db_info, self.parameter_dict) # deletes and add masktyes
                self.fill_patches_button.setVisible(True)
            if self.colony_type.currentText()=="colony":
                setup_masks(self.db,self.db_info, self.parameter_dict)
                self.fill_patches_button.setVisible(False)
            self.cp.reloadMaskTypes()  # reloading mask to display in clickpoints window
        else:
            self.colony_type.setCurrentIndex(self.colony_type_old_id) #reverting to old index


    def fill_patches(self):


        apply_to_frames(self.db, self.parameter_dict, analysis_function=fill_patches_for_cell_layer,
                        res_dict=self.res_dict, frames=self.all_frames, db_info=self.db_info)

        self.cp.reloadMask()
        self.cp.save()


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


        if self.check_box_def.isChecked():     
            self.calculate_deformation(frames)
        if self.check_box_tra.isChecked():
            self.calculate_traction(frames)
        self.calculate_general_properties(frames)
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




