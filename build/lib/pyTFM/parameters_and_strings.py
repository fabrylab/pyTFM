# contains the default parameters, parameters for plotting, messages that are printed while the programming is executed
# and tooltips for the tfm addon
import re
from collections import defaultdict

import numpy as np
from pyTFM.utilities_TFM import make_iterable, squeeze_list, convert_none_str, try_float_convert

# default parameters for the analysis
default_parameters = {
    "sigma": 0.49,  # poisson ratio of the substrate
    "young": 49000,  # young's modulus
    "pixelsize": 0.201,  # pixel size of the image with beads in  µm/pixel
    "window_size": 20,  # window size for particle image velocimetry in µm
    "overlap": 18,  # overlap  size for particle image velocimetry in µm. This should be at least window_size/2.
    "std_factor": 15,  # additional filter for extreme values in deformation field
    "h": 300,  # height of the substrate in µm
    "edge_padding": 0.1,  # fraction of the image close to the borders that is ignored for any analyzed value
    "padding_cell_layer": 0.2,
    # additional region ignored for stress analysis in "cell layer" mode. Average stresses and line
    # tension is only calculated on the area that is "edge_padding"+"padding_cell_layer" away from the image edge
    "sigma_cells":0.5, # Poisson's ratio of the Cell Sheet in the Monolayer Stress Microscopy. This parameter has little
    # influence on the resulting stresses and should therefore not be changed.
    "TFM_mode": "finite_thickness",  # mode of traction force microscopy ("finite_thickness" or "infinite_thcikness")
    "FEM_mode": "colony",  # mode for FEM type. Either perform FEM on a single colony (mode: "colony") or on the whole
    # filed of view (mode: "cell layer"). In "cell layer you select two areas and calculate stress and
    # contractile energy on them. In "colony" you select the area of a colony and draw cell borders. These
    # borders are used to analyze stresses along these borders.
    "filter_size": 3,  # size (sigma) of the gauss filter applied to the traction field in µm.
    "filter_type": "gaussian",  # size (sigma) of the gauss filter applied to the traction field in µm.
    "min_obj_size": 1500,  # all objects (cell patches/ isolated cells) below this size (in pixel) will be ignored
    "min_line_length": 0,
    # small lines below this value are filtered out. This is mainly for cosmetic reasons, as small
    # lines won't contribute much to the average line tension anyway
    "cv_pad": 0,  # padding when calculating the coefficient of variation in µm// only necessary if the
    # mask for the FEM area fits very closely to the mask for membranes
    "mask_properties": {
        "cell type1": {"use": ["defo", "forces", "area_layer", "FEM_layer", "stress_layer"], "FEM_mode": ["cell layer"],
                       "color": "#1322ff", "index": 1, "label": "cell type 1", "name": "cell type1"},
        "cell type2": {"use": ["defo", "forces", "area_layer", "FEM_layer", "stress_layer"], "FEM_mode": ["cell layer"],
                       "color": "#ebff05", "index": 3, "label": "cell type 2", "name": "cell type2"},
        "Cell Boundary": {"use": ["area_colony", "borders", "FEM_layer", "stress_colony",], # "forces"
                     "FEM_mode": ["cell layer", "colony"], "color": "#30ff0c", "index": 2, "label": "Cell Area",
                     "name": "Cell Boundary"},
        "Tractions": {"use": ["defo", "forces", "FEM_colony"], "FEM_mode": ["colony"], "color": "#ff0b23",
                           "index": 1,
                           "label": " Traction Area", "name": "Tractions"}}
}
# "FEM_area": {"use": ["FEM_colony"], "FEM_mode": ["colony"], "color": "#30ff0c", "index": 2, "label": " of FEM area",
#                     "name": "FEM_area"}}

#force measures

    # use: calculations this mask is used in
# FEM_mode: mode the mask is used in
# color: color of the mask in clickpoints
# index: index of the mask in clickpoints; FEM_mode and index mst be a unique pair.
# name: name of the mask in the clickpoints database
# label: string used for the output file


# plotting parameters
# available plots are ["deformation","traction","FEM_borders","stress_map","energy_points"]
default_fig_parameters = {
    "cmap": "rainbow",  # colormap for displaying magnitudes in deformation and traction fields
    "vmin": None,  # minimal value displayed in the colormap
    "vmax": None,  # maximal value displayed in the colormap
    "cbar_width": "2%",  # width of the color bar in % of the main image. Must be string with % at the end.
    "cbar_height": "50%",  # height of the color bar in % of the main image. Must be string with % at the end.
    "cbar_borderpad": 0.2,  # distance between the edge of the image and the color bar
    "scale_ratio": 0.2,  # scale arrows so that the longest arrow is "maximum image dimension" * "scale ratio" long
    "cbar_title_pad": 10,  # padding of the
    "headwidth": 3,  # width of the arrow heads (in pixels?)
    "headlength": 3,  # length of the arrow heads (in pixels?)
    "headaxislength": None,
    "width": 0.002,  # width of the arrow shaft (what unit?)
    "plot_t_vecs": False,  # plotting the stress vectors on the cell border stresses image
    "plot_n_arrows": False,  # plotting normal vectors on the cell border stresses image
    "linewidth": 4,  # line width when plotting the cell border stresses
    "border_arrow_filter": 1,  # plot only every n'th arrow for on the cell border stresses image
    "outer_cb_color": "grey",  # color of the line plotting the cell colony edge
    "outer_cb_style": "-",  # linestyle of the line plotting the cell colony edge. "--" for dashed line
    "boundary_resolution": 6,  # resolution when plotting the line tension. Highest is 1. Increase for lower resolution,
    # label of the color bar,
    "cbar_style": "clickpoints",  # if "clickpoints" the color bar is plottetd inside of the figure
    "plot_style": "clickpoints",
    "filter_factor": 1,  # this factor defines how many arrows are shown in deformation and traction images.
    # low number results in  many arrows, high number results in few arrows
    "background_color": "#330033",
    # set a color for background values. "cmap_0" fill read the zero color of the colormap. "white" would make the background white...
    # this doesn't effect images of deformation and traction
    "cbar_tick_label_size": 15,  # size of the tick labels on the color bar
    "cbar_axes_fraction": 0.2,
    # fraction of the axes in horrizontal direction, that the colorbar takes up, when colorbar is plotted outside
    # of the graph
    "plot_cbar": True,
    "cbar_str": {"deformation": "deformation\n[pixel]", "traction": "traction\n[Pa]",
                 "FEM_borders": "line tension\n[N/m]",
                 "stress_map": "avg. normal stress\nin N/m", "energy_points": "strain energy\nJ/pixel\n"
                 },
    "resolution": 200,  # dpi when saving plots
    "file_names": {"deformation": "deformation.png", "traction": "traction.png"  # filenames under wich plots are saved
        , "FEM_borders": "border_stress.png", "stress_map": "mean_normal_stress.png",
                   "energy_points": "energy_distribution.png"},
    # defining which plots are produced
    "plots": {"cell layer": ["deformation", "traction", "FEM_borders", "energy_points", "stress_map"],
              "colony": ["deformation", "traction", "FEM_borders", "stress_map"]},

    # dictionary specifying the name of the layer (in the cdb database) a plot is written to
    "plots_layers": {"deformation": "deformation", "traction": "traction", "FEM_borders": "FEM_borders"
        , "stress_map": "stress_map", "energy_points": "energy_points"},  #

}

# transform every entry with only one value to
for key, value in default_fig_parameters.items():
    if not isinstance(value, (dict, defaultdict)):
        default_fig_parameters[key] = defaultdict(lambda v=value: v)  # need special binding to reference to the value

default_parameters["fig_parameters"] = default_fig_parameters


# default parameters for plotting
def set_fig_parameters(shape, fig_shape, fig_parameters, figtype):
    # extracting value for this specific fig type
    fp = {}
    for key, value in fig_parameters.items():
        if isinstance(value, (defaultdict, dict)):  # read default value
            try:
                fp[key] = value[figtype]
            except:
                pass
        # in case someone passed it as a dict; will cause error if figtype is missing as a key
        if isinstance(value, dict):
            if key in value.keys():
                fp[key] = value[figtype]
        if not isinstance(value, (dict, defaultdict)):
            fp[key] = value
    # filtering: 1.minimal length of arrow, 2. draw only every n'th arrow (in x and y direction)
    fp["filter"] = [0, int(int(np.ceil(shape[0] / 50)) * fp["filter_factor"])]
    # figsize, so that saving the figure with dpi=dpi, gives an image of the shape fig_shape[0]
    # used to match the other images in the database
    fp["figsize"] = (fig_shape[1] / fp["resolution"], fig_shape[0] / fp["resolution"])
    return fp


def get_masks_by_key(default_parameters, key, prop, return_key=True):
    '''extracting lists of mask with certain properties'''
    if return_key:
        return [m for m, mdict in default_parameters["mask_properties"].items() if prop in mdict[key]]
    else:
        keys = [m for m, mdict in default_parameters["mask_properties"].items() if prop in mdict[key]]
        return [default_parameters["mask_properties"][key]["name"] for key in keys]


def get_properties_masks(default_parameters, masks, props):
    props, masks = make_iterable(props), make_iterable(masks)
    props_list = []
    for p in props:
        props_list.append([default_parameters["mask_properties"][m][p] for m in masks])
    return squeeze_list(props_list)


tooltips = defaultdict(lambda: "")
tooltips["button_start"] = "Start the calculation"
tooltips["image selection"] = "Select images,the output folder and name of the database"
tooltips["correct drift"] = "Correct a drift between images before and after relaxation"
tooltips["check_box_def"] = "Calculation of the deformation field"
tooltips["check_box_tra"] = "Calculation of the traction force field"
tooltips["check_box_FEM"] = "Calculation of cellular stresses and cell-cell force transfer. " \
                            "The mask 'Cell Boundary' is required."
tooltips["check_box_contract"] = "Calculation of strain energy and contractility. The mask 'Tractions' is required."
tooltips["apply to"] = "Apply the selected analysis to the current field of view or all field of views in the database"
tooltips["switch mode"] = "Switch between 'cell layer' mode (you can select two areas to calculate " \
                          "stresses and forces) and 'colony mode' " \
                          "(more accurate calculation of stresses especially for small cell patches)"
tooltips["segmentation"] = "A number of tools to help you mark the cell boundary or different cell types"
tooltips["young"] = "Set the Young's modulus of the substrate"
tooltips["sigma"] = "Set the Poisson's ratio of the substrate"
tooltips["pixelsize"] = "Set pixel size of your images"
tooltips["window_size"] = "Set size of correlation windows for the calculation of the deformation field with PIV. " \
                          "A good first guess is 7 times the radius of a bead. Increase or decrease " \
                          "the window size until you obtain a good result. You can use the ClickPoints add-on " \
                          "'Measure Tool' to directly track the the displacement of individual beads."
tooltips["overlap"] = "Set the overlap between correlation windows for the calculation of the deformation " \
                      "field with PIV. This parameters controls the resolution of the deformation field and should " \
                      "be over 95% of the window size when calculating cellular stresses. However, if you are " \
                      "just testing out parameters you can save a lot of calculation time by " \
                    "reducing the overlap."
tooltips["h"] = "Set the height of the substrate. You can type in 'infinite', if the substrate is very thick. " \
                "The assumption of infinite substrate height should should generally be valid above 300 µm. " \
                "As a simple test, you can varey the height and see if the strain energy is influenced. " \
                "If this is not the case the substrate height can be assumed as infinite."

tooltips["collect_button"] = "Generate a database with the selected images."
tooltips["folder_out"] = "Select the folder where the database file, output images and files will be stored."
tooltips["folder_after"] = "Select a folder with images after substrate relaxation."
tooltips["folder_before"] = "Select a folder with images before substrate relaxation."
tooltips["folder_cells"] = "Optional: Select a folder with images of the cells or type 'none' to skip this step."
tooltips["db_name_text"] = "Set the name of your database. It will be saved once you clicked the " \
                           "'collect images' button."
tooltips["folder_out_button"] = "Click to browse folders."
tooltips["folder_after_button"] = "Click to browse folders."
tooltips["folder_before_button"] = "Click to browse folders."
tooltips["folder_cells_button"] = "Click to browse folders"
tooltips["after"] = "Provide a regular expression to select images after substrate relaxation. " \
                    "If you leave this blank all images in the selected folder are selected."
tooltips["before"] = "Provide a regular expression to select images before substrate relaxation. " \
                     "If you leave this blank all images in the selected folder are selected."
tooltips["cells"] = "Optional: Provide a regular expression to select images of the cells. If you leave this blank " \
                    "all images in the selected folder are selected. Type 'none' to skip this."
tooltips["frames"] = "Provide a regular expression to identify the field of view/frame of each image. " \
                     "The identifier can be numbers, characters or any other sign, but must be present " \
                     "in all three (or two) images. You need to enclose the identifier with brackets '()'."

tooltips["segmentation_area"] = "Mark the entire image as cell type1 or cell type2. Use the slider to select " \
                                "an appropriate threshold and press the button to apply the threshold to all " \
                                "frames. The area marked as 'Cell Boundary' will be conserved. " \
                                "This only works in 'cell layer' mode."
tooltips["segmentation_membrane"] = "Mark the Cell Boundaries. Use the slider to select an appropriate threshold " \
                                "and press the button to apply the threshold to all frames."
tooltips["mark entire image"] = "Cover the entire image with a mask. Note that a stretch at the " \
                                "image border will always be omitted for the analysis."
tooltips["select mask segmentation"] = "Select the mask type with which to cover the entire image."
tooltips["fill_mode"] = "You can either fill out the image or draw an outline around it."

# automatically introducing some line breaks:
line_break = 65
for key, text in tooltips.items():
    try:
        breaks = int(len(text) / line_break)
        for lb in range(breaks):
            beak_point = (lb + 1) * line_break
            working = True
            i = 0
            while working:
                il = beak_point - i
                ir = beak_point + i
                if ir >= len(text) or il < 0:
                    break
                # replacing blank space
                if text[il] == " ":
                    text = text[:il] + "\n" + text[il+1:]
                    break
                if text[ir] == " ":
                    text = text[:ir] + "\n" + text[ir + 1:]
                    break
                i += 1
        tooltips[key] = text
    except Exception as e:
        raise e





# some message to be printed
calculation_messages = defaultdict(lambda: "%s")
calculation_messages["deformation"] = "calculating deformation on frames %s"
calculation_messages["traction_force"] = "calculating traction_forces on frames %s"
calculation_messages["FEM_full_analysis"] = "FEM_analysis on frames %s"
calculation_messages["get_contractillity_contractile_energy"] = "contractility/strain energy on frames %s"
calculation_messages["general_properties"] = "getting colony properties on frames %s"
calculation_messages["simple_shift_correction"] = "correct frame shift on frames %s"

# units of the returned stress and energy measures:
units = defaultdict(lambda: "")
#
units["sum deformations"] = "pixels"
units["sum deformations image"] = "pixels"
units["sum deformations Traction Area"] = "pixels"
units["sum traction forces Traction Area"] = "" # this is a weird unit

# areas and cell number
units["area"] = "m2"
units["area Cell Area"] = "m2"
units["area Traction Area"] = "m2"
units["cell numbe"] = ""
units["center of object"] = ""

# traction force measures
units["contractility"] = "N"
units["strain energy"] = "J"


# stresses
units["mean normal stress Cell Area"] = "N/m"
units["max normal stress Cell Area"] = "N/m"
units["max shear stress Cell Area"] = "N/m"
units["cv mean normal stress Cell Area"] = ""
units["cv max normal stress Cell Area"] = ""
units["cv max shear stress Cell Area"] = ""

# line tension
units["average magnitude line tension"] = "N/m"
units["average normal line tension"] = "N/m"
units["average shear line tension"] = "N/m"
units["std magnitude line tension"] = ""
units["std normal line tension"] = ""
units["std shear line tension"] = ""

# derivatives of the line tension
units["average cell force"] = "N/m"
units["average cell pressure"] = "N/m"
units["average cell shear"] = "N/m"
units["std cell force"] = ""
units["std cell pressure"] = ""
units["std cell shear"] = ""




def convert_config_input(x, type):
    x = convert_none_str(x)
    if type == "background_color":
        x = try_float_convert(x)  # try convert to float
        if not isinstance(x, float):  # else check tuple
            if re.search("\(.*\)", x):
                x = x.replace("(", "")
                x = x.replace(")", "")
                x = x.replace(" ", "")
                x = tuple([float(y) for y in x.split(",")])
    return x


default_layer_names = ["images_after", "images_before", "membranes"]

default_search_keys = {"after": "\d{1,4}after", "before": "\d{1,4}before",
                    "cells": "\d{1,4}bf_before",
                    "frames": "^(\d{1,4})"}

# could make options to read external config files in the addon and in normal applications.)
