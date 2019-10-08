# contains the default parameters, parameters for plotting, messages that are printed while the programming is executed
# and tooltips for the tfm addon

from collections import defaultdict
from matplotlib import cm  # making list from colormap
import numpy as np




# default parameters for the analysis
default_parameters={
    "sigma":0.49, # poisson ratio
    "young":49000, # young's modulus
    "pixelsize":0.201, # pixel size of the image with beads in  µm/pixel
    "window_size":20,  # window size for particle image velocimetry in µm
    "overlapp":19, # overlap  size for particle image velocimetry in µm. This should be at least window_size/2.
    "std_factor":15,  # additional filter for extreme values in deformation field
    "h":300, # height of the substrate in µm
    "edge_padding":0.1, # fraction of the image close to the borders that is ignored for any analyzed value
    "TFM_mode":"finite_thickness",  # mode of traction force microscopy ("finite_thickness" or "infinite_thcikness")
    "FEM_mode":"colony",  # mode for FEM type. Either perform FEM on a single colony (mode: "colony") or on the whole
                        # filed of view (mode: "cell layer"). In "cell layer you select two areas and calculate stress and
                        # contractile energy on them. In "colony" you select the area of a colony and draw cell borders. These
                        # borders are used to analyze stresses along these borders.
    "area masks":["cell type1","cell type2","membrane"] # use these mask to evaluate quantities only on this area.
                        # (mostly use full for sum of deformation...)
}

# dictionary setting a sting label for some calcualtions on mask. (summing deormations, finding the areas,
# averaging stresses
mask_label_dict={"cell type1":"cell type 1",
"cell type2":"cell type 2",
"membrane": "colony",
"contractillity_colony":"colony"
                }

# adding to the parameters dict
default_parameters["mask_labels"]=mask_label_dict




# setting all parameters for plotting. Parameters can be set for each plotting instance separately by setting the value
# in the first dictionary as another dictionary. The keys of this dictionary identify the plot they are applied to.
# "deformation": plotting the deformation field, "traction": plotting the traction field,"FEM" plotting the line
# stresses on cell-cell borders,"FEM_cell_layer" plotting the average normal stress in the field of view
# (this plot is currently not produced), "energy_points" plotting the contractile energy in the cell layer mode.
# If you don't use a dictionary the parameter is passed to all plots. If you use a dictionary and leave out a plot type
# nothing is passed to this plot.

default_fig_parameters={
    # list of which plots to generate depending on which analysis mode is chosen
    # available plots are ["deformation","traction","FEM_borders","stress_map","energy_points"]
    "plots":{"cell layer":["deformation","traction","energy_points","stress_map"],
            "colony":["deformation","traction","FEM_borders","stress_map"]},
    # dictionary specifying the name of the layer (in the cdb database) a lot is written to // you dont need to change this
    "plots_layers":{"deformation":"deformation","traction":"traction","FEM_borders":"FEM_borders"
        ,"stress_map":"stress_map","energy_points":"energy_points"}, #


    "cbar_str": {"deformation":"deformation\n[pixel]","traction":"traction\n[Pa]","FEM_borders":"line stress\n[N/µm]",
                 "stress_map":"avg. normal stress\nin N/m","energy_points":"contractile energy\nJ/pixel\n"
                 },  # label of the color bar
    "cmap": "rainbow",  # colormap for displaying magnitudes in deformation and traction fields
    "vmin": {"deformation":None,"traction":None,"FEM_borders":None,"stress_map":None,"energy_points":None},  # minimal value displayed in the colormap
    "vmax": {"deformation":None,"traction":None,"FEM_borders":None,"stress_map":None,"energy_points":None},  # maximal value displayed in the colormap
    "cbar_width": "2%",  # width of the color bar in % of the main image. Must be string with % at the end.
    "cbar_height": "50%",  # height of the color bar in % of the main image. Must be string with % at the end.
    "cbar_borderpad": 4,  # distance between the edge of the image and the color bar (in pixels???)
    "scale_ratio": 0.2,  # scale arrows so that the longest arrow is "maximum image dimension" * "scale ratio" long
    "headwidth": 3,  # width of the arrow heads (in pixels?)
    "headlength": 3,  # length of the arrow heads (in pixels?)
    "width": 0.002,  # width of the arrow shaft (what unit?)
    "plot_t_vecs":{"FEM_borders":False}, # plotting the stress vectors on the cell border stresses image
    "plot_n_arrows":{"FEM_borders":False}, # plotting normal vectors on the cell border stresses image
    "linewidth":{"FEM_borders":4}, # line width when plotting the cell border stresses
    "cm_cmap":{"FEM_borders":cm.jet}, # color map for plotting the cell border stresses. Needs a color maps object.
    "border_arrow_filter":{"FEM":1}, # plot only every n'th arrow for on the cell border stresses image
    #"filter":[0,4],
    #"figsize":(10,10),
    "file_names":{"deformation":"deformation.png","traction":"traction.png"   # filenames under wich plots are saved
        ,"FEM_borders":"border_stress_img.png","stress_map":"avg_normal.png","energy_points":"energy_distribution.png"}
}


# default parameters for plotting
def set_fig_parameters(shape, fig_shape,dpi, default_fig_parameters,figtype):
    fig_parameters = {
        # filtering: 1.minimal length of arrow, 2. draw only every n'th arrow (in x and y direction)
        "filter": [0, int(np.ceil(shape[0] / 50))],
        # figsize, so that saving the figure with dpi=dpi, gives an image of the shape fig_shape[0]
        # used to match the other images in the database
        "figsize": (fig_shape[1] / dpi, fig_shape[0] / dpi),
    }
    for key,value in default_fig_parameters.items(): #adding, or potentially overwriting other paramters
        if isinstance(value,dict): # check if values for multiple plot types exist
            if figtype in value.keys():
                fig_parameters[key] = value[figtype]
        else:
            fig_parameters[key]=value


    return fig_parameters




tooltips=defaultdict(lambda:"")
tooltips["check_box_def"]="Add the calculation of a deformation field"
tooltips["check_box_tra"]="Add the calculation of a traction force deformation field"
tooltips["check_box_FEM"]="Add stress analysis on the colony area. You need  to mark them membranes inside the cell colony."
tooltips["check_box_contract"]="Add the analysis of contractillity and contractile energy. You need to mark the area to be used for thsi calcualtion seperately"
tooltips["colony_type"]=""
tooltips["button_start"]="Start the calculation"
tooltips["apply to"]="Run the selected analysis on the current frame or on all frames in the database. THe later option will save the reuslts to a textfile"
tooltips["use_h_correction"]="Enable the use of height corrected tracktion force microscopy"
tooltips["young"]="Set the youngs modulus of the substrate that the cells are growing on"
tooltips["sigma"]="Set the possion ratio of the substrate that the cells are growing on"
tooltips["pixelsize"]="Set pixel size of your images"
tooltips["overlapp"]="Set the overlapp for calcualting the deformation fiel by PIV. THis should be about 95% of the windowsize if you include" \
            "the stress analysis. "
tooltips["window_size"]="Set the windowsize for calcualting the deormation field by PIV. A could guess is 7 times the radius of a bead."
tooltips["h"]="set the height of the substrate your cells are growing on"



# some message to be printed
calculation_messages=defaultdict(lambda:"%s")
calculation_messages["deformation"]="calculating deformation on frames %s"
calculation_messages["traction_force"]="calculating traction_forces on frames %s"
calculation_messages["FEM_full_analysis"]="FEM_analysis on frames %s"
calculation_messages["get_contractillity_contractile_energy"]="contractility/contractile energy on frames %s"
calculation_messages["general_properties"]="getting colony properties on frames %s"


# units of the returned stress and energy measures:
units=defaultdict(lambda: "")
units["avarage line stress"]="N/m"
units["avarage cell force"]="N"
units["avarage cell pressure"]="N/m"
units["avarage cell shear"]="N/m"
units["std line stress"]="N/m"
units["std cell force"]="N"
units["std cell pressure"]="N/m"
units["std cell shear"]="N/m"
units["contractillity on cell colony"]="N"
units["contractile energy on cell colony"]="J"
units["contractillity on cell type 1"]="N"
units["contractile energy on cell type 1"]="J"
units["contractillity on cell type 2"]="N"
units["contractile energy cell type 2"]="J"
units["avarage normal stress colony"]="N/m"
units["avarage shear stress colony"]="N/m"
units["avarage normal stress on cell type 1"]="N/m"
units["avarage shear stress on cell type 1"]="N/m"
units["avarage normal stress on cell type 2"]="N/m"
units["avarage shear stress on cell type 2"]="N/m"
units["area"]="m2"
units["area of colony"]="m2"
units["colony n_cells"]=""
units["sum deformations"]="pixels"
units["sum traction forces"]="N/m2"
units["sum deformations on cell type 1 "]="pixels"
units["sum traction forces on cell type 1"]="N/m2"
units["sum deformations on cell type 2"]="pixels"
units["sum traction forces on cell type 2"]="N/m2"
units["sum deformations on cell colony"]="pixels"

