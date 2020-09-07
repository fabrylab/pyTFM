Using Config Files, Hidden Parameters and Plotting Behavior
==========================================================================

Hidden Parameters
-------------------

The most important parameters (Young's modulus, Poisson's ratio, gel height, window size and overlap for Particle
Image Velocimetry) can be set in the clickpoints addons window. However, some parameters can only
be set by using config files. Config files also allow you to set a wide range of plotting options, most notably:
Color maps, the length of arrows, color bar positioning and dimensions and maximum and minimum
displayed values for all frames. You can find a complete list of parameters in `Overview of Analysis Parameters`_
and `Overview of Plotting Parameters`_.

Using Config Files
-------------------
You can provide a config file by placing a file "config.yaml" in the same folder as your database
(.cdb) file. YAML is a standard format for configuration files. In YAML files parameters are defined by
a series of indented key words like this:

.. code-block:: yaml

    key1:
      key2:
        key3: value


You can find a config.yaml file that reproduces all default parameters in the
`example data <https://github.com/fabrylab/example_data_for_pyTFM/archive/master.zip>`__.
The config file distinguishes two classes of parameters at the first level:
The "analysis_parameters" and the "fig_parameters". "analysis_parameters" set all parameters that are used
by the program not related to plotting.
They can only be numbers or strings. For example to change the minimum object size and change the
initial value for the Poisson's ration, your config
file should contain this:

.. code-block:: yaml

    analysis_parameters: # setting parameters for the analysis
      sigma: 0.4 # changing the Poisson's ratio set at the script start up
      min_obj_size: 1500 # changing the minimal object size (cannot be changed in the addon window)


"fig_parameters" control the plotting. Most parameters can be
specified as a single value that will apply to all plots. Alternatively you can specify a second key that
defines a type of plot first and then a value. This value will only apply to the plot type you specified. All other
plot types will use a default value. Possible plot types are:

| "deformation"    - Quiver plot of the deformation field
| "traction" - Quiver plot of the traction forces
| "FEM borders" - Plot of the line tension along cell-cell borders
| "stress map"  -  Plot of the mean normal stress
| "energy_points". - Plot of the strain energy density

Let's say you want to

| A): Change the color map to "jet" (all matplotlib color maps are accepted) for all plots.
| B): Set the maximum value of the colormap only in the deformation plots to 10 pixels.
| C): Generate a map of strain energy, when analysing cell colonies.

Then your yaml file should look like this:

.. code-block:: yaml

    fig_parameters:   # setting plotting parameters

      cmap: "jet" # changing the colormap for the deformation plot

      vmin: None  # minimal value displayed in the colormap, none means an automatic value for each frame
              # this applies to all plot types if none is specified
      vmax:       # maximum value displayed in the colormap
        deformation: 10 # applies only to deformation plots. All other plots use a default value (None)

     plots:  # defining which plots are produced
         colony: # you need to specify whether this is for "colony" or "cell layer"
          - "deformation" # this needs a list of values (marked by "-")
          - "traction"
          - "FEM_borders"
          - "stress_map"
          - "energy_points"

.. hint::
    It is best to mark strings with quotation marks (""). You can use None, True or False
    without quotation marks to set None or boolean values.


Overview of Analysis Parameters
---------------------------------

<table>
    <tr>
        <td>Foo</td>
    </tr>
</table>




+---------------------+--------------------+--------------------+----------------------------------------------------------+
|Parameter            |    Default Value   |   Type             |          Meaning                                         |
+=====================+====================+====================+==========================================================+
| **Main Parameters**                                                                                                      |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| sigma               | 0.49               | int,float          | Poisson's ration of the substrate.                       |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| young               | 49000              | int,float          | Young's modulus of the substrate in Pa.                  |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| pixelsize           | 0.201              | int,float          | Pixel size of the images of the beads.                   |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| window_size         | 20                 | int,float          |Size of the windows for PIV                               |
|                     |                    |                    |                                                          |
|                     |                    |                    |(Particle Image Velocimetry) in µm.                       |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| overlap             | 19                 | int,float          | Size of the overlap for PIV in µm.                       |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| FEM_mode            | "colony"           | string             | Analyzing colonies or cell layer. This changes the       |
|                     |                    |                    |                                                          |
|                     |                    |                    | behavior, concerning which masks are used,               |
|                     |                    |                    |                                                          |
|                     |                    |                    | which plots are generated and what area                  |
|                     |                    |                    |                                                          |
|                     |                    |                    | is used for stress measurements.                         |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| **Hidden Parameters**                                                                                                    |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| std_factor          | 15                 | int,float          | Additional filter for the deformation field.             |
|                     |                    |                    |                                                          |
|                     |                    |                    | Deformations greater then                                |
|                     |                    |                    | :math:`\mu+\sigma \times 15`                             |
|                     |                    |                    |                                                          |
|                     |                    |                    | (:math:`µ` and :math:`\sigma`:                           |
|                     |                    |                    | mean  and  standard deviation of the norm of             |
|                     |                    |                    |                                                          |
|                     |                    |                    | deformations) are replaced by the local mean             |
|                     |                    |                    |                                                          |
|                     |                    |                    | deformation.                                             |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| edge_padding        | 0.1                | float              | All masks are cut of close to the image edge, i.e. if    |
|                     |                    |                    |                                                          |
|                     |                    |                    | they are closer then edge_padding*axis_length. For FEM   |
|                     |                    |                    |                                                          |
|                     |                    |                    | analysis, all pixels at this edge are fixed so that      |
|                     |                    |                    |                                                          |
|                     |                    |                    | no displacement perpendicular to the axis is allowed.    |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| padding_cell_layer  | 0.2                | float              | If you are analyzing cell layers, and additional         |
|                     |                    |                    |                                                          |
|                     |                    |                    | region close to the image edge is ignored when           |
|                     |                    |                    |                                                          |
|                     |                    |                    | analyzing stresses, to avoid boundary effects.           |
|                     |                    |                    |                                                          |
|                     |                    |                    | The effectively ignored region for cell layers is        |
|                     |                    |                    |                                                          |
|                     |                    |                    | edge_padding + padding_cell_layer.                       |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| min_obj_size        | 1500               | int                | Minimum size of an object (cell or cell colony).         |
|                     |                    |                    |                                                          |
|                     |                    |                    | All masks are added up and all encircled areas are       |
|                     |                    |                    |                                                          |
|                     |                    |                    | filled to determine the object size.                     |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| cv_pad              | 0                  | int,float          | File names. Include the ending (e.g. ".png")             |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| TFM_mode            | "finite_thickness" | string             | Using a TFM algorithm assuming either                    |
|                     |                    |                    |                                                          |
|                     |                    |                    | finite substrate thickness ("finite_thickness")          |
|                     |                    |                    |                                                          |
|                     |                    |                    | for infinite substrate thickness ("infinte_thickness").  |
|                     |                    |                    |                                                          |
|                     |                    |                    | Always use "finite_thickness".                           |
+---------------------+--------------------+--------------------+----------------------------------------------------------+



.. _OverviewofPlottingParameters:

Overview of Plotting Parameters
---------------------------------

+---------------------+--------------------+--------------------+----------------------------------------------------------+
|Parameter            |    Default Value   |   Type             |          Meaning                                         |
+=====================+====================+====================+==========================================================+
| file_names          |     specific       | string             | File names. Include the ending (e.g. ".png")             |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| cmap                |     "rainbow"      | string             | Color maps. All matplotlib color maps                    |
|                     |                    |                    |                                                          |
|                     |                    |                    | are accepted.                                            |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| vmin                |     None           | float, int, None   | Minimal value of the color bar. None                     |
|                     |                    |                    |                                                          |
|                     |                    |                    | for automatic selection.                                 |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| vmax                |     None           | float, int, None   | Maximal value of the color bar. None                     |
|                     |                    |                    |                                                          |
|                     |                    |                    | for automatic selection.                                 |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| **Color bar Parameters**                                                                                                 |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| cbar_style          |    "clickpoints"   | "clickpoints" or   | Specifies whether the color bar is plotted               |
|                     |                    |                    |                                                          |
|                     |                    | "outside"          | inside or outside of the image.                          |
|                     |                    |                    |                                                          |
|                     |                    |                    | Plotting the color bar outside will lead                 |
|                     |                    |                    |                                                          |
|                     |                    |                    | to misaligned images in clickpoints.                     |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| cbar_axes_fraction  |    0.2             | float <1           | Height of the color bar when using cbar_style            |
|                     |                    |                    |                                                          |
|                     |                    |                    | "outside". This number signifies the fraction            |
|                     |                    |                    |                                                          |
|                     |                    |                    | of the length of the original image axis.                |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| cbar_width          |    "2%"            | string             | Width of the color bar when using cbar_style             |
|                     |                    |                    |                                                          |
|                     |                    |                    | "clickpoints". Has to be a string                        |
|                     |                    |                    |                                                          |
|                     |                    |                    | signifying the percentage of                             |
|                     |                    |                    |                                                          |
|                     |                    |                    | of the original image axis.                              |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| cbar_height         |    "50%"           | string             | Height of the color bar when using cbar_style            |
|                     |                    |                    |                                                          |
|                     |                    |                    | "clickpoints". Has to be a string                        |
|                     |                    |                    |                                                          |
|                     |                    |                    | signifying the percentage of                             |
|                     |                    |                    |                                                          |
|                     |                    |                    | of the original image axis.                              |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| cbar_borderpad      |    6               | int                | Distance between the color bar and                       |
|                     |                    |                    |                                                          |
|                     |                    |                    | the right image edge.                                    |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| cbar_str            |    specific        | string             | Title of the color bar.                                  |
|                     |                    |                    |                                                          |
|                     |                    |                    | Use quotation marks ("") in the config file.             |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| cbar_title_pad      |    10              | int                | Distance between the color bar and the                   |
|                     |                    |                    |                                                          |
|                     |                    |                    | color bar title.                                         |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| cbar_tick_label_size|    15              | int                | Size of the color bar tick labels.                       |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| **Arrows in Deformation and Traction Fields**                                                                            |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| filter_factor       |    1               | float,int > 0      | Factor that defines how many arrows are                  |
|                     |                    |                    |                                                          |
|                     |                    |                    | filtered out for plotting (traction and                  |
|                     |                    |                    |                                                          |
|                     |                    |                    | deformation fields). A high filter_factor                |
|                     |                    |                    |                                                          |
|                     |                    |                    | means less arrows are plotted.                           |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| scale_ratio         |    0.2             | float (0,1]        | Length of the arrows (deformation and                    |
|                     |                    |                    |                                                          |
|                     |                    |                    | traction fields). Arrows are scaled so that the          |
|                     |                    |                    |                                                          |
|                     |                    |                    | longest arrow has the length scale_ratio * longest       |
|                     |                    |                    |                                                          |
|                     |                    |                    | image axis.                                              |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| width               |    0.002           | float              | Width of the arrow shaft (traction and                   |
|                     |                    |                    |                                                          |
|                     |                    |                    | deformation fields).                                     |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| headlength          |    3               | float,int          | Length of the arrow heads (traction and                  |
|                     |                    |                    |                                                          |
|                     |                    |                    | deformation fields).                                     |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| headwidth           |    3               | float,int          | Width of the arrow head (traction and                    |
|                     |                    |                    |                                                          |
|                     |                    |                    | deformation fields)                                      |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| **Plotting the Line Tensions**                                                                                           |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| background_color    |    "#330033"       | string, tuple      | Color of the background. Can be any color                |
|                     |                    |                    |                                                          |
|                     |                    |                    | format accepted by matplotlib. You can use               |
|                     |                    |                    |                                                          |
|                     |                    |                    | "cmap_0" to use the color of zero in the                 |
|                     |                    |                    |                                                          |
|                     |                    |                    | colormap used for the plot.                              |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| plot_t_vecs         |    False           | bool               | Plotting the line tension vectors.                       |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| plot_n_arrows       |    False           | bool               | Plotting the normal vectors of the                       |
|                     |                    |                    |                                                          |
|                     |                    |                    | cell boundary lines.                                     |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| linewidth           |    4               | int, float         | Width of the lines representing the                      |
|                     |                    |                    |                                                          |
|                     |                    |                    | cell boundary lines.                                     |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| border_arrow_filter |    1               | int                | Filter defining how many arrows are                      |
|                     |                    |                    |                                                          |
|                     |                    |                    | plotted along the cell boundary lines.                   |
|                     |                    |                    |                                                          |
|                     |                    |                    | Only every n-th arrow is plotted, where                  |
|                     |                    |                    |                                                          |
|                     |                    |                    | n is the border_arrow_filter.                            |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| boundary_resolution |    6               | int                | Smoothness of the lines representing the                 |
|                     |                    |                    |                                                          |
|                     |                    |                    | cell boundary lines. A high boundary_resolution          |
|                     |                    |                    |                                                          |
|                     |                    |                    | means less smooth plotting. Very low values will cost    |
|                     |                    |                    |                                                          |
|                     |                    |                    | a considerable amount of computation time.               |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| **Choosing which Plots are generated**                                                                                   |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
| plots               | \-"deformation"    | list               | List of plots that are produced in "colony" or           |
|                     |                    |                    |                                                          |
| colony              | \-"traction"       |                    | "cell layer" mode.                                       |
|                     |                    |                    |                                                          |
|                     | \-"FEM_borders"    |                    |                                                          |
|                     |                    |                    |                                                          |
|                     | \-"stress map"     |                    |                                                          |
+---------------------+--------------------+--------------------+                                                          |
| plots               | \-"deformation"    | list               |                                                          |
|                     |                    |                    |                                                          |
| cell layer          | \-"traction"       |                    |                                                          |
|                     |                    |                    |                                                          |
|                     | \-"FEM_borders"    |                    |                                                          |
|                     |                    |                    |                                                          |
|                     | \-"stress map"     |                    |                                                          |
|                     |                    |                    |                                                          |
|                     | \-"energy points"  |                    |                                                          |
+---------------------+--------------------+--------------------+----------------------------------------------------------+
