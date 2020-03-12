Using pyTFM in Python
=========================

pyTFM makes it easy to perform Traction Force Microscopy and Monolayer Stress Microscopy in python. In this
tutorial we will calculate strain energy, contractillity, mean normal stress and average line tension for in a cell
colony. As always we need two images of fluorescent beads: One image before cell removal and one image after.
Additionally, since we are not using the clickpoints addon, we need to provide 3 masks:

.. TODO: make a wiki

A mask for the area of force calculation (encircling all deformations and forces originating from the
cell colony), a mask for the Finite Elements Analysis (encircling all forces originating from the cell colony) and
a mask of the cell boundaries. Of course the easies way to generate these masks is to use clickpoints, which also
makes it easy to export the masks.

In the `example data <https://github.com/fabrylab/example_data_for_pyTFM/archive/master.zip>`__ in the subfolder
"python_tutorial" you can find images and mask for the cell colony, that we are going to analyze here. I drew the masks
in clickpoints and saved them as .png files so that you can look at them with a standard image display tool.
You can use any other format, as long as you end up with a boolean (True and False) array at the end.

The "python_tutorial" folder also contains the complete code of the tutorial in a single python script "pyTFM_tutorial",
that works immediately and prints the same quantities and produces the same figures, that we are going to see in
this tutorial.

Calculating Deformation Fields
----------------------------------

First, let's import the function to calculate the deformation field and a function to plot it:

.. code-block:: python

    from pyTFM.TFM_functions import calculate_deformation
    from pyTFM.plotting import show_quiver

Calculating the deformation field requires the images of the beads before and after bead removal.
You can provide either file paths, or arrays (must have the datatype int32). The deformation field
is calculated with Particle Image Velocimetry, using a cross correlation algorithm. You need to
find appropriate values for the window_size and overlap, that produce a smooth deformation field.

.. TODO .. explain somewhere how to choose piv windowsize

For these images you an use:

.. code-block:: python

    # paths to the images
    im_path1 = r"/home/user/Software/example_data_for_pyTFM/python_tutorial/04after.tif"
    im_path2 = r"/home/user/Software/example_data_for_pyTFM/python_tutorial/04before.tif"
    # calculating the deformation
    u, v, mask_val, mask_std = calculate_deformation(im_path1, im_path2, window_size = 100, overlap = 60)
    # The unit of window size and overlap is pixels.

The overlap of 60 is a bit to small, especially for the FEM analysis later. If you want to be more accurate
use an overlap of up to 95. This will increase the calculation time from a few seconds to 5 minutes.

Let's plot the deformation field:

.. code-block:: python

    # plotting the deformation field
    fig1, ax = show_quiver(u, v, cbar_str="deformations\n[pixels]")

show_quiver accepts most of the plotting parameters listed in :ref:`OverviewofPlottingParameters`.

.. the reference works

Calculating Traction Fields
----------------------------------

Next, we are going to calculate the traction forces, that were necessary to cause the deformations, that
we have calculated in the previous section.
This is done with the "TFM_tractions" function:

.. code-block:: python

    from pyTFM.TFM_functions import calculate_deformation, TFM_tractions

Now we have to set the elastic parameters of the substrate (Young's modulus, Poisson's ratio and
the height of the substrate). You can also set the substrate height to "infinite". If you are not sure,
whether your substrate thickness is large enough to neglect it, you can simply check if you get different results
with h = "infinite".

We also need the pixelsize of the image of the beads and the pixelsize of the deformation field. The later
can be calculate if you know the dimensions and pixelsize of the image of the beads:

.. code-block:: python

    ps1 = 0.201 # pixel size of the image of the beads
    im1_shape = (1991, 2033) # dimensions of the image of the beads
    ps2 = ps1*np.mean(np.array(im1_shape) / np.array(u.shape)) # pixel size of of the deformation field
    young = 49000 # Young's modulus of the substrate in Pa
    sigma = 0.49 # Poisson's ratio of the substrate
    h = 300 # height of the substrate in µm, "infinite" is also accepted

Finally the traction field can be calculated by:

.. code-block:: python

    tx, ty = TFM_tractions(u, v, pixelsize1=ps1, pixelsize2=ps2, h=h, young=young, sigma=sigma)

We can plot it in the same way as we plotted the deformation field:

.. code-block:: python

    fig2, ax = show_quiver(tx, ty, cbar_str="tractions\n[Pa]")


Quantifying the Force Generation
-------------------------------------

In order to quantify the force generation of the cell colony, we have to select the area, where deformations
and tractions, that are generated by the colony are located. This selection requires a mask, which is a boolean
array, that has the value True in the area that we want to use. I produced the appropriate mask in clickpoints and
saved it as a grey scale image as "force_measurement.png". After loading the mask, there are two more things to be done:
First, we need to fill all holes in the mask in order to produce a continous area. Second, we need to resize the
mask to the dimensions of the deformation and traction fields:


.. code-block:: python

    import matplotlib.pyplot as plt
    from scipy.ndimage.morphology import binary_fill_holes
    from pyTFM.grid_setup_solids_py import interpolation # a simple function to resize the mask


    # loading a mask, marking the are that is used for measuring the force generation
    mask = plt.imread(r"/home/user/Software/example_data_for_pyTFM/python_tutorial/force_measurement.png").astype(bool)
    mask = binary_fill_holes(mask) # the mask should be a single patch without holes
    # changing the masks dimensions to fit to the deformation and traction field:
    mask = interpolation(mask, dims=u.shape)

This mask can now be used to calculate the contractillity and the strain energy:

.. code-block:: python

    from pyTFM.TFM_functions import strain_energy_points, contractillity

    # Strain energy:
    # First we calculate a map of strain energy
    energy_points = strain_energy_points(u, v, tx, ty, ps1, ps2) # J/pixel
    # Then we sum all energy points in the area defined by mask
    strain_energy = np.sum(energy_points[mask]) # 1.95*10**-13 J

    # Contractillity
    contractile_force, proj_x, proj_y, center = contractillity(tx, ty, ps2, mask) # 2.03*10**-6 N



Measuring stresses in Cell Colonies
-------------------------------------

Stresses are calculated with Finite Elements Methods, modelling the colony as a 2 dimensional sheet and
applying force opposite to the traction forces to it. We need select two areas for this:
First, the are used for the Finite Elements Analysis. This area should contain all traction
forces that are produced by the cell colony. Due to inaccuracies in the calculation of
traction forces this area is typically larger then the actual cell colony. All measures for stress
are however evaluated only on the area of the cell colony, which is the second area that we need
to provide for this step. Once again, I produced the appropriate masks with clickpoints. We can load
them like this:

.. code-block:: python

    # first mask: The area used for Finite Elements Methods.
    # should encircle all forces generated by the cell colony
    mask_FEM = plt.imread(r"/home/user/Software/example_data_for_pyTFM/python_tutorial/FEM_area.png").astype(bool)
    mask_FEM = binary_fill_holes(mask_FEM) # the mask should be a single patch without holes
    # changing the masks dimensions to fit to the deformation and traction field:
    mask_FEM = interpolation(mask_FEM, dims=tx.shape)

    # second mask: The area of the cells. Average stresses and other values are calculated only
    # on the actual area of the cell, represented by this mask.
    mask_cells = plt.imread(r"/home/user/Software/example_data_for_pyTFM/python_tutorial/cell_borders.png").astype(bool)
    mask_cells = binary_fill_holes(mask_cells)
    mask_cells = interpolation(mask_cells, dims=tx.shape)

The traction forces in the FEM area are typically slightly unbalanced, leading to a net force and torque acting on
the cell colony. We need to correct this:

.. code-block:: python

    from pyTFM.grid_setup_solids_py import prepare_forces

    # converting tractions (forces per surface area) to actual forces
    # and correcting imbalanced forces and torques

    # tx->traction forces in x direction, ty->traction forces in y direction
    # ps2->pixel size of the traction field, mask_FEM-> mask for FEM
    fx, fy = prepare_forces(tx, ty, ps2, mask_FEM)

Now we are ready to perform a Finite Elements Analysis. This is split into two steps:
First, a FEM grid is setup. The grid consists of nodes and for each node
the connectivity to other nodes, constraints on the displacements and forces acting on the node are defined .
Then, the FEM system is solved by calculating first the deformations and eventually the stress
that the applied forces would cause in the FEM grid.

.. code-block:: python

    from pyTFM.grid_setup_solids_py import grid_setup, FEM_simulation

    # constructing FEM grid
    nodes, elements, loads, mats = grid_setup(mask_FEM, -fx, -fy, sigma=0.5)
    # performing FEM analysis
    # Verbose prints the progress of numerically solving the FEM system of equations.
    UG_sol, stress_tensor = FEM_simulation(nodes, elements, loads, mats, mask_FEM, verbose=True)
    # UG_sol is a list of deformations for each node. We don't need it here.

The stress tensor completely defines the forces in the cell colony. We can vor example extract the
average mean normal stress and the coefficient of variation of the mean normal stress
(quantifying how much the stress varies in the colony). We will use the mask "mask_cells" which
marks the actual area of the cell colony for these measurements.

.. code-block:: python

    # mean normal stress
    ms_map = ((stress_tensor[:, :, 0, 0] + stress_tensor[:, :, 1, 1]) / 2) / (ps2 * 10**-6)
    # average on the area of the cell colony.
    ms = np.mean(ms_map[mask_cells]) # 0.0044 N/m

    # coefficient of variation
    cv = np.nanstd(ms_map[mask_cells]) / np.abs(np.nanmean(ms_map[mask_cells])) # 0.35 no unit



Calculating the Line Tension
-------------------------------------


A particularly interesting question is how much forces are transmitted across cell-cell boundaries.
This is quantified by the line tension. First we need to load a mask, marking all cell borders.
Note that this is the same mask that is used to get the area of the cell colony, only this time
we are not going to fill any holes or resize the mask.


 .. code-block:: python

    # loading a mask of the cell borders
    mask_borders = plt.imread(r"/home/user/Software/example_data_for_pyTFM/python_tutorial/cell_borders.png").astype(bool)

The cell-cell borders are stored in an object "borders", which among others contains a spline interpolation of
each border, assigns each border to a cell and contains a list of borders located at the edge of the cell borders:

.. code-block:: python

    from pyTFM.grid_setup_solids_py import find_borders

    # identifying borders, counting cells, performing spline interpolation to smooth the borders
    borders = find_borders(mask_borders, tx.shape)
    n_cells = borders.n_cells # For example you can get the number of cells from the "borders" object

We can use the cell-cell borders together with the stress tensor. To calculate the line tension.
The line tension is a force vector acting on a small slice of a cell border. We are going to
calculate the average length of this vector ("avg_line_tension") and the average normal
component of the line tension ("avg_normal_line_tension"):

.. code-block:: python

    from pyTFM.stress_functions import lineTension

    #calculating the line tension along the cell borders.
    lt, min_v, max_v = lineTension(borders.lines_splines, borders.line_lengths, stress_tensor, pixel_length=ps2)
    # lt is a nested dictionary. The first key is the id of a cell border. For each cell border
    # the line tension vectors ("t_vecs"), the normal and shear component of the line tension ("t_shear") and
    # the normal vectors of the cell border ("n_vecs") are calculated at a large number of points.

    # average norm of the line tension only borders not at colony edge are used.
    lt_vecs = np.concatenate([lt[l_id]["t_vecs"] for l_id in lt.keys() if l_id not in borders.edge_lines])
    avg_line_tension = np.mean(np.linalg.norm(lt_vecs, axis=1)) # 0.00555 N/m

    # average normal component of the line tension
    lt_normal = np.concatenate([lt[l_id]["t_normal"] for l_id in lt.keys() if l_id not in borders.edge_lines])
    avg_normal_line_tension = np.mean(np.abs(lt_normal)) # 0.00552 N/m,
    # here you can see that almost all line tension acts perpendicular to the cell borders.

Finally let's produce a plot of the line tension:

.. code-block:: python

    from pyTFM.plotting import plot_continuous_boundary_stresses

    # plotting the line tension
    fig3, ax = plot_continuous_boundary_stresses([borders.inter_shape, borders.edge_lines, lt, min_v, max_v], cbar_style="outside")


