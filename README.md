## Traction Force Microscopy
Many cellular functions dependent on the mechanical interactions of cells with their substrate. Traction force microscopy measures the forces that cells exert on their substrate. Additionally the menachnical properties of whole cell colonies, in particular the distribution of stresses accros the colony surface, as well as the forces transmitted accros cell-cell contacts, can be analyzed using a Finite Elements approach (Monolayer Stress Microscopy). Here both methods are implemented in a clickpoints addon, making it easy to analyzed whole datasets.

In traction force microscopy cells are seeded on a linearly elastic substrated and allowed to adhere to its surface. The substrate, e.g. Polyacrylamide Gel, is filled with micrometer sized flourescent beads. When the cell exerts forces on its substrate it cause deformations. These deformations can be trackt by imaging the beads. First an image of the beads is taken. Then the cells are detached from their substrated, e.g. by trypsinization. Lastly another image of the beads is taken. The two images of the beads are used to calculate a deformation field using Particle Image Velocimetry. The traction forces i.e. forces applied from the cells to the substrate surface are then calculated using Fourrier Transform Traction Force Microscopy (FTTC). Earlier implementations of FTTC relied on asuming inifinte substrate thikness. Here a correction for finite substrate thikness is included.
Forces that are exerted from a cell to its substrate must be belanced by the cell internally or at contact points to other cells. The internal stress state of a cell patch are calculate by Monolayer Stress Microscopy. In brief the cell sheat is modeled as a 2D surfaces. The traction forces calculated from FTTC are applied to this cell sheet. Internal streses are the recoverd by using standard 2D Finite Elements approache. The user can mark cell borders,along wich line stresses i.e. the force transmitted per line segment, are calculated. Additionally stress measures accros the whole colony area, such as the average shear and normal stress are calculated. 


# Literature

Standard Fourrier Transform Traction Force Microcopy:

**Traction fields, moments, and strain energy that cells exert on their surroundings**<br>
James P. Butler, Iva Marija Tolić-Norrelykke, Ben Fabry, and Jeffrey J. Fredberg<br>
[*Am J Physiol Cell Physiol 282: C595–C605, (2002)*](https://www.physiology.org/doi/pdf/10.1152/ajpcell.00270.2001)


Fourrier Transform Traction Force Microcopy with finite substrate thikness:

**Physical forces during collective cell migration**<br>
Xavier Trepat, Michael R. Wasserman, Thomas E. Angelini, Emil Millet, David A. Weitz,
James P. Butler and Jeffrey J. Fredberg<br>
[*Nature Physics volume 5, pages 426–430 (2009)*](https://www.nature.com/articles/nphys1269)

Monolayer stress microcopy:

**Monolayer Stress Microscopy: Limitations, Artifacts, and Accuracy of Recovered Intercellular Stresses**<br>
Dhananjay T. Tambe, Ugo Croutelle, Xavier Trepat, Chan Young Park, Jae Hun Kim, Emil Millet,
James P. Butler, Jeffrey J. Fredberg<br>
[*PLOS ONE 8(2): e55172 (2013)*](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0055172)


# Installation

It is recomanded to use this package with anaconda. 
If you are on windows you need the Microsoft Visual C++ build tools. Download and install them from [here] (https://visualstudio.microsoft.com/de/thank-you-downloading-visual-studio/?sku=BuildTools&rel=16).
Additionally you need to install clickpoints. See [here](https://clickpoints.readthedocs.io/en/latest/installation.html#windows-installer) for instructions.

Next install the traction force microsoft package. You can do so directly from git hub. First you need to install git.
For windows you can use
```
conda install git
```
then install this package by
```
pip install git+https://github.com/fabrylab/tracktion_force_microscopy.git
```

To activate the clickpoints addon download or clone this repository. Copy and paste the folder TFM_addon to the clickpoints addons subdirectory.

# Performing traction force microscopy with clickpoints. 

## Generating a clickpoints database from images

Prior to using the addon, you need to generate a clickpoints data base from images. Use the script build_cdb_database_TFM.py. It will use images in its current directory and sort them into frames and layers. A frame is identified by a number at the beginning of the filename (see the example folder). For each frame you need to provide two images of beads befor and after removing the cells. And another image showing the cells e.g. as a bright field image or with membrane staining. The images of the beads must start with a number, followed by eiher "before" or "after". The third image must start with a number followed by either "membrane" or "bf_before". For more details and to change the way images are sorted check the function TFM_functions_for_clickpoints.setup_database_for_tfm.
Note that the images of beads should idealy already be corrected for drift and rotation. You cna use the script function correcting_frame_shift.py for drift correction.


## Using the clickpoints addon

Open the clickpoints database. Addons can be activated by pressing on the central right most button. A window listing all available addons will open. Select "TFM_addon" and press activate. Youw should get a massage that the addon has beeen succesfully activated. Note that you can also activate other usefull addons here. One example is the "Measure Tool", used to measure distances.

![Analysis plot](images/opening_addon.png?raw=true "Optional Title")

After you have activated the Addon a new but appears on the right. Press this button to open the addon window.

![Analysis plot](images/opening_addon2.png?raw=true "Optional Title")

In this window you can set most paramters for your analysis. In the top right you can tick wich part of analysis you want to run. The raw output from these anlysis, such as the deformation field, is stored as an array in the same folder. That way it can be accesed later by other analysis steps. "Deformation", "traction force", and depending on the analysis mode "FEM analysis" or the "contractility analysis" will produce an image and added it to the clickpoints database. You can view the images by simply 
changing layers in your frame. During the analysis several measures (area of cells, contractile energy ...) are calculated. They are all stored in an text file called out.txt.
The field "apply to" allows you to run the anlysis on just the current frame or all frames at once. Note that the output file is only generated if you analyze all frames. If another outputfile exists already, it will be overwritten.
To start your analysis bress the start button on the top left.

![Analysis plot](images/main_window.png?raw=true "Optional Title")



You can choose between two diffrent analysis modes: "cell layer" and "colony". "Cell layer" assumes that the whole field of view is coverd in cells. Finite Elements analysis is performed with nodes at the edge of the field of view fixed, so that they can't move orthogonal to the image edge. You are supposed to marke two diffrent areas ("cell type1" and "cell type2"). On these areas you average stesses and average contractile energy are calcualted. "contractility_measures" will add a map of the contractile energy to the database.
A fixed area close to the image edges (10% of the image axis length) is excluded for the calcualtion of al measured. When you select "cell layer" a new button "fill cell area"
appears. This area helps you by trying to fill all areas encircled with the mask for cell type1 on all frames. Alternatively clickpoints also provides a tool to fill areas individually.
has to masks for two cell types.  

![Analysis plot](images/mode1.png?raw=true "Optional Title")


In the "colony" mode you can analyze an isolated cell colony. The finite elements analysis will not fix any node (provided you don't mark an area at the image edge). Instead the equilibrium of the system is guaranteed by correcting nodal loads for unbalanced Forces and torque. You are supposed to mark the colony edges and internal cell cell borders within the colony. 
Additionally you can circle an area around the colony to calculate contractility and contractile energy. "FEM analysis" in this mode will produce stress mesasures along the cell borders and on the wohle colony area. It will also add an image of the cell border stresses to the database. Note that the FEM analysis is only accurate if you have a high resolution in deformation and traction field. This can be achieved by choosing the "piv overlapp" close to the "piv window size". Evantually you should set this diffrence as low as 5 pixels or less (1 µm at 0.2 µm pixelsize . Unfortunately this will increase your calcualtion time. 

![Analysis plot](images/mode2.png?raw=true "Optional Title")



##  Further analyzing the results

All measured quantities are saved to an output file named "out.txt". This file is only generated when you perform the analysis
on all frames of the database. All parameters used for the analysis are written as header. All measured qunatities follow with 
the frame, name of the quantity, it's value, a unit and an additional wrning message, each in one column.

![Analysis plot](images/output_file.png?raw=true "Optional Title")

This package provdies some functions to read the output file, perform statistical test comapring two output files and plotting these results.
Check out the analysis in [output_data_analysis](/analysis_and_testing/output_data_analysis.py) for a detailed example on how to
compare the analysis results results for two cell types.
You can find a minimal example of a completed analysis in the example_analysis folder. It comapres the forces exerted by cell colonys of wild type and plektin knock out cells. 
Note that I change the path entry in the database so that they work immediately when you download this folder. 





# Common problems and things to look out for

Make sure you draw closed circles around the cell colony and the area selected for contractillity. The programm tries to fill the encircled areas, but will fail if there is a gap. If you get a warning, that your mask is tool small, this could be the reason.
Also any cell border in the cell colony mask that is not connected at both ends with another border is removed.

If there is a problem with the finite elements analysis, try incresing the resolution of deformation and traction field. 

The images in your are linked to the database with an absolute path. If you move an image or change the name of the 
folder the images are in, they will no longer be displayed. You can of cause eddit the path entry of the database.










