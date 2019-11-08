## Traction Force Microscopy
Many cellular functions depended on the mechanical interactions of cells with their substrate. Traction force microscopy measures the forces that cells exert on their substrate. Additionally the mechanical properties of whole cell colonies, can be analyzed using a Finite Elements approach (Monolayer Stress Microscopy). These properties are foe example  the distribution of stresses across the colony surface, as well as the forces transmitted across cell-cell contacts. Here both methods are implemented in a clickpoints addon, making it easy to quickly analyzed whole datasets.

In traction force microscopy cells are seeded on a linearly elastic substrate and allowed to adhere to its surface. The substrate, e.g. Polyacrylamide Gel, is filled with micrometer sized fluorescent beads. When the cell exerts forces on its substrate it causes deformations. These deformations can be tracked by imaging the beads. First an image of the beads is taken. Then the cells are detached from their substrates, e.g. by trypsinization. Now another image of the beads is taken. The two images of the beads are used to calculate a deformation field with Particle Image Velocimetry. The traction forces i.e. forces applied from the cells to the substrate surface are then calculated using Fourrier Transform Traction Force Microscopy (FTTC). Earlier implementations of FTTC relied on assuming infinite substrate thickness. Here a correction for finite substrate thickness is included.
Forces that are exerted from a cell to its substrate must be balanced by the cell internally or at contact points to other cells. The internal stress state of a cell patch is calculate by Monolayer Stress Microscopy. In brief the cell sheet is modeled as a 2D surfaces. The traction forces calculated from FTTC are applied to this cell sheet. Internal stresses are then recovered by using a standard 2D Finite Elements approach. The user can mark cell borders,along which line stresses i.e. the force transmitted per line segment, are calculated. Additionally stress measures across the whole colony area, such as the average shear and normal stress are also calculated. 


# Literature

Standard Fourrier Transform Traction Force Microcopy:

**Traction fields, moments, and strain energy that cells exert on their surroundings**<br>
James P. Butler, Iva Marija Tolić-Norrelykke, Ben Fabry, and Jeffrey J. Fredberg<br>
[*Am J Physiol Cell Physiol 282: C595–C605, (2002)*](https://www.physiology.org/doi/pdf/10.1152/ajpcell.00270.2001)


Fourrier Transform Traction Force Microcopy with finite substrate thickness:

**Physical forces during collective cell migration**<br>
Xavier Trepat, Michael R. Wasserman, Thomas E. Angelini, Emil Millet, David A. Weitz,
James P. Butler and Jeffrey J. Fredberg<br>
[*Nature Physics volume 5, pages 426–430 (2009)*](https://www.nature.com/articles/nphys1269)

Monolayer Stress Microcopy:

**Monolayer Stress Microscopy: Limitations, Artifacts, and Accuracy of Recovered Intercellular Stresses**<br>
Dhananjay T. Tambe, Ugo Croutelle, Xavier Trepat, Chan Young Park, Jae Hun Kim, Emil Millet,
James P. Butler, Jeffrey J. Fredberg<br>
[*PLOS ONE 8(2): e55172 (2013)*](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0055172)


# Installation

It is recommended to use this package with anaconda. 
If you are on Windows you need the Microsoft Visual C++ build tools. Download and install them from [here] (https://visualstudio.microsoft.com/de/thank-you-downloading-visual-studio/?sku=BuildTools&rel=16).
Additionally you need to install clickpoints. See [here](https://clickpoints.readthedocs.io/en/latest/installation.html#windows-installer) for instructions.

Next install the traction force microscopy package. I recommend to make an editable installation with pip.
Simply download or clone the repository. Unizp the files, then open a terminal and navigate to the folder tracktion_force_mircoscopy. Now perform a local installation with the command

```
pip install -e .
```

You can also install this package directly from github. First you need to install git.
For windows you can use
```
conda install git
```
then install this package by
```
pip install git+https://github.com/fabrylab/tracktion_force_microscopy.git
```

To activate the clickpoints addon download or clone this repository. Copy and paste the folder "TFM_addon" to the clickpoints addons subdirectory.If you don' know where clickpoints is installed, simply open a terminal and go into the python interpreter by typing the command "python". Then import clickpoints and get the location with:
```
import clickpoints
clickpoints.__file__
```
The addons folder is located directly in this folder.

# Performing traction force microscopy with clickpoints.

## Correction of frame shift in images:
Slight movements of the microscopes field of view while removing cells will reuslt in non-zero drift in the deformation field.
You can use the script "correcting_frame_shift.py" to correct for any shift between the images of the beads before and after the removal of cells. You can find it in the analysing_and_testing folder. This script searches a folder tree and finds pairs of images before and after cell removal. Then it uses image registration to find a shift. The common field of view of the image pair is determinied and they are cropped accordingly. If the script finds additional images of the cell, it will also cut these images to the same field of view. 

The basic folder structure of input and output looks like this:
![Analysis plot](images/frame_shift_folder_structure.png?raw=true "Optional Title"). 

On the right you can see variables that you can set.These variables define how images are named and identified (see below). The folder and filenames in this image would work for default values of these variables.

In general the script searches through the full tree of an input folder. Once it identifies a folder for the images after bead removal and before bead removal in the same subdirectory it stops and enters these folders. There is searches for the corresponding images and identifies their frame. Note that it searches for bright field or other images of the cells in both folders.The output images are saved in directories named the same as the directories the "before" and "after" folder were been found in. These directories are located directly in the input folder.
If you want to use the script with default settings open a terminal, navigate to the folder that contains your data. 
Then execute the script with. 

```
python path_to_the_script/correcting_frame_shift.py
```

Alternatively and if you want to change how files are recognized and named check out the correcting_frame_shift.py script.
There you will also find extensive explanations on all parameters.
You can for example open it in pycharm and play around with the various parameters 

You can basically change 3 parameters:
1. How the script finds the correct folders for "before" and "after" images. 
2. How the script identifies the image files, and the frame an image belongs to.
3. The filenames for the drift corrected images.
Inside the script you will find more detailed explanations on how to set the parameters.


## Generating a clickpoints database from images

Prior to using the addon, you need to generate a clickpoints data base from images. The data base will sort images into layers and frames. For each layer you need to provide 3 images: an image of beads before the cells were removed, an image of beads after the cells where removed and an image that should show the cells themselves. The default way to identify layers and frames is the following:
The frame is identified by a number at the beginning of the file name. The number can be up to 4 digits long and can be zeropadded. You can leave out frames, e.g provide images for frame 1 and 3, but leave out frame 2. Each frame needs 3 images: one that starts with the frame number followed by "after", one that starts with the frame number followed by "before" and one that starts with the frame number followed by "membrane"."bf_before" is also accepted. These are the images for beads after and before cells wher removed and an image of the cells themselves. You can check out the example_analysis for an example of how to name your image files. You can also provide your own keys to sort the images. See [generating a data base with the python interpreter](#Generating-a-data-base-with-the-python-interpreter) for more details.

### Using the build_cdb_database_TFM.py script
Building the data base can be done in two ways.
The first possibility is to use the script build_cdb_database_TFM.py. Open a terminal, navigate to the folder your images are located in and call the script with 
```
python build_cdb_database_TFM.py
```
If the script build_cdb_database_TFM.py is not located in the same folder, your need to provide the full path to the script.

### Generating a data base with the python interpreter: 
Alternativly you can build the database from a python interpreter directly. Either open python in the terminal or use the interpreter for example in PyCharm.
Import the function to setup a data base:
```
from andreas_TFM_package.database_functions import setup_database_for_tfm
```
set the folder that your images are located in. It could look like this:
```
folder = r"C:\Users\Andy\Desktop\exmaple_analysis\KO_shift_parts"
```
call the setup_database_for_tfm function with default parameters on this folder. You have to provide the folder and a name for your data base. The name has to end with ".cdb".
```
setup_database_for_tfm(folder,"database.cdb")
```
You can also pass your own keys to sort images as parameters to the function. For each parameter you need to provide a regular expression. The regular expression to identfy a frame needs to contain the actual frame number as a group marked with "()". Don't include the file ending. ".png" or ".tif" and so on is added automatically. 
```
# using custom search keys for the data base setup:

key1="\d{1,4}after" # finds any file name with 1 to 4 numbers followed by "after"
key2="\d{1,4}before" # finds any file name with 1 to 4 numbers followed by "before"
key3="\d{1,4}membrane" # finds any file name with 1 to 4 numbers followed by "membrane"
frame_key= "(\b\d{1,4})" # finds the first 4 leading numbers of the filename

# build the data base:
setup_database_for_tfm(folder,"database.cdb",key1=key1,key2=key2,key3=key3,frame_key=frame_key)
```

Note that the images of beads should ideally be corrected for drift and rotation. You can use the script correcting_frame_shift.py for drift correction.


## Using the clickpoints addon
For general instructions on how to use clickpoints go [here](https://clickpoints.readthedocs.io/en/latest/).

Open the clickpoints database. Addons can be activated by pressing on the central right most button. A window listing all available addons will open. Select "TFM_addon" and press activate. You should get a massage that the addon has been succesfully activated. Note that you can also activate other useful addons here. One example is the "Measure Tool", used to measure distances.

![Analysis plot](images/opening_addon.png?raw=true "Optional Title")

After you have activated the Addon a new button appears on the right. Press this button to open the addon window.

![Analysis plot](images/opening_addon2.png?raw=true "Optional Title")

In this window you can set most paramters for your analysis. In the top right you can tick which part of analysis you want to run. The raw output from these analysis, such as the deformation field, is stored as an array in the same folder. That way it can be accessed later by other analysis steps. Deformation and traction fields will also be plotted and added to your database in a new layer. Depending on the analysis mode "FEM analysis" or the "contractillity analysis" will also produce an image and added it to the clickpoints database. You can view the images by simply changing layers in your current frame in clickpoints.
During the analysis several measures (area of cells, contractile energy and so on) are calculated. They are all stored in an text file called out.txt.
The field "apply to" allows you to run the analysis on just the current frame or all frames at once. Note that the output file is only generated if you analyze all frames. If another output file exists already, it will be overwritten.
To start your analysis press the start button on the top left.

![Analysis plot](images/main_window.png?raw=true "Optional Title")



You can choose between two different analysis modes: "cell layer" and "colony". "Cell layer" assumes that the whole field of view is covered in cells. Finite Elements analysis is performed with nodes at the edge of the field of view fixed, so that they can't move orthogonal to the image edge. You are supposed to mark two different areas ("cell type1" and "cell type2"). On these areas you average stresses and average contractile energy are calculated. "contractillity_measures" will add a map of the contractile energy to the database.
A fixed area close to the image edges ( default 10% of the image axis length) is excluded for the calculation of all measures. When you select "cell layer" a new button "fill cell area"
appears. This button helps you by trying to fill all areas encircled with the mask for cell type1 on all frames. Alternatively clickpoints also provides a tool to fill areas individually.
has to masks for two cell types.  

![Analysis plot](images/mode1.png?raw=true "Optional Title")


In the "colony" mode you can analyze an isolated cell colony. The finite elements analysis will not fix any node (provided you don't mark an area at the image edge). Instead the equilibrium of the system is guaranteed by correcting nodal loads for unbalanced forces and torque. You are supposed to mark the colony edges and internal cell cell borders within the colony. 
Additionally you can circle an area around the colony to calculate contractillity and contractile energy. "FEM analysis" in this mode will produce stress measures along the cell borders and on the whole colony area. It will also add an image of the cell border stresses to the database. Note that the FEM analysis is only accurate if you have a high resolution in deformation and traction field. This can be achieved by choosing the "piv overlapp" close to the "piv window size". Eventually you should set this difference as low as 5 pixels or less (1 µm at 0.2 µm pixelsize) . Unfortunately this will increase your calculation time. 


![Analysis plot](images/mode2.png?raw=true "Optional Title")



##  Further analyzing the results

All measured quantities are saved to an output file named "out.txt". This file is only generated when you perform the analysis
on all frames of the database. The parameters used for the analysis are written as header. All measured quantities follow with 
the frame, name of the quantity, it's value, a unit and an additional warning message, each in one column. You should reconsider the results if you see a warning.

![Analysis plot](images/output_file.png?raw=true "Optional Title")

This package provides some functions to read the output file, perform statistical test comparing two output files and plotting these results.
Check out the analysis in [output_data_analysis](/analysis_and_testing/output_data_analysis.py) for a detailed example on how to
compare the analysis results results for two cell types.
You can find a minimal example of a completed analysis in the example_analysis folder. It compares the forces exerted by cell colonies of wild type and plectin knock out cells. 
Note that I change the path entry in the database so that they work immediately when you download this folder. 

## Changing plotting options and default parameters

Default parameters and plotting options can be found in the [file parameters_and_string] (/andreas_TFM_package/parameters_and_strings.py). Go to this file for more details. The parameters are mostly stored in dictionaries. You can freely change entries in these dictionaries, if you installed your package as editable. Note that currently reinstalling the package will overwrite changes you make there.

# Common problems and things to look out for


Make sure you draw closed circles around the cell colony and the area selected for contractillity. The program tries to fill the encircled areas, but will fail if there is a gap. If you get a warning, that your mask is tool small, this could be the reason.
Also any cell border in the cell colony mask that is not connected at both ends with another border is removed.

If there is a problem with the finite elements analysis, try increasing the resolution of deformation and traction field. 

The images in your are linked to the database with an absolute path. If you move an image or change the name of the 
folder the images are in, they will no longer be displayed. You can of cause edit the path entry of the database.









