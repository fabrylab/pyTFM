## tracktion_force_microscopy
Many cellular functions dependent on the mechanical interactions of cells with their substrate. Traktion force microscopy measures the forces that cells exert on their substrate. Additionally the menachnical properties of whole cellcolonies ,in particular the distribution of stresses acros the colony surface, as well as the forces transmitted accros cell-cell contacts, can be analyzed using a Finite Elements approach. 

In tracktion force microscopy cells are seeded on a linearly elastic substrated and allowed to adhere to its surface. The substrate, e.g. Polyacrylamide Gel, is filled with micrometer sized flourescent beads. When the cell exerts forces on its substrate it cause deformations. These deformations can be trackt by imaging the beads. First an image of the beads is taken. Then the cells are detached from their substrated, e.g. by trypsinization. Lastly another image of the beads is taken. The two images of the beads are used to calculate a deformation field using Particle Image Velocimetry. The tracktion forces i.e. forces applied from the cells to the substrate surface are then calculated using Fourrier Transform Tracktion Force Microscopy (FTTC). Earlier implementations of 
FTTC relied on asuming inifinte substrate thikness. Here a correction for finite substrate thikness is included.
Forces that are exerted from a cell to its substrate must be belanced by the cell internally or at contact points to other cells. The internal stress state of a cell patch are calculate by Monolayer Stress Microscopy. In brief the cell sheat is modeled as a 2D surfaces. The traction forces calculated from FTTC are applied to this cell sheet. Internal streses are the recoverd by using standard 2D Finite Elements approache. The user can mark cell borders,along wich line stresses i.e. the force transmitted per line segment, are calculated. Additionally stress measures accros the whole colony area, such as the average shear and normal stress are calculated. 






# literature

Standard Fourrier Transform Tracktion Force Microcopy:

**Traction fields, moments, and strain energy that cells exert on their surroundings**<br>
James P. Butler, Iva Marija Tolić-Norrelykke, Ben Fabry, and Jeffrey J. Fredberg<br>
[*Am J Physiol Cell Physiol 282: C595–C605, (2002)*](https://www.physiology.org/doi/pdf/10.1152/ajpcell.00270.2001)


Fourrier Transform Tracktion Force Microcopy with finite substrate thikness:

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

Next install the Tracktion force microsoft package. You can do so directly from git hub. First you need to install git.
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


You can choose between two diffrent analysis modes: "cell layer" and "colony". "Cell layer" assumes that the whole field of view is coverd in cells. Finite Elements analysis is performed with nodes at the edge of the field of view fixed, so that they can't move orthogonal to the image edge. You are supposed to marke two diffrent areas ("cell type1" and "cell type2"). On these areas you can calculate average stesses and average contractile energy. When you select "cell layer" an new button "fill cell area"







has to masks for two cell types.  
![Analysis plot](images/main_window.png?raw=true "Optional Title")
![Analysis plot](images/mode1.png?raw=true "Optional Title")
![Analysis plot](images/mode2.png?raw=true "Optional Title")





