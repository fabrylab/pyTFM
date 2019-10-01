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

# Using the clickpoints addon. 
Prior to using the addon, you need to generate a clickpoints data base from images. Use the script build_cdb_database_TFM.py. It will use images in its current directory and sort them into frames and layers. A frame is identified by a number at the beginning of the filename (see the example folder). For each frame you need to provide two images of beads befor and after removing the cells. And another image showing the cells e.g. as a bright field image or with membrane staining. The images of the beads must start with a number, followed by eiher "before" or "after". The third image must start with a number followed by either "membrane" or "bf_before". For more details and to change the way images are sorted check   











![Analysis plot](images/opening_addon.png?raw=true "Optional Title")
![Analysis plot](images/opening_addon2.png?raw=true "Optional Title")
![Analysis plot](images/main_window.png?raw=true "Optional Title")
![Analysis plot](images/mode1.png?raw=true "Optional Title")
![Analysis plot](images/mode2.png?raw=true "Optional Title")





