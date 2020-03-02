Tutorial
=================================


Analyzing Force Generation and Stresses in Cell-Colonies with Clickpoints
---------------------------------------------------------------------------

Prequisites
-------------
Using the clickpoints-addon in pyTFM obviously requires clickpoints. If you set up clickpoints correctly, you can
open images by right clicking on the image files and select "open with clickpoints".


The Data
-----------
The pyTFM package that you can find on github contains an example dataset in the /example_analysis folder.
It contains raw data for 2 types of cell colonies: In one group critical cytoskeletal protein has been knocked out.
We will compare these cell colonies to a set of wildtype colonies.
As you can see in Figure... there are 6 images in for each colony type. This corresponds to two field of views that
for each colonies. For each field of view you need 3 images. One image (03after.tif) shows the the colony and the
boundaries between cells. In this case the image shows fluorescently stained cell membranes.
The other two images show beads that are embedded in substrate that the cell lie on. One image was recorded before
the cells were removed (03before.tif) and the other was recorded after the cells were removed (03after.tif).
The number in front of the filename ("03","10" and so on) indicates which field of view that image belongs to.

.. figure:: data.png
  :width: 400
  :alt: Data structure in the example data set


Opening Clickpoints and sorting Images
------------------------------------------

The first step to analyze the data is to create a clickpoints database, in which the images are identfied correctly,
concerning their type (whether its an image of the cells or an image of the beads before or after the cell removal)
and the field of view (formul) they belong to.

To open a database simple right click on an image and select "open with clickpoints".

.. figure:: open_with_clickpoints.png
  :width: 550
  :alt: Data structure in the example data set

Clickpoints sorts images in two dimensions: Frames and layers. The frames are displayed in the bar at the bottom.
You can skip from frame to frame  using the left and right buttons on your key bord (formul). Layers can be changed
with the "Page Up" and "Page Down" buttons. When you opend the database, you will notice that there is only one layer
and every image is sorted into a new frame. Our goal is to sort each field of view into one frame, with three layers
each representing one type of image. In order to do this you need to open the pyTFM addon and open the "select image"
menu. Follow the steps described in Figure ...

.. figure:: open_select_images.png
  :width: 750
  :alt: Data structure in the example data set

  A: Open the Addon-browser in clickpoints. A new window
  , with all available addons will open. B: Activate the pyTFM
  addon by selecting pyTFM and clicking the "Activate"
  buton. A notification that the addon has been loaded
  successfully will appear. After you press "OK" a new Icon
  will appear in the clickpoints main window, to the right of the
  Addon-browser button. C: Click on this button to open the
  pyTFM Addon. D: Finally, open the menu to select images
  by pressing the "select images" button.


The "file selection" manu allows you to do three things: First you can select where images are located  how they are classified. You can also set an outputfolder, where the database file and all analysis
results will be saved and you can choose a name for the output database.
As mentioned above, the analysis requires three types of images. For each type you can select a folder
(left hand side) and a regular expression that identifies the image type from the image filename.

.. note::
   Regular expressions are the standart way to find patterns in text. For example, it allows you to
   identifie numbers of certain length, groups of characters or the beginning and end of a text. You
   find more inforamtion on regular expresions `here <https://docs.python.org/3/library/re.html>`_.
   Some usefull expressions are listed in the table below:


   ==============    ==============================================================
   search pattern     meaning
   ==============    ==============================================================
   after              all files with "after" in the filename
   ^after             all files with "after" at the beginning of the filename
   after$             all files with "after" at the end of the filename
   \*                 all files
   ^(\d{1-4}) 	      up to 4 numbers at beginning of the filename
   (\d{1-4}) 	      up to 4 consecutive numbers anywhere in the filename
   (\d{1-4})$ 	      up to 4 numbers at end of the filename
   ==============    ==============================================================



The "'after' image identifier" identifies images of beads after cell removal, the "'before' image identifier"
identifies images of beads before cell removal and the "cell image identifier" identifies images that
show the cells or cell membranes. Finally, there is a sperate regular expression, the
"frame identifier" that identifies the field view each image belongs to. This must point to a
number (e.g."02" or "2" and so on) in the image filename. This number must be specially marked as a group
(formul..) by using brackets "()". Note that the extension (".png",".tiff", ".jpeg" ...) must not
be included in the identifiers. The default identifiers fit to the example data set, meaning that
if you are using this data set or the same naming scheme for your images, you don't need to change
anything.

Once you have entered identifiers for image types, frames, the output folder and the database name
press the "collect image" button. You should see something like this:


.. figure:: output_select_images.png
  :width: 600
  :alt: Data structure in the example data set

Make sure your database didn't contain any masks that you dont't want to delete. If you just opend the
databse from new images, you can press OK. The path to the images that are sorted into the databse,
the type of the images (layer) and the field of view of the images (frame) are printed to the console.
Make sure all images are sorted correctly. The program generates a clickpoints database and sorts
images into layers and frames. Your clickpoints window updates automatically.

.. TODO: mention correct Drift

.. TODO: paramters seting and recomondation



Setting Parameters
--------------------------------------------
Lets continue with calculating the deformation and traction field. Go to the pyTFM addon window
(Figure...).


.. figure:: main.png
   :width: 600
   :alt: Main addon window



In this window you have to set the mechanical parameters of the substrate ("Youngs modulus" and "Poisson's
ratio"), the height of the substrate ("gel height") and the pixel size ("pixel size"). Then you have to
parameters for the calculation of the deformation field.

The deformation field is calculated with
particle image velocimetry. This method essentially cuts out square shaped patches from image of
the beads before cell removal, places them on the image of beads after cell removal
and checks how well they fit together. This is done for many positions and the vector
form the original position of the patch and the
position where the patch fits best to the image of beads after cell removal is the displacement vector.
You can control two things: the size of the patch that is cut out of the image of the beads
after cell  removal (with the parameter "PIV window size") and the resolution of the
resulting displacement field (with the parameter "PIV overlap"). I recommand to choose the window size
to be roughly 7 times the beads diameter, you should however try a few values an check which
window size yields a smooth yet accurate deformation field.

.. Note::
   You can measure the beads diameter directly in the clickpoints, using another addon:
   The Measrue Tool

The "PIV overlapp" mainly controlls the resolution of the resulting displacement field and must be
smaller then the "PIV window size" but at least half of the "PIV window size". You need
a high resoultion for analyzing stress. In this step the area of cells should at least contain 1000
pixels. However, if you are not calculating stresses, you can save a lot of calculation time by choosing a
"PIV overlapp" closer to half of the "PIV window size".

For this tutorial you can keep all parameters at their default value.

Calculating Traction and Deformation Fields
--------------------------------------------
Once you have set all parameters you can start the calculation: Use the tick boxes in the upper right to select
which part of the analysis you want to perfor. For now we are gonna select "deformation" and "traction". Then
use the "apply to" option to choose whether all frames should be analyzed or only the frame that you are currently
viewing. Finally press "start" in the upper left to begin the analysis. With the default parameters this takes
about 5 minutes per frame. "calculation completed" is printed to the console, once all frames have been analyzed.

The traction and deformation fields are added to the database as new layers. Switch to these layers using the "page up"
key on your keyboard. Traction and deformation for the first frame in the wild type data should look like this:


.. figure:: def_trac_res.png
   :width: 750
   :alt: Main addon window

If you do not see the display tool and mask names ("membrane", "force measures", "FEM area" on the right press F2


Quantifying Force Generation
-------------------------------

Force generation is quantified with the strain energy and the contractillity. You have to select an area on
which these measures are to be calculated. You can do this by drawing a mask in clickpoints. In the top right
of the clickpoints window you can see a set of tools to draw mask and three preset types of masks. If you
don't see these tools, press F2.

.. hint:: **Tips for masks in clickpoints.**
   Select a mask and use the brush tool |brush| to draw it. You can
   in crease and decrease the size of the prush with the plus and minus keys. If you want to
   erase a part of a mask use the eraser tool |rubber|. Additionally you can fill holes in your mask with
   the bucket tool |bucket|. Mask types cannot overlap, which means that you erase one mask type when you
   paint over it with another type.

  .. |brush| image:: brush.png
  .. |rubber| image:: rubber.png
  .. |bucket| image:: bucket.png

The mask type used to calculate strain energy and contractility is called "foce measures". Select this mask and
draw a circle around all deformations and force that you think belong to the cell colony. The area you encircle
is typically large then the cell colony itself. You don't need to fill the area you have encircle. This is done
automatically. However, if you get the "no mask found in frame .." message, you should first make sure that there
is no gap in the circle that you drew. I drew the mask like this:


.. figure:: mask_force_measures.png
   :width: 600
   :alt: Main addon window




Its no big deal if your selection is a bit to big, but you should make sure not to include deformations and
force that do not originate from the cell colony.

You could no press start again, and the programm would generate a text file with contractillity and strain energy
for all frames. In order be a bit more organized and get all results one text file, we will first prepare
to analyze stresses in the cell sheet


Measuring Stresses
-------------------------------

The tractions that we have just
calculated are are not exactly localized to the cell colony. You can easily see that some tractions are
predicted to originate from outside of the cell area to