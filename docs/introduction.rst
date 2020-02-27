Introduction to Traction Force Microscopy and Monolayer Stress Microscopy
================================================================================================

.. just a short intro for details refer to master thesis or to publication

Traction Force Microscopy and Monolayer Stress Microscopy are methods to measure the force generation and
internal forces of cells and cell colonies.
.. formul "force transfer"/ "internal forces"

A typical experiment starts by seeding cells on a substrate containing small flourescence labled beads. The cells adhere
to the substrate, start generating forces and consequently cause deformations in the substrate. These deformations can be
measured by tracking the labeled beads. In practice the beads are imaged two times: once when cells are attached to the
substrate, and once when the cells are removed from the substrate. In the first image the substrate is strained by the
cells, while the second image shows the completely relaxed substrate.

.. mention PIV

The deformations in the substrate can be used to calculate the traction forces that acted on the substrate surface.
This is done by Fourier Transformed Traction Force Microscopy (TFM)[1]_.
Any force that is applied from the cells to the substrate surface must be balanced by counteracting force
inside of the cells. Forces inside of materials are described by stress, that is forces
acting accross small surfaces. ( In a single point this surface can be oriented in different directions,
each orientation havig a diffrent force acting accross it.) Thinking about cells this, means that a cell
that generates two opposing forces at its ends, experiences a high stress between the points
of force generation.

The stress inside of cells and cell colonies is recovered by Monolayer Stress Microscopy (MSM) [2]_. In
MSM the cells are modeled as 2-Dimensional sheet. Following the argument of force-balance,
the traction forces calculated from TFM, with their sign inversed, are applied to the cell sheet. Then,
the stresses in the cell sheet are calculated with 2-Dimensional Finite Elements Methods








Tracking the labeled beads
allows to measure these deformation

These beads
allow to measure the deformations that the cells cause


.. [1] **Traction fields, moments, and strain energy that cells exert on their surroundings**
    James P. Butler, Iva Marija Tolić-Norrelykke, Ben Fabry, and Jeffrey J. Fredberg
    `Am J Physiol Cell Physiol 282: C595–C605, (2002) <https://www.physiology.org/doi/pdf/10.1152/ajpcell.00270.2001>`_
.. [2] **Monolayer Stress Microscopy: Limitations, Artifacts, and Accuracy of Recovered Intercellular Stresses**
    Dhananjay T. Tambe, Ugo Croutelle, Xavier Trepat, Chan Young Park, Jae Hun Kim, Emil Millet,
    James P. Butler, Jeffrey J. Fredberg
    `PLOS ONE 8(2): e55172 (2013) <https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0055172>`_
