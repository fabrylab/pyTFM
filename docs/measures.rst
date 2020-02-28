Typical Measures for Force Generation and Stresses in Cells
===============================================================
There is a number of diffrent ways to measure force generation and stresses. Here you can find
an overview of all quantities that can be calculated with this program.

Deformations in the Substrate
-------------------------------
The simplest way is to sum up all deformations on the substrate surface that the cells are attached to.
The deformations depend on the mechanical properties of the substrate, which means that this is not
possible to compare results, when diffrent substrates have been used.

Strain Energy
-----------------
The strain energy is the total work that cells have put in (fromul) to deform their substrate. It is
defined as :math:`\frac{1}{2} \int \vec{d} \times \vec{f}`, where :math:`\vec{d}` and :math:`\vec{f}`
are the deformation and the traction force vectors. That means that a high strain Energy only needs a
local alignement of both vectors.

.. actually force vector...


Contractillity
---------------
The contractillity is defined as the sum of the projection of all traction forces towards one point,
called the force epicenter. As such, the contractillity is high if all force are already orientated towards
one central point. Locally opposing forces and force not pointing in the general direction of force do
not contribute to the contractillity. A cell or cell colony wich is able to coordinate
its force generation in such a way that the forces seem to originate from a single point, can achieve
a high contractillity, while expending a comparably small aount of Strain Energy.

This is further illustrated in the Figure below. Case A represents a cell colony with two cells
with high coordination of force generation and Case B represents a cell with coordination. In case B
each cell generates contractile forces on its own. Accordingly case B has a lower contractility if we
assume equal strain energy in Case A and B.

.. figure:: contractillity.png
  :width: 400
  :alt: Illustration for contractility


To sum up: the strain energy is a measure for the total force generation, while the contractility is a
measure for the coordinated force generation.

Average Normal and Shear Stress
--------------------------------

.. mention maximum principal stress...
    TODO: implement maximum pricipal stress

Stress describes the forces that are transferred inside of a cell or cell sheet. For any given point
in the cell sheet the stress is defined by a tensor with 4 components. Each component represents forces
that would act on the edges of a square cut out of the cell sheet as shown bellow.

.. figure:: stress_tensor.png
  :width: 400
  :alt: Illustration for contractility

  :math:`\sigma_{xx}` and :math:`\sigma_{yy}`: Normal Stresses  in x and y direction.
  :math:`\sigma_{xy}` and :math:`\sigma_{yx}`: Shear stresses.


We can distinguish between two types stress: the shear stress, wich implies forces that act parallel
to the edges of this square, and normal stress, that act perpendicular to the edges of this square.
Due to geometric reasons both the shear components of the stress tensor must be identical. This
is not the case for the normal components. Since it is not important whether normal stress
comes predominately from the x or y directions, it is more usefull to calculate the mean of
both components. This leaves two stresses: the shear stress and the mean normal stress. These stresses
can then be averaged over the whole cell colony area.

.. note::
    The mean normal stress can be either negative, indicating a compressive stress, or positive
    ,indicating a tensile stress. A typical cell experiences tensile stress, as it is pulls
    from its edges on the underlying substrate. THe shear stress can also have a positive or a
    negative signe. However there is no usefull interpretation in this case.
.. TODO: decide weather to report the mean stress or the mean mean abs. value of the stress

Distribution of Stresses
--------------------------------

The distribution of Stresses can be described by the Coefficient of Variation (CV), that is
the standard deviation normalized with the mean, of mean normal or shear stress.

Forces acting accros Cell Boundaries
--------------------------------------
As hinted above, the stress tensor can be used to calculate forces that act accross edges in the
cell colony. This force is called the line tension, and has a straight forward interpretation:
Imagine you were to actually cut the cell sheet along the boundary of two cells. If the cells
continue to generate force the edges of this cut would drift appart or start overalpping
as you have just cut the material holding both edges together. In order to hold both edges in place
as they were before you cut them, you need to apply a force to them. This force, normalized
by the length of the cut is the line tension.

The line tension a vector with x and y components. Similar to stresses it can be split in a
shear component (the force acting parallel to the cut) and a normal component (the force acting
perpendicular to the cut). Both contribute to the magnitude (length of the line tension vector)
of the line tension.

.. note::
    Similar to normal stresses, the normal component of the line tension can be negative or
    positive, indicating that the two sides of the edge, along wich the line tension
    was calculated, are pushed together or pulled appart. The shear component of the line
    tension lacks such an interpretation. And of course the magnitude of the line tension
    can only be positive.

.. TODO: decide weather to return magnitude or not

..  Imagine the cell
    sheet to be as a sheet of paper
    If you where
    to cut a small
    figure:: line_tension.png
    width: 400
      :alt: Illustration for contractility



