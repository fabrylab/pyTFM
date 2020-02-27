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


To sum up: the strain energy is a measure for the total force generation, while the contractility is a
measure for the coordinated force generation.

.. image:: contractillity.png
  :width: 400
  :alt: Illustration for contractility



