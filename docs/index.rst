.. pyTFM documentation master file, created by
   sphinx-quickstart on Wed Feb 26 10:45:41 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to pyTFM's documentation!
======================================
pyTFM is a python package that allows you to analyze force generation and stresses in cells,
cell colonies and confluent cell layers growing on a 2 dimensional surface. This package implements the procedures of
`Traction Force Microscopy <https://www.ncbi.nlm.nih.gov/pubmed/11832345>`_ and
`Monolayer Stress Microscopy <https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0055172>`_.
In addition to the standard measures for stress and force generation, it
also includes the line tension, a measure for the force transfer exclusively across cell-cell boundaries.
pyTFM includes an addon for the image annotation tool `clickpoints <https://clickpoints.readthedocs.io/en/latest/>`_
allowing you to quickly analyze and vizualize large datasets.


.. toctree::
   :caption: Contents:
   :maxdepth: 2

   installation
   introduction
   measures
   tutorials
   configuration




Note
--------
If you encounter any bugs or other problems please report them `here <https://github.com/fabrylab/pyTFM/issues>`_
or contact me at andreas.b.bauer@fau.de.


.. TODO: drift correction