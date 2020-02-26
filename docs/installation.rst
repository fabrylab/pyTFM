Installation
============

It is recommended to use this package with the `Anaconda Distribution <https://www.anaconda.com/distribution/>`_.
If you are using Windows you additionally need to install `Microsoft Visual C++ build tools
<https://visualstudio.microsoft.com/de/thank-you-downloading-visual-studio/?sku=BuildTools&rel=16>`_.

This package Includes an Addon for the image display and annotation tool clickpoints. You can find detailed
information about installation and usage of clickpoints
`here<https://clickpoints.readthedocs.io/en/latest/installation.html>`_.
.. improve

pyTFM can be installed in several ways. The most straight forward way is to download the package from
`github<https://github.com/fabrylab/traction_force_microscopy>`_. Unzip the files, open a terminal and navigate
into the "pyTFM" folder. Next, install the package with pip:

    ``pip install -e .``

Pip is included in the Anaconda Distribution and should be available as a commannd line programm if Anaconda was
installed correctly.
The -e setting will enable you to edit package files, which you can use to manipulating options for plotting.
.. formul



Simply download or clone the repository. Unizp the files, then open a terminal and navigate to the folder tracktion_force_mircoscopy. Now perform a local installation with the command



You can also install this package directly from github. First you need to install git.
For Windows you can use

    ``conda install git``

then install this package by

    ``pip install git+https://github.com/fabrylab/tracktion_force_microscopy.git``

This will automatically add the an addon to clickpoints if you have clickpoints 1.9.0 or higher.
