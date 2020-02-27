Installation
============

It is recommended to use this package with the `Anaconda Distribution <https://www.anaconda.com/distribution/>`_.
If you are using Windows you additionally need to install `Microsoft Visual C++ build tools
<https://visualstudio.microsoft.com/de/thank-you-downloading-visual-studio/?sku=BuildTools&rel=16>`_.

This package Includes an Addon for the image display and annotation tool clickpoints. You can find detailed
information about installation and usage of clickpoints
`here <https://clickpoints.readthedocs.io/en/latest/installation.html>`_.

.. improve

pyTFM can be installed in several ways. The most straight forward way is to download the package from
`github <https://github.com/fabrylab/traction_force_microscopy>`_. Unzip the files, open a terminal and navigate
into the "pyTFM" folder. Next, install the package with pip:

    ``pip install -e .``

Pip is included in the Anaconda Distribution and should be available as a command line program if Anaconda was
installed correctly.  The -e setting will enable you to edit package files,
which you can use to manipulating options for plotting.

.. formul

You can also download and install the package in one step using the command line program git.
If you use the Anaconda Distribution, git can be installed with the command

    ``conda install git``

Now you can install pyTFM directly with the command

    ``pip install git+https://github.com/fabrylab/pyTFM.git``

This will download all files to into the your Anaconda subdirectly, for example to
"anaconda3/lib/python3.7/site-packages/pyTFM"


Dependencies
---------------------
The following packages will be installed automatically if necessary:
numpy, cython, scipy , scikit-image, matplotlib, tqdm, `solidspy<https://pypi.org/project/solidspy/>`_,
clickpoints higher then 1.9.0, `openpiv <http://www.openpiv.net/openpiv-python/>`_
version 0.20.8. Note that Clickpoints versions below 1.9.0 will fail to identify the addon. Also note that Openpiv
is still in deveolpement, meaning that more recent versions of Openpiv might not be compatible.


