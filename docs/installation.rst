Installation
============

It is recommended to use this package with the `Anaconda Distribution <https://www.anaconda.com/distribution/>`_.

.. improve

pyTFM can be installed in several ways. The most straight forward way is to download the package from
`github <https://github.com/fabrylab/traction_force_microscopy>`_, unzip the files, open a terminal and navigate
into the "pyTFM-master". Depending on how you unzip there you might need
to go to "pyTFM-master/pyTFM-master". In this folder you should see a file called "setup.py". Next, install the package with pip:

    ``pip install -e .``

Pip is included in the Anaconda Distribution and should be available as a command line program out of the box.
The -e setting will enable you to edit package files, and the "." signifes that the "setup.py" file is
located in the folder that you are currently in.

.. formul

You can also download and install the package in one step using the command line program "git".
If you use the Anaconda Distribution, "git" can be installed with the command

    ``conda install git``

Now you can install pyTFM directly with the command

    ``pip install git+https://github.com/fabrylab/pyTFM.git``


This will download the script file into your Anaconda subdirectory, for example to
"anaconda3/lib/python3.7/site-packages/pyTFM"
This package includes an addon for the image display and annotation tool clickpoints.
Clickpoints is installed automatically if you follow the steps outline below. However, if
you want to access clickpoints via the "open with" program option for image files, you have
to use the command:

    ``clickpoints register``

in the terminal.


.. hint::
    Execute the command "clickpoints register" from the terminal
    to add clickpoints to the "open width" menu for image files. You can find detailed
    information the usage of clickpoints
    `here <https://clickpoints.readthedocs.io/en/latest/installation.html>`_.




Dependencies
---------------------
The following packages will be installed automatically if necessary:
numpy, cython, scipy , scikit-image, matplotlib, tqdm, `solidspy <https://pypi.org/project/solidspy/>`_,
clickpoints with a version higher then 1.9.0, `OpenPIV <http://www.openpiv.net/openpiv-python/>`_
version 0.20.8. Note that Clickpoints versions below 1.9.0 will fail to identify the pyTFM addon. Also note that OpenPIV
is still in developement, meaning that more recent versions of OpenPIV might not be compatible with pyTFM. Currently, if you are on
windows, openpiv is installed from a compiled .whl file included in the pyTFM packge. If you want to install
openpiv on your own you are likely to need the `Microsoft Visual C++ build tools
<https://visualstudio.microsoft.com/de/thank-you-downloading-visual-studio/?sku=BuildTools&rel=16>`_.


