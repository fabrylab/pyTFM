#!/usr/bin/env python

from setuptools import setup
import os
import platform



install_requires=['numpy' ,'cython','openpiv == 0.22.0', 'pyyaml', 
		      'scipy', 'scikit-image >= 0.18.1', 'matplotlib >= 2.1.2', 'tqdm', 'solidspy','clickpoints >= 1.9.6', "natsort"]

version='1.3.4' # adding a file version automatically
file_path = os.path.join(os.getcwd(),os.path.join("pyTFM","_version.py"))
with open(file_path,"w") as f:
	f.write("__version__ = '%s'"%version)

# installing openpiv from local .whl if we are on windows
local_openpiv = os.path.join(os.getcwd(),"local_dependencies","OpenPIV-0.20.8-cp37-cp37m-win_amd64.whl")
print('\033[92m'+"Identified operating system: ",platform.system() + '\033[0m')

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


setup(
    name='pyTFM',
    packages=['pyTFM'],
    version=version,
    description='traction force microscopy and monolayer stress microscopy of cell sheets and isolated patches',
    long_description_content_type="text/markdown",
    long_description=long_description,
    url='https://pytfm.readthedocs.io/',
    download_url = 'https://github.com/fabrylab/pyTFM.git',
    author='Andreas Bauer',
    author_email='andreas.b.bauer@fau.de',
    license= 'GPLv3',
    install_requires=install_requires,
    keywords = ['traction force microscopy', 'finite elements'],
    classifiers = [],
    include_package_data=True,
    )


	
