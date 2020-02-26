#!/usr/bin/env python

from setuptools import setup
import os

version='1.1' # adding a version file automatically
file_path=os.path.join(os.getcwd(),os.path.join("andreas_TFM_package","_version.py"))
with open(file_path,"w") as f:
	f.write("__version__ = '%s'"%version)

setup(
    name='andreas_TFM_package',
    packages=['andreas_TFM_package'],
    version=version,
    description='traction force microscopy and FEM analysis of cell sheets',
    url='',
    download_url = '',
    author='Andreas Bauer',
    author_email='andreas.b.bauer@fau.de',
    license='',# ??
    install_requires=['numpy' ,'cython','openpiv==0.20.8','scipy', 'scikit-image', 'matplotlib', 'tqdm', 'solidspy','clickpoints >= 1.9.0'], #[clickpoints 18.3] could be rather problematic
    keywords = ['traction force microscopy','finite elements'],
    classifiers = [],
    include_package_data=True,
    )


	