#!/usr/bin/env python

from setuptools import setup
setup(
    name='andreas_TFM_package',
    packages=['andreas_TFM_package'],
    version='0.1',
    description='traction force microscopy and FEM analysis of cell sheets',
    url='',
    download_url = '',
    author='Andreas Bauer',
    author_email='andreas.b.bauer@fau.de',
    license='',# ??
    install_requires=['numpy' ,'cython','openpiv==0.20.8  ','scipy', 'scikit-image', 'matplotlib', 'tqdm', 'solidspy'], #[clickpoints 18.3] could be rather problematic
    keywords = ['traction force microscopy','finite elements'],
    classifiers = [],
    )
