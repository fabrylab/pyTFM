[![Python package](https://github.com/martinschatz-cz/pyTFM/actions/workflows/python-package.yml/badge.svg?event=push)](https://github.com/martinschatz-cz/pyTFM/actions/workflows/python-package.yml)

[![Python 3.7](https://img.shields.io/badge/python-3.7-green.svg)]() [![Python 3.8](https://img.shields.io/badge/python-3.8-red.svg)]() [![Python 3.9](https://img.shields.io/badge/python-3.9-red.svg)]()

## Readme

pyTFM is a python package that allows you to analyze force generation and stresses in cell colonies and confluent cell layers growing on a 2 dimensional surface. This package implements the procedures of [Traction Force Microscopy](https://www.ncbi.nlm.nih.gov/pubmed/11832345) and [Monolayer Stress Microscopy](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0055172). In addition to the standard measures for stress and force generation, it
also includes the line tension, a measure for the force transfer exclusively across cell-cell boundaries. 
pyTFM includes an addon for the image annotation tool [clickpoints](https://clickpoints.readthedocs.io/en/latest/) allowing you to quickly analyze and vizualize large datasets.

Please refer to the [Documentation](https://pytfm.readthedocs.io/en/latest/) of pyTFM for detailed instructions on installation and usage.

## Conda enviroment creation
You need to create python 3.6 enviroment. 
```
conda create --name pyTFM python=3.6
conda activate pyTFM
```

Note: python 3.7 works too, but clickpoints need python 3.6.

pyTFM package install
```
pip install git+https://github.com/fabrylab/pyTFM.git
```

If you want to use jupyterlab
```
pip install jupyterlab
jupyter-lab
```
Sometimes it is necessary to reinstall pyTFM again through the JupyterLab because of scikit-image package.
