import numpy as np
import pytest
import matplotlib.pyplot as plt
import numpy as np
import openpiv.filters
import scipy.fft

from openpiv.pyprocess import extended_search_area_piv
import openpiv.scaling
import openpiv.tools
import openpiv.validation
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pyTFM.utilities_TFM import suppress_warnings
from scipy.ndimage.filters import median_filter, gaussian_filter
from scipy.ndimage.filters import uniform_filter
from pyTFM.TFM_functions import ffttc_traction
from pyTFM.TFM_functions import get_xy_for_quiver
from pyTFM.TFM_functions import strain_energy_points
from pyTFM.TFM_functions import contractillity

def test_ffttc_traction_output_shape():
    # create test inputs
    u = np.zeros((10, 10))
    v = np.zeros((10, 10))
    pixelsize1 = 0.1
    pixelsize2 = 0.1
    young = 1e6

    # call the function
    tx_filter, ty_filter = ffttc_traction(u, v, pixelsize1, pixelsize2, young)

    # assert that the output shape is as expected
    assert tx_filter.shape == (10, 10)
    assert ty_filter.shape == (10, 10)

def test_ffttc_traction_output_type():
    # create test inputs
    u = np.zeros((10, 10))
    v = np.zeros((10, 10))
    pixelsize1 = 0.1
    pixelsize2 = 0.1
    young = 1e6

    # call the function
    tx_filter, ty_filter = ffttc_traction(u, v, pixelsize1, pixelsize2, young)

    # assert that the output is of the correct type
    assert isinstance(tx_filter, np.ndarray)
    assert isinstance(ty_filter, np.ndarray)

def test_ffttc_traction_output_values():
    # create test inputs
    u = np.zeros((10, 10))
    v = np.zeros((10, 10))
    pixelsize1 = 0.1
    pixelsize2 = 0.1
    young = 1e6

    # call the function
    tx_filter, ty_filter = ffttc_traction(u, v, pixelsize1, pixelsize2, young)

    # assert that the output is as expected
    assert np.allclose(tx_filter, 0)
    assert np.allclose(ty_filter, 0)

def test_ffttc_traction_input_shape():
    # create test inputs with incorrect shape
    u = np.zeros((10, 10, 10))
    v = np.zeros((10, 10))
    pixelsize1 = 0.1
    pixelsize2 = 0.1
    young = 1e6

    # assert that a ValueError is raised for incorrect input shape
    with pytest.raises(ValueError):
        ffttc_traction(u, v, pixelsize1, pixelsize2, young)

def test_ffttc_traction_input_type():
    # create test inputs with incorrect type
    u = np.zeros((10, 10))
    v = "not an array"
    pixelsize1 = 0.1
    pixelsize2 = 0.1
    young = 1e6

    # assert that a TypeError is raised for incorrect input type
    with pytest.raises(TypeError):
        ffttc_traction(u, v, pixelsize1, pixelsize2, young)
        
########get_xy_for_quiver
def test_get_xy_for_quiver():
    u = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    xs, ys = get_xy_for_quiver(u)
    
    # Test the shape of xs and ys
    assert np.shape(xs) == np.shape(u)
    assert np.shape(ys) == np.shape(u)
    
    # Test the values of xs
    for i in range(np.shape(u)[0]):
        assert np.array_equal(xs[i, :], np.arange(0, np.shape(u)[1], 1))
        
    # Test the values of ys
    for j in range(np.shape(u)[1]):
        assert np.array_equal(ys[:, j], np.arange(0, np.shape(u)[0], 1))
        
        
def test_contractility():
    tx = np.array([1, 2, 3, 4, 5])
    ty = np.array([1, 2, 3, 4, 5])
    pixelsize = 10
    mask = np.array([True, True, True, True, True])

    expected_output = (30.0, np.array([5.5, 11.0, 16.5, 22.0, 27.5]), np.array([5.5, 11.0, 16.5, 22.0, 27.5]), np.array([2.5, 2.5]))
    output = contractillity(tx, ty, pixelsize, mask)

    assert output == expected_output, f"Expected {expected_output}, but got {output}"

def test_strain_energy_points():
    u = np.array([1, 2, 3, 4, 5])
    v = np.array([1, 2, 3, 4, 5])
    tx = np.array([1, 2, 3, 4, 5])
    ty = np.array([1, 2, 3, 4, 5])
    pixelsize1 = 10
    pixelsize2 = 10

    expected_output = np.array([12.5, 25.0, 37.5, 50.0, 62.5])
    output = strain_energy_points(u, v, tx, ty, pixelsize1, pixelsize2)

    assert np.array_equal(output, expected_output), f"Expected {expected_output}, but got {output}"
