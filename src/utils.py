import pandas as pd
import numpy as np

def search_string(s, search):
    return search in str(s)

def get_max_height(coverage: np.ndarray) -> int:
    """
    get max height from coverage window for calculation of factor for arrow_line_loc
    DO NOT make this peak specific, as this software revolves around 
    """
    max_cov = coverage.max()

    return max_cov

def get_max_coordinates(coverage: np.ndarray):
    pass
