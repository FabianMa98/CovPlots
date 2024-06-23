import pandas as pd
import numpy as np

def search_string(s, search):
    """
    search specific string s in search string
    to be applied for in map function to each column of gff
    """
    return search in str(s)

def get_max_height(coverage: np.ndarray) -> int:
    """
    get max height from coverage window for calculation of factor for arrow_line_loc
    DO NOT make this peak specific, as this software revolves around 
    """
    max_cov = coverage.max()

    return max_cov

def get_max_coordinates(coverage: np.ndarray):
    """
    get coordinates from max_cov for automatic "peak"
    
    Parameters
    ----------
    coverage: np.ndarray
        chromosome filtered coverage
    """
    # Set coverage pos to index
    max_id = coverage["depth"].idxmax()
    peak_interval = [max_id - 100, max_id + 100]

    return peak_interval
