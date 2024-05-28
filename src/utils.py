import pandas as pd

def search_string(s, search):
    return search in str(s)

def get_max_height(peak_start: int, peak_end: int, coverage: pd.DataFrame) -> int:
    """
    get max height from coverage window for calculation of factor for arrow_line_loc
    DO NOT make this peak specific, as this software revolves around 
    """
    max_cov = coverage.iloc[peak_start: peak_end]["depth"].max()

    return max_cov