"""
@author: Fabian Matten
""" 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from constants import WINDOW, MUTATION_SCALE, ALPHA, LINEWIDTH, FACTOR
from utils import search_string, get_max_height

from typing import List, Union
from pathlib import Path

class CovPlot:
    """Initialize coverage Plots and spefic methods for 
    
    TODO:
        - Implement some check chromsome localization of locus
        - Implement peak colorization
        - Implement coverage cleanup, only coverages for needed chromosomes
            should be loaded into memory
        - Should also load coverage as n-dimensional nparray (more memory efficient?)
        - Figure out scaling factor line location in final plot
            - DONE: Around 6.5 (Factor from DAP-seq manuscript)
    
    """

    def __init__(self, coverage: pd.DataFrame, loci: List, gff_coordinates: pd.DataFrame):
        self._coverage = coverage
        self._lociList = loci
        self._gff = gff_coordinates

    @property
    def coverage(self):
        return self._coverage

    @coverage.setter
    def coverage(self, other_coverage):
        self._coverage = other_coverage

    @property
    def loci(self):
        return self._lociList

    @loci.setter
    def loci(self, other_loci):
        self._lociList = other_loci

    @property
    def peaks(self):
        return self._peaks
    
    @peaks.setter
    def peaks(self, other_peaks):
        self._peaks = other_peaks

    @property
    def gff_coordinates(self):
        return self._gff

    @gff_coordinates.setter
    def gff_coordinates(self, other_gff):
        self._gff = other_gff
    

    @staticmethod
    def load_gff(gff: Union[str, Path]) -> pd.DataFrame:
        """
        Generalized gff processing for this purpose
        Unoptimized right now, as this would be called in a loop,
        and having to read a large table everytime it does it
        """
        # Do not infer headers, we filter them out anyway and dont need them
        # Headers in gff files rely too much on the person which created it
        # Remove comments, as some gffs have coments in headers
        if not isinstance(gff, Union[str, Path]):
            raise TypeError("Gff file should be given as a string or path to file location")

        df = pd.read_table(gff, header = None, comment= "#")
        return df
    
    def filter_gff(self, locus) -> pd.DataFrame:
        """
        Moved loading to another function as this can be called once
        Should i keep it as a static method?
        Need 
        """
        gene_df = self._gff.apply(lambda x: x.map(lambda s: search_string(s, "gene")))
        filtered_df = gene_df.loc[gene_df.any(axis=1)]
        # from filtered gff get start end columns
        df = filtered_df[filtered_df[len(filtered_df.columns - 1)]].str.contains(locus).select_dtype("int")
        df_start_end = df.rename(columns={df.columns[0]: "start", df.columns[1]: "end"})
        df_start_end.insert(0, "gene_id", locus)

        return df_start_end

    def plot_locus(self, locus, output_path, arrow_color = "black", line_color = "grey", peak_color = "deeppink",
                   line_label = "non peak region", peak_label = "peak region"):    
        # arrow_line_loc shoulde calculated dynamically based on maximum peak height
        # arrow_line_loc = -get_maximum_peak_height(window(locus)) / SOME_FACTOR 
        # dont write text on arrows, this should be done manually
        # text_loc = - 4
        x = np.array(self.coverage.index)
        y = np.array(self.coverage.iloc[locus["start"] -  WINDOW : locus["end"] + WINDOW]["depth"])
        arrow_line_loc = -(round(get_max_height(y) / FACTOR)) 

        #x_upper = np.ma.masked_where(x > peak_coords["peak_start"], x)
        #x_lower = np.ma.masked_where(x < peak_coords["peak_end"], x)
        #x_middle = np.ma.masked_where((x > peak_coords["peak_end"]) | (x < peak_coords["peak_start"]), x)

        arrow_line = np.array([arrow_line_loc] * len(x))
        
        arrow = mpatches.FancyArrow((locus["start"], arrow_line_loc), (locus["end"], arrow_line_loc),
                                    mutation_scale = MUTATION_SCALE, alpha = ALPHA, color = arrow_color)
        
        # Initialize plot
        fig, ax = plt.subplots()
        #ax.plot(x_lower, y, color = line_color, label = line_label)
        #ax.plot(x_middle, y, color = peak_color, label = peak_label)
        #ax.plot(x_upper, y, color = line_color)
        ax.plot(x, y, color = line_color)
        ax.plot(x, arrow_line, color = arrow_color, linewidth = LINEWIDTH)
        ax.add_patch(arrow)
        ax.set_ylabel("read coverage")
        ax.set_xlabel("genomic coordinates")
        yticks = ax.yaxis.get_major_ticks()
        yticks[-1].set_visible(False)

        fig.savefig(output_path)
        plt.close(fig)

    def plot_all_loci(self, output_path):
        """
        Move plot logic to other function
        """
        for locus in self._lociList:
            # filter locus from gff
            self.filtered_gff = self.filter_gff(locus)

            # Create window around locus coordinates 
            # self.plot_locus(locus)
            self.plot_locus(self.filtered_gff, output_path = output_path)
        