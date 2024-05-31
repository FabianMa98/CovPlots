"""
@author: Fabian Matten
""" 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from .constants import WINDOW, MUTATION_SCALE, ALPHA, LINEWIDTH
from .utils import search_string
from .gff import gff
from typing import List, Union
from pathlib import Path

class CovPlot:
    """Initialize coverage Plots and spefic methods for 
    
    TODO:
        - Implement some check for genome length
        - Implement peak colorization
        - Implement coverage cleanup, only coverages for needed chromosomes
            should be loaded into memory
        - Should also load coverage as n-dimensional nparray (more memory efficient?)
        - Figure out scaling factor line location in final plot
            - DONE: Around 6.5 (Factor for DAP-seq manuscript)
    
    """

    def __init__(self, coverage: np.ndarray, loci: pd.DataFrame, gff_coordinates):
        self.coverage = coverage
        self.loci = loci
        self.gff = gff_coordinates

    @property
    def coverage(self):
        return self._coverage

    @coverage.setter
    def coverage(self, other_coverage):
        self._coverage = other_coverage

    @property
    def loci(self):
        return self._loci

    @loci.setter
    def loci(self, other_loci):
        self.loci = other_loci

    @property
    def peaks(self):
        return self.peaks
    
    @peaks.setter
    def peaks(self, other_peaks):
        self.peaks = other_peaks

    @property
    def gff_coordinates(self):
        return self.gff

    @gff_coordinates.setter
    def gff_coordinates(self, other_gff):
        self.gff = other_gff


    @staticmethod
    def load_loci(loci: Union[str, Path]) -> pd.DataFrame:
        """
        Loci should be parsed as a tsv
        """
        if not isinstance(loci, Union[str, Path]):
            raise TypeError("Table of loci needs to be path to file")
        df = pd.read_table(loci)

        return df
    
    @staticmethod
    def process_loci(loci: pd.DataFrame) -> List:
        pass

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
    
    @staticmethod
    def filter_gff(gff_df: pd.DataFrame, locus: str) -> pd.DataFrame:
        """
        Moved loading to another function as this can be called once
        Should i keep it as a static method?
        Need 
        """
        gene_df = gff_df.apply(lambda x: x.map(lambda s: search_string(s, "gene")))
        filtered_df = gene_df.loc[gene_df.any(axis=1)]
        # from filtered gff get start end columns
        df = filtered_df[filtered_df[len(filtered_df.columns - 1)]].str.contains(locus).select_dtype("int")
        df_start_end = df.rename(columns={df.columns[0]: "start", df.columns[1]: "end"})
        df_start_end.insert(0, "gene_id", locus)

        return df_start_end

    def standardise_chromosomes(gff: pd.DataFrame) -> pd.DataFrame:
        """
        Standardize chromosomes (especially relevant for bacteria 
        ToDo: Implement own MAP file for chromosomess
        This might have to moved to utils or standalone implementation of 
        gff file
        
        parameters
        ----------
            gff: pd.DataFrame
                already read and processed gff from specific organism

        returns
        -------
            standardised_df: pd.DataFrame


        """
        unique = gff[0].unique()
        chr_dict = {}
        for i in range(len(unique)):
            curr = i + 1
            chr_dict[unique[i]] = "Chr" + str(curr)

        standardised_df = gff.replace(chr_dict)

        return standardised_df 

    def plot_locus(self, locus, output_path, peak_coords, arrow_color = "black", line_color = "grey", peak_color = "deeppink",
                   line_label = "non peak region", peak_label = "peak region"):    
        # arrow_line_loc shoulde calculated dynamically based on maximum peak height
        # arrow_line_loc = -get_maximum_peak_height(window(locus)) / SOME_FACTOR 
        arrow_line_loc = -5
        # dont write text on arrows, this should be done manually
        # text_loc = - 4
        x = np.array(self.coverage.index)
        y = np.array(self.coverage.iloc[locus["start"] -  WINDOW : locus["end"] + WINDOW]["depth"])

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
        #ax.plot(x, arrow_line, color = arrow_color, linewidth = LINEWIDTH)
        ax.add_patch(arrow)
        ax.set_ylabel("read coverage")
        ax.set_xlabel("genomic coordinates")
        yticks = ax[0,0].yaxis.get_major_ticks()
        yticks[-1].set_visible(False)

        fig.savefig(output_path)
        plt.close(fig)


    def plot_all_loci(self):
        """
        Move plot logic to other function
        """
        gff = CovPlot.load_gff(self.gff)
        loci = CovPlot.load_loci(self.loci)
        #Pseudocode for locis (peak files from GEM)
        # peak_files = self.read_loci(self._loci)
        for locus in loci:
            # filter locus from gff
            # Create window around locus coordinates 
            # self.plot_locus(locus)
            self.plot(locus)


    def plot_by_coordinates(self):
        pass