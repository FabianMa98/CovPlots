"""
@author: Fabian Matten
""" 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from constants import WINDOW, MUTATION_SCALE, ALPHA, LINEWIDTH, FACTOR
from utils import search_string, get_max_height
from gff import gff
from loci import loci

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
        - Implement coerage
    
    """

    def __init__(self, coverage: pd.DataFrame, loci: List, gff_coordinates: gff):
        self._coverage = coverage
        self._lociList = loci
        self._gff = gff_coordinates

        self.loaded_coverage = None

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

    def load_coverage(self):
        coverage = pd.read_table(self.coverage, header = None)
        self.loaded_coverage = coverage.rename(columns= {0: "chromosome", 1: "pos", 2: "depth"})

        return self.loaded_coverage
    
    def filter_gff(self, locus) -> pd.DataFrame:
        """
        Moved loading to another function as this can be called once
        Should i keep it as a static method?
        Need 
        """
        gene_df = self._gff.standardised_df.apply(lambda x: x.map(lambda s: search_string(s, "gene")))
        filtered_df = self._gff.standardised_df.loc[gene_df.any(axis=1)]
        # from filtered gff get start end columns
        df = filtered_df[filtered_df[len(filtered_df.columns)-1].str.contains(locus)]
        chromosome = (df
              .drop_duplicates(subset = 0)
              .iloc[:, 0]
              .iloc[0]
             )
        
        start_end = (df
             .select_dtypes("int")
             .drop_duplicates()
            )
        self.filtered_gff = start_end.rename(columns={start_end.columns[0]: "start", start_end.columns[1]: "end"})
        self.filtered_gff.insert(0, "gene_id", locus)
        self.filtered_gff.insert(0, "chromosome", chromosome)

        return self.filtered_gff

    def plot_locus(self, locus, output_path, arrow_color = "black", line_color = "grey", peak_color = "deeppink",
                   line_label = "non peak region", peak_label = "peak region"):    
        self.loaded_coverage = self.load_coverage()
        chromosome_cov = (self.loaded_coverage[self.loaded_coverage["chromosome"] == locus["chromosome"].iloc[0]]
                          .reset_index()
                          .drop(columns = ["index"])
                          )
        print(chromosome_cov)
        genomic_coordinates = chromosome_cov.loc[int(locus.iloc[:,2].iloc[0] - WINDOW) : int(locus.iloc[:, 3].iloc[0]) + WINDOW]["depth"]

        y = np.array(genomic_coordinates)
        x = genomic_coordinates.index
        arrow_line_loc = -(round(get_max_height(y) / FACTOR)) 
        print(arrow_line_loc)

        arrow_line = np.array([arrow_line_loc] * len(x))
        print(locus["start"].iloc[0], locus["end"].iloc[0])
            
        arrow_1 = mpatches.FancyArrowPatch((locus.iloc[0]["start"], arrow_line_loc), (locus.iloc[0]["end"], arrow_line_loc), mutation_scale = MUTATION_SCALE, alpha = 0.8,
                                    color = "black")
        # Initialize plot
        fig, ax = plt.subplots()
        ax.plot(x, y, color = line_color)
        ax.plot(x, arrow_line, color = arrow_color, linewidth = LINEWIDTH)
        ax.add_patch(arrow_1)
        ax.set_ylabel("read coverage")
        ax.set_xlabel(f"genomic coordinates for {locus.iloc[0]["gene_id"]} on chromosome {locus.iloc[0]["chromosome"]}")
        ax.set_xlim(left = locus.iloc[0]["start"] - WINDOW , right = locus.iloc[0]["end"] + WINDOW)
        yticks = ax.yaxis.get_major_ticks()
        yticks[-1].set_visible(False)

        print(output_path + str(locus["gene_id"].iloc[0]) + ".png")
        plt.savefig(output_path + str(locus["gene_id"].iloc[0]) + ".png")

        return

    def plot_all_loci(self, output_path):
        """
        Move plot logic to other function
        """
        for locus in self._lociList:
            # filter locus from gff
            filtered = self.filter_gff(locus)

            # Create window around locus coordinates 
            # self.plot_locus(locus)
            self.plot_locus(filtered, output_path = output_path)
        