"""
@author: Fabian Matten
""" 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from constants import WINDOW, MUTATION_SCALE, ALPHA, LINEWIDTH, FACTOR
from utils import search_string, get_max_height, get_max_coordinates
from gff import gff
from loci import loci

from typing import List, Union
from pathlib import Path

class CovPlot:
    """Initialize coverage Plots and spefic methods for 
    
    TODO:
        - Peak Colorization
    
    """

    def __init__(self, coverage: Union[str, Path], loci: List, gff_coordinates: gff):
        """
        Constructor for Covplot class

        Parameters
        ----------
        coverage: pd.DataFrame
            coverage file from samtools output as DataFrame
            with cols: chromosome, pos, depth
        loci: List
            inherited from list class, but already converted
            to a List for ease of used
        gff_coordinates: gff
            Own gff implementation for proper filtering
            Can be standardized for chromosoome filtering
        """
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
        """
        """
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

    def plot_locus(self, locus, output_path, peak_mode = "auto", arrow_color = "black", line_color = "grey", peak_color = "deeppink",
                   line_label = "non peak region", peak_label = "peak region", plot_label = False):    
        self.loaded_coverage = self.load_coverage()
        chromosome_cov = (self.loaded_coverage[self.loaded_coverage["chromosome"] == locus["chromosome"].iloc[0]]
                          .reset_index()
                          .drop(columns = ["index"])
                          )
        genomic_coordinates = chromosome_cov.loc[int(locus.iloc[:,2].iloc[0] - WINDOW) : int(locus.iloc[:, 3].iloc[0]) + WINDOW]["depth"]
        full_window = chromosome_cov.loc[int(locus.iloc[:,2].iloc[0] - WINDOW) : int(locus.iloc[:, 3].iloc[0]) + WINDOW]
        y = np.array(genomic_coordinates)
        x = np.array(genomic_coordinates.index)
        arrow_line_loc = -(round(get_max_height(y) / FACTOR)) 

        arrow_line = np.array([arrow_line_loc] * len(x))
            
        arrow_1 = mpatches.FancyArrowPatch((locus.iloc[0]["start"], arrow_line_loc), (locus.iloc[0]["end"], arrow_line_loc), mutation_scale = MUTATION_SCALE, alpha = 0.8,
                                    color = "black")
        # Initialize plot
        if peak_mode == "auto":
            # If peak mode is set to auto, a peak will be defined around +/- 100 bp around the summit at maximum coverage value
            # within window
            peak_coordinates = get_max_coordinates(full_window)
            print(peak_coordinates)
            x_upper = np.ma.masked_where(x > peak_coordinates[0], x)
            x_lower = np.ma.masked_where(x < peak_coordinates[1], x)
            x_middle = np.ma.masked_where((x > peak_coordinates[1]) | (x < peak_coordinates[0]), x)
            fig, ax = plt.subplots()
            if not plot_label:
                print("Here")
                ax.plot(x_lower, y, color = line_color)
                ax.plot(x_middle, y, color = peak_color)
                ax.plot(x_upper, y, color = line_color)
            else:
                ax.plot(x_lower, y, color = line_color, label = line_label)
                ax.plot(x_middle, y, color = peak_color, label = peak_label)
                ax.plot(x_upper, y, color = line_color)
            ax.plot(x, arrow_line, color = arrow_color, linewidth = LINEWIDTH)
            ax.add_patch(arrow_1)
            ax.set_ylabel("read coverage")
            ax.set_xlabel(f"genomic coordinates for {locus.iloc[0]["gene_id"]}")
            ax.set_xlim(left = locus.iloc[0]["start"] - WINDOW , right = locus.iloc[0]["end"] + WINDOW)
            yticks = ax.yaxis.get_major_ticks()
            yticks[-1].set_visible(False)

#arrow_1 = mpatches.FancyArrowPatch((second_operon["start"], arrow_line_loc), (second_operon["end"], arrow_line_loc), mutation_scale = mutation_scale, alpha = 0.8,
#                                    color = "black")
#ax[0,1].plot(x_lower, y, color = "grey", label = "non peak region")
#ax[0,1].plot(x_middle, y, color = "deeppink", label = "peak region")
#ax[0,1].plot(x_upper, y, color = "grey")
#ax[0,1].plot(x, y_line, color = "black", linewidth = 1.5)
#ax[0,1].add_patch(arrow_1)
#ax[0,1].set_ylabel("read coverage")
#ax[0,1].set_xlabel("genomic coordinates")
#ax[0,1].text(x = ((second_operon["start"] + second_operon["end"])/2 - 150) , y = text_loc, s = second_operon["old_locus_tag"], color = "black", size = 10)
#ax[0,1].set_xlim(left = second_operon["start"] - 500, right = second_operon["end"] + 500)

            plt.savefig(output_path + str(locus["gene_id"].iloc[0]) + ".png")

        elif peak_mode == "manual":
            # parse peak file into 
            pass

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
        