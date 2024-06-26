"""
@author: Fabian Matten
""" 
from functools import cached_property
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import ticker

from constants import WINDOW, MUTATION_SCALE, ALPHA, LINEWIDTH, FACTOR, OFFSET
from utils import search_string, get_max_height, get_max_coordinates
from gff import gff

from typing import List, Union
from pathlib import Path

# Define plt.Params, should be moved to constants
plt.rcParams["figure.figsize"] = (10,5)
plt.rcParams["figure.dpi"] = 300

class CovPlot:
    """Initialize coverage Plots for  
    
    TODO:
        - Major Refactoring needed:
            - Remove redundancies 
            - Refactor for memory efficiency 
        - Right now coverage is loaded everytime in a loop, which takes up a lot of memoery in mode_loci
        - Do not allow for negative genomic coordinates at window calculation
        - Add strand support 
        - Figure out better way for centralization of loci name and implement it in to mode_loci
    """

    def __init__(self, coverage: Union[str, Path], gff_coordinates: gff, loci: Union[List, None] = None, peaks: Union[Path, str, None] = None):
        """
        Constructor for Covplot class

        Parameters
        ----------
        coverage: pd.DataFrame
            coverage file from samtools output as DataFrame
            with cols: chromosome, pos, depth
        loci: List, None
            inherited from list class, but already converted
            to a List for ease of used
        gff_coordinates: gff
            Own gff implementation for proper filtering
            Can be standardized for chromosoome filtering
        peaks: Path, str, None
            depends on mode
        """
        self._coverage = coverage
        self._lociList = loci
        self._gff = gff_coordinates
        self._peaks = peaks             

        self.load_coverage()

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
    def load_peaks(peaks):
        # Right now support for MACS file
        peak_df = pd.read_table(peaks, header = None)
        chr_peak_dict = {}

        for _ , row in peak_df.iterrows():
            chr_peak_dict[row[1]] = row[0]

        return chr_peak_dict
    
    def get_peaks(self):
        if self._peaks:
            df = pd.read_table(self._peaks, header = None, comment = '#')
            self.loaded_peaks = df.rename(columns = {0: "chromosome", 1: "start", 2: "end", 3: "peak_name"})

            return self.loaded_peaks
    
    def load_coverage(self):
        """
        """
        coverage = pd.read_table(self._coverage, header = None)
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
    
    def filter_peaks_by_loci(self, peak_start, locus, window):
        # Find peaks associated with given loci
        peak_loci_association = {}
        # check if peak start is inside WINDOW of a given loci and chromosome
        if peak_start - window.index[0] >= 0 and peak_start - window.index[-1] <= 0:
            peak_loci_association[peak_start] = locus

            return peak_loci_association
        
    def get_window(self, locus):
        chromosome_cov = (self.loaded_coverage[self.loaded_coverage["chromosome"] == locus["chromosome"].iloc[0]]
                          .reset_index()
                          .drop(columns = ["index"])
                          )
        full_window = chromosome_cov.loc[int(locus.iloc[:,2].iloc[0] - WINDOW) : int(locus.iloc[:, 3].iloc[0]) + WINDOW]

        return full_window

    def get_loci_inside_peak_window(self, peak, gff):
        """Get all loci inside predefined peak window for a speficic peak

        TODO:
            - Change conditions for which loci are considered
                Right now if the end of loci is not inside peak window it wont be drawn
                but start up until window end could and should be made possible

        Parameters
        ----------
        peak: DataType still to be decided
            peak_start and end?
            Right now just peak start, but needs to be adjusted 
        gff: pd.DataFrame
            standardised dataframe from gff class

        Returns
        -------
            loci: pd.DataFrame
                DataFrame with all loci with following columns
        """
        loci = []
        peak_window = [peak["start"] - WINDOW, peak["end"] + WINDOW]
        gff_chr = gff[gff[0] == peak["chromosome"]]
        for _, row in gff_chr.iterrows():
            # At worst O(n) if peak near "last" gene
            # Changed condition to peak_start
            if row[3] - peak_window[0] >= 0 and row[3] <= peak_window[1]:
                # Chromosome, start, end, strand, gene_id
                candidate = [row[0], row[3], row[4], row[6], row.iloc[-1]]
                loci.append(candidate)
                # Change to start 
                if row[3] + WINDOW + 1 >= peak_window[1]:
                    break
        # Iterate through both loops, but think of something for more computationally feasible
        if loci:
            df = pd.DataFrame(loci)
            last_column = df.iloc[:, -1]
            names = last_column.str.split(";")
            replace_list = last_column.tolist()
            final = []
            for i, row in enumerate(names):
                for element in row:
                    if element.startswith("Name="):
                        name = element.split("=")[1]
                        final.append(name)
            final = pd.Series(final)
            # Create dict with values this replace wont work
            final_df = df.replace(to_replace = replace_list, value = final)
            final_df.rename(columns = {0: "chromosome", 1: "start", 2: "end", 3: "strand", 4: "name"}, inplace = True)
            #print(final_df)
    
            return final_df
        
        else:
            pass

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
                print("Printing line label")
                ax.plot(x_lower, y, color = line_color, label = line_label)
                ax.plot(x_middle, y, color = peak_color, label = peak_label)
                ax.plot(x_upper, y, color = line_color)
                fig.legend(("non peak region", "peak region"), loc="upper right")
            ax.plot(x, arrow_line, color = arrow_color, linewidth = LINEWIDTH)
            ax.add_patch(arrow_1)
            ax.set_ylabel("read coverage")
            ax.set_xlabel(f"genomic coordinates for {locus.iloc[0]["gene_id"]}")
            ax.set_xlim(left = locus.iloc[0]["start"] - WINDOW , right = locus.iloc[0]["end"] + WINDOW)
            yticks = ax.yaxis.get_major_ticks()
            yticks[-1].set_visible(False)
            ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))

            plt.savefig(output_path + str(locus["gene_id"].iloc[0]) + ".png")

        return
    
    def plot_peak(self, output_path, peak, loaded_coverage, arrow_color = "black", line_color = "grey", peak_color = "deeppink",
                   line_label = "non peak region", peak_label = "peak region", plot_label = False):
        """
        Plot coverage around spefic peak region, with loci being detected inside that peak window

        Parameters
        ----------
            output_path:
                output for plt images
            peak: pd.DataFrame
                peak dataframe: to be more spefific a pd.DataFrame row through df.iterrows()
        """
        df = self._gff.standardised_df
        self.gff_gene = df[df[2] == "gene"]
        self.peak_loci = self.get_loci_inside_peak_window(peak, self.gff_gene)#self._gff.standardise_chromosomes())
        
        peak_start = peak["start"]
        peak_end = peak["end"]
        chromosome_cov = (loaded_coverage[loaded_coverage["chromosome"] == peak["chromosome"]]
                          .reset_index()
                          .drop(columns = ["index"])
                          )
        genomic_coordinates = chromosome_cov.loc[peak_start - WINDOW : peak_end + WINDOW]["depth"]
        x = np.array(genomic_coordinates.index)
        y = np.array(genomic_coordinates)
        x_upper = np.ma.masked_where(x > peak_start, x)
        x_lower = np.ma.masked_where(x < peak_end, x)
        x_middle = np.ma.masked_where((x > peak_end) | (x < peak_start), x)

        arrow_line_loc = -(round(get_max_height(y) / FACTOR)) 
        arrow_line = np.array([arrow_line_loc] * len(x))
        fig, ax = plt.subplots()
        if not plot_label:
            ax.plot(x_lower, y, color = line_color)
            ax.plot(x_middle, y, color = peak_color)
            ax.plot(x_upper, y, color = line_color)
        else:
            ax.plot(x_lower, y, color = line_color, label = line_label)
            ax.plot(x_middle, y, color = peak_color, label = peak_label)
            ax.plot(x_upper, y, color = line_color)
            fig.legend(("non peak region", "peak region"), loc="upper right")
        if self.peak_loci is not None:
            print(self.peak_loci)
            for _, locus in self.peak_loci.iterrows():
                arrow = mpatches.FancyArrowPatch((locus["start"], arrow_line_loc), (locus["end"], arrow_line_loc), mutation_scale = MUTATION_SCALE, alpha = ALPHA,
                                        color = "black")
                text_loc = arrow_line_loc / 1.5 # Move this to constants again

                ax.plot(x, arrow_line, color = arrow_color, linewidth = LINEWIDTH)
                ax.add_patch(arrow)
                ax.text(x = ((locus["start"] + locus["end"]) / 2) - 150, y = text_loc, s = locus["name"], size = 8)
        
        ax.set_ylabel("read coverage")
        ax.set_xlabel(f"genomic coordinates for {peak["peak_name"]}")
        ax.set_xlim(left = peak["start"] - WINDOW , right = peak["end"] + WINDOW)
        yticks = ax.yaxis.get_major_ticks()
        yticks[-1].set_visible(False)
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))

        plt.savefig(output_path + str(peak_start) + ".png")
        
        return

    def plot_all_loci(self, output_path, loci):
        """
        Move plot logic to other function
        """
        for i, locus in enumerate(loci):
            # filter locus from gff
            filtered = self.filter_gff(locus)
            if i == 0:
                self.plot_locus(filtered, output_path = output_path, plot_label=True)
            else:
                self.plot_locus(filtered, output_path = output_path)
        
    def plot_all_peaks(self, output_path: Union[str, Path]):
        """
        Plot all peaks, get coverage and peaks

        Parameters
        ----------
        output_path: Union[str, Path]

        """
        loaded_coverage = self.load_coverage()
        self.peaks = self.get_peaks()
        i = 0
        for _, peak in self.peaks.iterrows():
            print(peak)
            if i == 0:
                self.plot_peak(output_path = output_path, peak = peak, loaded_coverage = loaded_coverage, plot_label = True)
                i += 1
            else:
                self.plot_peak(output_path = output_path, peak = peak, loaded_coverage = loaded_coverage, plot_label = False)


    def plotter(self, mode, output_path, loci = None):
        """
        plotter handler
        Check if peak exist for loci in loci List
        If not generate plot anyway but in manual mode
        """
        # If args.mode = loci or peaks
        # Plot locus/peaks
        if mode not in ["peak", "loci"]:
            raise OSError("Please provide acceptable mode: peak / loci")
        
        if mode == "peak":
            if self._peaks == None:
                raise TypeError("Please provide an appropiate MACS peak file")
            self.plot_all_peaks(output_path = output_path)
        elif mode == "loci":
            if not loci:
                self.plot_all_loci(output_path = output_path, loci = loci)
            else:
                raise FileNotFoundError("No Loci File was given although")
        else:
            raise Exception("Unkown mode detected. Please provide proper modes")