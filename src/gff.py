"" 
import pandas as pd
import numpy as np

from typing import Dict

class gff:
    """
    Own implementation for specific processing of gff files
    """

    # Constructor
    def __init__(self, path):
        self._path = path
        self._gff = self.load_gff(self._path)
        self._chr_dict = self.chromosome_maps()
        self._standardised_df = self.standardise_chromosomes()

    @property
    def path(self):
        return self._path
    
    @path.setter
    def path(self, other_path):
        self._path = other_path 

    @property
    def gff(self):
        return self._gff

    @gff.setter
    def gff(self, other_gff):
        self._gff = self.load_gff(other_gff)

    @property
    def standardised_df(self):
        return self._standardised_df
    
    @standardised_df.setter
    def standardised_df(self, other_df):
        self._standardised_df = other_df

    @property
    def chr_dict(self):
        return self._chr_dict

    @chr_dict.setter
    def chr_dict(self, other_dict):
        self._chr_dict = other_dict

    def load_gff(self, path):
        self._gff = pd.read_table(path, header = None, comment = '#')

        return self._gff

    def chromosome_maps(self) -> Dict[str, str]:
        """
        Compute chromsome maps

        Returns
        -------
            chr_dict: Dict[str, str]
                Each specific "omosome" as a key and Chr + i as its value
        """
        unique = self._gff[0].unique()
        self.chr_dict = {}
        for i in range(len(unique)):
            curr = i + 1
            self.chr_dict[unique[i]] = "Chr_" + str(curr)

        return self.chr_dict

    def standardise_chromosomes(self) -> pd.DataFrame:
        """
        Standardize chromosomes (especially relevant for bacteria 
        ToDo: Implement own MAP file for chromosomess
        This might have to moved to utils or standalone implementation of 
        gff file
        
        parameters
        ----------
            gff: pd.DataFrame
                GFF as pd DataFrame with column[0] as "omosomes" (for bacteria genome + plasmids)
            chr_map: Dict[str, str]
                chromsome map, for unique "omosomes" 

        returns
        -------
            standardised_df: pd.DataFrame


        """
        self._standardised_df = self.gff.replace(self.chr_dict)

        return self._standardised_df
