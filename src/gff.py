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
        self.path = path

    @property
    def path(self):
        return self.path
    
    @path.setter
    def path(self, other_path):
        self.path = other_path

    def load_gff(self):
        self.gff = pd.read_table(self.path, header = None, comment = '#')

    def standardise_chromosomes(self, chr_map: Dict[str, str]) -> pd.DataFrame:
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
        self.standardised_df = self.gff.replace(chr_map)

        return self.standardised_df
    
    def chromosome_maps(self) -> Dict[str, str]:
        """
        Compute chromsome maps

        Returns
        -------
            chr_dict: Dict[str, str]
                Each specific "omosome" as a key and Chr + i as its value
        """
        unique = self.gff[0].unique()
        chr_dict = {}
        for i in range(len(unique)):
            curr = i + 1
            chr_dict[unique[i]] = "Chr" + str(curr)

        return chr_dict