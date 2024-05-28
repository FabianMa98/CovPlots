"" 
import pandas as pd

class loci:
    """
    Custom loci class
    """

    def __init__(self, path):
        self.path = path
        self.loci = None

    @property
    def path(self):
        return self.path
    
    @path.setter
    def path(self, other_path):
        self.path = other_path

    @property
    def loci(self):
        return self.loci
    
    @loci.setter
    def loci(self):
        """
        Load loci, which are assumed to be parsed in as tsv
        """
        if not isinstance(self.path, str):
            raise TypeError("Path should be a string")
        self.loci = pd.read_table(self.path, header = None, comment = '#')