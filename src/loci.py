"" 
import pandas as pd

class loci:
    """
    Custom loci class
    """

    def __init__(self, path):
        self._path = path
        self._loci = None
        self.list = None

    @property
    def path(self):
        return self._path
    
    @path.setter
    def path(self, other_path):
        self._path = other_path

    #@property
    #def loci(self):
    #    return self._loci
    
    #@loci.setter
    #def loci(self, other_loci):
    #    """
    #    Load loci, which are assumed to be parsed in as tsv
    #    """
    #    if not isinstance(self.path, str):
    #        raise TypeError("Path should be a string")
    #    self._loci = pd.read_table(self.path, header = None, comment = '#')

    def load_loci(self):
        self._loci = pd.read_table(self._path, header = None, comment = '#')

    def to_list(self):
        locus_list = []
        loci = self.load_loci()
        for i in range(len(loci)):
            for j in range(len(loci[i])):
                locus_list.append(loci[i][j])

        self.list = locus_list

        return 