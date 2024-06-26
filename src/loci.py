"" 
import pandas as pd

class loci:
    """
    Custom loci class
    """

    def __init__(self, path):
        self._path = path
        self._loci = self.load_loci()
        self.list = None

    @property
    def path(self):
        return self._path
    
    @path.setter
    def path(self, other_path):
        self._path = other_path

    @property
    def loci(self):
        return self._loci
    
    @loci.setter
    def loci(self, other_loci):
        self._loci = other_loci

    def load_loci(self):
        self._loci = pd.read_table(self._path, header = None, comment = '#')
        
        return self._loci

    def to_list(self):
        locus_list = []
        loci = self._loci.values.tolist()
        for i in range(len(loci)):
            for j in range(len(loci[i])):
                locus_list.append(loci[i][j])

        self._list = locus_list

        return self._list
    
    def to_df(self):
        pass