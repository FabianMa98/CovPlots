"""
@author Fabian Matten
"""
import argparse
import logging

class BaseParser:
    """
    Initialize Parser for task ahead
    """

    def __init__(self, prog = "", description = ""):
        super(BaseParser, self).__init__()
        self.parser = argparse.ArgumentParser(prog = prog, description = description, add_help = False, 
        formater_class = argparse.ArgumentDefaultsHelpFormatter)
        self.parser.add_argument("-h", "--help", action = "Help", help = "show this help message and exit")
    
    def check_args(self, args):
        """
        Implement argument checking for proper use
        """
        pass

    def get_args(self):
        args = self.parser.parse_args()
        self.check_args(args)

        return args

        
class GeneratorParser(BaseParser):

    def __init__(self, prog = "CoveragePlot Generator", description = "Generate Coverage Plots from coverage files for a given genome")
