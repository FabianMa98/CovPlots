"""
@author Fabian Matten
"""
import argparse
import os
import logging

"""
TODO:

"""

class BaseParser:
    """
    Initialize Parser for task ahead
    """

    def __init__(self, prog = "", description = ""):
        super(BaseParser, self).__init__()
        self.parser = argparse.ArgumentParser(prog = prog, description = description, add_help = False, 
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
        self.parser.add_argument("-h", "--help", action = "help", help = "show this help message and exit")
    
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

    def __init__(self):
        super().__init__(prog = "CoveragePlot Generator", description = "Generate Coverage Plots from coverage files for a given genome")
        self.io_group = self.parser.add_argument_group("IO for data stream")
        self.io_group.add_argument("-c", "--coverage", help = "Coverage file for each position on genome")
        self.io_group.add_argument("-f", "--gff", type = str, help = "Gff input", dest = "gff")
        self.io_group.add_argument("-l", "--loci", help = "Loci file", required = False)
        self.io_group.add_argument("-o", "--output_dir", help = "Output directory for plot files")
        self.io_group.add_argument("-p", "--peaks", help = "MACS file with all peaks", required = False)

        # Flag arguments
        self.mode_group = self.parser.add_mutually_exclusive_group(required = True)
        self.mode_group.add_argument("-mp", "--mode_peak", help = "This will require MACS peaks file for final image",
                                     action = "store_true")
        self.mode_group.add_argument("-ml", "--mode_loci", help = "This will require loci tsv for final images",
                                     action = "store_true")
        
    def check_args(self, args):
        if not os.path.exists(args.coverage):
            raise FileNotFoundError(f"{args.coverage} file was not found")
        if not isinstance(args.coverage, str):
            raise TypeError("Path to coverage file must be a string")
        
        if not os.path.exists(args.gff):
            raise FileNotFoundError(f"{args.gff} could not be found")
        
        #if not os.path.exists(args.loci):
        #    raise FileNotFoundError(f"{args.loci} file was not found")
        
        if not os.path.exists(args.output_dir):
            os.mkdir(args.output_dir)
            if not os.path.isdir(args.output_dir):
                raise OSError(f"{args.output_dir} is not a directory")
            
        if not os.path.isdir(args.output_dir):
            raise OSError(f"{args.output_dir} is not a directory")
