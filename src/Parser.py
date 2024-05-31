"""
@author Fabian Matten
"""
import argparse
import os
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

    def __init__(self):
        super().__init__(prog = "CoveragePlot Generator", description = "Generate Coverage Plots from coverage files for a given genome")
        self.io_group = self.parser.add_argument_group("IO for data stream")
        self.io_group.add_argument("-g", "--genome", help = "Input genome as fasta")
        self.io_group.add_argument("-c", "--coverage", help = "Coverage file for each position on genome")
        self.io_group.add_argument("-f", "-gff", type = str, help = "Gff input")
        self.io_group.add_argument("-l", "--loci", help = "Loci file")
        self.io_group.add_argument("-o", "--output_dir", help = "Output directory for plot files")

    def check_args(self, args):
        if not isinstance(args.genome, str):
            raise TypeError("Path the genome file must be a string")
        if not os.path.exists(args.genome):
            raise FileNotFoundError(f"{args.genome} file was not found")
        if args.genome.lower().endswith((".fa", ".faa", ".fna", ".fasta")):
            raise TypeError(f"{args.genome.lower} is not a fasta file")

        if not os.path.exists(args.coverage):
            raise FileNotFoundError(f"{args.coverage} file was not found")
        if not isinstance(args.coverage, str):
            raise TypeError("Path to coverage file must be a string")

        for path in [args.output_dir]:
            if not os.path.exists(path):
                raise OSError(f"{path} does not exist")