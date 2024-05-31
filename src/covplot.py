from .Generator import CovPlot
from .Parser import GeneratorParser
from .gff import gff
from .loci import loci

def main():
    parser = GeneratorParser()
    args = parser.get_args()
    cov_gff = gff(args.gff).standardise_chromosomes()
    cov_loci = loci(args.loci).to_list()
    CovPlots = CovPlot(args.coverage, cov_loci, cov_gff)
    CovPlots.plot_all_loci()

if __name__ == "__main__":
    main()