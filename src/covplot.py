from Generator import CovPlot
from Parser import GeneratorParser
from gff import gff
from loci import loci

def main():
    parser = GeneratorParser()
    args = parser.get_args()
    cov_gff = gff(args.gff)
    cov_gff_standard = cov_gff.standardised_df
    cov_loci = loci(args.loci).to_list()
    CovPlots = CovPlot(args.coverage, cov_loci, cov_gff_standard)
    CovPlots.plot_all_loci(args.output_dir)

if __name__ == "__main__":
    main()