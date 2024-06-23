from Generator import CovPlot
from Parser import GeneratorParser
from gff import gff
from loci import loci

def main():
    parser = GeneratorParser()
    args = parser.get_args()
    cov_gff = gff(args.gff)
    if args.mode_loci:
        if args.loci:
            cov_loci = loci(args.loci).to_list()
            CovPlots = CovPlot(args.coverage, cov_gff, loci = cov_loci)
            CovPlots.plotter(mode = "loci", output_path = args.output_dir, loci = cov_loci)
        else:
            raise FileNotFoundError("Loci mode detected, but no loci file was found")
    elif args.mode_peak:
        if args.peaks:
            CovPlots = CovPlot(args.coverage, cov_gff, peaks = args.peaks)
            CovPlots.plotter(mode = "peak", output_path = args.output_dir)
        else:
            raise FileNotFoundError("Peak mode detected, but no loci file was found")
if __name__ == "__main__":
    main()