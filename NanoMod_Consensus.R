###Script to merge results from NanoMod and generate a consensus putative modified positions list:

##Import libraries:
suppressMessages(library('plyr'))
suppressMessages(library('dplyr'))
suppressMessages(library('VennDiagram'))
suppressMessages(library('argparse'))
suppressMessages(library('ggplot2'))
suppressMessages(library('GenomicRanges'))
suppressMessages(library('stringr'))
suppressMessages(library('scales'))
suppressMessages(library('ggnewscale'))

##Import accessory functions:
source('Accessory_functions_consensusNanoMod.R')

##Argument parser:
#Create parser object
parser <- ArgumentParser()

#Define desired outputs:
#GLOBAL FEATURES:
parser$add_argument("-output", "--Output_name", type="character", help="Output(s) filenames.")
parser$add_argument("-fasta", "--Fasta_file", type="character", help="Genome fasta file.")
parser$add_argument("-ini_pos", "--Initial_position", type="integer", default=50, help="Initial position [default %(default)].")
parser$add_argument("-fin_pos", "--Final_position", type="integer", help="Final position.")
parser$add_argument("-plot", "--Plotting", action="store_true", help="Plot significant positions for all methods.")
parser$add_argument("-chr", "--Chr", type="character", help="Character to study.")
parser$add_argument("--MZS_thr", default=2.5, type="double", 
                    help="Modified Z-Score threshold for all results [default %(default)]")
parser$add_argument("-subset", "--Subsetting", action="store_true", help="Input data contain multiple chr.")

#EPINANO:
parser$add_argument("-Epi_Sample", "--Epinano_Sample", nargs=1, type="character", help="Path to Epinano features sample results.")
parser$add_argument("-Epi_IVT", "--Epinano_IVT", nargs=1, type="character", help="Path to Epinano features IVT results.")

#NANOPOLISH:
parser$add_argument("-NP_Sample", "--Nanopolish_Sample", nargs=1, type="character", help="Path to Nanopolish mean per position sample results.")
parser$add_argument("-NP_IVT", "--Nanopolish_IVT", nargs=1, type="character", help="Path to Nanopolish mean per position IVT results.")

#TOMBO:
parser$add_argument("-Tombo", "--Tombo_Sample", nargs=1, type="character", help="Path to Tombo pairwise comparison results.")

#NANOCOMPORE:
parser$add_argument("-Nanocomp", "--Nanocomp_Sample", nargs=1, type="character", help="Path to Nanocompore pairwise comparison results.")

parser$add_argument("--nanocomp_metric", default="GMM_logit_pvalue_context_4", type="character", 
                    help="Metric to use for Nanocompore analysis [default %(default)]")


#Get command line options, if help option encountered - print help and exit:
args <- parser$parse_args()

print('Processing data')
##EPINANO processing: 
epinano_data <- epinano_processing(args$Epinano_Sample, args$Epinano_IVT, args$Initial_position, args$Final_position, args$MZS_thr)

##NANOPOLISH processing: 
nanopolish_data <- nanopolish_processing(args$Nanopolish_Sample, args$Nanopolish_IVT, args$Initial_position, args$Final_position, args$MZS_thr)

##TOMBO processing: 
tombo_data <- tombo_processing(args$Tombo_Sample, args$thr_tombo_pos, args$thr_tombo_kmer, args$Initial_position, args$Final_position, args$MZS_thr)

##NANOCOMPORE processing:
nanocompore_data <- nanocomp_processing(args$Nanocomp_Sample, args$nanocomp_metric, args$thr_nanocomp, args$Initial_position, args$Final_position, args$MZS_thr)

##Plotting significant positions across all methods:
list_plotting <- list(epinano_data[[1]], nanopolish_data[[1]], tombo_data[[1]], nanocompore_data[[1]])
list_significant <- list(epinano_data[[2]], nanopolish_data[[2]], tombo_data[[2]], nanocompore_data[[2]])

##Only include positions from a specific chromosome - indicated by the user:
if (args$Subsetting==TRUE) {
  list_plotting <- subset_by_chr(list_plotting, args$Chr)
  list_significant <- subset_by_chr(list_significant, args$Chr)
}

if(args$Plotting==TRUE){
  print('Plotting across all methods')
  barplot_plotting(list_plotting, list_significant, args$Output_name)
}

print('Overlap analysis and Venn Diagram')
#Analysis of significant positions across methods: 
analysis_significant_positions(list_significant, list_plotting, args$Fasta_file, args$Output_name,  args$Initial_position, args$Final_position)
