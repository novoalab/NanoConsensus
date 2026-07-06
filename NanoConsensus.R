#!/usr/local/bin/Rscript --vanilla

##Import accessory functions:
source(Sys.which("Accessory_functions_consensusNanoMod.R"))

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
suppressMessages(library('ggrepel'))
suppressMessages(library('gtable'))
suppressMessages(library('xfun'))


##Argument parser:
#Create parser object
parser <- ArgumentParser()

#Define desired outputs:
#GLOBAL FEATURES:
parser$add_argument("-output", "--Output_name", type="character", help="Output(s) filenames.")
parser$add_argument("-fasta", "--Fasta_file", type="character", help="Genome fasta file.")
parser$add_argument("-ini_pos", "--Initial_position", type="integer", default=50, help="Initial position [default %(default)].")
parser$add_argument("-fin_pos", "--Final_position", type="integer", help="Final position.")
parser$add_argument("-chr", "--Chr", type="character", help="Character to study.")
parser$add_argument("--MZS_thr", default=5, type="double",
                    help="Modified Z-Score threshold for all results [default %(default)]")
parser$add_argument("--NC_thr", default=5, type="double",
                    help="NanoConsensus score threshold for all results [default %(default)]")
parser$add_argument("-exclude", "--Exclude", nargs='+', type="integer", help="Exclude these positions from the analysis (SNPs) - it will exclude the 17-mer.")
parser$add_argument("--model_score", default="global", type="character",
                    help="Model used to calculate NanoConsensus score [default %(default)]")
parser$add_argument("--coverage", default=1, type="integer",
                    help="Minimum coverage per position to be included in the analysis [default %(default)]")
parser$add_argument("--bed", help="Path to RNA modification annotation (*.bed)")
parser$add_argument("--ablines", action='store_true', help="Plot reported modified sites from the bed file.")
parser$add_argument("--extended_outputs", action='store_true', help="Generate extended outputs (kmer_tracks and Raw_kmers file).")

#bedRmod format metadata:
parser$add_argument("--organism", default="", type="character", help="Organism name for bedRmod header")
parser$add_argument("--assembly", default="", type="character", help="Genome assembly for bedRmod header")
parser$add_argument("--annotation_source", default="", type="character", help="Annotation source for bedRmod header")
parser$add_argument("--annotation_version", default="", type="character", help="Annotation version for bedRmod header")
parser$add_argument("--basecalling", default="", type="character", help="Basecalling method for bedRmod header")
parser$add_argument("--bioinformatics_workflow", default="", type="character", help="Bioinformatics workflow for bedRmod header")
parser$add_argument("--experiment", default="", type="character", help="Experiment description for bedRmod header")
parser$add_argument("--strand", default="", type="character", help="Strand information for bedRmod data")

#EPINANO:
parser$add_argument("-Epi_Sample", nargs=1,  default="", type="character", help="Path to Epinano features sample results.")
parser$add_argument("-Epi_IVT", nargs=1, default="", type="character", help="Path to Epinano features IVT results.")

##Bedmethyl inputs:
#BaseQ:
parser$add_argument("-BaseQ", nargs=1, default="", type="character", help="Path to baseQ pairwise comparison results.")

#nanoRMS SI:
parser$add_argument("-nanoRMS_SI", nargs=1, default="", type="character", help="Path to nanoRMS (SI) pairwise comparison results.")

#nanoRMS DT:
parser$add_argument("-nanoRMS_DT", nargs=1, default="", type="character", help="Path to nanoRMS (DT) pairwise comparison results.")

#nanoRMS SD:
parser$add_argument("-nanoRMS_SD", nargs=1, default="", type="character", help="Path to nanoRMS (SD) pairwise comparison results.")

#Get command line options, if help option encountered - print help and exit:
args <- parser$parse_args()

##Create and update log file:
write('NanoConsensus - v 2.0', file = paste("NanoConsensus_", args$Output_name,".log", sep=""))
write(paste('Analysing sample: ',args$Output_name, sep=""), file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)
write(paste('Minimum coverage: ',args$coverage, sep=""), file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)
write(paste('Transcript analysed: ',args$Chr, " - from ", args$Initial_position, " to ", args$Final_position, sep=""), file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)
write(paste('Z score threshold: ',args$MZS_thr, sep=""), file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)
write(paste('NanoConsensus score threshold: ',args$NC_thr, "*median(NanoConsensus Score)", sep=""), file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)
write('Step 1: Processing data from individual softwares', file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)

##EPINANO processing:
if (file.exists(args$Epi_Sample) && (file.exists(args$Epi_IVT))) {
  epinano_data <- epinano_processing(args$Epi_Sample, args$Epi_IVT, args$Initial_position, args$Final_position, args$MZS_thr, args$Chr, args$Exclude, args$coverage)

} else {
  epinano_data <- list(data.frame(Reference= character(), Position=integer(), Difference=double(), Feature=character()),data.frame(Reference= character(), Position=integer(), Difference=double(), Feature=character()))

}

##Bedmethyl processing:
files_to_parse <- c(args$BaseQ, args$nanoRMS_SI, args$nanoRMS_DT, args$nanoRMS_SD)
features_to_parse <- c('baseQ', 'nanoRMS_SI', 'nanoRMS_DT', 'nanoRMS_SD')

for (i in seq(1,length(files_to_parse))){
  print(paste(features_to_parse[i],'data', sep="_"))
  if (file.exists(files_to_parse[i])) {
    assign(paste(features_to_parse[i],'data', sep="_"), bedmethyl_processing(files_to_parse[i], args$Initial_position, args$Final_position, args$MZS_thr, args$Chr, args$Exclude, args$coverage, features_to_parse[i]))

  } else {
    assign(paste(features_to_parse[i],'data', sep="_"), list(data.frame(Reference= character(), Position=integer(), Difference=double(), Feature=character()),data.frame(Reference= character(), Position=integer(), Difference=double(), Feature=character())))

  }

}

##DATA PROCESSING:
#Generate list with all positions and significant positions respectively:
list_plotting <- list(epinano_data[[1]], baseQ_data[[1]], nanoRMS_SI_data[[1]], nanoRMS_DT_data[[1]], nanoRMS_SD_data[[1]])
list_significant <- list(epinano_data[[2]], baseQ_data[[2]], nanoRMS_SI_data[[2]], nanoRMS_DT_data[[2]], nanoRMS_SD_data[[2]])

#Extract first available coverage data (epinano, baseQ, or nanoRMS)
coverage_data <- NULL
coverage_sources <- list(
  if (length(epinano_data) >= 3) epinano_data[[3]] else NULL,
  if (length(baseQ_data) >= 3) baseQ_data[[3]] else NULL,
  if (length(nanoRMS_SI_data) >= 3) nanoRMS_SI_data[[3]] else NULL,
  if (length(nanoRMS_DT_data) >= 3) nanoRMS_DT_data[[3]] else NULL,
  if (length(nanoRMS_SD_data) >= 3) nanoRMS_SD_data[[3]] else NULL
)

for (cov in coverage_sources) {
  if (!is.null(cov) && nrow(cov) > 0) {
    coverage_data <- cov
    break
  }
}

#Exit if no coverage data is available
if (is.null(coverage_data)) {
  write("ERROR: No coverage data found in any input files", file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)
  stop("No coverage data available - cannot proceed with analysis")
}

#If there is annotation, process it:
if (length(args$bed)!=0){
  annotation <- process_bed(args$bed, args$Chr)
} else {
  annotation <- c()
}

#Create Z-Scores plotting object:
write('Step 2: Plotting ZScores from individual softwares', file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)
barplot_4soft <- barplot_plotting(list_plotting, list_significant, args$Output_name, args$MZS_thr, args$Initial_position, args$Final_position, annotation, args$ablines)

##Analysis of SIGNIFICANT POSITIONS across methods:
write('Step 3: Overlapping analysis', file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)
analysis_significant_positions(list_significant, list_plotting, args$Fasta_file, paste(args$Output_name, args$Chr, sep='_'),  args$Initial_position, args$Final_position, args$MZS_thr, args$NC_thr, args$model_score, barplot_4soft, annotation, args$ablines, args$Chr, coverage_data, args$organism, args$assembly, args$annotation_source, args$annotation_version, args$basecalling, args$bioinformatics_workflow, args$experiment, args$extended_outputs, args$strand)
