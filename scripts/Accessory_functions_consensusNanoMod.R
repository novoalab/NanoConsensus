#!/usr/bin/env Rscript

###Script which contains multiple R functions used to generate consensus putative modified positions from NanoMod results.

#Read gzipped or flat files
read_gzipped <- function(input_file) {
	input <- input_file
	if (file_ext(input_file)=="gz") {
		input = gzfile(input_file)
	}
	return(input)
}

read_tab_file <- function(input_file) {
    file_content <- read.delim(read_gzipped(input_file))
	return(file_content)
}

read_csv_file <- function(input_file) {
    file_content <- read.csv(read_gzipped(input_file), stringsAsFactors = FALSE)
	return(file_content)
}

#Processing Epinano results:
epinano_processing <- function(sample_file, ivt_file, initial_position, final_position, MZS_thr, chr, exclude_SNP, Coverage) {

  #Import and clean data:
  sample <- read_csv_file(sample_file)
  sample <- subset(sample, cov>=Coverage)
  sample <- subset(sample, pos>=initial_position)
  sample <- subset(sample, pos<=final_position)
  sample$reference <- paste(sample$X.Ref, sample$pos, sep='_')
  sample$Difference <- as.numeric(sample$mis)+as.numeric(sample$ins)+as.numeric(sample$del)

  #Extract coverage information before column filtering
  sample_coverage <- data.frame(Position = sample$pos, Coverage = sample$cov, stringsAsFactors = FALSE)

  sample <- sample[,c(1,2,14,13)]
  colnames(sample) <- c('Reference', 'Position', 'Difference_sample', 'Merge')

  ivt <- read_csv_file(ivt_file)
  ivt <- subset(ivt, cov>=Coverage)
  ivt <- subset(ivt, pos>=initial_position)
  ivt <- subset(ivt, pos<=final_position)
  ivt$reference <- paste(ivt$X.Ref, ivt$pos, sep='_')
  ivt$Difference <- as.numeric(ivt$mis)+as.numeric(ivt$ins)+as.numeric(ivt$del)

  #Extract coverage information from ivt before column filtering
  ivt_coverage <- data.frame(Position = ivt$pos, Coverage = ivt$cov, stringsAsFactors = FALSE)

  ivt <- ivt[,c(1,2,14,13)]
  colnames(ivt) <- c('Reference', 'Position', 'Difference_IVT', 'Merge')

  if (nrow(sample)!=0 && nrow(ivt)!=0) {
    #Combine coverage from sample and ivt by calculating median per position
    combined_coverage <- rbind(sample_coverage, ivt_coverage)
    coverage_data <- aggregate(Coverage ~ Position, data = combined_coverage, FUN = median, na.rm = TRUE)

    #Join both dataframes and clean unecessary columns:
    plotting_positions <- join(sample, ivt, by="Merge")
    plotting_positions <- subset(plotting_positions, Reference == chr)

    #Exclude SNPs and 10 positions before and after (21mer):
    if (length(exclude_SNP)!=0) {
      excluded_positions <- c()

      for (single_position in exclude_SNP){
        excluded_positions <- c(excluded_positions, seq(single_position-10,single_position+10))
      }

      plotting_positions <- subset(plotting_positions, !Position %in% unique(excluded_positions))
    }

    plotting_positions$Difference <- abs(as.numeric(plotting_positions$Difference_sample) - as.numeric(plotting_positions$Difference_IVT))
    plotting_positions$Feature <- "Epinano"
    plotting_positions <- plotting_positions[,c(4,2,8,9)]

    #Calculate the threshold:
    threshold <- median(plotting_positions$Difference, na.rm = TRUE)

    #Calculate fold change and re-order:
    plotting_positions$Score <- plotting_positions$Difference/threshold
    plotting_positions$Modified_ZScore <- (plotting_positions$Score-median(plotting_positions$Score, na.rm = TRUE))/sd(plotting_positions$Score, na.rm = TRUE)

    plotting_positions <- plotting_positions[,c(1,2,5,4,6)]
    colnames(plotting_positions) <- c('Reference', 'Position', 'Score', 'Feature', 'Modified_ZScore')

    #Extract significant positions based on the specific threshold:
    significant_positions <- subset(plotting_positions, Modified_ZScore>MZS_thr)

    #Filter coverage to match plotting_positions
    coverage_data <- coverage_data[coverage_data$Position %in% plotting_positions$Position, ]

  } else {
    plotting_positions <- data.frame(Reference= character(), Position=integer(), Difference=double(), Feature=character())
    significant_positions <- data.frame(Reference= character(), Position=integer(), Difference=double(), Feature=character())
    
  }

  return(list(plotting_positions, significant_positions, coverage_data))
}

bedmethyl_processing <- function(sample_file, initial_position, final_position, MZS_thr, chr, exclude_SNP, Coverage, feature) {
  #Import data:
  sample <- read_tab_file(sample_file)

    if (nrow(sample)>0) {
      #Apply some filters and labels:
      sample$Feature <- feature
      sample <- sample[sample[,10]>=Coverage,c(1,2,10,5,12)]
      colnames(sample) <- c('Chr', 'Position', 'Coverage', 'KS', 'Feature')

    #If some positions arent reported by baseQ/nanoRMS, create data for all positions with threshold values:
    all_positions <- data.frame(
      Chr = chr,
      Position = seq(initial_position, final_position)

    )

    # Merge full positions with sample data:
    merged_sample <- merge(all_positions, sample, by = c("Chr", "Position"), all.x = TRUE)

    #Fill NA values with threshold values:
    merged_sample$Coverage[is.na(merged_sample$Coverage)] <- Coverage
    min_value <- min(abs(merged_sample$KS), na.rm=TRUE)
    merged_sample$KS[is.na(merged_sample$KS)] <- min_value

    merged_sample$Feature[is.na(merged_sample$Feature)] <- feature

    #Subset:
    merged_sample <- subset(merged_sample, Chr == chr)
    merged_sample <- subset(merged_sample, Position >= initial_position)
    merged_sample <- subset(merged_sample, Position <= final_position)
    merged_sample$Reference <- paste(merged_sample$Chr, merged_sample$Position, sep='_')

    #Extract coverage information before column filtering
    coverage_data <- data.frame(Position = merged_sample$Position, Coverage = merged_sample$Coverage, stringsAsFactors = FALSE)

    #Exclude SNPs and 10 positions before and after (21mer):
    if (length(exclude_SNP)!=0) {
      excluded_positions <- c()

      for (single_position in exclude_SNP){
        excluded_positions <- c(excluded_positions, seq(single_position-10,single_position+10))
      }

      merged_sample <- subset(merged_sample, !Position %in% unique(excluded_positions))
    }

    #Transform metric:
    merged_sample$KS <- abs(merged_sample$KS)
    threshold_position <- median(merged_sample$KS, na.rm = TRUE)

    #Calculate fold change:
    merged_sample$Score <- merged_sample$KS/threshold_position
    #sample$Modified_ZScore <- (sample$pvalue_KS-median(sample$pvalue_KS, na.rm = TRUE))/sd(sample$pvalue_KS, na.rm = TRUE)
    merged_sample$Modified_ZScore <- (merged_sample$Score-median(merged_sample$Score, na.rm = TRUE))/sd(merged_sample$Score, na.rm = TRUE)

    #Due to the inputted values, transform all Zscores = 0 into NA not to affect NanoConsensus score calculation:
    merged_sample$Modified_ZScore[merged_sample$Modified_ZScore == 0 & merged_sample$Score == 1] <- NA

    #Filter columns to get data in plotting format:
    plotting_positions <- merged_sample[,c(6,2,7,5,8)]

    #Extract significant positions:
    significant_positions <- subset(plotting_positions, Modified_ZScore>MZS_thr)

    #Filter coverage to match plotting_positions
    coverage_data <- coverage_data[coverage_data$Position %in% plotting_positions$Position, ]

  } else {
    plotting_positions <- data.frame(Reference= character(), Position=integer(), Difference=double(), Feature=character())
    significant_positions <- data.frame(Reference= character(), Position=integer(), Difference=double(), Feature=character())
    
  }

  return(list(plotting_positions, significant_positions, coverage_data))
}

process_bed <- function(bed_file, chr) {
  whole_bed <- read.delim(bed_file, header=FALSE)
  chr_bed <- subset(whole_bed, V1==chr)

  return(chr_bed)
}

barplot_plotting <- function (list_plotting, list_significant, output_name, MZS_thr, initial_pos, final_pos, annotation, ablines){

  #Rbind all data - already in long format:
  initial_join <- TRUE

  for (i in 2:length(list_plotting)){
    if (initial_join==TRUE){
      initial_df <- rbind(list_plotting[[i-1]], list_plotting[[i]])
      putative_positions <- rbind(list_significant[[i-1]], list_significant[[i]])
      initial_join <- FALSE
    } else {
      initial_df <- rbind(initial_df, list_plotting[[i]])
      putative_positions <- rbind(putative_positions, list_significant[[i]])
    }
  }

  #Set Feature into a factor for plotting purposes:
  initial_df$sample_f <- factor(initial_df$Feature, levels = c('Epinano', 'baseQ', 'nanoRMS_SI', 'nanoRMS_DT', 'nanoRMS_SD'))
  putative_positions$sample_f <- factor(putative_positions$Feature, levels = c('Epinano', 'baseQ', 'nanoRMS_SI', 'nanoRMS_DT', 'nanoRMS_SD'))

  ##Plotting:
  #If there are annotated positions:
  if (nrow(annotation)!=0 && ablines){
    barplot_4soft <- ggplot(initial_df, aes(x=Position, y=Modified_ZScore, fill=sample_f)) + ggtitle(output_name) +
      geom_bar(data=subset(initial_df, Modified_ZScore < MZS_thr), stat= "identity", width=4, fill = "#dcdcdd") +
      new_scale_color() + xlim(initial_pos, final_pos) + ylab('Z-Score ((x-median)/sd)') + xlab("") +
      geom_bar(data=subset(initial_df, Modified_ZScore >= MZS_thr), stat = "identity", width=4) +
      scale_fill_manual(values = c("#00A651", "#00AEEF", "#F59364", "#C48240", "#AB7A62"), breaks = c('Epinano', 'baseQ', 'nanoRMS_SI', 'nanoRMS_DT', 'nanoRMS_SD')) +
      geom_vline(xintercept=as.numeric(annotation$V3), linetype="dashed") +
      theme_bw() +theme(plot.title = element_text(face = "bold", hjust = 0.5), text = element_text(size=25),
                        axis.text = element_text(size = 25), strip.text.y = element_text(size = 20),
                        legend.text=element_text(size=22), legend.position = "none") +
      facet_grid(sample_f ~ . , scales="fixed", drop = FALSE)

  } else {
    barplot_4soft <- ggplot(initial_df, aes(x=Position, y=Modified_ZScore, fill=sample_f)) + ggtitle(output_name) +
      geom_bar(data=subset(initial_df, Modified_ZScore < MZS_thr), stat= "identity", width=4, fill = "#dcdcdd") +
      new_scale_color() + xlim(initial_pos, final_pos) + ylab('Z-Score ((x-median)/sd)') + xlab("") +
      geom_bar(data=subset(initial_df, Modified_ZScore >= MZS_thr), stat = "identity", width=4) +
      scale_fill_manual(values = c("#00A651", "#00AEEF", "#F59364", "#C48240", "#AB7A62"), breaks = c('Epinano', 'baseQ', 'nanoRMS_SI', 'nanoRMS_DT', 'nanoRMS_SD')) +
      theme_bw() +theme(plot.title = element_text(face = "bold", hjust = 0.5), text = element_text(size=25),
                        axis.text = element_text(size = 25), strip.text.y = element_text(size = 20),
                        legend.text=element_text(size=22), legend.position = "none") +
      facet_grid(sample_f ~ . , scales="fixed", drop = FALSE)
  }

  return(barplot_4soft)
}

Nanoconsensus_plotting <- function(data, supported_kmers, output_name, barplot_4soft, initial_pos, final_pos, annotation, ablines) {

  #Extracting supported kmers:
  supported_positions <- c()
  kmers_limits <- c()
  barplot_4soft <- barplot_4soft

  #Format data:
  data$Position <- data$Start+2
  data$Feature <- "NanoConsensus"

  #Generate NanoConsensus track:
  if (nrow(supported_kmers)!=0) {
    #If supported kmers, include them in the plot object:
    for (i in seq(1, nrow(supported_kmers))) {
      supported_positions <- c(supported_positions, seq(supported_kmers[i,2], supported_kmers[i,3]))
      kmers_limits <- c(kmers_limits, supported_kmers[i,2], supported_kmers[i,3])

    }

    #Retrieve borders of supported kmers:
    limits_supp_kmers <- subset(data[,c(-3,-2,-1)], Position %in% kmers_limits)

    #If there are annotated positions:
    if (nrow(annotation)!=0 && ablines){
      #Create plot object:
      nanoconsensus_plot <- ggplot(data, aes(x=Position, y=Merged_Score)) +
             geom_bar(stat= "identity", width=4, fill = "#dcdcdd") + ylim(0,1) +
             geom_bar(data=subset(data, Position %in% supported_positions), stat= "identity", width=4, fill = "#BE1E2D") +
             geom_label_repel(data=limits_supp_kmers,aes(label = Position, x=Position, y = Merged_Score), size = 8, label.size = 0.75) +
             ylab('NanoConsensus Score') +  xlim(initial_pos, final_pos) +
             geom_vline(xintercept=as.numeric(annotation$V3), linetype="dashed") +
             theme_bw() +theme(plot.title = element_text(face = "bold", hjust = 0.5), text = element_text(size=25),
                               axis.text = element_text(size = 25), strip.text.y = element_text(size = 20),
                               legend.text=element_text(size=22)) +
             facet_grid(Feature ~ . , scales="fixed")
    } else {
      nanoconsensus_plot <- ggplot(data, aes(x=Position, y=Merged_Score)) +
        geom_bar(stat= "identity", width=4, fill = "#dcdcdd") + ylim(0,1) +
        geom_bar(data=subset(data, Position %in% supported_positions), stat= "identity", width=4, fill = "#BE1E2D") +
        geom_label_repel(data=limits_supp_kmers,aes(label = Position, x=Position, y = Merged_Score), size = 8, label.size = 0.75) +
        ylab('NanoConsensus Score') + xlim(initial_pos, final_pos) +
        theme_bw() +theme(plot.title = element_text(face = "bold", hjust = 0.5), text = element_text(size=20),
                          axis.text = element_text(size = 25), strip.text.y = element_text(size = 20),
                          legend.text=element_text(size=22)) +
        facet_grid(Feature ~ . , scales="fixed")
    }

  } else {
    #Create plot object if there arent any supported kmers:
    if (nrow(annotation)!=0 && ablines){
      nanoconsensus_plot <- ggplot(data, aes(x=Position, y=Merged_Score)) +  geom_bar(stat= "identity", width=4, fill = "#dcdcdd") + ylim(0,1) +
             ylab('NanoConsensus Score') + xlim(initial_pos, final_pos) +
             geom_vline(xintercept=as.numeric(annotation$V3), linetype="dashed") +
             theme_bw() +theme(plot.title = element_text(face = "bold", hjust = 0.5), text = element_text(size=25),
                               axis.text = element_text(size = 25), strip.text.y = element_text(size = 25),
                               legend.text=element_text(size=22)) + facet_grid(Feature ~ . , scales="fixed")
    } else {
      nanoconsensus_plot <- ggplot(data, aes(x=Position, y=Merged_Score)) +  geom_bar(stat= "identity", width=4, fill = "#dcdcdd") + ylim(0,1) +
        ylab('NanoConsensus Score') + xlim(initial_pos, final_pos) +
        theme_bw() +theme(plot.title = element_text(face = "bold", hjust = 0.5), text = element_text(size=25),
                          axis.text = element_text(size = 25), strip.text.y = element_text(size = 20), 
                          legend.text=element_text(size=22)) + facet_grid(Feature ~ . , scales="fixed")
    }

  }

  #Plot both plots in the same pdf file:
  pdf(file=paste(output_name,"NanoConsensus_Scores.pdf", sep = "-"), bg = "transparent", width = 26, height = 18.50 )
  g2 <- ggplotGrob(barplot_4soft)
  g3 <- ggplotGrob(nanoconsensus_plot)
  g <- rbind(g2, g3, size = "last")
  g$widths <- unit.pmax(g2$widths, g3$widths)
  #grid.newpage()
  grid.draw(g)

  dev.off()

}

extract_length_from_GRobjects <- function(GRange_object) {

  if (length(GRange_object)==0){
    n_length <- 0
  } else {
    n_length <- length(GRange_object)
  }

  return(n_length)

}

overlapping_GRobjects <- function(GRange_object_1, GRange_object_2, length_object1, length_object2) {

  #Only perform the intersection if both GRange are valid:
  if (is.null(GRange_object_1)==FALSE & is.null(GRange_object_2)==FALSE) {

    #Perform the intersections:
    if (length_object1 >= length_object2) {
      intersect_object <- subsetByOverlaps(GRange_object_2, GRange_object_1, minoverlap=1)
    } else {
      intersect_object <- subsetByOverlaps(GRange_object_1, GRange_object_2, minoverlap=1)
    }

  } else {

    intersect_object <- GRanges()
  }

  return(intersect_object)

}

extract_kmers <- function (bedfile, fasta) {
  #Create Temp files:
  a.file=tempfile()
  out=tempfile()

  #Format the ranges to obtain 9-mers, centered in the 5mer identified by NanoMod - REMEMBER: bedtools understands bed files as 0-based!
  bedfile$Start <- bedfile$Start-3
  bedfile$End <- bedfile$End+2

  #Write formatted dataframes to tempfile
  write.table(bedfile,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)

  #Create the command for the bedtools command and execute it:
  command=paste("bedtools getfasta -fi", fasta, "-bed",a.file,"-tab >",out,sep=" ")
  try(system(command))

  #Save results into a dataframe:
  res=read.table(out,header=F)

  #Check if there is the RRACH motif using regular expressions:
  motif <- c()
  pattern <- "([A|G]{2})AC([A|C|T]{1})"

  for (i in 1:nrow(res)) {
    motif <- c(motif,str_detect(res$V2[i], pattern))
  }

  #Add new columns:
  colnames(res) <- c('Data', 'Kmer')
  res$RRACH_motif <- motif

  return(list(res$Kmer, res$RRACH_motif))

}

overwrite_NaNs <- function (input) {

  if (is.nan(input) == TRUE || is.infinite(input) == TRUE || length(input) == 0) {
    out_value <- NA
  } else {
    out_value <- input
  }

  return(out_value)
}

extracting_status <- function (positions_df, list_number, summit, MZS_thr) {

  ##Declaring initial variables:
  soft_rawScore <- c()
  soft_modifiedScore <- c()
  soft_status <- c()

  #Loop across kmers:
  for (i in seq(1:length(positions_df$Start))){
    initial_position <- positions_df$Start[i]
    final_position <- positions_df$End[i]

    ##Searching for the highest value - summit:
    if (summit == TRUE) {
      #Looping within the kmer to find the highest value - summit:
      for (x in seq(initial_position, final_position)){

        if (x == initial_position){
          highest_rawScore <- overwrite_NaNs(list_plotting[[list_number]][which(list_plotting[[list_number]]$Position == x), 3])
          highest_modifiedScore <- overwrite_NaNs(list_plotting[[list_number]][which(list_plotting[[list_number]]$Position == x), 5])

        } else {
          new_rawScore <- overwrite_NaNs(list_plotting[[list_number]][which(list_plotting[[list_number]]$Position == x), 3])
          new_modifiedScore <- overwrite_NaNs(list_plotting[[list_number]][which(list_plotting[[list_number]]$Position == x), 5])

          ##Checking for higher score - rawScore:
          if (is.na(highest_rawScore) == TRUE || is.nan(highest_rawScore) == TRUE){
            if (is.na(new_rawScore) == FALSE & is.nan(new_rawScore) == FALSE) {
              highest_rawScore <- new_rawScore
            }
          } else if (is.na(new_rawScore) == TRUE || is.nan(new_rawScore) == TRUE) {
            next

          } else {
            if (new_rawScore > highest_rawScore){
              highest_rawScore <- new_rawScore
            }
          }


          ##Checking for higher score - modified score:
          if (is.na(highest_modifiedScore) == TRUE || is.nan(highest_modifiedScore) == TRUE){
            if (is.na(new_modifiedScore) == FALSE & is.nan(new_modifiedScore) == FALSE) {
              highest_modifiedScore <- new_modifiedScore
            }
          } else if (is.na(new_modifiedScore) == TRUE || is.nan(new_modifiedScore) == TRUE) {
            next

          } else {
            if (new_modifiedScore > highest_modifiedScore){
              highest_modifiedScore <- new_modifiedScore
            }
          }

        }

      }

      ##Adding high score to final output:
      soft_rawScore[i] <- highest_rawScore
      soft_modifiedScore[i] <- highest_modifiedScore

      #Check if a specific software identified it:
      if (is.na(highest_modifiedScore) == TRUE || is.nan(highest_modifiedScore) == TRUE) {
        soft_status <- c(soft_status, 'NO')

      } else {
        if(soft_modifiedScore[i] >= MZS_thr){
          soft_status <- c(soft_status, 'YES')
        } else {
          soft_status <- c(soft_status, 'NO')
        }
      }

    } else {
      ##Searching for position 0 value:
      position <- initial_position + 2

      ##Software - extract values:
      rawScore <- list_plotting[[list_number]][which(list_plotting[[list_number]]$Position == position), 3]
      if (length(rawScore)==0 || is.infinite(rawScore) == TRUE){
        soft_rawScore <- c(soft_rawScore, NA)
      } else {
        soft_rawScore <- c(soft_rawScore, rawScore)
      }

      modifiedScore <- list_plotting[[list_number]][which(list_plotting[[list_number]]$Position == position), 5]

      if (length(modifiedScore)==0 || is.infinite(modifiedScore) == TRUE){
        soft_modifiedScore <- c(soft_modifiedScore, NA)
      } else {
        soft_modifiedScore <- c(soft_modifiedScore, modifiedScore)
      }

      #Define values and overwrite NaNs:
      single_pos_0 <- overwrite_NaNs(list_plotting[[list_number]][which(list_plotting[[list_number]]$Position == initial_position), 5])
      single_pos_1 <- overwrite_NaNs(list_plotting[[list_number]][which(list_plotting[[list_number]]$Position == initial_position+1), 5])
      single_pos_2 <- overwrite_NaNs(list_plotting[[list_number]][which(list_plotting[[list_number]]$Position == initial_position+2), 5])
      single_pos_3 <- overwrite_NaNs(list_plotting[[list_number]][which(list_plotting[[list_number]]$Position == initial_position+3), 5])
      single_pos_4 <- overwrite_NaNs(list_plotting[[list_number]][which(list_plotting[[list_number]]$Position == initial_position+4), 5])

      #Loop over the kmer to find if specific softwares identified it:
      kmer_positions <- c(single_pos_0, single_pos_1, single_pos_2, single_pos_3, single_pos_4)
      if(any(kmer_positions >= MZS_thr, na.rm = TRUE)){
        soft_status <- c(soft_status, 'YES')
      } else {
        soft_status <- c(soft_status, 'NO')
      }
    }

  }

  #Create final dataframe:
  final <- data.frame(soft_rawScore, soft_modifiedScore, soft_status)
  colnames(final) <- c('rawScore', 'modifiedScore', 'status')

  return(final)
}

calcNanoConsensusScore <- function(data, type) {
  
  if (type=="m66A") {
    w <- c(0.36,0.1,0.21,0.33)
    processed_data <- sweep(data, MARGIN=2, w, "*")
    return(apply(processed_data,1,sum,na.rm = TRUE))

  } else if (type=="m5C") {
    w <- c(0.25,0.05,0.70,0)
    processed_data <- sweep(data, MARGIN=2, w, "*")
    return(apply(processed_data,1,sum,na.rm = TRUE))

  } else if (type=="m7G") {
    w <- c(0.33,0.25,0.42,0)
    processed_data <- sweep(data, MARGIN=2, w, "*")
    return(apply(processed_data,1,sum,na.rm = TRUE))

  } else if (type=="Am") {
    w <- c(0.23,0.07,0.55,0.15)
    processed_data <- sweep(data, MARGIN=2, w, "*")
    return(apply(processed_data,1,sum,na.rm = TRUE))

  } else if (type=="Um") {
    w <- c(0.27,0.12,0.5,0.11)
    processed_data <- sweep(data, MARGIN=2, w, "*")
    return(apply(processed_data,1,sum,na.rm = TRUE))

  } else if (type=="pU") {
    w <- c(0.22,0.19,0.55,0.04)
    processed_data <- sweep(data, MARGIN=2, w, "*")
    return(apply(processed_data,1,sum,na.rm = TRUE))

  } else {
    print(data)
    return(apply(data,1,median,na.rm = TRUE))
  }
}

extracting_modified_ZScores <- function (GRange_supported_kmers, MZS_thr, summit, Consensus_score, model_score) {

  #If there aren't any supported kmers:
  if(is.null(GRange_supported_kmers)==FALSE){

    #Parse data into a data frame:
    positions_df <- data.frame(start(GRange_supported_kmers), end(GRange_supported_kmers))
    colnames(positions_df) <- c('Start', 'End')
    positions_df$Chr <- seqlevels(GRange_supported_kmers)
    positions_df <- positions_df[,c(3,1,2)]
    
    #Extracting scores and software status:
    epinano_data <- extracting_status(positions_df, 1, summit, MZS_thr)
    baseQ_data <- extracting_status(positions_df, 2, summit, MZS_thr)
    nanoRMS_SI_data <- extracting_status(positions_df, 3, summit, MZS_thr)
    nanoRMS_DT_data <- extracting_status(positions_df, 4, summit, MZS_thr)
    nanoRMS_SD_data <- extracting_status(positions_df, 5, summit, MZS_thr)
    
    #Add data to the final dataframe:
    positions_df$Epinano_RawScore <- epinano_data$rawScore
    positions_df$baseQ_RawScore <- baseQ_data$rawScore
    positions_df$NanoRMS_SI_RawScore <- nanoRMS_SI_data$rawScore
    positions_df$NanoRMS_DT_RawScore <- nanoRMS_DT_data$rawScore
    positions_df$NanoRMS_SD_RawScore <- nanoRMS_SD_data$rawScore

    positions_df$Epinano_Score <- epinano_data$modifiedScore
    positions_df$baseQ_Score <- baseQ_data$modifiedScore
    positions_df$NanoRMS_SI_Score <- nanoRMS_SI_data$modifiedScore
    positions_df$NanoRMS_DT_Score <- nanoRMS_DT_data$modifiedScore
    positions_df$NanoRMS_SD_Score <- nanoRMS_SD_data$modifiedScore

    positions_df$Epinano_Status <- epinano_data$status
    positions_df$baseQ_Status <- baseQ_data$status
    positions_df$NanoRMS_SI_Status <- nanoRMS_SI_data$status
    positions_df$NanoRMS_DT_Status <- nanoRMS_DT_data$status
    positions_df$NanoRMS_SD_Status <- nanoRMS_SD_data$status
    
    positions_NanoConsensus <- c()
    
    ##Calculate the merged_score:
    #Re-scaling:
    if (summit == F){
      data <- data.frame(positions_df$Epinano_Score, positions_df$baseQ_Score, positions_df$NanoRMS_SI_Score, positions_df$NanoRMS_DT_Score, positions_df$NanoRMS_SD_Score)

      #Re-scale Modified Z-Scores between 0 and 1:
      for (i in seq(1:length(data))) {
        data[,i] <- rescale(unlist(data[i]), to=c(0,1), na.rm=TRUE)

      }
      
      #Rescale outputs 0.5 when the software gives the same MZS for all positions - correcting for it if needed:
      #If all values are NA, keep NA; if all non-NA values are the same, replace with 0
      if (!all(is.na(data$positions_df.Epinano_Score)) && length(unique(na.omit(data$positions_df.Epinano_Score))) == 1) {
        data$positions_df.Epinano_Score <- 0
      }

      if (!all(is.na(data$positions_df.baseQ_Score)) && length(unique(na.omit(data$positions_df.baseQ_Score))) == 1) {
        data$positions_df.baseQ_Score <- 0
      }

      if (!all(is.na(data$positions_df.NanoRMS_SI_Score)) && length(unique(na.omit(data$positions_df.NanoRMS_SI_Score))) == 1) {
        data$positions_df.NanoRMS_SI_Score <- 0
      }

      if (!all(is.na(data$positions_df.NanoRMS_DT_Score)) && length(unique(na.omit(data$positions_df.NanoRMS_DT_Score))) == 1) {
        data$positions_df.NanoRMS_DT_Score <- 0
      }

      if (!all(is.na(data$positions_df.NanoRMS_SD_Score)) && length(unique(na.omit(data$positions_df.NanoRMS_SD_Score))) == 1) {
        data$positions_df.NanoRMS_SD_Score <- 0
      }

      #Calculate NanoConsensus score:
      write(paste("Step 4: Calculating NanoConsensus scores with model: ", model_score, sep = ""), file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)
      positions_df$Merged_Score <- calcNanoConsensusScore(data, model_score)
      
      threshold <- Consensus_score * median(positions_df$Merged_Score, na.rm = TRUE)
      print(threshold)

      if (!is.na(threshold) && threshold == 0) {
        positions_NanoConsensus <- subset(positions_df, Merged_Score > threshold)
      } else {
        positions_NanoConsensus <- subset(positions_df, Merged_Score >= threshold)
      }

    }
  } else {
    positions_df <- ""
    positions_NanoConsensus <- ""

  }

  return(list(positions_df, positions_NanoConsensus))
}

bedRmod_tracks <- function (data, output_name, color, methods, organism = "", assembly = "", annotation_source = "", annotation_version = "", basecalling = "", bioinformatics_workflow = "", experiment = "", strand = "") {

  #Check if output directory exists - if not, create it:
  if (!dir.exists("./BedRmod_tracks")){
    dir.create("BedRmod_tracks", showWarnings = FALSE)
  }
  
  for (i in 4:(ncol(data))){
    #Sliced dataset:
    subset_data <- data.frame(Chr = data[,1], Start = data[,2], End = data[,3], Score = as.numeric(data[,i]))

    #From kmers to individual positions - center on single position:
    subset_data$Start <- subset_data$Start + 2
    subset_data$End <- subset_data$Start + 1

    #Replace NAs with 0:
    subset_data$Score[is.na(subset_data$Score)] <- 0

    #Format data for bedRmod:
    bedRmod_data <- data.frame(
      Chr = subset_data$Chr,
      Start = subset_data$Start,
      End = subset_data$End,
      Name = 'xX',
      Score = sprintf("%.3f", subset_data$Score),
      Strand = strand,
      ThickStart = data[,2],
      ThickEnd = data[,3],
      ItemRgb = "0,0,0",
      Coverage = "",
      Frequency = "",
      stringsAsFactors = FALSE
    )

    #Generate bedRmod headers:
    headers <- create_bedRmod_headers(organism = organism, modification_names = 'xX', assembly = assembly,
                                      annotation_source = annotation_source, annotation_version = annotation_version,
                                      basecalling = basecalling, bioinformatics_workflow = methods[i-3],
                                      experiment = experiment)

    #Generate bedRmod output file:
    output_file <- paste("BedRmod_tracks/", methods[i-3], "-", str_split_fixed(output_name,"_Raw_kmers.txt",2)[1],'.bedrmod', sep='')
    writeLines(headers, con = output_file)

    #Append data to file:
    write.table(bedRmod_data, file = output_file, sep = '\t', row.names = FALSE,
                col.names = FALSE, quote = FALSE, append = TRUE)
  }
}

bed_tracks <- function (data, output_name, color, methods) {
  #Check if output directory exists - if not, create it:
  if (!dir.exists("./Kmer_tracks")){
    dir.create("Kmer_tracks", showWarnings = FALSE)
  }

  for (i in 1:length(data)){
    bed_data <- data.frame(seqnames(data[[i]]),
                     start(data[[i]]),
                     end(data[[i]]),
                     c(rep(".", length(data[[i]]))),
                     c(rep(0, length(data[[i]]))),
                     strand(data[[i]]), start(data[[i]]),end(data[[i]]), c(rep(color[i], length(data[[i]]))))


    #Prepare the header:
    header_track <- paste(" \'1i track name=", methods[i],"_kmers visibility=2 itemRgb=\"On\"\'", sep="")

    #Generate bed tracks:
    write.table(bed_data, file = paste("Kmer_tracks/", methods[i], "-", output_name, "-kmers.bed", sep=''),
                sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

    #Include the header to be able to load the track into IGV:
    command=paste("sed -i", header_track, paste(" ./Kmer_tracks/", methods[i], "-", output_name, "-kmers.bed", sep=''), sep="")
    try(system(command))

  }
}

nearest_distance_mod <- function(all_ranges, annotation) {
  distance <- c()
  mods <- c()

  #Loop through all the supported kmers:
  for (i in 1:nrow(all_ranges)){

    #Define variables to determine distance to nearest modified site:
    initial <- as.numeric(all_ranges[i,2])
    final <- as.numeric(all_ranges[i,3])
    single_distance <- c()
    single_mods <- c()
    within <- FALSE
    
    #Skip if annotation is empty for this specific chr:
    if (nrow(annotation) == 0) next
    
    #Loop through all the annotated positions:
    for (j in 1:nrow(annotation)){
      annotated_position <- as.numeric(annotation[j,3])

      #Annotated position within the modified kmer:
      if (annotated_position<=final && annotated_position>=initial && within==TRUE){
        single_distance <- c(single_distance, 0)
        single_mods <- c(single_mods, paste(annotation[j,4],annotation[j,3], sep="-"))

      } else if (annotated_position<=final && annotated_position>=initial) {
        single_distance <- 0
        single_mods <- paste(annotation[j,4],annotation[j,3], sep="-")
        within <- TRUE

      } else {
        #Annotated position outside the modified kmer:
        d <- min(abs(as.numeric(annotated_position)-as.numeric(initial)), abs(as.numeric(annotated_position)-as.numeric(final)))
        if (d<single_distance || length(single_distance)==0){
          single_distance <- d
          single_mods <- paste(annotation[j,4],annotation[j,3], sep="-")

        } else if (d==single_distance) {
          single_mods <- c(single_mods, paste(annotation[j,4],annotation[j,3], sep="-"))
        }
      }

    }

    distance <- c(distance, unique(single_distance))
    mods <- c(mods, str_c(single_mods, collapse = ","))

  }
  #Adding the data to the final results:
  all_ranges$Distance_nearest_mods <- distance
  all_ranges$Nearest_mods <- mods

  return(all_ranges)
}

create_bedRmod_headers <- function(organism = "", modification_names = "", assembly = "",
                                    annotation_source = "", annotation_version = "",
                                    basecalling = "", bioinformatics_workflow = "", experiment = "", strand = "") {
  #Check if mandatory fields are empty
  missing_fields <- c()
  if (organism == "") missing_fields <- c(missing_fields, "organism")
  if (assembly == "") missing_fields <- c(missing_fields, "assembly")
  if (annotation_source == "") missing_fields <- c(missing_fields, "annotation_source")
  if (annotation_version == "") missing_fields <- c(missing_fields, "annotation_version")

  if (length(missing_fields) > 0) {
    missing_str <- paste(missing_fields, collapse = ", ")
    print(paste("Mandatory files required by bedRmod format are left blank as the user didn't provide the information. Missing fields: ", missing_str, sep = ""))
  }

  c(
    "#fileformat=bedRModv2",
    paste0("#organism=", organism),
    "#modification_type=RNA",
    paste0("#modification_names=", modification_names),
    paste0("#assembly=", assembly),
    paste0("#annotation_source=", annotation_source),
    paste0("#annotation_version=", annotation_version),
    "#sequencing_platform=ONT",
    paste0("#basecalling=", basecalling),
    paste0("#bioinformatics_workflow=", bioinformatics_workflow),
    paste0("#experiment=", experiment),
    "#external_source=",
    "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tcoverage\tfrequency"
  )
}

write_bedRmod_output <- function(all_ranges, output_name, bed_data = NULL, annotation = NULL, coverage_data = NULL, organism = "", assembly = "", annotation_source = "", annotation_version = "", basecalling = "", bioinformatics_workflow = "", experiment = "", strand = "") {
  bedRmod_rows <- list()
  collected_mod_names <- c()

  # Process all ranges from all_ranges
  for (i in 1:nrow(all_ranges)) {
    range_chr <- all_ranges$Chr[i]
    range_start <- as.numeric(all_ranges$Start[i])
    range_end <- as.numeric(all_ranges$End[i])

    # Check if there's an overlapping annotation for this range
    overlapping_annotation <- NULL
    if (!is.null(annotation) && nrow(annotation) > 0) {
      overlapping_annotation <- annotation[annotation[,1] == range_chr &
                                          annotation[,3] >= range_start &
                                          annotation[,3] <= range_end, ]
    }

    if (!is.null(overlapping_annotation) && nrow(overlapping_annotation) > 0) {
      # Case: Range overlaps with annotation(s)
      # Process each overlapping annotation separately
      for (ann_row in 1:nrow(overlapping_annotation)) {
        ann_pos <- as.numeric(overlapping_annotation[ann_row, 2])
        ann_mod_name <- overlapping_annotation[ann_row, 4]
        collected_mod_names <- c(collected_mod_names, ann_mod_name)

        # Get max score around annotation position
        if (!is.null(bed_data)) {
          three_mer <- bed_data[bed_data$Chr == range_chr &
                               bed_data$Start >= (ann_pos - 3) &
                               bed_data$Start <= (ann_pos - 1), ]
          max_score <- if (nrow(three_mer) > 0) max(three_mer$Merged_Score, na.rm = TRUE) else NA
        } else {
          max_score <- NA
        }

        # Format score with 3 decimals
        score_formatted <- ifelse(!is.na(max_score), sprintf("%.3f", max_score), "")

        # Look up coverage for this position
        cov_value <- ""
        if (!is.null(coverage_data) && nrow(coverage_data) > 0) {
          cov_match <- coverage_data[coverage_data$Position == ann_pos, ]
          if (nrow(cov_match) > 0) {
            cov_value <- as.integer(cov_match$Coverage[1])
          }
        }

        bedRmod_rows[[length(bedRmod_rows) + 1]] <- data.frame(
          Chr = range_chr,
          Start = ann_pos,
          End = ann_pos + 1,
          Name = ann_mod_name,
          Score = score_formatted,
          Strand = strand,
          ThickStart = range_start,
          ThickEnd = range_end,
          ItemRgb = "0,0,0",
          Coverage = cov_value,
          Frequency = "",
          stringsAsFactors = FALSE
        )
      }
    } else {
      # Case: Range doesn't overlap with annotation (or no annotation)
      collected_mod_names <- c(collected_mod_names, "xX")

      # Find position within kmer with highest Merged_Score
      if (!is.null(bed_data)) {
        kmer_data <- bed_data[bed_data$Chr == range_chr &
                             bed_data$Start >= range_start &
                             bed_data$End <= range_end, ]

        if (nrow(kmer_data) > 0) {
          max_idx <- which.max(kmer_data$Merged_Score)
          max_score <- kmer_data$Merged_Score[max_idx]
          central_pos <- kmer_data$Start[max_idx] + 2
          max_pos_start <- central_pos
          max_pos_end <- central_pos + 1
        } else {
          max_score <- NA
          max_pos_start <- range_start
          max_pos_end <- range_end
        }
      } else {
        max_score <- NA
        max_pos_start <- range_start
        max_pos_end <- range_end
      }

      # Format score with 3 decimals
      score_formatted <- ifelse(!is.na(max_score), sprintf("%.3f", max_score), "")

      # Look up coverage for this position
      cov_value <- ""
      if (!is.null(coverage_data) && nrow(coverage_data) > 0) {
        cov_match <- coverage_data[coverage_data$Position == max_pos_start, ]
        if (nrow(cov_match) > 0) {
          cov_value <- as.integer(cov_match$Coverage[1])
        }
      }

      bedRmod_rows[[length(bedRmod_rows) + 1]] <- data.frame(
        Chr = range_chr,
        Start = max_pos_start,
        End = max_pos_end,
        Name = "xX",
        Score = score_formatted,
        Strand = strand,
        ThickStart = range_start,
        ThickEnd = range_end,
        ItemRgb = "0,0,0",
        Coverage = cov_value,
        Frequency = "",
        stringsAsFactors = FALSE
      )
    }
  }

  # Combine all rows into single dataframe
  if (length(bedRmod_rows) > 0) {
    bedRmod_data <- do.call(rbind, bedRmod_rows)
    rownames(bedRmod_data) <- NULL
  } else {
    bedRmod_data <- data.frame()
  }

  # Generate and write bedRmod headers with collected modification names
  modification_names <- paste(unique(collected_mod_names), collapse = ",")
  headers <- create_bedRmod_headers(organism = organism, modification_names = modification_names, assembly = assembly, annotation_source = annotation_source, annotation_version = annotation_version, basecalling = basecalling, bioinformatics_workflow = bioinformatics_workflow, experiment = experiment)
  writeLines(headers, con = output_name)

  # Append data to file only if there are rows
  if (nrow(bedRmod_data) > 0) {
    write.table(bedRmod_data, file = output_name, sep = '\t', row.names = FALSE,
                col.names = FALSE, quote = FALSE, append = TRUE)
  }
}

kmer_analysis <- function (all_ranges, fasta_file, output_name, tracks, annotation, sup_kmers, color_beds, methods_name, bedRmod = FALSE, bed_data = NULL, coverage_data = NULL, organism = "", assembly = "", annotation_source = "", annotation_version = "", basecalling = "", bioinformatics_workflow = "", experiment = "", extended_outputs = NULL, strand = "") {
  print('Kmer analysis')
  kmer_data <- extract_kmers(all_ranges, fasta_file)
  all_ranges$Kmer <- kmer_data[[1]]
  all_ranges$RRACH_motif <- kmer_data[[2]]
  all_ranges <- all_ranges[order(all_ranges$Start, decreasing = FALSE),]
  
  #If needed, generate bedRmod tracks:
  if (tracks){
    bedRmod_tracks(all_ranges[,c(1,2,3,9,10,11,12,13,19)], output_name, c(color_beds,"190,30,45"), c(methods_name, 'NanoConsensus'), organism, assembly, annotation_source, annotation_version, basecalling, bioinformatics_workflow, experiment, strand)
  }

  #If annotation file is provided, calculate distance to the nearest + annotated site:
  if (sup_kmers && length(annotation)!=0){
    all_ranges <- nearest_distance_mod(all_ranges, annotation)
  }
  
  #Output formatting:
  if (bedRmod) {
    write_bedRmod_output(all_ranges, output_name, bed_data, annotation, coverage_data, organism, assembly, annotation_source, annotation_version, basecalling, bioinformatics_workflow, experiment, strand)
  } else if (extended_outputs) {
    write.table(all_ranges, file = output_name, sep = '\t', row.names = FALSE, quote = FALSE)
  }

}

analysis_significant_positions <- function (list_significant, list_plotting, fasta_file, output_name, initial_position, final_position, MZS_thr, Consensus_score, model_score, barplot_4soft, annotation, ablines, chr, coverage_data = NULL, organism = "", assembly = "", annotation_source = "", annotation_version = "", basecalling = "", bioinformatics_workflow = "", experiment = "", extended_outputs = FALSE, strand = "") {
  epinano <- list_significant[[1]]
  baseQ <- list_significant[[2]]
  nanoRMS_SI <- list_significant[[3]]
  nanoRMS_DT <- list_significant[[4]]
  nanoRMS_SD <- list_significant[[5]]

  methods <- list(epinano, baseQ, nanoRMS_SI, nanoRMS_DT, nanoRMS_SD)
  methods_name <- c('Epinano', 'baseQ', 'nanoRMS_SI', 'nanoRMS_DT', 'nanoRMS_SD')

  #Create grRange objects with kmers per each method:
  #print('Transforming data into GRange objects')
  for (j in 1:length(methods_name)) {

    #Transform the positions in kmers using GRanges library:
    grList <- c()
    sliced <- methods[[j]]

    if (length(sliced[[1]])>0) {
      for (i in 1:nrow(sliced)){
        features <- strsplit(as.character(sliced[[1]][i]), "[_]")

        #Create GR objects:
        if (length(features[[1]])==2){
          chr <- features[[1]][1]
        } else {
          elements <- c()
          for (i in 1:length(features[[1]])-1){
            elements <- c(elements, features[[1]][i])
          }
          chr <- paste(elements, collapse="_")

        }
        grNew <- GRanges(seqnames=chr,ranges=IRanges(as.integer(features[[1]][length(features[[1]])])-2, end = as.integer(features[[1]][length(features[[1]])])+2))

        if(is.null(grList)==TRUE){
          grList <- grNew
        } else {
          grList <- reduce(c(grList,grNew))
        }

      }

      assign(paste('gr',methods_name[j],sep=""), reduce(unique(grList)))

    } else {
      assign(paste('gr',methods_name[j],sep=""), unique(grList))
    }

  }

  ##Perform intersections:
  #Check how many elements are in each GRange object and if it is null, create an empty one:
  GRanges_list <- list(grEpinano, grbaseQ, grnanoRMS_SI, grnanoRMS_DT, grnanoRMS_SD)
  supported_kmers_per_software <- c()

  for (count in seq(1, length(GRanges_list))){

    if (is.null(GRanges_list[[count]])==TRUE) {
      GRanges_list[[count]] <- GRanges()
      supported_kmers_per_software <- c(supported_kmers_per_software, 0)

    } else {
      supported_kmers_per_software <- c(supported_kmers_per_software, length(GRanges_list[[count]]))
    }

    ##Update log file:
    write(paste('-Positions identified by', methods_name[count], ':', supported_kmers_per_software[count], sep = " "), file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)

  }
  
  ##Generate bed files:
  color_beds <- c("0,166,81", "0,174,239","242,101,34", "196,130,64", "171,122,98")
  if (extended_outputs) {
    bed_tracks(GRanges_list, output_name, color_beds, methods_name)
  }

  ##Overlappings of supported kmers:
  write('Kmers supported by multiple softwares:', file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)

  #Merge all GRanges into a non-redundant set
  all_ranges <- do.call(c, GRanges_list)
  unique_regions <- reduce(all_ranges)
  
  #Create a logical matrix: rows = significant kmers, columns = original GRanges object
  overlap_matrix <- do.call(cbind, lapply(GRanges_list, function(gr) {
    overlaps <- findOverlaps(unique_regions, gr)
    hits <- logical(length(unique_regions))
    hits[queryHits(overlaps)] <- TRUE
    hits
  }))

  #Ensure overlap_matrix is always a 2D matrix
  if (is.null(dim(overlap_matrix))) {
    overlap_matrix <- as.matrix(overlap_matrix)
  }

  #Report the regions supported by 2, 3, 4, 5 and 6 softwares:
  initial <- TRUE
  
  for (i in seq(2,5)){
    intersection <- unique_regions[rowSums(overlap_matrix) >= i]
    length_intersection <- extract_length_from_GRobjects(intersection)

    #Update log file:
    write(paste('-Positions identified by', i, 'softwares:', length_intersection, sep = " "), file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)

    #Update supported_kmers:
    if (length_intersection > 0 & initial) {
      supported_kmers_list <- c(intersection)
      initial <- FALSE

    } else if (length_intersection == 0 & initial) {
      next

    } else {
      supported_kmers_list <- c(supported_kmers_list, intersection)
    }

  }
  
  #If no supported kmers are found, create NULL object:
  if (!exists('supported_kmers_list')){
    supported_kmers <- NULL
  } else {
    supported_kmers <- reduce(supported_kmers_list)
  }

  ##Kmer analysis:
  #Analysis of all kmers across the chromosome:
  all_kmers_raw <- GRanges(seqnames = chr, ranges = IRanges(initial_position:(final_position-4), end = (initial_position+4):final_position))
  all_kmers <- extracting_modified_ZScores(all_kmers_raw, MZS_thr, FALSE, Consensus_score, model_score)
  kmer_analysis(all_kmers[[1]], fasta_file, paste(output_name,'Raw_kmers.txt', sep='_'), TRUE, annotation, FALSE, color_beds, methods_name, extended_outputs = extended_outputs, strand = strand)

  #Analyse the supported kmers - only if they are present:
  if (is.null(supported_kmers)==FALSE) {
    filtered_supported_kmers <- overlapping_GRobjects(reduce(supported_kmers), GRanges(seqnames=all_kmers[[2]][,c('Chr')],ranges=IRanges(all_kmers[[2]][,c('Start')], end = all_kmers[[2]][,c('End')])),1,2)

    if(extract_length_from_GRobjects(filtered_supported_kmers)!=0){
      all_ranges <- extracting_modified_ZScores(filtered_supported_kmers, MZS_thr, TRUE, Consensus_score, model_score)
      kmer_analysis(all_ranges[[1]], fasta_file, paste(output_name,'Supported_kmers.bedrmod', sep='_'), FALSE, annotation, TRUE, color_beds, methods_name, bedRmod = TRUE, bed_data = all_kmers[[1]], coverage_data = coverage_data, organism = organism, assembly = assembly, annotation_source = annotation_source, annotation_version = annotation_version, basecalling = basecalling, bioinformatics_workflow = bioinformatics_workflow, experiment = experiment, strand = strand)

      #Plot NanoConsensus score across transcripts:
      Nanoconsensus_plotting(all_kmers[[1]], all_ranges[[1]], output_name, barplot_4soft, initial_position, final_position, annotation, ablines)
      write("Step 5: Plotting NanoConsensus scores across the transcript", file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)

    } else {
      all_ranges <- data.frame()
      #Plot NanoConsensus score across transcripts:
      Nanoconsensus_plotting(all_kmers[[1]], all_ranges, output_name, barplot_4soft, initial_position, final_position, annotation, ablines)

      write("Step 5: Plotting NanoConsensus scores across the transcript", file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)

    }

  } else {

    all_ranges <- data.frame()
    #Plot NanoConsensus score across transcripts:
    Nanoconsensus_plotting(all_kmers[[1]], all_ranges, output_name, barplot_4soft, initial_position, final_position, annotation, ablines)
    write("Step 5: Plotting NanoConsensus scores across the transcript", file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)

  }
  write("ANALYSIS COMPLETED SUCCESSFULLY", file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)

}
