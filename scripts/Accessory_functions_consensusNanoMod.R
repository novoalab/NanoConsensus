###Script which contains multiple R functions used to generate consensus putative modified positions from NanoMod results. 

#Processing Epinano results:
epinano_processing <- function(sample_file, ivt_file, initial_position, final_position, MZS_thr, chr, exclude_SNP) {
  
  #Import and clean data:
  sample <- read.csv(sample_file,stringsAsFactors = FALSE)
  sample <- subset(sample, cov>30)
  sample <- subset(sample, pos>=initial_position) 
  sample <- subset(sample, pos<=final_position)
  sample$reference <- paste(sample$X.Ref, sample$pos, sep='_')
  sample$Difference <- as.numeric(sample$mis)+as.numeric(sample$ins)+as.numeric(sample$del)
  sample <- sample[,c(1,2,12,11)]
  colnames(sample) <- c('Reference', 'Position', 'Difference_sample', 'Merge')
  
  ivt <- read.csv(ivt_file,stringsAsFactors = FALSE)
  ivt <- subset(ivt, cov>30)
  ivt <- subset(ivt, pos>=initial_position) 
  ivt <- subset(ivt, pos<=final_position)
  ivt$reference <- paste(ivt$X.Ref, ivt$pos, sep='_')
  ivt$Difference <- as.numeric(ivt$mis)+as.numeric(ivt$ins)+as.numeric(ivt$del)
  ivt <- ivt[,c(1,2,12,11)]
  colnames(ivt) <- c('Reference', 'Position', 'Difference_IVT', 'Merge')
  
  if (nrow(sample)!=0 && nrow(ivt)!=0) {
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
  
  } else {
    plotting_positions <- data.frame(Reference= character(), Position=integer(), Difference=double(), Feature=character())
    significant_positions <- data.frame(Reference= character(), Position=integer(), Difference=double(), Feature=character())
  }

  return(list(plotting_positions,significant_positions))
}

nanopolish_processing <- function(sample_file, ivt_file, initial_position, final_position, MZS_thr, chr, exclude_SNP) {
  #Import data:
  sample <- read.delim(sample_file)
  
  #Add sample information:
  sample$read_name <- 'Nanopolish'
  sample <- subset(sample, coverage>30)
  colnames(sample)<- c("contig_wt","position","reference_kmer_wt", "feature_wt", "event_level_median_wt", 'coverage')
  sample<- subset(sample, contig_wt == chr)
  sample$reference <- paste(sample$contig_wt, sample$position, sep='_')
  
  #Import KO: 
  raw_data_ivt <-read.delim(ivt_file)
  raw_data_ivt$read_name <- 'IVT'
  raw_data_ivt <- subset(raw_data_ivt, coverage>30)
  colnames(raw_data_ivt)<- c("contig_ko","position","reference_kmer_ko", "feature", "event_level_median_ko", 'coverage')
  raw_data_ivt <- subset(raw_data_ivt, contig_ko == chr)
  raw_data_ivt$reference <- paste(raw_data_ivt$contig_ko, raw_data_ivt$position, sep='_')
  
  #Join tables, calculate differences between means/medians:
  plotting_data <- join(sample, raw_data_ivt, by="reference", type='inner')
  plotting_data$diff <- abs(plotting_data$event_level_median_ko-plotting_data$event_level_median_wt)
  plotting_positions <- data.frame(plotting_data$reference, plotting_data$position, plotting_data$diff, plotting_data$feature_wt)
  colnames(plotting_positions) <- c('Reference', 'Position', 'Difference', 'Feature')

  plotting_positions <- subset(plotting_positions, Position>=initial_position)
  plotting_positions <- subset(plotting_positions, Position<=final_position)
  
  #Exclude SNPs and 10 positions before and after (21mer):
  if (length(exclude_SNP)!=0) {
    excluded_positions <- c()
    
    for (single_position in exclude_SNP){
      excluded_positions <- c(excluded_positions, seq(single_position-10,single_position+10))
    }
    
    plotting_positions <- subset(plotting_positions, !Position %in% unique(excluded_positions))
  } 
  
  #Calculate the threshold:
  threshold <- median(plotting_positions$Difference, na.rm = TRUE)

  #Calculate fold change:
  plotting_positions$Score <- plotting_positions$Difference/threshold
  plotting_positions$Modified_ZScore <- (plotting_positions$Score-median(plotting_positions$Score, na.rm = TRUE))/sd(plotting_positions$Score, na.rm = TRUE)
  
  #Format data for plotting:
  plotting_positions <- plotting_positions[,c(1,2,5,4,6)]
  
  #Extract significant positions:
  significant_positions <- subset(plotting_positions, Modified_ZScore>MZS_thr)
  
  return(list(plotting_positions,significant_positions))
  
}

tombo_processing <- function(sample_file, t_position, t_kmer, initial_position, final_position, MZS_thr, chr, exclude_SNP) {
  #Import data:
  sample <- read.delim(sample_file)

  if (nrow(sample)>0) {
    #Apply some filters and labels:
    sample$Feature <- 'Tombo'
    sample <- subset(sample, Coverage_Sample>30 & Coverage_IVT>30)
    colnames(sample) <- c('Reference', 'Chr', 'Position', 'Difference', 'Coverage_Sample', 'Coverage_IVT', 'Statistic_kmer',
                             'Feature')
    
    sample <- subset(sample, Chr == chr)
    sample <- subset(sample, Position >= initial_position)
    sample <- subset(sample, Position <= final_position)
    
    #Exclude SNPs and 10 positions before and after (21mer):
    if (length(exclude_SNP)!=0) {
      excluded_positions <- c()
      
      for (single_position in exclude_SNP){
        excluded_positions <- c(excluded_positions, seq(single_position-10,single_position+10))
      }
      
      sample <- subset(sample, !Position %in% unique(excluded_positions))
    } 
    
    #Calculate the thresholds:
    threshold_position <- median(sample$Difference, na.rm = TRUE)
    threshold_kmer <- median(sample$Statistic_kmer, na.rm = TRUE)
    
    #Calculate fold change:
    sample$Score <- sample$Difference/threshold_position
    sample$Score_kmer <- sample$Statistic_kmer/threshold_kmer
    sample$Modified_ZScore <- (sample$Score-median(sample$Score, na.rm = TRUE))/sd(sample$Score, na.rm = TRUE)
    sample$Modified_ZScore_kmer <- (sample$Score_kmer-median(sample$Score_kmer, na.rm = TRUE))/sd(sample$Score_kmer, na.rm = TRUE)
  
    #Filter columns to get data in plotting format: 
    plotting_positions <- sample[,c(1,3,9,8,11)]
    
    #Extract significant positions and kmers and then perform the intersection:
    positions <- subset(sample, Modified_ZScore > MZS_thr)
    kmer <- subset(sample, Modified_ZScore_kmer > MZS_thr)
    
    significant_positions <- join(kmer, positions, by = 'Reference', type = "inner")
    significant_positions <- significant_positions[,c(1,3,9,8,11)]
  
  } else {
    plotting_positions <- data.frame(Reference= character(), Position=integer(), Difference=double(), Feature=character())
    significant_positions <- data.frame(Reference= character(), Position=integer(), Difference=double(), Feature=character())
  }
  
  return(list(plotting_positions, significant_positions))
}

nanocomp_processing <- function(sample_file, nanocomp_metric, t_nanocomp, initial_position, final_position, MZS_thr, chr, exclude_SNP){
  #Import data:
  sample <- read.delim(sample_file)
  if (nrow(sample)>0) {
    
    #Transform metric:
    sample$stat <- log(sample$GMM_logit_pvalue_context_4)
    sample$log_stat <- (sample$stat)*(-1)
    
    #Prepare plotting data:
    sample$reference <- paste(sample$ref_id, sample$pos, sep='_')
    sample$Feature <- 'Nanocompore'
    
    sample <- sample[which(sample$ref_id==chr),]
    plotting_data <- sample[,c(27, 1, 26, 28)]
    colnames(plotting_data) <- c('Reference', 'Position', 'Difference', 'Feature')
    plotting_data <- subset(plotting_data, Position>=initial_position)
    plotting_data <- subset(plotting_data, Position <= final_position)
    
    #Exclude SNPs and 10 positions before and after (21mer):
    if (length(exclude_SNP)!=0) {
      excluded_positions <- c()
      
      for (single_position in exclude_SNP){
        excluded_positions <- c(excluded_positions, seq(single_position-10,single_position+10))
      }
      
      plotting_data <- subset(plotting_data, !Position %in% unique(excluded_positions))
    } 
    
    #Calculate the thresholds:
    threshold <- median(plotting_data$Difference, na.rm = TRUE)

    #Calculate fold change:
    plotting_data$Score <- plotting_data$Difference/threshold
    plotting_data$Modified_ZScore <- (plotting_data$Score-median(plotting_data$Score, na.rm = TRUE))/sd(plotting_data$Score, na.rm = TRUE)

    #Format data for plotting:
    plotting_data <- plotting_data[,c(1,2,5,4,6)]

    #Extract significant positions:
    significant_positions <- subset(plotting_data, Modified_ZScore > MZS_thr)
    
  } else {
    plotting_data <- data.frame(Reference= character(), Position=integer(), Difference=double(), Feature=character())
    significant_positions <- data.frame(Reference= character(), Position=integer(), Difference=double(), Feature=character())
}
  
  
  return(list(plotting_data, significant_positions))
}

barplot_plotting <- function (list_plotting, list_significant, output_name, MZS_thr, autoscaling, initial_pos, final_pos){
  
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
  initial_df$sample_f <- factor(initial_df$Feature)
  putative_positions$sample_f <- factor(putative_positions$Feature)
  
  #Plotting:
  pdf(file=paste(output_name,".pdf", sep = ""), bg = "transparent", width = 26, height = 16 )
  plot(ggplot(initial_df, aes(x=Position, y=Modified_ZScore, fill=Modified_ZScore)) + ggtitle(output_name) +
          geom_bar(data=subset(initial_df, Modified_ZScore < MZS_thr), stat= "identity", width=4, fill = "#dcdcdd") +
          new_scale_color() + xlim(initial_pos, final_pos) +
          geom_bar(data=subset(initial_df, Modified_ZScore >= MZS_thr), stat = "identity", width=4) + 
          scale_fill_gradient(low="#ff7f7f", 
                              high="#ff0000") + 
          theme_bw() +theme(plot.title = element_text(face = "bold", hjust = 0.5), text = element_text(size=25),
                            axis.text = element_text(size = 25), strip.text.y = element_text(size = 25),
                            legend.text=element_text(size=22)) + 
         facet_grid(sample_f ~ . , scales="fixed") )
  dev.off()
  
  if (autoscaling == TRUE) {
    #Plotting:
    pdf(file=paste(output_name,"_AUTOSCALE.pdf", sep = ""),bg = "transparent", width = 26, height = 16 )
    plot(ggplot(initial_df, aes(x=Position, y=Modified_ZScore, fill=Modified_ZScore)) + ggtitle(output_name) +
           geom_bar(data=subset(initial_df, Modified_ZScore < MZS_thr), stat= "identity", width=1.5, fill = "#dcdcdd") +
           new_scale_color() +
           geom_bar(data=subset(initial_df, Modified_ZScore >= MZS_thr), stat = "identity", width=1.5) + 
           scale_fill_gradient(low="#ff7f7f", 
                               high="#ff0000") + 
           theme_bw() +theme(plot.title = element_text(face = "bold", hjust = 0.5), text = element_text(size=25),
                             axis.text = element_text(size = 25), strip.text.y = element_text(size = 25),
                             legend.text=element_text(size=22)) + 
           facet_grid(sample_f ~ . , scales="free_y") )
    dev.off()
  }
  
}

Nanoconsensus_plotting <- function(data, supported_kmers, output_name) {
  #Extracting supported positions:
  supported_positions <- c()
  kmers_limits <- c()
  
  if (nrow(supported_kmers)!=0) {
    for (i in seq(1, nrow(supported_kmers))) {
      supported_positions <- c(supported_positions, seq(supported_kmers[i,2], supported_kmers[i,3]))
      kmers_limits <- c(kmers_limits, supported_kmers[i,2], supported_kmers[i,3])
    }
    
    #Plotting:
    data$Position <- data$Start+2
    pdf(file=paste(output_name,"NanoConsensus_Score.pdf", sep = "-"), bg = "transparent", width = 26, height = 16 )
    plot(ggplot(data, aes(x=Position, y=Merged_Score)) + ggtitle(output_name) + geom_bar(stat= "identity", width=2) + ylim(0,1) +
           geom_bar(data=subset(data, Position %in% supported_positions), stat= "identity", width=2, fill = "red") + 
           geom_label_repel(data=subset(data, Position %in% kmers_limits),aes(label = Position, y = Merged_Score), size = 8, label.size = 0.75) +
           ylab('NanoConsensus Score') +
           theme_bw() +theme(plot.title = element_text(face = "bold", hjust = 0.5), text = element_text(size=25),
                             axis.text = element_text(size = 25), strip.text.y = element_text(size = 25),
                             legend.text=element_text(size=22)))
    dev.off()
    
  } else {
    data$Position <- data$Start+2
    pdf(file=paste(output_name,"NanoConsensus_Score.pdf", sep = "-"), bg = "transparent", width = 26, height = 16 )
    plot(ggplot(data, aes(x=Position, y=Merged_Score)) + ggtitle(output_name) + geom_bar(stat= "identity", width=2) + ylim(0,1) +
           ylab('NanoConsensus Score') +
           theme_bw() +theme(plot.title = element_text(face = "bold", hjust = 0.5), text = element_text(size=25),
                             axis.text = element_text(size = 25), strip.text.y = element_text(size = 25),
                             legend.text=element_text(size=22)))
    dev.off()
  }
 
  
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
  
  if (length_object1 >= length_object2) {
    intersect_object <- subsetByOverlaps(GRange_object_2, GRange_object_1, minoverlap=1)
  } else {
    intersect_object <- subsetByOverlaps(GRange_object_1, GRange_object_2, minoverlap=1)
  }
  
  return(intersect_object)
  
}

draw_pairwise_venn_diagram <- function (group_1, group_2, intersect_12, groups, output_name){
  
  #Draw Venn Diagram:
  grid.newpage()
  venn.plot <- draw.pairwise.venn(group_1, group_2, intersect_12, 
                                  category = groups, fill = c("darksalmon", "dodgerblue"), cat.pos = c(0, 0), alpha = 0.5
  )
  
  # Writing to file
  png(filename = paste(output_name,'VennDiagram.png', sep="_"))
  grid.draw(venn.plot)
  dev.off()
  
}

draw_triple_venn_diagram <- function (group_1, group_2, group_3, intersect_12, intersect_13, intersect_23, intersect_123, groups, output_name){
  
  #Draw Venn Diagram:
  grid.newpage()
  venn.plot <- draw.triple.venn(group_1, group_2, group_3, intersect_12, intersect_23, intersect_13, 
                              intersect_123, category = groups, fill = c("darksalmon", "dodgerblue", "lightseagreen"), cat.pos = c(-45, 0, 45), alpha = 0.5
  )
  
  # Writing to file
  png(filename = paste(output_name,'VennDiagram.png', sep="_"))
  grid.draw(venn.plot)
  dev.off()
  
}

draw_venn_diagram <- function (group_1, group_2, group_3, group_4, intersect_12, intersect_13, intersect_14, intersect_23, intersect_24,
                  intersect_34, intersect_123, intersect_124, intersect_134, intersect_234, intersect_1234, groups, output_name){
  
  #Draw Venn Diagram:
  grid.newpage()
  venn.plot <- draw.quad.venn(group_1, group_2, group_3, group_4, intersect_12, intersect_13, intersect_14, intersect_23, intersect_24,
                              intersect_34, intersect_123, intersect_124, intersect_134, intersect_234, intersect_1234, 
                              category = groups, fill = c("darksalmon", "dodgerblue", "lightseagreen", "darkorange"), cat.pos = c(0, 0, 0, 0), alpha = 0.5
  )
  
  # Writing to file
  png(filename = paste(output_name,'VennDiagram.png', sep="_"))
  grid.draw(venn.plot)
  dev.off()
  
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
      
      #Loop over the kmer to find if Epinano identified it:
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

extracting_modified_ZScores <- function (GRange_supported_kmers, list_plotting, MZS_thr, summit, Consensus_score) {
  
  #Create vectors to store software data:
  epinano_rawScore <- c()
  nanopolish_rawScore <- c()
  tombo_rawScore <- c()
  nanocompore_rawScore <- c()
  
  epinano_modifiedScore <- c()
  nanopolish_modifiedScore <- c()
  tombo_modifiedScore <- c()
  nanocompore_modifiedScore <- c()
  
  epinano_status <- c()
  nanopolish_status <- c()
  tombo_status <- c()
  nanocompore_status <- c()
  
  #Parse data into a data frame:
  positions_df <- data.frame(start(GRange_supported_kmers), end(GRange_supported_kmers))
  colnames(positions_df) <- c('Start', 'End')
  positions_df$Chr <- seqlevels(GRange_supported_kmers)
  positions_df <- positions_df[,c(3,1,2)]
  
  #Extracting scores and software status:
  epinano_data <- extracting_status(positions_df, 1, summit, MZS_thr)
  nanopolish_data <- extracting_status(positions_df, 2, summit, MZS_thr)
  tombo_data <- extracting_status(positions_df, 3, summit, MZS_thr)
  nanocompore_data <- extracting_status(positions_df, 4, summit, MZS_thr)
  
  #Add data to the final dataframe:
  positions_df$Epinano_RawScore <- epinano_data$rawScore
  positions_df$Nanopolish_RawScore <- nanopolish_data$rawScore
  positions_df$Tombo_RawScore <- tombo_data$rawScore
  positions_df$Nanocompore_RawScore <- nanocompore_data$rawScore
  
  positions_df$Epinano_Score <- epinano_data$modifiedScore
  positions_df$Nanopolish_Score <- nanopolish_data$modifiedScore
  positions_df$Tombo_Score <- tombo_data$modifiedScore
  positions_df$Nanocompore_Score <- nanocompore_data$modifiedScore
  
  positions_df$Epinano_Status <- epinano_data$status
  positions_df$Nanopolish_Status <- nanopolish_data$status
  positions_df$Tombo_Status <- tombo_data$status
  positions_df$Nanocompore_Status <- nanocompore_data$status
  
  positions_NanoConsensus <- c()
  ##Calculate the merged_score:
  #Re-scaling:
  if (summit == F){
    data <- data.frame(positions_df$Epinano_Score, positions_df$Nanopolish_Score, positions_df$Tombo_Score, positions_df$Nanocompore_Score)
    

    #Re-scale Modified Z-Scores between 0 and 1:
    for (i in seq(1:length(data))) {
      data[,i] <- rescale(unlist(data[i]), to=c(0,1))
    
    }
    
    data[is.na(data)] <- 0
    
    #Rescale outputs 0.5 when the software gives the same MZS for all positions - correcting for it if needed:
    if (length(unique(data$positions_df.Epinano_Score)) == 1 || all(is.na(unique(data$positions_df.Epinano_Score)))) {
      data$positions_df.Epinano_Score <- 0 
    } 
    
    if (length(unique(data$positions_df.Nanopolish_Score)) == 1 || all(is.na(unique(data$positions_df.Nanopolish_Score)))) {
      data$positions_df.Nanopolish_Score <- 0 
    }

    if (length(unique(data$positions_df.Tombo_Score)) == 1 || all(is.na(unique(data$positions_df.Tombo_Score)))) {
      data$positions_df.Tombo_Score <- 0 
    } 
    
    if (length(unique(data$positions_df.Nanocompore_Score)) == 1 || all(is.na(unique(data$positions_df.Nanocompore_Score)))) {
      data$positions_df.Nanocompore_Score <- 0 
    }
    
    print(head(data))
    #Calculate NanoConsensus score:
    positions_df$Merged_Score <- apply(data,1,median,na.rm = TRUE)
    threshold <- Consensus_score*median(positions_df$Merged_Score,  na.rm = TRUE)
    print(threshold)
    positions_NanoConsensus <- subset(positions_df, Merged_Score >= threshold)
  }

  return(list(positions_df, positions_NanoConsensus))
}

kmer_analysis <- function (all_ranges, fasta_file, output_name) {
  print('Kmer analysis')
  kmer_data <- extract_kmers(all_ranges, fasta_file)
  all_ranges$Kmer <- kmer_data[[1]]
  all_ranges$RRACH_motif <- kmer_data[[2]]
  all_ranges <- all_ranges[order(all_ranges$Start, decreasing = FALSE),]
  
  #Merging data per kmer: 
  write.table(all_ranges, file = output_name, sep = '\t', row.names = FALSE)
}

analysis_significant_positions <- function (list_significant, list_plotting, fasta_file, output_name, initial_position, final_position, MZS_thr, Consensus_score) {
  epinano <- list_significant[[1]]
  nanopolish <- list_significant[[2]]
  tombo <- list_significant[[3]]
  nanocompore <- list_significant[[4]]
  
  methods <- list(epinano, nanopolish, tombo, nanocompore)
  methods_name <- c('Epinano', 'Nanopolish', 'Tombo', 'Nanocompore')
  
  #Create grRange objects with kmers per each method: 
  print('Transforming data into GRange objects')
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
        grList <- pc(grList,grNew)

      }
      
      if(elementNROWS(grList) == 1) {
        grList <- grList
      } else {
        grList <- do.call(c, grList)
      }
        
      assign(paste('gr',methods_name[j],sep=""), reduce(unique(grList)))
    
    } else {
      assign(paste('gr',methods_name[j],sep=""), unique(grList))
    }
    
  }
  
  ##Perform intersections:
  #Check how many elements are in each GRange object and if it is null, create an empty one:
  print('Perform intersections')
  
  if (is.null(grEpinano)==TRUE){
    grEpinano <- GRanges()
    n1 <- 0
  } else {
    n1 <- length(grEpinano)
  }
  
  if (is.null(grNanopolish)==TRUE){
    grNanopolish <- GRanges()
    n2 <- 0
  } else {
    n2 <- length(grNanopolish)
  }
  
  if (is.null(grTombo)==TRUE){
    grTombo <- GRanges()
    n3 <- 0
  } else {
    n3 <- length(grTombo)
  }
  
  if (is.null(grNanocompore)==TRUE){
    grNanocompore <- GRanges()
    n4 <- 0
  } else {
    n4 <- length(grNanocompore)
  }
  
  print(c(n1,n2,n3,n4))
  
  if (n1 != 0 & n2 != 0 & n3 != 0 & n4 == 0 ) {
    #Overlappings: checking which software has identified less significant positions and then it uses it as query
    intersect_12 <- overlapping_GRobjects(grEpinano, grNanopolish, n1, n2)
    length_intersect_12 <- extract_length_from_GRobjects(intersect_12)
    
    intersect_13 <- overlapping_GRobjects(grEpinano, grTombo, n1, n3)
    length_intersect_13 <- extract_length_from_GRobjects(intersect_13)
    
    intersect_23 <- overlapping_GRobjects(grNanopolish, grTombo, n2, n3)
    length_intersect_23 <- extract_length_from_GRobjects(intersect_23)

    intersect_123 <- overlapping_GRobjects(intersect_12, grTombo, length_intersect_12, n3)
    length_intersect_123 <- extract_length_from_GRobjects(intersect_123)

    #Venn Diagram:  
    print(c(n1, n2, n3, length_intersect_12, length_intersect_13, length_intersect_23, length_intersect_123))
    methods_name <-  c('Epinano', 'Nanopolish', 'Tombo')
    draw_triple_venn_diagram(n1, n2, n3, length_intersect_12, length_intersect_13, length_intersect_23, length_intersect_123, methods_name, output_name)
  
    #Extract kmers supported by two or more softwares: 
    supported_kmers <- unique(c(unlist(intersect_12), unlist(intersect_13), unlist(intersect_23), unlist(intersect_123)))
    
  } else if (n1 != 0 & n2 == 0 & n3 != 0 & n4 != 0 ) {
    #Overlappings: checking which software has identified less significant positions and then it uses it as query
    intersect_13 <- overlapping_GRobjects(grEpinano, grTombo, n1, n3)
    length_intersect_13 <- extract_length_from_GRobjects(intersect_13)
                                                         
    intersect_14 <- overlapping_GRobjects(grEpinano, grNanocompore, n1, n4)
    length_intersect_14 <- extract_length_from_GRobjects(intersect_14)                                                     
    
    intersect_34 <- overlapping_GRobjects(grTombo, grNanocompore, n3, n4)
    length_intersect_34 <- extract_length_from_GRobjects(intersect_34) 
    
    intersect_134 <- overlapping_GRobjects(intersect_13, grNanocompore, length_intersect_13, n4)
    length_intersect_134 <- extract_length_from_GRobjects(intersect_134)
    
    #Venn Diagram:  
    print(c(n1, n3, n4, length_intersect_13, length_intersect_14, length_intersect_34, length_intersect_134))
    methods_name <-  c('Epinano', 'Tombo', 'Nanocompore')
    draw_triple_venn_diagram(n1, n3, n4, length_intersect_13, length_intersect_14, length_intersect_34, length_intersect_134, methods_name, output_name)
    
    #Extract kmers supported by two or more softwares: 
    supported_kmers <- unique(c(unlist(intersect_13), unlist(intersect_14), unlist(intersect_34), unlist(intersect_134)))
    
  } else if (n1 == 0 & n2 != 0 & n3 != 0 & n4 != 0 ) {
    #Overlappings: checking which software has identified less significant positions and then it uses it as query
    intersect_23 <- overlapping_GRobjects(grNanopolish, grTombo, n2, n3)
    length_intersect_23 <- extract_length_from_GRobjects(intersect_23)

    intersect_24 <- overlapping_GRobjects(grNanopolish, grNanocompore, n2, n4)
    length_intersect_24 <- extract_length_from_GRobjects(intersect_24)
    
    intersect_34 <- overlapping_GRobjects(grTombo, grNanocompore, n3, n4)
    length_intersect_34 <- extract_length_from_GRobjects(intersect_34)
    
    intersect_234 <- overlapping_GRobjects(intersect_23, grNanocompore, length_intersect_23, n4)
    length_intersect_234 <- extract_length_from_GRobjects(intersect_234)
    
    #Venn Diagram:  
    print(c(n2, n3, n4, length_intersect_23, length_intersect_24, length_intersect_34, length_intersect_234))
    methods_name <-  c('Nanopolish', 'Tombo', 'Nanocompore')
    draw_triple_venn_diagram(n2, n3, n4, length_intersect_23, length_intersect_24, length_intersect_34, length_intersect_234, methods_name, output_name)
    
    #Extract kmers supported by two or more softwares: 
    supported_kmers <- unique(c(unlist(intersect_23), unlist(intersect_24), unlist(intersect_34), unlist(intersect_234)))
    
  } else if (n1 != 0 & n2 != 0 & n3 == 0 & n4 != 0 ) {
    #Overlappings: checking which software has identified less significant positions and then it uses it as query
    intersect_12 <- overlapping_GRobjects(grEpinano, grNanopolish, n1, n2)
    length_intersect_12 <- extract_length_from_GRobjects(intersect_12)
    
    intersect_14 <- overlapping_GRobjects(grEpinano, grNanocompore, n1, n4)
    length_intersect_14 <- extract_length_from_GRobjects(intersect_14)

    intersect_24 <- overlapping_GRobjects(grNanopolish, grNanocompore, n2, n4)
    length_intersect_24 <- extract_length_from_GRobjects(intersect_24)

    intersect_124 <- overlapping_GRobjects(intersect_12, grNanocompore, length_intersect_12, n4)
    length_intersect_124 <- extract_length_from_GRobjects(intersect_124)
    
    #Venn Diagram:  
    print(c(n1, n2, n4, length_intersect_12, length_intersect_14, length_intersect_24, length_intersect_124))
    methods_name <-  c('Epinano', 'Nanopolish', 'Nanocompore')
    draw_triple_venn_diagram(n1, n2, n4, length_intersect_12, length_intersect_14, length_intersect_24, length_intersect_124, methods_name, output_name)
    
    #Extract kmers supported by two or more softwares: 
    supported_kmers <- unique(c(unlist(intersect_12), unlist(intersect_14), unlist(intersect_24), unlist(intersect_124)))

  } else if (n1 == 0 & n2 == 0) {
    #Overlappings: checking which software has identified less significant positions and then it uses it as query
    intersect_34 <- overlapping_GRobjects(grTombo, grNanocompore, n3, n4)
    length_intersect_34 <- extract_length_from_GRobjects(intersect_34)
    
    #Venn Diagram:  
    print(c(n3, n4, length_intersect_34))
    methods_name <-  c('Tombo', 'Nanocompore')
    draw_pairwise_venn_diagram(n3, n4, length_intersect_34, methods_name, output_name)
    
    #Extract kmers supported by two or more softwares: 
    supported_kmers <- unlist(intersect_34)
    
  } else if (n3 == 0 & n4 == 0) {
    #Overlappings: checking which software has identified less significant positions and then it uses it as query
    intersect_12 <- overlapping_GRobjects(grEpinano, grNanopolish, n1, n2)
    length_intersect_12 <- extract_length_from_GRobjects(intersect_12)
    
    #Venn Diagram:  
    print(c(n1, n2, length_intersect_12))
    methods_name <-  c('Epinano', 'Nanopolish')
    draw_pairwise_venn_diagram(n1, n2, length_intersect_12, methods_name, output_name)
    
    #Extract kmers supported by two or more softwares: 
    supported_kmers <- unlist(intersect_12)
    
  } else if (n2 == 0 & n3 == 0) { 
    #Overlappings: checking which software has identified less significant positions and then it uses it as query
    intersect_14 <- overlapping_GRobjects(grEpinano, grNanocompore, n1, n4)
    length_intersect_14 <- extract_length_from_GRobjects(intersect_14)
    
    #Venn Diagram:  
    print(c(n1, n4, length_intersect_14))
    methods_name <-  c('Epinano', 'Nanocompore')
    draw_pairwise_venn_diagram(n1, n4, length_intersect_14, methods_name, output_name)
    
    #Extract kmers supported by two or more softwares: 
    supported_kmers <- unlist(intersect_14)
    
  } else if (n1 == 0 & n3 == 0) { 
    #Overlappings: checking which software has identified less significant positions and then it uses it as query
    intersect_24 <- overlapping_GRobjects(grNanopolish, grNanocompore, n2, n4)
    length_intersect_24 <- extract_length_from_GRobjects(intersect_24)
    
    #Venn Diagram:  
    print(c(n2, n4, length_intersect_24))
    methods_name <-  c('Nanopolish', 'Nanocompore')
    draw_pairwise_venn_diagram(n2, n4, length_intersect_24, methods_name, output_name)
    
    #Extract kmers supported by two or more softwares: 
    supported_kmers <- unlist(intersect_24)
    
  } else {
    #Overlappings: checking which software has identified less significant positions and then it uses it as query
    intersect_12 <- overlapping_GRobjects(grEpinano, grNanopolish, n1, n2)
    length_intersect_12 <- extract_length_from_GRobjects(intersect_12)
    
    intersect_13 <- overlapping_GRobjects(grEpinano, grTombo, n1, n3)
    length_intersect_13 <- extract_length_from_GRobjects(intersect_13)
    
    intersect_14 <- overlapping_GRobjects(grEpinano, grNanocompore, n1, n4)
    length_intersect_14 <- extract_length_from_GRobjects(intersect_14)
    
    intersect_23 <- overlapping_GRobjects(grNanopolish, grTombo, n2, n3)
    length_intersect_23 <- extract_length_from_GRobjects(intersect_23)
    
    intersect_24 <- overlapping_GRobjects(grNanopolish, grNanocompore, n2, n4)
    length_intersect_24 <- extract_length_from_GRobjects(intersect_24)
    
    intersect_34 <- overlapping_GRobjects(grTombo, grNanocompore, n3, n4)
    length_intersect_34 <- extract_length_from_GRobjects(intersect_34)
    
    intersect_123 <- overlapping_GRobjects(intersect_12, grTombo, length_intersect_12, n3)
    length_intersect_123 <- extract_length_from_GRobjects(intersect_123)
    
    intersect_124 <- overlapping_GRobjects(intersect_12, grNanocompore, length_intersect_12, n4)
    length_intersect_124 <- extract_length_from_GRobjects(intersect_124)

    intersect_134 <- overlapping_GRobjects(intersect_13, grNanocompore, length_intersect_13, n4)
    length_intersect_134 <- extract_length_from_GRobjects(intersect_134)

    intersect_234 <- overlapping_GRobjects(intersect_23, grNanocompore, length_intersect_23, n4)
    length_intersect_234 <- extract_length_from_GRobjects(intersect_234)
    
    intersect_1234 <- overlapping_GRobjects(intersect_12, intersect_34, length_intersect_12, length_intersect_34)
    length_intersect_1234 <- extract_length_from_GRobjects(intersect_1234)

    #Venn Diagram:  
    print(c(n1, n2, n3, n4, length_intersect_12, length_intersect_13, length_intersect_14, length_intersect_23, length_intersect_24,
            length_intersect_34, length_intersect_123, length_intersect_124, length_intersect_134, length_intersect_234, length_intersect_1234))
    draw_venn_diagram(n1, n2, n3, n4, length_intersect_12, length_intersect_13, length_intersect_14, length_intersect_23, length_intersect_24,
                      length_intersect_34, length_intersect_123, length_intersect_124, length_intersect_134, length_intersect_234, length_intersect_1234, methods_name, output_name)
    
    #Extract kmers supported by two or more softwares: 
    supported_kmers <- unique(c(unlist(intersect_12), unlist(intersect_13), unlist(intersect_14), unlist(intersect_23), unlist(intersect_24),
                                unlist(intersect_34), unlist(intersect_123), unlist(intersect_124), unlist(intersect_134), unlist(intersect_234), unlist(intersect_1234)))
  
  }
  
  ##Kmer analysis: 
  #Analysis of all kmers across the chromosome:
  all_kmers_raw <- GRanges(seqnames = chr, ranges = IRanges(initial_position:(final_position-4), end = (initial_position+4):final_position))
  all_kmers <- extracting_modified_ZScores(all_kmers_raw, list_plotting, MZS_thr, FALSE, Consensus_score)
  kmer_analysis(all_kmers[[1]], fasta_file, paste(output_name,'Raw_kmers.txt', sep='_'))
  
  #Analyse the supported kmers - only if they are present:
  filtered_supported_kmers <- overlapping_GRobjects(reduce(supported_kmers), GRanges(seqnames=all_kmers[[2]][,c('Chr')],ranges=IRanges(all_kmers[[2]][,c('Start')], end = all_kmers[[2]][,c('End')])),1,2)
  if(extract_length_from_GRobjects(filtered_supported_kmers)!=0){
    all_ranges <- extracting_modified_ZScores(filtered_supported_kmers, list_plotting, MZS_thr, TRUE, Consensus_score)
    kmer_analysis(all_ranges[[1]], fasta_file, paste(output_name,'Supported_kmers.txt', sep='_'))
    
    #Plot NanoConsensus score across transcripts:
    Nanoconsensus_plotting(all_kmers[[1]], all_ranges[[1]], output_name)
    
  } else {
    all_ranges <- data.frame()
    #Plot NanoConsensus score across transcripts:
    Nanoconsensus_plotting(all_kmers[[1]], all_ranges, output_name)
  }

  
  
}
