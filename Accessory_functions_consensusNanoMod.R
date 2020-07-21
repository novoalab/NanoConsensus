###Script which contains multiple R functions used to generate consensus putative modified positions from NanoMod results. 

#Processing Epinano results:
epinano_processing <- function(sample_file, ivt_file, initial_position, final_position, MZS_thr) {
  
  #Import and clean data:
  sample <- read.csv(sample_file,stringsAsFactors = FALSE)
  sample <- subset(sample, cov>50)
  sample <- subset(sample, pos>=initial_position) 
  sample <- subset(sample, pos<=final_position)
  sample$reference <- paste(sample$X.Ref, sample$pos, sep='_')
  sample$Difference <- as.numeric(sample$mis)+as.numeric(sample$ins)+as.numeric(sample$del)
  sample <- sample[,c(1,2,12,11)]
  colnames(sample) <- c('Reference', 'Position', 'Difference_sample', 'Merge')
  
  ivt <- read.csv(ivt_file,stringsAsFactors = FALSE)
  ivt <- subset(ivt, cov>50)
  ivt <- subset(ivt, pos>=initial_position) 
  ivt <- subset(ivt, pos<=final_position)
  ivt$reference <- paste(ivt$X.Ref, ivt$pos, sep='_')
  ivt$Difference <- as.numeric(ivt$mis)+as.numeric(ivt$ins)+as.numeric(ivt$del)
  ivt <- ivt[,c(1,2,12,11)]
  colnames(ivt) <- c('Reference', 'Position', 'Difference_IVT', 'Merge')
  
  
  #Join both dataframes and clean unecessary columns:
  plotting_positions <- join(sample, ivt, by="Merge")
  plotting_positions$Difference <- abs(plotting_positions$Difference_sample - plotting_positions$Difference_IVT)
  plotting_positions$Feature <- "Epinano"
  plotting_positions <- plotting_positions[,c(4,2,8,9)]

  #Calculate the threshold:
  threshold <- median(plotting_positions$Difference, na.rm = TRUE)
  
  #Calculate fold change and re-order:
  plotting_positions$Score <- rescale((plotting_positions$Difference/threshold), to=c(0,1))
  plotting_positions$Modified_ZScore <- (plotting_positions$Score-median(plotting_positions$Score))/sd(plotting_positions$Score)
  
  plotting_positions <- plotting_positions[,c(1,2,5,4,6)]
  colnames(plotting_positions) <- c('Reference', 'Position', 'Score', 'Feature', 'Modified_ZScore')
  
  #Extract significant positions based on the specific threshold:
  significant_positions <- subset(plotting_positions, Modified_ZScore>MZS_thr)

  return(list(plotting_positions,significant_positions))
}

nanopolish_processing <- function(sample_file, ivt_file, initial_position, final_position, MZS_thr) {
  #Import data:
  sample <- read.delim(sample_file)
  
  #Add sample information:
  sample$read_name <- 'Nanopolish'
  sample <- subset(sample, coverage>50)
  colnames(sample)<- c("contig_wt","position","reference_kmer_wt", "feature_wt", "event_level_median_wt", 'coverage')
  sample$reference <- paste(sample$contig_wt, sample$position, sep='_')
  
  #Import KO: 
  raw_data_ivt <-read.delim(ivt_file)
  raw_data_ivt$read_name <- 'IVT'
  raw_data_ivt <- subset(raw_data_ivt, coverage>50)
  colnames(raw_data_ivt)<- c("contig_ko","position","reference_kmer_ko", "feature", "event_level_median_ko", 'coverage')
  raw_data_ivt$reference <- paste(raw_data_ivt$contig_ko, raw_data_ivt$position, sep='_')
  
  #Join tables, calculate differences between means/medians:
  plotting_data <- join(sample, raw_data_ivt, by="reference", type = "inner")
  plotting_data$diff <- abs(plotting_data$event_level_median_ko-plotting_data$event_level_median_wt)
  plotting_positions <- data.frame(plotting_data$reference, plotting_data$position, plotting_data$diff, plotting_data$feature_wt)
  colnames(plotting_positions) <- c('Reference', 'Position', 'Difference', 'Feature')
  plotting_positions <- subset(plotting_positions, Position>=initial_position)
  plotting_positions <- subset(plotting_positions, Position<=final_position)
  
  #Calculate the threshold:
  threshold <- median(plotting_positions$Difference, na.rm = TRUE)
  
  #Calculate fold change:
  plotting_positions$Score <- rescale((plotting_positions$Difference/threshold), to=c(0,1))
  plotting_positions$Modified_ZScore <- (plotting_positions$Score-median(plotting_positions$Score))/sd(plotting_positions$Score)
  
  #Format data for plotting:
  plotting_positions <- plotting_positions[,c(1,2,5,4,6)]
  
  #Extract significant positions:
  significant_positions <- subset(plotting_positions, Modified_ZScore>MZS_thr)
  
  return(list(plotting_positions,significant_positions))
  
}

tombo_processing <- function(sample_file, t_position, t_kmer, initial_position, final_position, MZS_thr) {
  #Import data:
  sample <- read.delim(sample_file)

  #Apply some filters and labels:
  sample$Feature <- 'Tombo'
  sample <- subset(sample, Coverage_Sample>50 & Coverage_IVT>50)
  colnames(sample) <- c('Reference', 'Chr', 'Position', 'Difference', 'Coverage_Sample', 'Coverage_IVT', 'Statistic_kmer',
                           'Feature')
  
  sample <- subset(sample, Position >= initial_position)
  sample <- subset(sample, Position <= final_position)
  
  #Calculate the thresholds:
  threshold_position <- median(sample$Difference)
  threshold_kmer <- median(sample$Statistic_kmer, na.rm = TRUE)
  
  #Calculate fold change:
  sample$Score <- rescale((sample$Difference/threshold_position), to=c(0,1))
  sample$Score_kmer <- rescale((sample$Statistic_kmer/threshold_kmer), to=c(0,1))
  sample$Modified_ZScore <- (sample$Score-median(sample$Score))/sd(sample$Score)
  sample$Modified_ZScore_kmer <- (sample$Score_kmer-median(sample$Score_kmer, na.rm = TRUE))/sd(sample$Score_kmer, na.rm = TRUE)

  #Filter columns to get data in plotting format: 
  plotting_positions <- sample[,c(1,3,9,8,11)]
  
  #Extract significant positions and kmers and then perform the intersection:
  positions <- subset(sample, Modified_ZScore > MZS_thr)
  kmer <- subset(sample, Modified_ZScore_kmer > (MZS_thr))
  
  significant_positions <- join(kmer, positions, by = 'Reference', type = "inner")
  significant_positions <- significant_positions[,c(1,3,9,8,11)]
  
  return(list(plotting_positions, significant_positions))
}


nanocomp_processing <- function(sample_file, nanocomp_metric, t_nanocomp, initial_position, final_position, MZS_thr){
  #Import data:
  sample <- read.delim(sample_file)
  if (nrow(sample)>0) {
    #Transform metric:
    sample$stat <- 1/(sample$GMM_logit_pvalue_context_4)
    sample$log_stat <- log(sample$stat)
    
    #Prepare plotting data:
    sample$reference <- paste(sample$ref_id, sample$pos, sep='_')
    sample$Feature <- 'Nanocompore'
    plotting_data <- sample[,c(27, 1, 26, 28)]
    colnames(plotting_data) <- c('Reference', 'Position', 'Difference', 'Feature')
    plotting_data <- subset(plotting_data, Position>=initial_position)
    plotting_data <- subset(plotting_data, Position <= final_position)
    
    #Calculate the thresholds:
    threshold <- median(plotting_data$Difference, na.rm = TRUE)

    #Calculate fold change:
    plotting_data$Score <- rescale((plotting_data$Difference/threshold), to=c(0,1))
    plotting_data$Modified_ZScore <- (plotting_data$Score-median(plotting_data$Score, na.rm = TRUE))/sd(plotting_data$Score, na.rm = TRUE)
    
    #Format data for plotting:
    plotting_data <- plotting_data[,c(1,2,5,4,6)]
    
    #Extract significant positions:
    significant_positions <- subset(plotting_data, Modified_ZScore > MZS_thr)
    
    return(list(plotting_data, significant_positions))
  }
}

subset_by_chr <- function(list_methods, chr){
 for (i in 1:length(list_methods)){
   list_methods[[i]] <- subset(list_methods[[i]], startsWith(as.character(Reference), chr))
 }
 return(list_methods)
}

barplot_plotting <- function (list_plotting, list_significant, output_name, MZS_thr){
  
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
  png(file=paste(output_name,".png", sep = ""),bg = "transparent", height=1480, width = 2640)
  plot(ggplot(initial_df, aes(x=Position, y=Modified_ZScore, fill=Modified_ZScore)) + ggtitle(output_name) +
          geom_bar(data=subset(initial_df, Modified_ZScore < MZS_thr), stat= "identity", width=1, fill = "#dcdcdd") +
          new_scale_color() +
          geom_bar(data=subset(initial_df, Modified_ZScore > MZS_thr), stat = "identity", width=1) + 
          scale_fill_gradient(low="#dcdcdd", 
                              high="#ff0000") + 
          theme_bw() +theme(plot.title = element_text(face = "bold", hjust = 0.5), text = element_text(size=20),
                            axis.text = element_text(size = 20), strip.text.y = element_text(size = 20)) + 
          facet_grid(sample_f ~ . , scales="fixed") )
  dev.off()
  
}

extract_length_from_GRobjects <- function(GRange_object) {

  if (length(GRange_object)==0){
    n_length <- 0
  } else {
    n_length <- elementNROWS(GRange_object)
  }
  
  return(n_length)
  
}

overlapping_GRobjects <- function(GRange_object_1, GRange_object_2, length_object1, length_object2) {
  if (length_object1 > length_object2) {
    intersect_object <- subsetByOverlaps(GRange_object_2, GRange_object_1, minoverlap = 1)
  } else {
    intersect_object <- subsetByOverlaps(GRange_object_1, GRange_object_2, minoverlap = 1)
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
                              intersect_123, category = groups, fill = c("darksalmon", "dodgerblue", "lightseagreen"), cat.pos = c(0, 0, 0), alpha = 0.5
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

#merging_data_per_kmer <- function (data) {
  
  #Merge data from the same kmer:
#  data$Ref <- paste(data$Chr, data$Start, sep="_")
#  final_data <- data.frame()
#  for (i in 1:nrow(data)){
#    if (i != nrow(data)) {
#      if (data$Ref[i] != data$Ref[i+1]) {
#        final_data <- rbind(final_data, data[i,])
#
#      }
#    } else {
#      final_data <- rbind(final_data, data[i,])
#    }
#  }
  
#  return(final_data[,c(1,2,3,4,5,6,7,8,9)])
#}

extracting_software_status <- function() {
  
}

extracting_modified_ZScores <- function (GRange_supported_kmers, list_plotting) {
  
  #Create vectors to store software data:
  epinano_rawScore <- c()
  nanopolish_rawScore <- c()
  tombo_rawScore <- c()
  nanocompore_rawScore <- c()
  
  epinano_modifiedScore <- c()
  nanopolish_modifiedScore <- c()
  tombo_modifiedScore <- c()
  nanocompore_modifiedScore <- c()
  
  #Parse data into a data frame:
  positions_df <- data.frame(start(GRange_supported_kmers), end(GRange_supported_kmers))
  colnames(positions_df) <- c('Start', 'End')
  positions_df$Chr <- seqlevels(GRange_supported_kmers)
  positions_df <- positions_df[,c(3,1,2)]

  for (single_pos in positions_df$Start){
    position <- single_pos+2
    epinano_rawScore <- c(epinano_rawScore, list_plotting[[1]][which(list_plotting[[1]]$Position == position), 3])
    epinano_modifiedScore <- c(epinano_modifiedScore, list_plotting[[1]][which(list_plotting[[1]]$Position == position), 5])
    
    nanopolish_rawScore <- c(nanopolish_rawScore, list_plotting[[2]][which(list_plotting[[2]]$Position == position), 3])
    nanopolish_modifiedScore <- c(nanopolish_modifiedScore, list_plotting[[2]][which(list_plotting[[2]]$Position == position), 5])
    
    
    if (length(list_plotting[[3]][which(list_plotting[[3]]$Position == position), 3]) == 0) {
      tombo_rawScore <- c(tombo_rawScore, 'NA')
      tombo_modifiedScore <- c(tombo_modifiedScore, 'NA')
      
    } else {
      tombo_rawScore <- c(tombo_rawScore, list_plotting[[3]][which(list_plotting[[3]]$Position == position), 3])
      tombo_modifiedScore <- c(tombo_modifiedScore, list_plotting[[3]][which(list_plotting[[3]]$Position == position), 5])
      
    }

    nanocompore_rawScore <- c(nanocompore_rawScore, list_plotting[[4]][which(list_plotting[[4]]$Position == position), 3])
    nanocompore_modifiedScore <- c(nanocompore_modifiedScore, list_plotting[[4]][which(list_plotting[[4]]$Position == position), 5])
    
  }
  
  #Add data to the final dataframe:
  positions_df$Epinano_RawScore <- epinano_rawScore
  positions_df$Nanopolish_RawScore <- nanopolish_rawScore
  positions_df$Tombo_RawScore <- tombo_rawScore
  positions_df$Nanocompore_RawScore <- nanocompore_rawScore
  
  positions_df$Epinano_Score <- epinano_modifiedScore
  positions_df$Nanopolish_Score <- nanopolish_modifiedScore
  positions_df$Tombo_Score <- tombo_modifiedScore
  positions_df$Nanocompore_Score <- nanocompore_modifiedScore
  
  return(positions_df)
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

analysis_significant_positions <- function (list_significant, list_plotting, fasta_file, output_name, initial_position, final_position) {
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
        grNew <- GRanges(seqnames=chr,ranges=IRanges(start=as.integer(features[[1]][length(features[[1]])])-2, width=5))
        grList <- pc(grList,grNew)
        
      }
    
      assign(paste('gr',methods_name[j],sep=""), unique(grList))
    
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
    n1 <- elementNROWS(grEpinano)
  }
  
  if (is.null(grNanopolish)==TRUE){
    grNanopolish <- GRanges()
    n2 <- 0
  } else {
    n2 <- elementNROWS(grNanopolish)
  }
  
  if (is.null(grTombo)==TRUE){
    grTombo <- GRanges()
    n3 <- 0
  } else {
    n3 <- elementNROWS(grTombo)
  }
  
  if (is.null(grNanocompore)==TRUE){
    grNanocompore <- GRanges()
    n4 <- 0
  } else {
    n4 <- elementNROWS(grNanocompore)
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
  all_kmers <- extracting_modified_ZScores(all_kmers_raw, list_plotting)
  kmer_analysis(all_kmers, fasta_file, paste(output_name,'Raw_kmers.txt', sep='_'))

  #Analyse the supported kmers:
  all_ranges <- extracting_modified_ZScores(supported_kmers, list_plotting) 
  kmer_analysis(all_ranges, fasta_file, paste(output_name,'Supported_kmers.txt', sep='_'))
  
}
