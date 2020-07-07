###Script which contains multiple R functions used to generate consensus putative modified positions from NanoMod results. 

#Processing Epinano results:
epinano_processing <- function(sample_file, ivt_file, initial_position, final_position, FC_thr) {
  
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
  threshold <- mean(plotting_positions$Difference, na.rm = TRUE)
  
  #Calculate fold change and re-order:
  plotting_positions$Fold_change <- log2((plotting_positions$Difference/threshold))
  plotting_positions <- plotting_positions[,c(1,2,5,4)]
  colnames(plotting_positions) <- c('Reference', 'Position', 'Fold_change', 'Feature')
  
  #Extract significant positions based on the specific threshold:
  significant_positions <- subset(plotting_positions, Fold_change>=FC_thr)

  return(list(plotting_positions,significant_positions))
}

nanopolish_processing <- function(sample_file, ivt_file, initial_position, final_position, FC_thr) {
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
  threshold <- mean(plotting_positions$Difference, na.rm = TRUE)
  
  #Calculate fold change:
  plotting_positions$Fold_change <- log2((plotting_positions$Difference/threshold))
  
  #Format data for plotting:
  plotting_positions <- plotting_positions[,c(1,2,5,4)]
  
  #Extract significant positions:
  significant_positions <- subset(plotting_positions, Fold_change>=FC_thr)
  
  return(list(plotting_positions,significant_positions))
  
}

tombo_processing <- function(sample_file, t_position, t_kmer, initial_position, final_position, FC_thr) {
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
  threshold_position <- mean(sample$Difference)
  threshold_kmer <- mean(sample$Statistic_kmer, na.rm = TRUE)
  
  #Calculate fold change:
  sample$Fold_change <- log2((sample$Difference/threshold_position))
  sample$Fold_change_kmer <- log2((sample$Statistic_kmer/threshold_kmer))
  
  #Filter columns to get data in plotting format: 
  plotting_positions <- sample[,c(1,3,9,8)]
  
  #Extract significant positions and kmers and then perform the intersection:
  positions <- subset(sample, Fold_change >= FC_thr)
  kmer <- subset(sample, Fold_change_kmer >= FC_thr)
  
  significant_positions <- join(kmer, positions, by = 'Reference', type = "inner")
  significant_positions <- significant_positions[,c(1,3,9,8)]
  
  return(list(plotting_positions, significant_positions))
}


nanocomp_processing <- function(sample_file, nanocomp_metric, t_nanocomp, initial_position, final_position, FC_thr){
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
    threshold <- mean(plotting_data$Difference, na.rm = TRUE)

    #Calculate fold change:
    plotting_data$Fold_change <- log2((plotting_data$Difference/threshold))
    
    #Format data for plotting:
    plotting_data <- plotting_data[,c(1,2,5,4)]
    
    #Extract significant positions:
    significant_positions <- subset(plotting_data, Fold_change>=FC_thr)
    
    return(list(plotting_data, significant_positions))
  }
}

subset_by_chr <- function(list_methods, chr){
 for (i in 1:length(list_methods)){
   list_methods[[i]] <- subset(list_methods[[i]], startsWith(as.character(Reference), chr))
 }
 return(list_methods)
}

barplot_plotting <- function (list_plotting, list_significant, output_name){
  
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
  plot(ggplot(initial_df, aes(x=Position, y=Fold_change)) + ggtitle(output_name) +
          geom_bar(stat = "identity", width=1, fill="#2a7886") +
          geom_bar(data=putative_positions, stat = "identity", width=1, fill="red") + 
          #geom_text_repel(data=putative_positions, aes(Position, Difference, label=Position,size=3), segment.size  = 0.4,segment.color = "grey50")+
          theme_bw() +theme(plot.title = element_text(face = "bold", hjust = 0.5), text = element_text(size=20),
                            axis.text = element_text(size = 20), strip.text.y = element_text(size = 20)) + 
          facet_grid(sample_f ~ . , scales="fixed") )
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

merging_data_per_kmer <- function (data) {
  
  #Merge data from the same kmer:
  data$Ref <- paste(data$Chr, data$Start, sep="_")
  final_data <- data.frame()
  for (i in 1:nrow(data)){
    if (i != nrow(data)) {
      if (data$Ref[i] != data$Ref[i+1]) {
        final_data <- rbind(final_data, data[i,])

      }
    } else {
      final_data <- rbind(final_data, data[i,])
    }
  }
  
  return(final_data[,c(1,2,3,4,5,6,7,8,9)])
}


analysis_significant_positions <- function (list_significant, fasta_file, output_name) {
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

  #Perform intersections:
  print('Perform intersections')
  n1 <- length(grEpinano[[1]])
  n2 <- length(grNanopolish[[1]])
  n3 <- length(grTombo[[1]])
  n4 <- length(grNanocompore[[1]])
  print(c(n1,n2,n3,n4))
  
  if (n1 != 0 & n2 != 0 & n3 != 0 & n4 == 0 ) {
    intersect_12 <- subsetByOverlaps(grNanopolish[[1]], grEpinano[[1]], maxgap = 4)
    intersect_13 <- subsetByOverlaps(grTombo[[1]], grEpinano[[1]], maxgap = 4)
    intersect_14 <- c()
    intersect_23 <- subsetByOverlaps(grNanopolish[[1]], grTombo[[1]], maxgap = 4)
    intersect_24 <- c()
    intersect_34 <- c()
    intersect_123 <- subsetByOverlaps(intersect_12, grTombo[[1]], maxgap = 4)
    intersect_124 <- c()
    intersect_134 <- c()
    intersect_234 <- c()
    intersect_1234 <- c()
    
  } else if (n1 != 0 & n2 == 0 & n3 != 0 & n4 != 0 ) {
    intersect_12 <- c()
    intersect_13 <- subsetByOverlaps(grEpinano[[1]], grTombo[[1]], maxgap = 4)
    intersect_14 <- subsetByOverlaps(grEpinano[[1]], grNanocompore[[1]], maxgap = 4)
    intersect_23 <- c()
    intersect_24 <- c()
    intersect_34 <- subsetByOverlaps(grNanocompore[[1]], grTombo[[1]], maxgap = 4)
    intersect_123 <- c()
    intersect_124 <- c()
    intersect_134 <- subsetByOverlaps(intersect_13, grNanocompore[[1]], maxgap = 4)
    intersect_234 <- c()
    intersect_1234 <- c()
    
  } else if (n1 == 0 & n2 != 0 & n3 != 0 & n4 != 0 ) {
    intersect_12 <- c()
    intersect_13 <- c()
    intersect_14 <- c()
    intersect_23 <- subsetByOverlaps(grNanopolish[[1]], grTombo[[1]], maxgap = 4)
    intersect_24 <- subsetByOverlaps(grNanopolish[[1]], grNanocompore[[1]], maxgap = 4)
    intersect_34 <- subsetByOverlaps(grTombo[[1]], grNanocompore[[1]], maxgap = 4)
    intersect_123 <- c()
    intersect_124 <- c()
    intersect_134 <- c()
    intersect_234 <- subsetByOverlaps(intersect_23, grNanocompore[[1]],  maxgap = 4)
    intersect_1234 <- c()
    
  } else if (n1 != 0 & n2 != 0 & n3 == 0 & n4 != 0 ) {
    intersect_12 <- subsetByOverlaps(grNanopolish[[1]], grEpinano[[1]], maxgap = 4)
    intersect_13 <- c()
    intersect_14 <- subsetByOverlaps(grEpinano[[1]], grNanocompore[[1]], maxgap = 4)
    intersect_23 <- c()
    intersect_24 <- subsetByOverlaps(grNanopolish[[1]], grNanocompore[[1]], maxgap = 4)
    intersect_34 <- c()
    intersect_123 <- c()
    intersect_124 <- subsetByOverlaps(intersect_12, grNanocompore[[1]], maxgap = 4)
    intersect_134 <- c()
    intersect_234 <- c()
    intersect_1234 <- c()
    
  } else if (n1 == 0 & n2 == 0) {
    intersect_12 <- c()
    intersect_13 <- c()
    intersect_14 <- c()
    intersect_23 <- c()
    intersect_24 <- c()
    intersect_34 <- subsetByOverlaps(grNanocompore[[1]], grTombo[[1]], maxgap = 4)
    intersect_123 <- c()
    intersect_124 <- c()
    intersect_134 <- c()
    intersect_234 <- c()
    intersect_1234 <- c()
    
  } else if (n3 == 0 & n4 == 0) {
    intersect_12 <- subsetByOverlaps(grNanopolish[[1]], grEpinano[[1]], maxgap = 4)
    intersect_13 <- c()
    intersect_14 <- c()
    intersect_23 <- c()
    intersect_24 <- c()
    intersect_34 <- c()
    intersect_123 <- c()
    intersect_124 <- c()
    intersect_134 <- c()
    intersect_234 <- c()
    intersect_1234 <- c()
    
  } else if (n2 == 0 & n3 == 0) { 
    intersect_12 <-  c()
    intersect_13 <-  c()
    intersect_14 <- subsetByOverlaps(grEpinano[[1]], grNanocompore[[1]], maxgap = 4)
    intersect_23 <-  c()
    intersect_24 <-  c()
    intersect_34 <-  c()
    intersect_123 <-  c()
    intersect_124 <-  c()
    intersect_134 <-  c()
    intersect_234 <-  c()
    intersect_1234 <-  c()
      
  } else {
    intersect_12 <- subsetByOverlaps(grEpinano[[1]], grNanopolish[[1]], maxgap = 4)
    intersect_13 <- subsetByOverlaps(grTombo[[1]], grEpinano[[1]], maxgap = 4)
    intersect_14 <- subsetByOverlaps(grEpinano[[1]], grNanocompore[[1]], maxgap = 4)
    intersect_23 <- subsetByOverlaps(grTombo[[1]], grNanopolish[[1]], maxgap = 4)
    intersect_24 <- subsetByOverlaps(grNanocompore[[1]], grNanopolish[[1]], maxgap = 4)
    intersect_34 <- subsetByOverlaps(grTombo[[1]], grNanocompore[[1]], maxgap = 4)
    intersect_123 <- subsetByOverlaps(intersect_12, grTombo[[1]], maxgap = 4)
    intersect_124 <- subsetByOverlaps(intersect_12, grNanocompore[[1]], maxgap = 4)
    intersect_134 <- subsetByOverlaps(intersect_13, grNanocompore[[1]], maxgap = 4)
    intersect_234 <- subsetByOverlaps(intersect_23, grNanocompore[[1]],  maxgap = 4)
    intersect_1234 <- subsetByOverlaps(intersect_134, intersect_234, maxgap = 4)
  
  }
    
  #Venn Diagram:  
  print(c(n1, n2, n3, n4, length(intersect_12), length(intersect_13), length(intersect_14), length(intersect_23), length(intersect_24),
          length(intersect_34), length(intersect_123), length(intersect_124), length(intersect_134), length(intersect_234), length(intersect_1234)))
  draw_venn_diagram(n1, n2, n3, n4, length(intersect_12), length(intersect_13), length(intersect_14), length(intersect_23), length(intersect_24),
                    length(intersect_34), length(intersect_123), length(intersect_124), length(intersect_134), length(intersect_234), length(intersect_1234), methods_name, output_name)
  
  #Output plain text final result: 
  supported_by_2 <- list(intersect_12, intersect_13, intersect_14, intersect_23, intersect_24, intersect_34)
  supported_by_3 <- list(intersect_123, intersect_124, intersect_134, intersect_234)
  
  #Parse data supported by two algorithms:
  initial <- TRUE
  for (i in 1:length(supported_by_2)){
     
    #Create vectors to store software data:
    epinano_algorithm <- c()
    nanopolish_algorithm <- c()
    tombo_algorithm <- c()
    nanocompore_algorithm <- c()
    
    #Parse data into a data frame:
    if (length(supported_by_2[[i]])>0) {
      positions_df <- data.frame(start(supported_by_2[[i]]), end(supported_by_2[[i]]))
      colnames(positions_df) <- c('Start', 'End')
      positions_df$Chr <- seqlevels(supported_by_2[[i]])
      positions_df <- positions_df[,c(3,1,2)]
      
      if (i==1){
        epinano_algorithm <- c(epinano_algorithm, 'YES')
        nanopolish_algorithm <- c(nanopolish_algorithm, 'YES')
        tombo_algorithm <- c(tombo_algorithm, 'NO')
        nanocompore_algorithm <- c(nanocompore_algorithm, 'NO')
      } else if (i==2) {
        epinano_algorithm <- c(epinano_algorithm, 'YES')
        nanopolish_algorithm <- c(nanopolish_algorithm, 'NO')
        tombo_algorithm <- c(tombo_algorithm, 'YES')
        nanocompore_algorithm <- c(nanocompore_algorithm, 'NO')
      } else if (i==3) {
        epinano_algorithm <- c(epinano_algorithm, 'YES')
        nanopolish_algorithm <- c(nanopolish_algorithm, 'NO')
        tombo_algorithm <- c(tombo_algorithm, 'NO')
        nanocompore_algorithm <- c(nanocompore_algorithm, 'YES')
      } else if (i==4) {
        epinano_algorithm <- c(epinano_algorithm, 'NO')
        nanopolish_algorithm <- c(nanopolish_algorithm, 'YES')
        tombo_algorithm <- c(tombo_algorithm, 'YES')
        nanocompore_algorithm <- c(nanocompore_algorithm, 'NO')
      } else if (i==5) {
        epinano_algorithm <- c(epinano_algorithm, 'NO')
        nanopolish_algorithm <- c(nanopolish_algorithm, 'YES')
        tombo_algorithm <- c(tombo_algorithm, 'NO')
        nanocompore_algorithm <- c(nanocompore_algorithm, 'YES')
      } else {
        epinano_algorithm <- c(epinano_algorithm, 'NO')
        nanopolish_algorithm <- c(nanopolish_algorithm, 'NO')
        tombo_algorithm <- c(tombo_algorithm, 'YES')
        nanocompore_algorithm <- c(nanocompore_algorithm, 'YES')
      }
      
      #Add data to the final dataframe:
      positions_df$Epinano <- epinano_algorithm
      positions_df$Nanopolish <- nanopolish_algorithm
      positions_df$Tombo <- tombo_algorithm
      positions_df$Nanocompore <- nanocompore_algorithm
      
      if (initial == TRUE){
        final_df2 <- positions_df
        initial <- FALSE
      } else {
        final_df2 <- rbind(final_df2, positions_df)
      }
    }
  }
  
  #Parse data supported by three algorithms:
  if (sum(c(n1,n2,n3,n4)==0)<2) {
    initial <- TRUE
    
    for (i in 1:length(supported_by_3)){
      #Create vectors to store software data:
      epinano_algorithm <- c()
      nanopolish_algorithm <- c()
      tombo_algorithm <- c()
      nanocompore_algorithm <- c()
      
      #Parse data into a data frame: 
      if (length(supported_by_3[[i]])>0) {
        positions_df <- data.frame(start(supported_by_3[[i]]), end(supported_by_3[[i]]))
        colnames(positions_df) <- c('Start', 'End')
        positions_df$Chr <- seqlevels(supported_by_3[[i]])
        positions_df <- positions_df[,c(3,1,2)]
        
        if (i==1){
          epinano_algorithm <- c(epinano_algorithm, 'YES')
          nanopolish_algorithm <- c(nanopolish_algorithm, 'YES')
          tombo_algorithm <- c(tombo_algorithm, 'YES')
          nanocompore_algorithm <- c(nanocompore_algorithm, 'NO')
        } else if (i==2) {
          epinano_algorithm <- c(epinano_algorithm, 'YES')
          nanopolish_algorithm <- c(nanopolish_algorithm, 'YES')
          tombo_algorithm <- c(tombo_algorithm, 'NO')
          nanocompore_algorithm <- c(nanocompore_algorithm, 'YES')
        } else if (i==3) {
          epinano_algorithm <- c(epinano_algorithm, 'YES')
          nanopolish_algorithm <- c(nanopolish_algorithm, 'NO')
          tombo_algorithm <- c(tombo_algorithm, 'YES')
          nanocompore_algorithm <- c(nanocompore_algorithm, 'YES')
        }  else {
          epinano_algorithm <- c(epinano_algorithm, 'NO')
          nanopolish_algorithm <- c(nanopolish_algorithm, 'YES')
          tombo_algorithm <- c(tombo_algorithm, 'YES')
          nanocompore_algorithm <- c(nanocompore_algorithm, 'YES')
        }
        
        #Add data to the final dataframe:
        positions_df$Epinano <- epinano_algorithm
        positions_df$Nanopolish <- nanopolish_algorithm
        positions_df$Tombo <- tombo_algorithm
        positions_df$Nanocompore <- nanocompore_algorithm
        
        if (initial == TRUE){
          final_df3 <- positions_df
          initial <- FALSE
        } else {
          final_df3 <- rbind(final_df3, positions_df)
        }
      }
    }
  } else {
    final_df3 <- data.frame()
  }
  
  #Parse data supported by four algorithms:
  if (length(intersect_1234)>0) {
    final_df4 <- data.frame(start(intersect_1234), end(intersect_1234))
    colnames(final_df4) <- c('Start', 'End')
    final_df4$Chr <- seqlevels(intersect_1234)
    positions_df <- positions_df[,c(3,1,2)]
    
    #Add data to the final dataframe:
    positions_df$Epinano <- 'YES'
    positions_df$Nanopolish <- 'YES'
    positions_df$Tombo <- 'YES'
    positions_df$Nanocompore <- 'YES'
    final_df4 <- positions_df
    
  } else {
    final_df4 <- data.frame()
  }
  
  #Merge and remove duplicates:
  all_ranges <- unique(rbind(final_df2, final_df3, final_df4))
  
  #Kmer analysis: 
  print('Kmer analysis')
  kmer_data <- extract_kmers(all_ranges, fasta_file)
  all_ranges$Kmer <- kmer_data[[1]]
  all_ranges$RRACH_motif <- kmer_data[[2]]
  all_ranges <- all_ranges[order(all_ranges$Start, decreasing = FALSE),]

  #Merging data per kmer: 
  final <- merging_data_per_kmer(all_ranges)
  write.table(final, file = paste(output_name,'FinalOverlap.txt', sep="_"), sep = '\t', row.names = FALSE)

}