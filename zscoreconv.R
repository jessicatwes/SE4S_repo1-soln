##### Import libraries
library(ComplexHeatmap)
library(GenomicRanges)

##### Set file paths for original files
mdfile = "/Users/jessicawestfall/Documents/Dowell_lab/Experiment_data/Exp73_normalize_barcode/SRR1810071_MDS.tsv"

##### Read in org files and convert to GRanges 
# read the .tsv file of motif distribution
md = read.csv(mdfile, sep='\t')
md_skip = read.csv(mdfile, sep='\t', skip=659)
md_skip = md_skip[1:641,]
colnames(md_skip) <- c('ID','NON','TSS','COMBINE')
keep <- c('ID','COMBINE')
md_skip = md_skip[keep]
combine = data.frame(do.call("rbind", strsplit(as.character(md_skip$COMBINE), ",", fixed = TRUE)))
md <-cbind(md_skip$ID, combine)
# transform data so each column is a single TF
t_md <-t(md)
# select a single TF
graph1 = as.numeric(t_md[2:3000])
graph1 = graph1 + 0.00001
ComplexHeatmap::Heatmap(graph1, cluster_rows = FALSE)


##### z-score bin conversion
zscore_conv <- function(tf_col) {
  bin_num <- cut(tf_col, seq(-1, max(tf_col), 1)) # divide d into intervals defined by seq
  bin_num_df <- as.data.frame(cbind(bin_num, tf_col))
  bin_num_df$index <- 1:nrow(bin_num_df)
  #bin_md_df <- merge(graph1, bin_num_df, by.x="graph1", by.y="bin_num")
  
  n_occur_dists <- data.frame(table(bin_num_df$bin_num)) # Var1 is bin number
  n_occur_dists$Var1 <- as.integer(n_occur_dists$Var1) #bin as integer
  n_occur_dists <- n_occur_dists[order(n_occur_dists$Freq),] #ordered based on freq
  q1_bin = as.integer((nrow(n_occur_dists))/4) + 1
  q2_median_bin = as.integer(((nrow(n_occur_dists))/4)*2) + 1
  q3_bin = as.integer(((nrow(n_occur_dists))/4)*3) + 1
  # what is frequency of bins that fall in quantile
  q1 = n_occur_dists[q1_bin,]$Freq # frequency at q1
  q2_median = n_occur_dists[q2_median_bin,]$Freq
  q3 = n_occur_dists[q3_bin,]$Freq
  IQR = q3 - q1
  print(IQR)
  
  n_occur_dists$Freq_IQR <- ifelse((n_occur_dists$Freq > q3 | 
                                      n_occur_dists$Freq < q1), 
                                   NA, n_occur_dists$Freq)
  print(head(n_occur_dists))
  n_occur_dists_ordered_by_bin_position <-n_occur_dists[order(n_occur_dists$Var1),]
  
  n_occur_dists_zscore <- n_occur_dists %>% 
    mutate(zscore = (Freq-mean(Freq_IQR, na.rm=TRUE)) / sd(Freq_IQR, na.rm=TRUE))
  #upper2_zscore = n_occur_dists_zscore[q2_median_bin,]$zscore + ((n_occur_dists_zscore[q3_bin,]$zscore-n_occur_dists_zscore[q1_bin,]$zscore)*1.5)
  upper2_zscore = n_occur_dists_zscore[q2_median_bin,]$zscore + n_occur_dists_zscore[IQR,]$zscore*1.5
  
  lower2_zscore= n_occur_dists_zscore[q2_median_bin,]$zscore - n_occur_dists_zscore[IQR,]$zscore*1.5
  q1_zscore = n_occur_dists_zscore[q1_bin,]$zscore
  q2_zscore = n_occur_dists_zscore[q2_median_bin,]$zscore
  q3_zscore = n_occur_dists_zscore[q3_bin,]$zscore
  upper3_zscore = n_occur_dists_zscore[q2_median_bin,]$zscore + n_occur_dists_zscore[IQR*3,]$zscore
  lower3_zscore= n_occur_dists_zscore[q2_median_bin,]$zscore - n_occur_dists_zscore[IQR*3,]$zscore
  col2_fun = colorRamp2(c(lower3_zscore, lower2_zscore, q1_zscore, q2_zscore, q3_zscore, upper2_zscore, upper3_zscore), 
                        c("red", "#ffdaf0","#fff3fa", "white", "#c9F5C0","#A1F8A0" ,"green"))
  n_occur_dists_zscore_ordered_by_bin_position <-n_occur_dists_zscore[order(n_occur_dists_zscore$Var1),]
  
  ##### z-score link back to org table
  position_zscore <- merge(bin_num_df, n_occur_dists_zscore, by.x="bin_num", by.y="Var1")
  sort_position <- position_zscore[order(position_zscore$index),]
  print(head(sort_position))
  }

zscore_conv(graph1)

bin_num <- cut(graph1, seq(-1, max(graph1), 1)) # divide d into intervals defined by seq
bin_num_df <- as.data.frame(cbind(bin_num, graph1))
bin_num_df$index <- 1:nrow(bin_num_df)
#bin_md_df <- merge(graph1, bin_num_df, by.x="graph1", by.y="bin_num")

n_occur_dists <- data.frame(table(bin_num_df$bin_num)) # Var1 is bin number
n_occur_dists$Var1 <- as.integer(n_occur_dists$Var1) #bin as integer
n_occur_dists <- n_occur_dists[order(n_occur_dists$Freq),] #ordered based on freq
q1_bin = as.integer((nrow(n_occur_dists))/4) 
q2_median_bin = as.integer(((nrow(n_occur_dists))/4)*2) 
q3_bin = as.integer(((nrow(n_occur_dists))/4)*3) 
# what is frequency of bins that fall in quantile
q1 = n_occur_dists[q1_bin,]$Freq # frequency at q1
q2_median = n_occur_dists[q2_median_bin,]$Freq
q3 = n_occur_dists[q3_bin,]$Freq
IQR = q3 - q1

n_occur_dists$Freq_IQR <- ifelse((n_occur_dists$Freq > q3 | 
                                    n_occur_dists$Freq < q1), 
                                 NA, n_occur_dists$Freq)

n_occur_dists_ordered_by_bin_position <-n_occur_dists[order(n_occur_dists$Var1),]

n_occur_dists_zscore <- n_occur_dists %>% 
  mutate(zscore = (Freq-mean(Freq_IQR, na.rm=TRUE)) / sd(Freq_IQR, na.rm=TRUE))
#upper2_zscore = n_occur_dists_zscore[q2_median_bin,]$zscore + ((n_occur_dists_zscore[q3_bin,]$zscore-n_occur_dists_zscore[q1_bin,]$zscore)*1.5)
upper2_zscore = n_occur_dists_zscore[q2_median_bin,]$zscore + n_occur_dists_zscore[IQR,]$zscore*1.5

lower2_zscore= n_occur_dists_zscore[q2_median_bin,]$zscore - n_occur_dists_zscore[IQR,]$zscore*1.5
q1_zscore = n_occur_dists_zscore[q1_bin,]$zscore
q2_zscore = n_occur_dists_zscore[q2_median_bin,]$zscore
q3_zscore = n_occur_dists_zscore[q3_bin,]$zscore
upper3_zscore = n_occur_dists_zscore[q2_median_bin,]$zscore + n_occur_dists_zscore[IQR*3,]$zscore
lower3_zscore= n_occur_dists_zscore[q2_median_bin,]$zscore - n_occur_dists_zscore[IQR*3,]$zscore
col2_fun = colorRamp2(c(lower3_zscore, lower2_zscore, q1_zscore, q2_zscore, q3_zscore, upper2_zscore, upper3_zscore), 
                      c("red", "#ffdaf0","#fff3fa", "white", "#c9F5C0","#A1F8A0" ,"green"))
n_occur_dists_zscore_ordered_by_bin_position <-n_occur_dists_zscore[order(n_occur_dists_zscore$Var1),]

##### z-score link back to org table
position_zscore <- merge(bin_num_df, n_occur_dists_zscore, by.x="bin_num", by.y="Var1")
sort_position <- position_zscore[order(position_zscore$index),]
write.csv(sort_position,"/Users/jessicawestfall/Documents/Dowell_lab/Experiment_data/Exp73_normalize_barcode/SRR1810071_E4F1.csv", row.names = FALSE)

ComplexHeatmap::Heatmap(sort_position$graph1, col=col2_fun, cluster_rows = FALSE)
ComplexHeatmap::Heatmap(sort_position$zscore, col=col2_fun, cluster_rows = FALSE)


