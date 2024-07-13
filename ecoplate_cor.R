one_swab<-read.delim('Bullfrog-nutrient-enrichment/one_swab_ecoplate_data.txt', header=T)
five_swab<-read.delim('Bullfrog-nutrient-enrichment/five_swab.txt', header=T)

one_split<-split(one_swab, one_swab$IndivID)
five_split<-split(five_swab, five_swab$Pond)

calculate_pairwise_spearman <- function(data_list, column_name) {
  # Number of data frames
  n <- length(data_list)
  
  # Initialize matrix to store correlations
  corr_matrix <- matrix(NA, n, n)
  rownames(corr_matrix) <- colnames(corr_matrix) <- paste0("DF", 1:n)
  
  # Calculate pairwise Spearman correlations
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      # Extract the column of interest from both data frames
      values_i <- data_list[[i]][[column_name]]
      values_j <- data_list[[j]][[column_name]]
      
      # Calculate Spearman correlation
      corr <- cor(values_i, values_j, method = "spearman", use = "complete.obs")
      
      # Store the correlation in the matrix
      corr_matrix[i, j] <- corr
      corr_matrix[j, i] <- corr
    }
  }
  
  return(corr_matrix)
}

correlation_matrix1 <- calculate_pairwise_spearman(one_split, "OD590")
correlation_matrix5 <- calculate_pairwise_spearman(five_split, "OD590")

library(reshape2)

cor_m1<-melt(correlation_matrix1)
cor_m5<-melt(correlation_matrix5)
cor_m5<-cor_m5[!is.na(cor_m5$value),]
cor_m1<-cor_m1[!is.na(cor_m1$value),]

library(ggplot2)
ggplot(cor_m1, aes(value))+
  geom_boxplot()+
  geom_boxplot(cor_m1, aes(value))+
  theme_bw()

ggplot(cor_m5, aes(value))+
  geom_boxplot(fill='red')+
  geom_boxplot(data=cor_m1, fill='grey')+
  theme_bw()+
  xlab('Spearman Rho')+
  ylab("")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())




