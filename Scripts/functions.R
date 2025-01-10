RLE_plot <- function(exprs_data, title) {
  row_medians_assayData <- rowMedians(as.matrix(exprs_data))
  RLE_data <- sweep(exprs_data, 1, row_medians_assayData)
  RLE_data <- as.data.frame(RLE_data)
  RLE_data_gathered <- tidyr::gather(RLE_data, samples, log2_expression_deviation)
  
  ggplot2::ggplot(RLE_data_gathered, aes(samples, log2_expression_deviation)) + 
    geom_boxplot(outlier.shape = NA) + ylim(c(-2, 2)) + 
    ggtitle(title) +
    theme(axis.text.x = element_text(colour = "aquamarine4", angle = 60, size = 6.5, hjust = 1, face = "bold"))
}


PCA_plot <- function(expr_data, phenotype_data, sample_names, title = "PCA Plot") {
  PCA <- prcomp(t(expr_data), scale = FALSE)
  percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  sd_ratio <- sqrt(percentVar[2] / percentVar[1])
  
  dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], Phenotype = phenotype_data, Sample = sample_names)
  ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(color = Phenotype)) +
    ggtitle(title) +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    coord_fixed(ratio = sd_ratio) + 
    theme_bw() +
    geom_text(aes(label = Sample), vjust = -0.5, size = 3, color = "black")
}
