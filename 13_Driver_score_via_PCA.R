############################################################
# 13 Driver score via PCA
############################################################

limmanetwork_df =
  read.csv("../results/limma_network_table.csv")


pca_testing =
  data.frame(
    logFC = limmanetwork_df$logFC,
    nlog10FDR = -log10(limmanetwork_df$adj.P.Val),
    degree = log(limmanetwork_df$degree+1),
    betweenness = log(limmanetwork_df$betweenness+1),
    row.names = limmanetwork_df$Symbol
  )


pca_testing_data =
  prcomp(
    pca_testing,
    scale=TRUE
  )


pca_testing_x =
  as.data.frame(
    pca_testing_data$x
  )


var =
  pca_testing_data$sdev^2 /
  sum(pca_testing_data$sdev^2)


alpha =
  var[2]/var[1]


DriverScore =
  pca_testing_x$PC1 +
  (alpha * pca_testing_x$PC2)


Driver_table =
  data.frame(
    UNIPROT = limmanetwork_df$UNIPROT,
    Symbol = limmanetwork_df$Symbol,
    DriverScore = DriverScore,
    ExprScore = pca_testing_x$PC1,
    NetScore = pca_testing_x$PC2
  )


write.csv(
  Driver_table,
  "../results/driver_scores.csv",
  row.names=FALSE
)
