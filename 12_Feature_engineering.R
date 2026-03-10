############################################################
# 12 Feature engineering
############################################################

exnetscore_df =
  limmanetwork_df[
    , c("UNIPROT","logFC","degree","betweenness")
  ]

exnetscore_df$ExprScore =
  as.vector(abs(scale(exnetscore_df$logFC)))

exnetscore_df$NetScore =
  as.vector(
    scale(log(exnetscore_df$degree + 1)) +
    scale(log(exnetscore_df$betweenness + 1))
  )
