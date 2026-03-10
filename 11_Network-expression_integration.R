############################################################
# 11 Network-expression integration
############################################################

library(biomaRt)

mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

mapping <- getBM(
  attributes = c("hgnc_symbol", "uniprotswissprot"),
  filters = "hgnc_symbol",
  values = network_summary$Symbol,
  mart = mart
)

colnames(mapping) <- c("SYMBOL","UNIPROT")

mapping <- mapping[mapping$UNIPROT != "", ]
mapping <- mapping[!duplicated(mapping$SYMBOL), ]

network_summary2 <- merge(
  network_summary,
  mapping,
  by.x="Symbol",
  by.y="SYMBOL"
)


############################################################
# Merge with limma results
############################################################

limmanetwork_df = merge(
  network_summary2,
  imputed_result,
  by = "UNIPROT"
)

library(dplyr)

limmanetwork_df = limmanetwork_df %>%
  relocate(Symbol, UNIPROT)

write.csv(
  limmanetwork_df,
  "../results/limma_network_table.csv"
)
