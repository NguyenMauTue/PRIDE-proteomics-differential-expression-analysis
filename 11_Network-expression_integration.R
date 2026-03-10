############################################################
# 11 Network-expression integration
############################################################

library(biomaRt)
library(dplyr)

network_summary =
  read.csv("../results/network_summary.csv")

imputed_result =
  read.csv("../results/imputed_result.csv")


mart =
  useMart("ensembl",
          dataset="hsapiens_gene_ensembl")

mapping =
  getBM(
    attributes=c(
      "hgnc_symbol",
      "uniprotswissprot"
    ),
    filters="hgnc_symbol",
    values=network_summary$Symbol,
    mart=mart
  )

colnames(mapping)=
  c("SYMBOL","UNIPROT")

mapping =
  mapping[mapping$UNIPROT!="", ]

mapping =
  mapping[!duplicated(mapping$SYMBOL), ]


network_summary2 =
  merge(
    network_summary,
    mapping,
    by.x="Symbol",
    by.y="SYMBOL"
  )


limmanetwork_df =
  merge(
    network_summary2,
    imputed_result,
    by="UNIPROT"
  )


write.csv(
  limmanetwork_df,
  "../results/limma_network_table.csv",
  row.names=FALSE
)
