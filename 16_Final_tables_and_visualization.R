############################################################
# 16 Final tables and visualization
############################################################

library(ggplot2)
library(ggrepel)
library(openxlsx)

filtered_themed =
  read.csv("../results/BiomarkerCandidates_themed.csv")

#Filter candidate
filtered_themed =
  filtered_themed[
    filtered_themed$Num_Themes_Match >= 1,
  ]

#Summary table
summarized_df =
  filtered_themed[
    ,
    c(
      "UNIPROT",
      "Symbol",
      "logFC",
      "adj.P.Val",
      "ExprScore",
      "NetScore",
      "DriverScore",
      "Localization",
      "PubMed_hits"
    )
  ]

#Split by localization
listLocal = unique(summarized_df$Localization)

get_Subset_Local = function(df,loc){

  clean_name = gsub("[^A-Za-z0-9]","_",loc)

  assign(
    paste0(clean_name,"_df"),
    subset(df,df$Localization==loc),
    envir=.GlobalEnv
  )
}

for(loc in listLocal){
  get_Subset_Local(summarized_df,loc)
}

#Excel export
wb = createWorkbook()

addWorksheet(wb,"All")
writeData(wb,"All",summarized_df)

addWorksheet(wb,"Membrane")
writeData(wb,"Membrane",membrane_df)

addWorksheet(wb,"Extracellular")
writeData(wb,"Extracellular",extracellular_df)

addWorksheet(wb,"Intracellular")
writeData(wb,"Intracellular",intracellular_df)

addWorksheet(wb,"Vesicle_exosome")
writeData(wb,"Vesicle_exosome",vesicle_exosome_df)

saveWorkbook(
  wb,
  "../results/Localization_tables.xlsx",
  overwrite=TRUE
)

#Plot driver landscape
top10_df =
  summarized_df[
    order(-summarized_df$DriverScore),
  ][1:10,]

p1 =
  ggplot(summarized_df,
         aes(x=ExprScore,y=NetScore)) +
  geom_point(color="grey50") +
  geom_point(data=top10_df,
             color="red",
             size=3) +
  geom_text_repel(
    data=top10_df,
    aes(label=Symbol)
  ) +
  theme_classic()

ggsave(
  "../results/DriverScore_network_plot.png",
  p1,
  width=6,
  height=5
)

#Novelty analysis
filtered_themed$novel =
  filtered_themed$PubMed_hits <
  quantile(filtered_themed$PubMed_hits,0.25) &
  filtered_themed$DriverScore > 1

p2 =
  ggplot(filtered_themed,
         aes(x=PubMed_hits,
             y=DriverScore)) +
  geom_point(aes(color=novel),
             size=3,
             alpha=0.8) +
  scale_x_log10() +
  geom_hline(yintercept=1,
             linetype="dashed") +
  geom_vline(
    xintercept=
      quantile(filtered_themed$PubMed_hits,0.25),
    linetype="dashed"
  ) +
  scale_color_manual(values=c("grey70","red")) +
  theme_classic() +
  geom_text_repel(
    data=subset(filtered_themed,novel),
    aes(label=Symbol),
    size=3
  )

ggsave(
  "../results/Literature_novelty_plot.png",
  p2,
  width=6,
  height=5
)
