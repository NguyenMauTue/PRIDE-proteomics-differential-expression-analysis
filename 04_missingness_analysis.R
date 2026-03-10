############################################################
# 04 Missingness mechanism analysis
############################################################

library(ggplot2)

log_lfq_matrix = readRDS("../data/log_lfq_matrix.rds")
metadata = read.csv("../data/sample_metadata.csv")

############################################################
# Detection matrix
############################################################

detection_matrix =
  ifelse(is.na(log_lfq_matrix), 0, 1)

############################################################
# Missingness per protein
############################################################

protein_missing =
  rowSums(is.na(log_lfq_matrix))

protein_detected =
  rowSums(!is.na(log_lfq_matrix))

missing_fraction =
  protein_missing / ncol(log_lfq_matrix)

missing_summary =
  data.frame(
    Protein = rownames(log_lfq_matrix),
    Detected = protein_detected,
    Missing = protein_missing,
    Missing_fraction = missing_fraction
  )

write.csv(
  missing_summary,
  "../results/protein_missingness_summary.csv",
  row.names = FALSE
)

############################################################
# Detection probability vs intensity
############################################################

intensity_vector =
  as.vector(log_lfq_matrix)

detected_vector =
  ifelse(is.na(intensity_vector), 0, 1)

missing_df =
  data.frame(
    Intensity = intensity_vector,
    Detected = detected_vector
  )

missing_df =
  missing_df[!is.na(missing_df$Intensity), ]

############################################################
# Logistic regression model
############################################################

missing_model =
  glm(
    Detected ~ Intensity,
    data = missing_df,
    family = binomial
  )

capture.output(
  summary(missing_model),
  file = "../results/missingness_logistic_model.txt"
)

############################################################
# Detection probability plot
############################################################

missing_plot =
  ggplot(missing_df,
         aes(x = Intensity, y = Detected)) +
  geom_jitter(height = 0.05, alpha = 0.2) +
  geom_smooth(
    method = "glm",
    method.args = list(family = "binomial"),
    color = "red"
  ) +
  theme_minimal() +
  labs(
    title = "Detection probability vs intensity",
    x = "log2 LFQ intensity",
    y = "Detection probability"
  )

ggsave(
  "../results/detection_probability_curve.png",
  missing_plot
)

############################################################
# 3v0 pattern detection
############################################################

group_vector =
  metadata$Group

presence_matrix =
  !is.na(log_lfq_matrix)

tumor_presence =
  rowSums(presence_matrix[, group_vector == "Tumor"])

normal_presence =
  rowSums(presence_matrix[, group_vector == "Normal"])

three_vs_zero =
  data.frame(
    Protein = rownames(log_lfq_matrix),
    Tumor_presence = tumor_presence,
    Normal_presence = normal_presence
  )

three_vs_zero =
  three_vs_zero[
    (three_vs_zero$Tumor_presence > 0 &
     three_vs_zero$Normal_presence == 0) |
    (three_vs_zero$Tumor_presence == 0 &
     three_vs_zero$Normal_presence > 0),
  ]

write.csv(
  three_vs_zero,
  "../results/three_vs_zero_proteins.csv",
  row.names = FALSE
)

############################################################
# Summary statistics
############################################################

summary_stats =
  data.frame(
    Total_proteins = nrow(log_lfq_matrix),
    Proteins_with_missing =
      sum(protein_missing > 0),
    Proteins_3v0 =
      nrow(three_vs_zero)
  )

write.csv(
  summary_stats,
  "../results/missingness_summary_stats.csv",
  row.names = FALSE
)
