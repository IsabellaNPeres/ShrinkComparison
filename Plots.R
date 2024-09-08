library(tibble)
library(tidyr)
library(dplyr)
library(readxl)
library(ggplot2)

setwd('/Data/ClusterGeneration/Normal')
cov_matrix <- as.matrix(read_excel("/Data/ClusterGeneration/cov_matriz_clusterGeneration.xlsx"))
max_cols <- ncol(cov_matrix)

m_names <- c("stein_org_", "disutility_", "frechet_", "frobenius_org_", "inverse_frobenius_",
             "inverse_quadratic_", "inverse_stein_", "min_var_", "quadratic_", "rmse_var_",
             "rmse_cov_", "symmetrized_stein_", "weighted_frobenius_")
for (measure_name in m_names) {
  measure <- data.frame()
  for (n_cols in c(10, 50, 100)) {
    for (n_obs in c(550, 1000, 2000)) {
      measure_parcial <- read.csv(paste0(measure_name, n_obs, "_", n_cols, ".csv"))[,-1]
      measure_parcial$n_obs <- n_obs
      measure_parcial$n_cols <- n_cols
      measure<- rbind(measure, measure_parcial)
    }
  }
  measure_longer <- pivot_longer(measure, cols = amos:nld, names_to = "estimador")  # amos:spikedSt
  ggplot(measure_longer) + 
    geom_violin(aes(x = estimador, y = value, fill = factor(n_obs))) +
    facet_grid(n_cols ~ ., scales = "free_y") +
    theme_bw()  +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          legend.position = "bottom") +
    labs(fill = 'Num. Observations') + 
    guides(x =  guide_axis(angle = 90))
  ggsave(paste0("boxplot_", measure_name, ".pdf"), width = 297, height = 210, units = "mm")
}











resultados <- measure_longer  |> group_by(n_cols, n_obs, estimador) |> 
  summarise(media = mean(value, na.rm = TRUE))

