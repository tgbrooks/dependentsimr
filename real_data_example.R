#  Run dependent_sim on a real RNA-seq dataset
library(ggplot2)

# Load the real data ---------------------------------
# from GEO: GSE77221
head(Weger18)

# Run dependent_sim on this data --------------------
rs <- get_random_structure(list(counts=Weger18), rank=2, type="DESeq2")
draws <- draw_from_multivariate_corr(rs, n_samples=100)$counts

# Generate simulated data without any dependence
rs_indep <- remove_dependence(rs)
indep_draws <- draw_from_multivariate_corr(rs_indep, n_samples=100)$counts

# Scale to Counts Per Million ------------------------
cpm <- function(x) { # counts per million
  return(t(t(x) / apply(x, 2, sum) * 1e6))
}
scaled_read_data <- cpm(Weger18)
scaled_draws <- cpm(draws)
scaled_indep_draws <- cpm(indep_draws)

# Plot marginal distributions --------------------------
simulated_mean <- apply(scaled_draws, 1, mean)
simulated_variance <- apply(scaled_draws, 1, var)
real_mean <- apply(scaled_read_data, 1, mean)
real_variance <- apply(scaled_read_data, 1, var)
ggplot() +
  geom_point(aes(x=log(real_mean+1),y=log(simulated_mean+1))) +
  geom_abline(slope=1, intercept=0)
ggplot() +
  geom_point(aes(x=log(real_variance+1),y=log(simulated_variance+1))) +
  geom_abline(slope=1, intercept=0)

# Plot gene-gene correlation ---------------------------
N_rows <- 3000
high_expr_rows <- which(apply(Weger18, 1, mean) > 100)
selected_rows <- sample(high_expr_rows, N_rows)
selected_samples <- sample(ncol(draws), ncol(Weger18))
real_corr <- cor(t(scaled_read_data[selected_rows,]))
sim_corr <- cor(t(scaled_draws[selected_rows,selected_samples]))
indep_sim_corr <- cor(t(scaled_indep_draws[selected_rows,selected_samples]))
probs <- (1:1000)/1000
real_corr_q <- quantile(real_corr, probs = probs, na.rm=TRUE)
sim_corr_q <- quantile(sim_corr, probs = probs, na.rm=TRUE)
indep_sim_corr_q <- quantile(indep_sim_corr, probs = probs, na.rm=TRUE)
ggplot(data = rbind(data.frame(
    real_corr_q = real_corr_q,
    sim_corr_q = sim_corr_q,
    type="dependent"
  ),
  data.frame(
    real_corr_q = real_corr_q,
    sim_corr_q = indep_sim_corr_q,
    type = "independent"
  ))) +
  geom_point(aes(x = real_corr_q, y = sim_corr_q, color=type)) +
  geom_abline(slope=1, intercept=0) +
  labs(
    x = "Correlation coefficient (real data)",
    y = "Correlation coefficient (simulated data)",
    title = "QQ of Simulated Correlation",
    subtitle = "In 2000 randomly chosen genes, simulated with or without dependence.")

# Plot the top principal components --------------------
constant_rows <- apply(Weger18, 1, function(x) all(x - mean(x) == 0))
pca <- prcomp(t(scaled_read_data[!constant_rows,]), rank.=2, scale.=TRUE)
projected_read_data <- predict(pca, t(scaled_read_data[!constant_rows, ]))
projected_draws <- predict(pca, t(scaled_draws[!constant_rows, ]))
projected_indep_draws <- predict(pca, t(scaled_indep_draws[!constant_rows, ]))
both_data <- data.frame(list(
  PC1 = c(projected_read_data[,"PC1"], projected_draws[,"PC1"], projected_indep_draws[,"PC1"]),
  PC2 = c(projected_read_data[,"PC2"], projected_draws[,"PC2"], projected_indep_draws[,"PC2"]),
  type = c(
    rep("real", nrow(projected_read_data)),
    rep("sim", nrow(projected_draws)),
    rep("indep_sim", nrow(projected_indep_draws))
  )
))
both_data$type <- factor(both_data$type, levels=c("real", "sim", "indep_sim"))
ggplot(data = both_data, aes(x=PC1, y=PC2,color=type)) +
   geom_point() +
  labs(
    title = "Top 2 PCA Components of the Real Data",
    subtitle = "Real and simulated data (with and without independence) projected onto the PCA components of the real"
  ) +
  scale_color_manual(values = list(real="blue", sim="orange", indep_sim="red"))
