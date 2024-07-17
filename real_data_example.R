#  Run dependent_sim on a real RNA-seq dataset
library(dependentsimr)
library(ggplot2)
library(ggpointdensity)
library(tidyr)
library(dplyr)

# Load the real data ---------------------------------
# from GEO: GSE77221
head(Weger18)
weger18_library_sizes <- apply(Weger18, 2, sum)

# Run dependent_sim on this data --------------------
N_SAMPLES <- 96
library_sizes <- (weger18_library_sizes + rep(0, N_SAMPLES)) / mean(weger18_library_sizes) # Recycle the existing library sizes
rs <- get_random_structure(list(counts=Weger18), method="pca", rank=2, type="DESeq2")
draws <- draw_from_multivariate_corr(rs, n_samples=N_SAMPLES, size_factors=library_sizes)$counts

rs_corpcor <- get_random_structure(list(counts=Weger18), method="corpcor", type="DESeq2")
draws_corpcor <- draw_from_multivariate_corr(rs_corpcor, n_samples=N_SAMPLES, size_factors=library_sizes)$counts


# Generate simulated data without any dependence
rs_indep <- remove_dependence(rs)
indep_draws <- draw_from_multivariate_corr(rs_indep, n_samples=N_SAMPLES, size_factors=library_sizes)$counts

# Scale to Counts Per Million ------------------------
cpm <- function(x) { # counts per million
  return(t(t(x) / apply(x, 2, sum) * 1e6))
}
scaled_read_data <- cpm(Weger18)
scaled_draws <- cpm(draws)
scaled_indep_draws <- cpm(indep_draws)

# Plot marginal distributions --------------------------
simulated_mean <- apply(scaled_draws[,1:ncol(Weger18)], 1, mean)
simulated_variance <- apply(scaled_draws[,1:ncol(Weger18)], 1, var)
real_mean <- apply(scaled_read_data, 1, mean)
real_variance <- apply(scaled_read_data, 1, var)
marg_dist <- data.frame(
  real_mean = real_mean,
  simulated_mean = simulated_mean,
  real_variance = real_variance,
  simulated_variance = simulated_variance
  )
ggplot(marg_dist|> filter(real_mean > 0.1))+
  geom_pointdensity(aes(x=log(real_mean+1),y=log(simulated_mean+1))) +
  geom_abline(slope=1, intercept=0)
ggplot(marg_dist|> filter(real_mean > 0.1))+
  geom_pointdensity(aes(x=log(real_variance+1),y=log(simulated_variance+1))) +
  geom_abline(slope=1, intercept=0)

# Plot gene-gene correlation ---------------------------
N_rows <- 3000
high_expr_rows <- which(apply(Weger18, 1, mean) > 300)
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
    subtitle = paste("In", N_rows, "randomly chosen genes"))

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
    subtitle = "Projections onto the PCA components of the real data"
  ) +
  scale_color_manual(values = list(real="blue", sim="#00BFC4", indep_sim="#F8766D"))

# Compare SVD with and without dependence ----------------
PCA <- function(data) {
  constant_rows <- apply(data, 1, function(x) all(x - mean(x) == 0))
  return(prcomp(t(data[!constant_rows,]), rank.=2, scale.=TRUE))
}
pca_sdevs <- tibble(
  type = "real",
  component = 1:12,
  sdev = PCA(scaled_read_data)$sdev,
)
print("Running PCAs")
for (i in 0:7) {
  pca_sdevs <- rbind(
    pca_sdevs,
    tibble(
      type = "sim dep",
      component = 1:12,
      sdev = PCA(scaled_draws[,(1+12*i):(12*(i+1))])$sdev,
    ),
    tibble(
      type = "sim indep",
      component = 1:12,
      sdev = PCA(scaled_indep_draws[,(1+12*i):(12*(i+1))])$sdev,
    ),
    deparse.level=1
  )
}
ggplot(data=pca_sdevs |> filter(component<6)) +
  facet_grid(cols=vars(component)) +
  geom_boxplot(aes(x=type, y=sdev)) +
  labs(title="Size of top PCA components")

# Run DESeq2 on it ---------------------------------------
# Do a 5v5 study. We expect no significant differences since all are drawn from the same samples
low_disp_rs <- rs
low_disp_rs$marginals$counts$dispersion <- rs$marginals$counts$dispersion/1
low_disp_draws <- draw_from_multivariate_corr(low_disp_rs, n_samples=30, size_factors=library_sizes[1:30])$counts
dds <- DESeq2::DESeqDataSetFromMatrix(draws[,1:12], tibble(group = rep(c("ctrl", "case"), 6)), ~group)
dds <- DESeq2::DESeq(dds)
res <- DESeq2::results(dds)
ggplot(tibble(pvalue=res$pvalue), aes(sample=pvalue)) + stat_qq(distribution=stats::qunif) + geom_abline(slope=1, intercept=0) + scale_x_log10() + scale_y_log10()
ggplot(tibble(
  deseq_dispersion = DESeq2::dispersions(dds),
  gene_dispersion = S4Vectors::mcols(dds)$dispGeneEst,
  true_dispersion = rs$marginals$counts$dispersion
)) + geom_pointdensity(aes(true_dispersion, gene_dispersion)) +
    scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope=1, intercept=0)


# Other comparisons of real to simulated data ------------
long_marg_dist <- marg_dist |>
  pivot_longer(cols=everything(), names_to=c("source", "var"), names_pattern="(.*)_(.*)", values_to="value")
ggplot(long_marg_dist |> filter(var == "variance")) +
  geom_freqpoly(aes(log10(value), color=source)) +
  labs(x="log10(variance)")
ggplot(long_marg_dist |> filter(var == "mean")) +
  geom_freqpoly(aes(log10(value), color=source)) +
  labs(x="log10(mean)")
ggplot(marg_dist, aes(x=log10(real_mean), y=log10(real_variance))) + geom_point()
ggplot(tibble(
  deseq_dispersion = DESeq2::dispersions(dds),
  est_dispersion = S4Vectors::mcols(dds)$dispGeneEst,
  simulated_dispersion = rs$marginals$counts$dispersion,
  simulated_mean = rs$marginals$counts$q,
  deseq_mean = res$baseMean
), aes(x=deseq_mean, y=deseq_dispersion)) +
  geom_point(aes(x=simulated_mean, y=simulated_dispersion), color="blue", size=3) +
  geom_point(aes(x=deseq_mean, y=est_dispersion), color="red", size=2) +
  geom_point(size=1) +
  scale_x_log10() + scale_y_log10()

