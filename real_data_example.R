#  Run dependent_sim on a real RNA-seq dataset
library(ggplot2)

# Load the real data ---------------------------------
liver <- read.delim("C:/Users/tgb/data/circadian_controls/results/Liver/num_reads_all_samples.txt")
# we choose just 12 samples to use from this large dataset
read_data <- as.matrix(liver[ ,3:14])
mode(read_data) <- "integer"
rownames(read_data) <- liver[ ,1]

# Run dependent_sim on this data --------------------
rs <- get_random_structure(list(counts=read_data), rank=2, type="DESeq2")
draws <- draw_from_multivariate_corr(rs, n_samples=100)$counts

rs_indep <- remove_dependence(rs)
indep_draws <- draw_from_multivariate_corr(rs_indep, n_samples=100)$counts

# Scale to Counts Per Million ------------------------
cpm <- function(x) { # counts per million
  return(t(t(x) / apply(x, 2, sum) * 1e6))
}
scaled_read_data <- cpm(read_data)
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
N_rows <- 1000
high_expr_rows <- which(apply(read_data, 1, mean) > 100)
sample_rows <- sample(high_expr_rows, N_rows)
real_corr <- cor(t(scaled_read_data[sample_rows,]))
sim_corr <- cor(t(scaled_draws[sample_rows,1:ncol(read_data)]))
indep_sim_corr <- cor(t(scaled_indep_draws[sample_rows,1:ncol(read_data)]))
probs <- (1:1000)/1000
real_corr_q <- quantile(real_corr, probs = probs, na.rm=TRUE)
sim_corr_q <- quantile(sim_corr, probs = probs, na.rm=TRUE)
indep_sim_corr_q <- quantile(indep_sim_corr, probs = probs, na.rm=TRUE)
ggplot(data = data.frame(
    real_corr_q = real_corr_q,
    sim_corr_q = sim_corr_q,
    indep_sim_corr_q = indep_sim_corr_q
  )) +
  geom_point(aes(x = real_corr_q, y = sim_corr_q), color="red") +
  geom_point(aes(x = real_corr_q, y = indep_sim_corr_q), color="blue") +
  geom_abline(slope=1, intercept=0)

# Plot the top principal components --------------------
constant_rows <- apply(read_data, 1, function(x) all(x - mean(x) == 0))
pca <- prcomp(t(scale(read_data[!constant_rows,])), rank.=2, scale.=TRUE)
projected_read_data <- predict(pca, t(scale(read_data[!constant_rows, ])))
projected_draws <- predict(pca, t(scale(draws[!constant_rows, ])))
projected_indep_draws <- predict(pca, t(scale(indep_draws[!constant_rows, ])))
both_data <- data.frame(list(
  PC1 = c(projected_read_data[,"PC1"], projected_draws[,"PC1"], projected_indep_draws[,"PC1"]),
  PC2 = c(projected_read_data[,"PC2"], projected_draws[,"PC2"], projected_indep_draws[,"PC2"]),
  type = c(
    rep("real", nrow(projected_read_data)),
    rep("sim", nrow(projected_draws)),
    rep("indep_sim", nrow(projected_indep_draws))
  )
))
ggplot(data = both_data, aes(x=PC1, y=PC2,color=type)) +
   geom_point()
