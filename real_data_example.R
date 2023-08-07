#  Run dependent_sim on a real RNA-seq dataset
library(ggplot2)

# Load the real data ---------------------------------
liver <- read.delim("C:/Users/tgb/data/circadian_controls/results/Liver/num_reads_all_samples.txt")
# we choose just 8 samples to use from this large dataset
read_data <- as.matrix(liver[ ,3:10])
mode(read_data) <- "integer"
rownames(read_data) <- liver[ ,1]

# Run dependent_sim on this data
rs <- get_random_structure(read_data, rank=2, type="DESeq2")
draws <- draw_from_multivariate_corr(rs, n_samples=100)

rs_indep <- remove_dependence(rs)
indep_draws <- draw_from_multivariate_corr(rs_indep, n_samples=100)

# Plot marginal distributions --------------------------
simulated_mean <- apply(draws, 1, mean)
simulated_variance <- apply(draws, 1, var)
real_mean <- apply(read_data, 1, mean)
real_variance <- apply(read_data, 1, var)
plot(log(real_mean+1), log(simulated_mean+1))
plot(log(real_variance+1), log(simulated_variance+1))

# Plot correlation -------------------------------------
# Randomly sample pairs of genes to assess correlation
N_pairs = 10000
real_corr <- rep(NA, N_pairs)
simulated_corr <- rep(NA, N_pairs)
for (i in 1:N_pairs) {
  while(TRUE) {
    gene1 <- sample(nrow(read_data), 1)
    if (var(read_data[gene1, ]) != 0
          && var(draws[gene1, ]) != 0) {
      break
    }
  }
  while(TRUE) {
    gene2 <- sample(nrow(read_data), 1)
    if (var(read_data[gene2, ]) != 0
        && var(draws[gene2, ]) != 0
        && gene1 != gene2) {
          break
    }
  }
  real_corr[i] <- cor(read_data[gene1, ], read_data[gene2, ])
  simulated_corr[i] <- cor(draws[gene1, ], draws[gene2, ])
}
plot(real_corr, simulated_corr)

# Plot the top principal components --------------------
constant_rows <- apply(read_data, 1, function(x) all(x - mean(x) == 0))

# using PCA
#pca <- prcomp(t(read_data[!constant_rows,]), rank.=2, scale.=TRUE)
#projected_read_data <- t(pca$rotation) %*% scale(read_data[!constant_rows, ])
#projected_draws <- t(pca$rotation) %*% scale(draws[!constant_rows, ])
#projected_indep_draws <- t(pca$rotation) %*% scale(indep_draws[!constant_rows, ])
# using the existing projection from rs$cov
projected_read_data <- t(rs$cov$u[!constant_rows,]) %*% scale(read_data[!constant_rows, ])
projected_draws <- t(rs$cov$u[!constant_rows,]) %*% scale(draws[!constant_rows, ])
projected_indep_draws <- t(rs$cov$u[!constant_rows,]) %*% scale(indep_draws[!constant_rows, ])
both_data <- data.frame(list(
  x = c(projected_read_data[1,], projected_draws[1,], projected_indep_draws[1,]),
  y = c(projected_read_data[2,], projected_draws[2,], projected_indep_draws[2,]),
  type = c(
    rep("real", ncol(projected_read_data)),
    rep("sim", ncol(projected_draws)),
    rep("indep_sim", ncol(projected_indep_draws))
  )
))
ggplot(data = both_data, aes(x=x, y=y,color=type)) +
   geom_point()
