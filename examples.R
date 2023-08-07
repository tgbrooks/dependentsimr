## MULTIVARIATE NORMAL DATA ----------------------------------
# Generate the example MV normal data
library(MASS)
Sigma = matrix(c(
  1,      0.8,    0,  0,  0,  0,
  0.8,    1,      0,  0,  0,  0,
  0,      0,      1,  0,  0,  0,
  0,      0,      0,  1,  0,  0,
  0,      0,      0,  0,  1,  0.3,
  0,      0,      0,  0,  0.3, 1
), nrow=6, ncol=6)
norm_data <- t(mvrnorm(n=20, mu=c(0,0,0,0,0,0), Sigma=Sigma))

# Simulate draws mimicking that data
rs_normal <- get_random_structure(norm_data, rank=2, type="normal")
draws_normal <- draw_from_multivariate_corr(rs_normal, n_samples=30)


## POISSON DATA -------------------------------------------
# Generate dependent Poisson data from the MV normal data 
pois_data <- qpois(pnorm(norm_data), lambda = c(50, 1000, 2, 10, 0.2, 0.5))

# Simulate draws mimicking that data
rs_poisson <- get_random_structure(pois_data, rank=1, type="poisson")
draws_poisson <- draw_from_multivariate_corr(rs_poisson, n_samples=30)

## DESEQ2 DATA -------------------------------------------------
# A small amount of RNA-seq read count data
d <- c(
         3563,      3839,      3407,      3745,      4063,      3089,      3219,      4577,      3769,
         0,         0,         0,         0,         0,         0,         0,         0,         0,
         0,         0,         0,         0,         0,         0,         0,         0,         0,
         3626,      2630,      5470,      4230,      4739,      4499,      5862,      6524,      3621,
         317,       310,       254,       449,       431,       463,       292,       325,       299,
         0,         0,         0,         0,         0,         0,         0,         0,         0,
         0,         0,         0,         0,         0,         0,         0,         0,         0,
         0,         0,         0,         0,         0,         0,         0,         0,         0,
         0,         0,         0,         0,         0,         0,         0,         0,         0,
         3,         2,         0,         1,         5,         3,         9,         1,         2,
         2,         1,         0,         2,         0,         0,         0,         0,         0,
         0,         0,         0,         0,         0,         0,         0,         0,         0,
         1,         0,         0,         0,         1,         0,         0,         0,         0,
         1,         1,         2,         0,         0,         1,         0,         0,         1,
         0,         0,         0,         0,         0,         0,         0,         0,         0,
         984,      1014,       777,      1211,      1174,      1096,       973,      1193,      1083,
         0,         0,         0,         0,         0,         0,         0,         0,         0,
         16,        15,        16,         9,        23,        10,        25,        10,         8,
         4473,      3954,      4476,      3965,      4606,      3788,      4084,      4316,      3658,
         0,         0,         0,         0,         0,         0,         0,         1,         0
)

read_counts <- matrix(d, nrow=20, ncol=9, byrow=TRUE)

# Simulate draws mimicking that data
rs_deseq <- get_random_structure(read_counts, rank=2, type="DESeq2")
draws_deseq <- draw_from_multivariate_corr(rs_deseq, n_samples=30)
