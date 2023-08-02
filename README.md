# Generating Correlated Data for Omics Simulation

Omics data is in the "p >> n" condition where there are fewer samples than measurements per sample.
This creates dual challenges in generating realistic simulated data for the purposes of benchmarking.
First, there isn't enough data to be able to compute a dependence structure (e.g., a full-rank correlation matrix).
Second, generating omics-scale data with a specified correlation matrix is slow due to the typical $O(p^3)$ nature of these datatypes.
These often mean that simulators assume independence of the measurements, which may not be appropriate.

Here we give simple solutions to both of these problems by using a low-rank correlation matrix to both approximate realistic dependencies in the data and efficiently generate simulated data mimicking a real dataset.
Using a NORTA (Normal to Anything) approach, the marginal (univariate) distributions can have realistic forms like negative binomial appropriate for omics datasets like RNA-seq read counts.

Using this, simulating data that matches a given `data` (with samples in columns and measurements/features in rows) with each feature having a negative binomial distribution is as easy as:

``` R
random_structure <- get_random_structure(data, rank = 3, type = "negative binomial")
simulated_data <- draw_from_multivariate_corr(random_structure, n_samples = 20)
```

## Covariance matrices
Suppose that we have a dataset with `n` samples and `p` measurements and the data stored in the `p x n` matrix `X`.
In omics, we assume that `n < p`.
Then the $p \times p$ correlation (or covariance) matrices $\Sigma$ computed from the experimental dataset will have rank at most $n$ and so is low rank.
If the measured $\Sigma$ is then used to generate data from a multivariate normal distribution, all the data will lie on an $n$-dimensional plane in the $p$-dimensional space of possible measurements.
This is clearly unrealistic for real data; if a new datapoint were collected experimentally we would expect it to not lie on that plane at all by experimental noise alone.
Therefore, this is inappropriate for simulating data for the purposes of benchmarking.

Instead, we implement a simple alternative.
Pick a rank $k < n$ and then approximate $\Sigma$ by a rank $k$ matrix $\Sigma_k$.
This is equivalent to taking the $k$-dimensional principal component analysis (PCA) projection of the dataset and using the correlation matrix of that.
However, $\Sigma_k$ does not capture the full correlation of the data.
To remedy this, we add in the 'missing' variance from the remainign $n-k$ dimensions by assuming independence.
This means that for $k=0$, this method reduces to the basic simulation where complete independence is assumed.
For $k = n$, the method reduces to using the entire empirical correlation matrix, causing data to lie on a $n$-dimensional plane.
Choosing a number in between these extremes avoids both of these problems and can be done analogously to choosing a dimension for PCA.

## Marginal distributions
Omics data is rarely normally distributed.
The NORTA approach to generating data with dependence is to first generate data with a multivariate normal distribution of a specified covariance matrix.
Then each variable is transformed to its desired marginal distribution via first the CDF of the normal distribution and then the inverse CDF of its marginal distribution.
This easily allows data to be simulated with arbitrary distributions, including discrete distributions such as Poisson or negative binomial.
However, the step of generating multivariate normal data from a covariance matrix is generally slow in existing packages as it requires a Cholesky decomposition of the covariance matrix, an operation that takes $O(p^3)$ time and is computationally expensive for omics-scale problems.

To work with our low-rank covariance, we modify the above approach to directly generate low-rank multivariate normal data, thereby avoiding the expensive and redundant decomposition.
Ontop of this, the remaining variance is added back as (independent) normal data.
Finally, the data is then transformed, as in NORTA, to have the desired marginal distributions.

The difficulty in both NORTA and this setup is obtaining the appropriate covariance matrix for generating multivariate normal data.
The objective is to match the provided data, but the covariance matrix of the transformed data will not equal that of the multivariate normal data used to generate it.
As the standard NORTA approach does not preserve the low-rank structure of correlation matrix, we take an alternative approach.
First, we transform the data to have each feature have a standard normal marginal distribution, via the CDF of its MLE-fit marginal distribution.
Then the correlation matrix of this normalized dataset is used; we take it's rank $k$ approximation instead of the correlation matrix of the original data.

## Literature Notes
- https://www.cse.unr.edu/~fredh/papers/conf/196-agpuantaafhdms/paper.pdf
  - GPU accelerated version of NORTA
- https://citeseerx.ist.psu.edu/doc/10.1.1.48.281
  - Original NORTA paper
- DOI: 10.1002/asmb.901
  - For poisson variables specifically, gives efficient approximation to find the right input correlation matrix
