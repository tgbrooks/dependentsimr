# Generating Correlated Data for Omics Simulation

Omics data is in the "p >> n" regime where there are fewer samples than measurements per sample.
This creates dual challenges in generating realistic simulated data for the purposes of benchmarking.
First, there isn't enough data to be able to compute a dependence structure (e.g., a full-rank correlation matrix).
Second, generating omics-scale data with a specified correlation matrix is slow due to the typical $O(p^3)$ nature of these algorithms.
These often mean that simulators assume independence of the measurements, which does not reflect reality.

Here, we give a simple solution to both of these problems by using a low-rank correlation matrix to both approximate realistic dependencies in a real dataset and generate simulated data mimicking the real data.
Using a NORTA (Normal to Anything) approach, the marginal (univariate) distributions can have realistic forms like the negative binomial appropriate for omics datasets like RNA-seq read counts.
Our implementation supports normal, Poisson, and [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)-based (negative binomial with sample-specific size factors) marginal distributions.
This makes it particularly suited to RNA-seq data but also widely applicable.

Using this, simulating data that matches a given `data` (with samples in columns and measurements/features in rows) with each feature having a negative binomial distribution (fit by DESeq2) is fast and simple:

``` R
random_structure <- get_random_structure(data, rank = 2, type = "DESeq2")
simulated_data <- draw_from_multivariate_corr(random_structure, n_samples = 20)
```

## Covariance matrices
Suppose that we have a dataset with $n$ samples and $p$ measurements and the data stored in the $p x n$ matrix $X$.
In omics, we assume that $n < p$.
The $p \times p$ correlation (or covariance) matrix $\Sigma$ computed from the experimental dataset will have rank at most $n-1$.
If the measured $\Sigma$ is then used to generate data from a multivariate normal distribution, all the data will lie on an $n-1$-dimensional plane in the $p$-dimensional space of possible measurements.
This is clearly unrealistic for real data; if a new datapoint were collected experimentally we would expect it to not lie on that plane at all due to experimental noise as well as unmeasured biological variation.
Therefore, this full sample correlation matrix is inappropriate for simulating data and is arguably worse than including no dependence.

We suggest a simple alternative.
Pick a rank $k < n$ and then approximate $\Sigma$ by a rank $k$ matrix $\Sigma_k$.
This is equivalent to taking the $k$-dimensional principal component analysis (PCA) projection of the dataset and using the correlation matrix of that.
However, $\Sigma_k$ does not capture the full correlation of the data.
To remedy this, we add in the 'missing' variance from the remaining $n-k$ dimensions by assuming independence (i.e., as a diagonal correlation matrix).
For $k=0$, this method reduces to the basic simulation where complete independence is assumed.
For $k = n$, this instead reduces to using the entire sample correlation matrix, causing generated data to lie on a $n$-dimensional plane.
Choosing a $k$ between these extremes avoids both of these problems and can be done analogously to choosing a dimension for PCA.

## Marginal distributions
Omics data is rarely normally distributed.
The NORTA approach to generating data from other distributions with dependence is to first generate data with a multivariate normal distribution with a specified covariance matrix.
Then each variable is transformed to its desired marginal distribution via first the CDF of the normal distribution and then the inverse CDF of its marginal distribution.
This easily allows data to be simulated with arbitrary distributions, including discrete distributions such as Poisson or negative binomial.
However, the step of generating multivariate normal data from a covariance matrix is generally slow in existing packages as it requires a Cholesky decomposition of the covariance matrix, an operation that takes $O(p^3)$ time and is computationally expensive for omics-scale problems.

Since we already want to use a low-rank covariance matrix, we adapt the NORTA approach to directly generate low-rank multivariate normal data, thereby avoiding the expensive and redundant decomposition.
Then, the remaining variance is added back as (independent) normal data.
Finally, the data is then transformed, as in NORTA, to have the desired marginal distributions.

The remaining difficulty is that in NORTA the covariance matrix used to generate the multivariate normal data does not equal the covariance matrix of the generated end results with non-normal marginal distribution.
Various methods have been developed to determine an approximate matrix to give the desired end results.
Moreover, the standard NORTA approach does not preserve the low-rank structure of correlation matrix.

We avoid these difficulties by noting that we wish to match an existing dataset rather than being able to specify an arbitrary correlation matrix.
Therefore, we first normalize the data so that each feature has a standard normal marginal distribution, via the CDF of its MLE-fit marginal distribution.
This is the inverse of the operation we will do when generating our data from multivariate normal draws.
Therefore, we use the correlation matrix of this normalized dataset as the input to our low-rank approximation described above.

## Algorithm
Start with a $p \times n$ data matrix $X$ with samples in the columns and features (e.g., genes) in the rows and a chosen rank $0 < k < n$.
We will generate an output $p$-vector $X'$ of a single sample whose distribution approximates that of $X$. 

1. Fit marginal distributions to each feature in $X$ to determine CDFs $F_{i}$ for each feature.
2. Transform $X$ to normalized values by $Z_{ij} = \Phi^{-1}(F_{i}(X_{ij}))$ where $\Phi$ is the CDF of the standard normal distribution.
3. Take the PCA of $Z_{ij}$, i.e., the $k$ highest singular values $\lambda_1, \ldots, \lambda_k$ and their left-singular vectors $U = \left[u_1, \ldots, u_k\right]$. Note that $U D^2 U^T$ is the rank $k$ approximation to the correlation matrix $Z Z^T$ of $Z_{ij}$ where $D$ is the diagonal matrix with $\lambda_1, \ldots, \lambda_k$ on its diagonal. This is because the normalization of $Z$ means that $Z$ is already centered and scaled.
4. Compute the remaining variance for each variance by $M_{i} = \sum_{j=1}^n Z_{ij}^2/(n-1) - (UD^2U^T)_{ii}/(n-1)$.
5. Generate $k$ i.i.d. standard normally distributed values $w^T = \left[w_1, \ldots, w_k\right]$.
6. Generate $p \times n$ independent normally distributed values $V_{ij}$ with mean 0 and standard deviation $\sqrt{M_i}$.
7. Set $Z' = UDw/\sqrt{n-1} + V$.
8. Output $X'_i = F_i^{-1}(\Phi(Z'))$.

Each $Z'_{i}$ is normally distributed with the same standard deviation as $Z_{i \cdot}$, which was constructed to approximately be from the standard normal distribution.
Therefore, the output $X'_i$ has approximately the fit distribution with CDF $F_i$.
Moreover, using the transformation from step 2, the $\Phi^{-1} \circ F_i$ to each component of $X'$ gives $Z'$ which has a multivariate normal distribution with covariance matrix $U D^2 U^T/ + diag(\sqrt(M))$.
This agrees with the rank $k$ approximation of the empirical covariance matrix of $Z_{ij}$ as desired.
Note that we use $n-1$ rather than $n$ to obtain an unbiased estimate of variance.

## Implementation

We provide an R package `dependent_sim` that implements this method for normal, Poisson, and negative binomial marginal distributions, using the popular DESeq2 package to fit negative binomial datasets with varying size factors.
One advantage of this approach is its simplicity: the implementation is under 200 lines of R.
Performance is fast and generating 100 samples with 40,000 features (genes) each requires about 15 seconds with most of that time being spent on the DESeq2 fit.

One drawback is that extreme outliers are not well-captured.
Indeed, DESeq2 discards outliers and so the generated data does not include outliers.
Choosing to not discard outliers instead gives poor marginal distributional fits for those features.
Instead, we recommend adding outliers in as a separate post-processing step if desired.

Another drawback is that this method intentionally truncates the correlation matrix at a certain number of principal components.
Real data no doubt has additional dependencies that this does not capture; however, the typical sample sizes of omics studies make those smaller-scale dependencies impractical to measure.
Lastly, the dependencies are all captured linearly and therefore this is not a good match for situations such as data with clustered samples.
Instead, those clusters should each be simulated separately.

## Alternatives

The R package [SPsimSeq](https://github.com/CenterForStatistics-UGent/SPsimSeq) provides a dedicated RNA-seq and scRNA-seq simulator using a NORTA approach to simulate gene dependence.
In contrast to this package, it uses WGCNA to determine the correlation matrix, which is a gene network approach.
However, this method takes significant time to run.
Indeed, the `SPsimSeq` paper generated data for just 5000 genes based on a randomly sampled 5000 gene subset of the RNA-seq data and trying to generate a full sample exhausted the memory of a 24GB computer.
In contrast, our method runs in seconds on to generate full 40,000 gene (Mus musculus) samples on the same computer.
`SPsimSeq` is more specialized and full-featured for RNA-seq simulation, providing, for example, native differential expression (DE) options.
Using our package requires manually setting marginal expression values to simulate DE, but also supports other marginal distributions for situations outside of RNA-seq.

Other NORTA-based R packages that may be applicable, at least for datasets with smaller numbers of features, include `bindata`, `GenOrd`, and `SimMultiCorrData`, the last of these being the most comprehensive.
The `bigsimr` package provides faster implementations of these methods to scale up to omics-level data.
Though even this is computationally demanding; their paper references generating 20,000-dimensional vectors in "under an hour" using 16 threads.
The `copula` package provides even more flexible dependence options through use of copulas (the NORTA approach is equivalent to using Gaussian copulas).
All of these packages provide more flexibility in specifying dependence than our package, which can only mimic an existing dataset, and therefore the longer run-times may be unavoidable for use cases where researchers need to parameterize the dependence structure.

The `corpcor` package provides a James-Stein style shrinkage estimator for the correlation matrix for sample-limited cases like omics.
For inference, this is superior to our simple rank $k$ approximation, however for simulation, shrinkage of variance is not desirable as generated data does not reflect a realistic level of variation.
Moreover, `corpcor` yields full-rank correlation matrices that are not efficient to sample from.
