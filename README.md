
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dependentsimr

## Generating Correlated Data for Omics Simulation

Omics data is in the “p \>\> n” regime where there are fewer samples
than measurements per sample. This creates dual challenges in generating
realistic simulated data for the purposes of benchmarking. First, there
isn’t enough data to be able to compute a dependence structure (e.g., a
full-rank correlation matrix). Second, generating omics-scale data with
a specified correlation matrix is slow due to the typical $O(p^3)$
nature of these algorithms. These often mean that simulators assume
independence of the measurements, which does not reflect reality.

Here, we give a simple solution to both of these problems by using a
low-rank correlation matrix to both approximate realistic dependencies
in a real dataset and generate simulated data mimicking the real data.
Using a NORTA (Normal to Anything) approach, the marginal (univariate)
distributions can have realistic forms like the negative binomial
appropriate for omics datasets like RNA-seq read counts. Our
implementation supports normal, Poisson,
[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)-based
(negative binomial with sample-specific size factors), and empirical
(for ordinal data) marginal distributions. This makes it particularly
suited to RNA-seq data but also widely applicable.

Using this, simulating data that matches a given `data` (with samples in
columns and measurements/features in rows) with each feature having a
negative binomial distribution (fit by DESeq2) is fast and simple:

``` r
library(dependentsimr)
head(Weger18) # An RNA-seq dataset
#>                    GSM2046160 GSM2046184 GSM2046157 GSM2046183 GSM2046155
#> ENSMUSG00000087193          2          0          1          2          0
#> ENSMUSG00000039545         11          5          3         10          8
#> ENSMUSG00000085599          1          1          2          1          7
#> ENSMUSG00000084908          0          0          0          0          0
#> ENSMUSG00000086484          1          0          0          0          0
#> ENSMUSG00000085087          4          5          5          5          4
#>                    GSM2046159 GSM2046182 GSM2046180 GSM2046156 GSM2046158
#> ENSMUSG00000087193          0          0          2          0          0
#> ENSMUSG00000039545          6          3          2          2          1
#> ENSMUSG00000085599          2          0          2          0          0
#> ENSMUSG00000084908          1          0          0          0          0
#> ENSMUSG00000086484          0          0          0          0          0
#> ENSMUSG00000085087          5          1          5          7          0
#>                    GSM2046181 GSM2046179
#> ENSMUSG00000087193          0          0
#> ENSMUSG00000039545          5          0
#> ENSMUSG00000085599          1          0
#> ENSMUSG00000084908          0          0
#> ENSMUSG00000086484          0          0
#> ENSMUSG00000085087          2          1
```

``` r
rs <- get_random_structure(list(counts=Weger18), method="pca", rank=2, type="DESeq2")
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
#> -- replacing outliers and refitting for 319 genes
#> -- DESeq argument 'minReplicatesForReplace' = 7 
#> -- original counts are preserved in counts(dds)
#> estimating dispersions
#> fitting model and testing
```

``` r
simulated_data <- draw_from_multivariate_corr(rs, n_samples=20)$counts
head(simulated_data)
#>                    [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
#> ENSMUSG00000087193    1    0    0    0    0    1    0    0    1     1     2
#> ENSMUSG00000039545    2    1    6    9    0    2    3    2    1     2     2
#> ENSMUSG00000085599    2    0    1    1    0    0    2    0    0     0     0
#> ENSMUSG00000084908    0    0    0    0    0    0    0    0    0     0     0
#> ENSMUSG00000086484    0    0    0    0    0    0    0    0    0     0     0
#> ENSMUSG00000085087    3    3    5    2    2    4    0    1    5     2     4
#>                    [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20]
#> ENSMUSG00000087193     1     1     0     1     1     0     0     0     0
#> ENSMUSG00000039545     0     1     2     3     7     6     1     9     5
#> ENSMUSG00000085599     0     2     1     0     0     1     0     1     1
#> ENSMUSG00000084908     0     0     0     0     0     0     0     0     0
#> ENSMUSG00000086484     0     0     0     0     0     1     0     0     0
#> ENSMUSG00000085087     1     0     5     5     9     3     4     2     1
```

Finally, this also supports simultaneously generating multiple ‘modes’
of data, such as happens in multi-omics, where each node can have a
distinct marginal distribution type. For example, proteomics might have
normal margins and RNA-seq the DESeq2 margins. This captures cross-mode
dependencies observed in the data as well as intra-mode.

## Installation

You can install the development version of dependentsimr from GitHub
with:

``` r
# install.packages("remotes")
remotes::install_github("tgbrooks/dependentsimr")
```

## Covariance matrices

Suppose that we have a dataset with $n$ samples and $p$ measurements and
the data stored in the $p \times n$ matrix $X$. In omics, we assume that
$n < p$. The $p \times p$ correlation (or covariance) matrix $\Sigma$
computed from the experimental dataset will have rank at most $n-1$. If
the measured $\Sigma$ is then used to generate data from a multivariate
normal distribution, all the data will lie on an $n-1$-dimensional plane
in the $p$-dimensional space of possible measurements. This is clearly
unrealistic for real data; if a new datapoint were collected
experimentally we would expect it to not lie on that plane at all due to
experimental noise as well as unmeasured biological variation.
Therefore, this full sample correlation matrix is inappropriate for
simulating data and is arguably worse than including no dependence.

We propose two alternatives: a principal component analysis (PCA)
method, and a shrinkage method based off the `corpcor` R-package. For
the PCA method, pick a rank $k < n$ and then approximate $\Sigma$ by a
rank $k$ matrix $\Sigma_k$. This is equivalent to taking the
$k$-dimensional PCA projection of the dataset and using the correlation
matrix of that. However, $\Sigma_k$ does not capture the full
correlation of the data. To remedy this, we add in the ‘missing’
variance from the remaining $n-k$ dimensions by assuming independence
(i.e., as a diagonal correlation matrix). For $k=0$, this method reduces
to the basic simulation where complete independence is assumed. For
$k = n$, this instead reduces to using the entire sample correlation
matrix, causing generated data to lie on a $n$-dimensional plane.
Choosing a $k$ between these extremes avoids both of these problems and
can be done analogously to choosing a dimension for PCA. This decomposes
the covariance matrix as a sum of a low-rank part (the $k$-dimensional
PCA) and a diagonal (independent) part.

For the corpcor method, we use the `corpcor` to determine a shrinkage
from the sample correlation matrix to the identity correlation matrix.
This is done in `corpcor` by simply linearly interpolating between the
two by an amount determined by the data. Since the sample correlation
matrix is has rank bounded by the number of samples (more precisely, at
most $n-1$), the covariance matrix is a sum of a diagonal matrix and a
low-rank matrix, matching the form of the PCA method above.

Generating multivariate normal data from a large covariance matrix is
generally slow in existing packages as it requires a Cholesky
decomposition of the covariance matrix, an operation that takes $O(p^3)$
time and is computationally expensive for omics-scale problems. However,
our decomposition into low-rank and diagonal parts allows for extremely
efficient generation of random values, see Algorithms.

## Marginal distributions

Omics data is rarely normally distributed. The NORTA approach to
generating data from other distributions with dependence is to first
generate data with a multivariate normal distribution with a specified
covariance matrix. Then each variable is transformed to its desired
marginal distribution via first the CDF of the normal distribution and
then the inverse CDF of its marginal distribution. This easily allows
data to be simulated with arbitrary distributions, including discrete
distributions such as Poisson or negative binomial and arbitrary
empirical distributions.

Since we already want to use a low-rank covariance matrix, we adapt the
NORTA approach to directly generate low-rank multivariate normal data,
thereby avoiding the expensive and redundant decomposition. Then, the
remaining variance is added back as (independent) normal data. Finally,
the data is then transformed, as in NORTA, to have the desired marginal
distributions.

The remaining difficulty is that in NORTA the covariance matrix used to
generate the multivariate normal data does not equal the covariance
matrix of the generated end results with non-normal marginal
distribution. Various methods have been developed to determine an
approximate matrix to give the desired end results, though these depend
upon computationally demanding steps that do not preserve our low-rank
decomposition of the covariance matrix.

We avoid these difficulties by noting that we wish to match an existing
dataset rather than being able to specify an arbitrary correlation
matrix. We first normalize the data so that each feature has a standard
normal marginal distribution, via the CDF of its MLE-fit marginal
distribution. This is the inverse of the operation we will do when
generating our data from multivariate normal draws. Therefore, we use
the correlation matrix of this normalized dataset as the input to our
low-rank approximation described above.

## Algorithms

We first describe the PCA method. Start with a $p \times n$ data matrix
$X$ with samples in the columns and features (e.g., genes) in the rows
and a chosen rank $0 < k < n$. We will generate an output $p$-vector
$X'$ of a single sample whose distribution approximates that of $X$.

1.  Fit marginal distributions to each feature in $X$ to determine CDFs
    $F_{i}$ for each feature.
2.  Transform $X$ to normalized values by
    $Z_{ij} = \Phi^{-1}(F_{i}(X_{ij}))$ where $\Phi$ is the CDF of the
    standard normal distribution.
3.  Take the PCA of $Z$, i.e., the $k$ highest singular values
    $\lambda_1, \ldots, \lambda_k$ and their left-singular vectors
    $U = \left[u_1, \ldots, u_k\right]$ of $Z$.
4.  Compute the size factors $w$ by solving $A w = B$ where
    $A_{ij} = \delta_{ij} - \sum_\ell U_{\ell,i}^2 U_{\ell,j}^2$ and
    $B_{i} = \lambda_i^2/(n-1)^2 - \sum_\ell U_{\ell,i}^2 V_{\ell\ell}$,
    where $\delta_{ij}$ is the Kronecker delta and $V = Z^T Z / (n-1)$
    is the covariance matrix of $Z$. Set $W$ to be the diagonal matrix
    with $w$ along its diagonal.
5.  Set $D$ to be the diagonal matrix with
    $D_{ii} = \sqrt{V_{ii} - (UW^2U^T)_{ii}}$, which is the remaining
    variance.
6.  Generate $k$ i.i.d. standard normally distributed values $u$ and $p$
    i.i.d standard normally distributed values in $v$.
7.  Set $Z' = UWu + D v$.
8.  Output the vector $X'$ where $X'_i = F_i^{-1}(\Phi(Z'))$.

Each $Z_i'$ is normally distributed with the same standard deviation as
$Z_{i \cdot}$, which was constructed to approximately be from the
standard normal distribution. Therefore, the output $X_i'$ has
approximately the fit distribution with CDF $F_i$. Moreover, using the
transformation from step 2, the $\Phi^{-1} \circ F_i$ to each component
of $X'$ gives $Z'$ which has a multivariate normal distribution with
covariance matrix $\Sigma = U W^2 U^T + D^2$. The size factors $W$ were
chosen so that $\Sigma$ matches the rank $k$ approximation of the
empirical covariance matrix $V$ of $Z_{ij}$ while also having the same
diagonal elements as $V$. That is,
$U_{\cdot i}^T \Sigma U_{\cdot i} = U_{\cdot i}^T V U_{\cdot i}$ and
also $\Sigma_{ii} = V_{ii}$.

The corpcor method modifes this steps 3-7 with steps 1, 2, and 8
unchanged:

3.  Compute the $\lambda_1$ and $\lambda_2$ values from
    `corpcor::estimate.lambda` and `corpcor::estimate.lambda.var`
    functions on $Z$, respectively.
4.  Set $D$ to be diagonal with
    $D_{ii} = \sqrt{\lambda_1} (\lambda_2 \sigma_{med} + (1 - \lambda_2) \sigma_i)$
    where $\sigma_i$ is the standard deviation of the $Z_{i\cdot}$ and
    $\sigma_{med}$ is the median of the $\sigma_i$.
5.  Set $U$ to be $\sqrt{1-\lambda} S Z / \sqrt{n-1}$ where $S$ is the
    diagonal matrix with
    $S_{ii} = \lambda_2 \sigma_{med} / \sigma_i + (1 - \lambda_2)$.
6.  Generate $n$ i.i.d. standard normally distributed values $u$ and $p$
    i.i.d. standard normally distributed values $v$.
7.  Set $Z' = U u + D v$.

## R package

We provide an R package `dependent_sim` that implements this method for
normal, Poisson, and negative binomial marginal distributions, using the
popular DESeq2 package to fit negative binomial datasets with varying
size factors. One advantage of this approach is its simplicity: the
implementation is under 300 lines of R. Performance is fast and
generating 100 samples with 40,000 features (genes) each requires about
15 seconds with most of that time being spent on the DESeq2 fit.

One drawback is that extreme outliers are not well-captured. Indeed,
DESeq2 discards outliers and so the generated data does not include
outliers. Choosing to not discard outliers instead gives poor marginal
distributional fits for those features. Instead, we recommend adding
outliers in as a separate post-processing step if desired.

Another drawback is that this method intentionally truncates the
correlation matrix at a certain number of principal components. Real
data no doubt has additional dependencies that this does not capture;
however, the typical sample sizes of omics studies make those
smaller-scale dependencies impractical to measure. Lastly, the
dependencies are all captured linearly and therefore this is not a good
match for situations such as data with clustered samples. Instead, those
clusters should each be simulated separately.

## Alternatives

The R package
[SPsimSeq](https://github.com/CenterForStatistics-UGent/SPsimSeq)
provides a dedicated RNA-seq and scRNA-seq simulator using a NORTA
approach to simulate gene dependence. In contrast to this package, it
uses WGCNA to determine the correlation matrix, which is a gene network
approach. However, this method takes significant computational
resources. Indeed, the `SPsimSeq` paper generated data for just 5000
genes based on a randomly sampled 5000 gene subset of the RNA-seq data
and our attempts to use `SPsimSeq` to generate a full sample exhausted
the memory of a 24GB computer. In contrast, our method runs in seconds
to generate full 40,000 gene samples on the same computer. `SPsimSeq` is
more specialized and full-featured for RNA-seq simulation, providing,
for example, native differential expression (DE) options. In contrast,
our `dependentsimr` package requires manually setting marginal
expression values to inject DE, but also supports other marginal
distributions for situations outside of RNA-seq.

The
[scDesign2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02367-2)
simulator for scRNA-seq uses NORTA and, like our method, uses the
approach of estimating the correlation matrix from the normalized
dataset. However, it limits this correlation matrix to top-expressed
genes.

Other NORTA-based R packages that may be applicable, at least for
datasets with smaller numbers of features, include `bindata`, `GenOrd`,
and `SimMultiCorrData`, the last of these being the most comprehensive.
The `bigsimr` package provides faster implementations of these methods
to scale up to omics-level data. Though even this is computationally
demanding; their paper references generating 20,000-dimensional vectors
in “under an hour” using 16 threads. The `copula` package provides even
more flexible dependence options through use of copulas (the NORTA
approach is equivalent to using Gaussian copulas). All of these packages
provide more flexibility in specifying dependence than our package,
which can only mimic existing datasets, and therefore the longer
run-times may be unavoidable for use cases where researchers need to
parameterize the dependence structure.
