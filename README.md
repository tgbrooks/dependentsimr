
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
Poisson distribution is fast and simple:

``` r
library(dependentsimr)
head(read_counts) # An RNA-seq dataset
#>              gene_id Arpp19_10_3 Arpp19_4_1 Arpp19_5_1 Arpp19_5_3 Arpp19_5_4
#> 1 ENSMUSG00000033653        2377       2162       2169       2650       2217
#> 2 ENSMUSG00000025192         908        799        730        941        761
#> 3 ENSMUSG00000069011         316        342        373        373        325
#> 4 ENSMUSG00000066892         584        650        466        687        551
#> 5 ENSMUSG00000018736        2739       2531       2594       2988       2614
#> 6 ENSMUSG00000029624         807        753        746        894        806
#>   Arpp19_6_1 Arpp19_7_2 Arpp19_9_1 Dach1_11_1 Dach1_2_2 Dach1_3_1 Dach1_5_3
#> 1       2552       1947       2249       2414      2253      1933      1978
#> 2        780        643        778        812       691       627       640
#> 3        405        300        369        364       316       320       264
#> 4        606        532        604        573       551       485       513
#> 5       3081       2584       2656       2961      2643      2298      2421
#> 6        819        569        783        763       687       688       597
```

``` r
rs <- get_random_structure(list(counts=as.matrix(read_counts[,-1])), method="pca", rank=2, type="DESeq2")
#> converting counts to integer mode
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
```

``` r
simulated_data <- draw_from_multivariate_corr(rs, n_samples=5)$counts
head(simulated_data)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,] 2242 2277 2254 2209 2385
#> [2,]  703  741  772  827  780
#> [3,]  328  325  318  365  317
#> [4,]  558  529  505  606  587
#> [5,] 2660 2706 2783 2665 2728
#> [6,]  751  736  714  753  810
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
remotes::install_github("tgbrooks/dependent_sim")
```

## Vignette

For an extended example of using this package to generate RNA-seq data,
please see [this
vignette](https://html-preview.github.io/?url=https://github.com/tgbrooks/dependentsimr/blob/main/doc/simulate_data.html).
