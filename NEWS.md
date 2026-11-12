# dependentsimr 2.0.0

## Breaking changes

* The simulation is now three steps, replacing `get_random_structure()` and
  `draw_from_multivariate_corr()`:

  ```r
  marginals  <- fit_marginals(data, method = "DESeq2")
  dependence <- fit_dependence(marginals, method = "pca", rank = 2)
  sim        <- simulate_data(marginals, dependence, n_samples = 10)
  ```

  Fitting the marginals (particularly with DESeq2) is the slow step, so it can now
  be reused while trying different dependence methods.
* `get_random_structure()` and `draw_from_multivariate_corr()` still work but are
  deprecated and give a warning. Random draws differ from version 1 even with the
  same seed.

## New features

* `fit_marginals()` accepts a `design` matrix of known sample covariates (e.g. groups,
  batches, continuous covariates). Each feature's marginal distribution is fit with
  that design (`lm` for normal, `glm` for Poisson, DESeq2's GLM for DESeq2) and the
  dependence is estimated from the residuals. `simulate_data()` then takes a
  `new_design` for the simulated samples. The `"empirical"` method does not support
  a design.
* `add_effect()` adds (or overwrites) a coefficient, e.g. log2 fold changes for
  simulating differential expression, which `simulate_data()` applies via the
  matching column of `new_design`.
* `fit_dependence(method = "independent")` for simulating with no dependence
  between features. `remove_dependence()` also works on the new objects.
* `fit_marginals()` accepts a single matrix or data frame (whose first column may
  hold feature names) instead of requiring a list, and `simulate_data()` then
  returns a single matrix.

## Bug fixes and improvements

* `"poisson"` marginals now use the numerically stable randomized quantile residuals,
  fixing infinite values for extreme observations (the stable version was computed
  but not used).
* Simulated values no longer become infinite for extreme latent normal draws: the
  quantile functions now work in log space on the nearer tail.
* `"normal"` marginals no longer give `NaN` for constant features.
* `"spiked Wishart"` draws now have the fitted covariance: variance `spiked_sd^2`
  along the spiked directions and `population_sd^2` in the space orthogonal to them.
  Previously the independent part used `sqrt(population_sd)` as its standard
  deviation instead of `population_sd` (inflating the independent variance whenever
  `population_sd < 1`) and was added in every direction, including the spiked ones.
* `sample_spiked_wishart()` and `sample_spiked_wishart_and_jac()` gave the wrong
  distribution when there were fewer variables than observations (the main diagonal
  used the wrong chi-squared degrees of freedom). The usual omics case of more
  variables than observations was unaffected.
* The dependence is estimated with the correct degrees of freedom when there is a
  design (n minus the number of design columns).

# dependentsimr v1.0.0

* Initial CRAN submission.
