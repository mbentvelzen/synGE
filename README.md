# synGE

**Synthetic data Generation and Evaluation** for tabular healthcare data.

`{synGE}` provides four synthetic-data generators and a collection of evaluators
for comparing and evaluating real and synthetic tabular datasets. It supports comparing 
synthetic data generation methods and evaluating synthetic data quality on multiple aspect.
The code accompanies the manuscript **“Generating and evaluating realistic synthetic healthcare data:
a comprehensive framework”** by Mieke L.M. Bentvelzen, Laura Jochim, Majed
Ramadan, Simon C. Mastbergen, Paco M.J. Welsing, and Said el Bouhaddani.

The manuscript evaluates five complementary pillars: **marginal fidelity,
bivariate correlations, multivariate fidelity, prediction-based
distinguishability, and privacy**. The manuscript studies quality, stability across
repetitions, and generalizability under cross-validation. **The current package
provides individual generators, evaluators, and plotting functions; it does not
yet provide an automated quality/stability/generalizability pipeline.**

> **Healthcare privacy:** Similarity scores and nearest-neighbour privacy metrics
> do not certify that synthetic patient-level data are anonymous or safe to
> publish. Assess disclosure risk, regulatory requirements and intended use
> separately before sharing any sensitive data.

## Installation

Install the package from GitHub:

```r
install.packages("remotes")           # once, if needed
remotes::install_github("mbentvelzen/synGE")
library(synGE)
```

Only packages used throughout the package are mandatory dependencies.
Most method-specific R packages are optional (`Suggests`) and must be installed
before using the corresponding functions. For example:

```r
install.packages(c("synthpop", "e1071", "FNN"))
```

For other methods, install the relevant dependencies, for example:

```r
install.packages(c(
  "rvinecopulib", "kde1d",  # vine copula
  "reticulate",             # CTGAN/TVAE Python interface
  "caret", "glmnet", "randomForest", "pROC",  # prediction
  "isotree", "DataSimilarity", "fasano.franceschini.test", # multivariate
  "ggnewscale", "ggrepel",  # optional plots
  "cluster", "proxy", "transport", "gower"  # distance metrics
))
```

The CTGAN and TVAE functions require the **Python** `ctgan` module as well
as the R package `reticulate`. In the Python environment used by reticulate:

```r
reticulate::py_install("ctgan")
reticulate::py_module_available("ctgan")
```

Configure the intended Python environment **before** reticulate first imports
Python in that R session. GPU availability and library versions can affect
reproducibility of neural-network generators.

## Quick start

The following example uses only numeric columns so that the moment comparison
has no categorical variables to omit:

```r
library(synGE)

real_df <- iris[, 1:4]
syn_df <- synth_sp(real_df, seed = 123)

# Compare marginal moments. Result: list with by_moment and meta.
moments <- eval_moments(real_df, syn_df)
head(moments$by_moment)

# Compare pairwise correlations. Result: list with result and cor_dif.
correlations <- eval_Cor(real_df, syn_df)
head(correlations$result)

# Optional figure objects:
eval_moments_calplot(moments)
eval_Cor_plot(correlations, plot = "barplot")
```

The functions do **not** automatically create a real-to-real reference.
To reproduce that part of the manuscript, bootstrap two samples from real data
and run the same evaluator on those samples repeatedly.

## Synthetic data generators

| Function | Method | Backend |
| --- | --- | --- |
| `synth_sp()` | Sequential conditional synthesis (CART by default) | R `synthpop` |
| `synth_vc()` | Vine copula with estimated marginal distributions | R `rvinecopulib` and `kde1d` |
| `synth_ctgan()` | Conditional tabular GAN | Python `ctgan` via `reticulate` |
| `synth_vae()` | Tabular variational autoencoder (TVAE) | Python `ctgan` via `reticulate` |

The default `seed` in all functions is `NULL`: pass a seed explicitly when reproducibility
matters. Do not pass a `discrete_columns` argument to `synth_vc()`;
factor/character columns are handled internally. CTGAN and TVAE *do* require a
`discrete_columns` argument listing their discrete columns.


Both `synth_ctgan()` and `synth_vae()` accept `tune = TRUE` and an optional
`tune_grid`; their `tune_eval` is either `"MWS"` (requires `eval_MWS()` and its
optional distance dependencies) or `"quality"` (a KS/TVD-based score). This
is the tuning approach **implemented in these R files**; it should not be
confused with the SDMetrics-based tuning described in the manuscript.

## Evaluators

| Pillar | Main functions | Purpose |
| --- | --- | --- |
| Marginal | `eval_moments()` and `eval_moments_heatmap()`, `eval_moments_calplot()` | Mean, variance, SD, skewness and kurtosis |
| Marginal | `eval_JSD()`, `eval_JSD_plot()` | Jensen–Shannon divergence |
| Marginal | `eval_KS()`, `eval_KS_plot()` | Two-sample Kolmogorov–Smirnov |
| Marginal | `eval_CvM()`, `eval_CvM_plot()` | Cramér–von Mises |
| Marginal | `eval_AD()`, `eval_AD_plot()` | Anderson–Darling |
| Marginal | `eval_pMse()` | Propensity-score mean squared error via `synthpop` |
| Bivariate | `eval_Cor()`, `eval_Cor_plot()` | Pairwise correlation differences |
| Multivariate | `eval_MKS()` | Multivariate KS / Fasano–Franceschini test |
| Multivariate | `eval_MMD()` | Maximum mean discrepancy |
| Multivariate | `eval_PCA()` | PCA loading cosine similarity |
| Multivariate | `eval_IF()`, `eval_IF_plot()` | Isolation Forest anomaly-score similarity |
| Multivariate (additional) | `eval_MWS()` | Optimal-transport / Wasserstein-style distance |
| Prediction-based | `eval_logreg()`, `eval_randomForest()` | Real-versus-synthetic classifier performance |
| Privacy | `DCR()`, `NNDR()` | Distance to closest record and nearest-neighbour distance ratio |


Some functions run permutations, cross-validation or nearest-neighbour/optimal-
transport computations. They can be expensive on large datasets. Interpretation
of results also depends on sample size, variable types, and the metric's target.

## Relationship to the manuscript

The paper compares synthpop, vine copulas, VAE and CTGAN across healthcare
datasets of different sizes and dimensions. It argues for interpreting multiple
metric families together, and for comparing synthetic-to-real discrepancies
against real-to-real sampling variability. It describes three experimental
layers: one-dataset quality, repeated-run stability, and held-out-fold
generalizability. These are **study procedures**, not public wrapper functions
in the current R package.

The code provides research components. It is not a validated clinical decision
support system or a standalone privacy guarantee.


## Citation

If you use the methods in research, cite the associated manuscript:

> Bentvelzen, M. L. M., Jochim, L., Ramadan, M., Mastbergen, S. C.,
> Welsing, P. M. J., and el Bouhaddani, S. *Generating and evaluating
> realistic synthetic healthcare data: a comprehensive framework*.
> The manuscript is currently under consideration. Journal and DOI will be updated.

## License

GNU GENERAL PUBLIC LICENSE. See [LICENSE](LICENSE) for the full text.
