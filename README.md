# Tree methods for length-biased survival data

We propose tree-based methods for length-biased right-censored (LBRC) data. 

Although existing tree methods for left-truncated right-censored (LTRC) data can be applied to LBRC data, they are inefficient because they ignore information from the truncation process.

In brief, the likelihood for LBRC data can be decomposed as:
```math
\underbrace{\mathcal{L}_F}_{\substack{\text{Full-likelihood of}\\ (T,A)}} = \underbrace{\mathcal{L}_C}_{\substack{\text{Conditional-likelihood of}\\ T \text{ given }A}} \times \underbrace{\mathcal{L}_M}_{\substack{\text{Marginal-likelihood of}\\ A}}
```
where $T$ is an observed failure time and $A$ is an observed truncation time. A unique feature of LBRC data is that the marginal distributions of truncation time $A$ and residual time $V=T-A$, are identical $f_A=f_V$.

Many statistical methods of LTRC data can be seen as utilizing only the conditional likelihood, not including the distributional information of observed truncation time, which leads to inefficiency.

## Our proposal

We extend the existing LTRC conditional inference tree and forest methods (LTRC-CIT/CIF [1,2]) to exploit LBRC-specific structures. This enhancement improves efficiency in both **tree construction** and **survival prediction**.

Specifically:

* **Tree construction** uses the score function of the LBRC full likelihood as the influence function.

* **Survival prediction** and **tree construction** both use one of two nonparametric survival estimators:

  1. **Full-likelihood NPMLE** [3] — most efficient but computationally heavier (EM-based).
      * **LBRC-CIT-F:** tree using the full-likelihood NPMLE for both construction and prediction.
      * **LBRC-CIF-F:** forest using the full-likelihood NPMLE for both construction and prediction.
  2. **Composite conditional-likelihood NPMLE** [4] — closed-form and faster, with some loss in efficiency.
      * **LBRC-CIT-C:** tree using the composite conditional-likelihood NPMLE for both construction and prediction.
      * **LBRC-CIF-C:** forest using the composite conditional-likelihood NPMLE for both construction and prediction.

Full algorithm details and simulation results are provided in our paper: *[link to paper]*.

## Fitting function

The function for LBRC-CITs and LBRC-CIFs are called `lbrccit` and `lbrccif` respectively, where they can be found in models_lbrc.R from methods folder.

The function works as follows:

```R
source("./methods/models_lbrc.R")		# functions for LBRC-CIT/CIFs
source("./methods/plot_lbrc.R")			# functions to plot the tree of LBRC-CITs
source("./methods/predictProb_lbrc.R")	# functions for prediction of fitted LBRC-CIT/CIFs
source("./methods/sbrier_lbrc.R")		# functions for calculating brier scores of LBRC-CIT/CIFs
source("./methods/tune.lbrccif.R")		# functions for tuning mtry parameter of LBRC-CIFs

obj <- lbrccit(formula, data, perm_test_est)
obj <- lbrccit(formula, data, perm_test_est)
```

On progress...

## References

[1] Wei Fu and Jeffrey S Simonoff. “Survival trees for left-truncated and right-censored data, with application to time-varying covariate data”. In: Biostatistics 18.2 (2017), pp. 352–369.

[2] Weichi Yao et al. “Ensemble methods for survival function estimation with time-varying covariates”. In: Statistical Methods in Medical Research 31.11 (2022), pp. 2217–2236

[3] Yehuda Vardi. “Multiplicative censoring, renewal processes, deconvolution and decreasing density: nonparametric estimation”. In: Biometrika 76.4 (1989), pp. 751–761.

[4] Yifan He and Yong Zhou. “Nonparametric and semiparametric estimators of restricted mean survival time under length-biased sampling”. In: Lifetime Data Analysis 26.4 (2020), pp. 761–788.
