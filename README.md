# Tree methods for length-biased survival data

We propose tree-based methods for length-biased right-censored (LBRC) data. 

Although existing tree methods for left-truncated right-censored (LTRC) data can be applied to LBRC data, they are inefficient because they ignore information from the truncation process.

In brief, the likelihood for LBRC data can be decomposed as:
$$
\underbrace{\mathcal{L}_F}_{\substack{\text{Full-likelihood of}\\ (T,A)}} = \underbrace{\mathcal{L}_C}_{\substack{\text{Conditional-likelihood of}\\ T \text{ given }A}} \times \underbrace{\mathcal{L}_M}_{\substack{\text{Marginal-likelihood of}\\ A}}
$$
where $T$ is an observed failure time and $A$ is an observed truncation time. A unique feature of LBRC data is that the marginal distributions of truncation time $A$ and residual time $V=T-A$, are identical $f_A=f_V$.

Many statistical methods of LTRC data can be seen as utilizing only the conditional likelihood, not including the distributional information of observed truncation time, which leads to inefficiency.

## Our proposal

We extend the existing LTRC conditional inference tree and forest methods (LTRC-CIT/CIF [1,2]) to exploit LBRC-specific structures. This enhancement improves efficiency in both **tree construction** and **survival prediction**.

Specifically:

* **Tree construction** uses the score function of the LBRC full likelihood as the influence function.

* **Survival prediction** and **tree construction** both use one of two nonparametric survival estimators:

  1. **Full-likelihood NPMLE** [3] — most efficient but computationally heavier (EM-based).
      **LBRC-CIT-F:** tree using the full-likelihood NPMLE for both construction and prediction.
      **LBRC-CIF-F:** forest using the full-likelihood NPMLE for both construction and prediction.

     **Composite conditional-likelihood NPMLE** [4] — closed-form and faster, with some loss in efficiency.
      **LBRC-CIT-C:** tree using the composite conditional-likelihood NPMLE for both construction and prediction.
      **LBRC-CIF-C:** forest using the composite conditional-likelihood NPMLE for both construction and prediction.

Full algorithm details and simulation results are provided in our paper: *[link to paper]*.
