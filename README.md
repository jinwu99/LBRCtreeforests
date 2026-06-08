# Tree methods for length-biased survival data
We propose tree-based methods for length-biased right-censored (LBRC) data. 
Although existing tree methods for left-truncated right-censored (LTRC) data can be applied to LBRC data, they are inefficient because they ignore information from the truncation process.
In brief, LBRC data are a special case of LTRC data in which disease onset follows a stationary Poisson process. The likelihood for LBRC data can be decomposed as:
```math
\underbrace{\mathcal{L}_F}_{\substack{\text{Full-likelihood of}\\ (T,A)}} = \underbrace{\mathcal{L}_C}_{\substack{\text{Conditional-likelihood of}\\ T \text{ given }A}} \times \underbrace{\mathcal{L}_M}_{\substack{\text{Marginal-likelihood of}\\ A}}
```
where $T$ is an observed failure time and $A$ is an observed truncation time. A unique feature of LBRC data is that the marginal distributions of truncation time $A$ and residual time $V=T-A$, are identical $f_A=f_V$​.
Many statistical methods of LTRC data can be seen as utilizing only the conditional likelihood, not including the distributional information of observed truncation time, which leads to inefficiency.

## Repository Structure & Navigation

This repository is organized into two main directories tailored to your specific goals:

* 📁 **[`pkg/`](./pkg)** **(For End-Users):** Contains the R implementation of the proposed LBRC-CIT and LBRC-CIF methods. **Navigate to this folder if you want to apply our algorithms to your own dataset.** It includes user-friendly functions for model fitting, hyperparameter tuning, and survival prediction.
* 📁 **[`analysis/`](./analysis)** **(For Reviewers & Reproducibility):** Contains all the code required to reproduce the simulation studies, real-data application, figures, and tables presented in our manuscript and supplementary material. **Navigate to this folder if you want to verify or reproduce our published results.**

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

## Publication
* Title: "Tree-Based Methods for Length-Biased Survival Data"
* Authors: Jinwoo Lee, Donghwan Lee, Hyunwoo Lee, and Jiyu Sun
* Code Author: Jinwoo Lee (josh_99@ncc.re.kr)
The code was written/evaluated in R with the following software versions:
```R
> sessionInfo()
R version 4.4.0 (2024-04-24 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Matrix products: default


locale:
[1] LC_COLLATE=Korean_Korea.utf8  LC_CTYPE=Korean_Korea.utf8    LC_MONETARY=Korean_Korea.utf8
[4] LC_NUMERIC=C                  LC_TIME=Korean_Korea.utf8    

time zone: Asia/Seoul
tzcode source: internal

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] pdftools_3.5.0     magick_2.8.5       purrr_1.0.2        scales_1.4.0       broom_1.0.5       
 [6] patchwork_1.3.0    emmeans_1.11.1     stringr_1.5.1      simplecolors_0.1.2 ggpubr_0.6.0      
[11] ggforce_0.4.2      ggplot2_4.0.3      dplyr_1.1.4        tidyr_1.3.1        Rcpp_1.0.12       
[16] rsample_1.3.0      survival_3.5-8     partykit_1.2-22    mvtnorm_1.2-4      libcoin_1.0-10    

loaded via a namespace (and not attached):
 [1] gtable_0.3.6       rstatix_0.7.2      lattice_0.22-6     vctrs_0.6.5        tools_4.4.0       
 [6] generics_0.1.3     sandwich_3.1-0     parallel_4.4.0     tibble_3.2.1       fansi_1.0.6       
[11] pkgconfig_2.0.3    Matrix_1.7-0       RColorBrewer_1.1-3 S7_0.2.2           lifecycle_1.0.4   
[16] compiler_4.4.0     farver_2.1.1       codetools_0.2-20   carData_3.0-5      Formula_1.2-5     
[21] pillar_1.9.0       furrr_0.3.1        car_3.1-2          MASS_7.3-60.2      multcomp_1.4-26   
[26] rpart_4.1.23       abind_1.4-5        parallelly_1.38.0  tidyselect_1.2.1   digest_0.6.35     
[31] inum_1.0-5         stringi_1.8.4      future_1.34.0      listenv_0.9.1      forcats_1.0.0     
[36] splines_4.4.0      polyclip_1.10-7    colorspace_2.1-0   cli_3.6.2          magrittr_2.0.3    
[41] utf8_1.2.4         TH.data_1.1-2      withr_3.0.0        backports_1.4.1    estimability_1.5.1
[46] globals_0.16.3     qpdf_1.3.4         ggsignif_0.6.4     askpass_1.2.0      zoo_1.8-12        
[51] coda_0.19-4.1      rlang_1.1.3        xtable_1.8-4       glue_1.7.0         tweenr_2.0.3      
[56] rstudioapi_0.16.0  R6_2.5.1
```

## References
[1] Wei Fu and Jeffrey S Simonoff. “Survival trees for left-truncated and right-censored data, with application to time-varying covariate data”. In: Biostatistics 18.2 (2017), pp. 352–369.

[2] Weichi Yao et al. “Ensemble methods for survival function estimation with time-varying covariates”. In: Statistical Methods in Medical Research 31.11 (2022), pp. 2217–2236

[3] Yehuda Vardi. “Multiplicative censoring, renewal processes, deconvolution and decreasing density: nonparametric estimation”. In: Biometrika 76.4 (1989), pp. 751–761.

[4] Yifan He and Yong Zhou. “Nonparametric and semiparametric estimators of restricted mean survival time under length-biased sampling”. In: Lifetime Data Analysis 26.4 (2020), pp. 761–788.
