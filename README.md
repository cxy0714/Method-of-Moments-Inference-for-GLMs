
<h1 align="center" style="margin-bottom:0px; border-bottom:0px; padding-bottom:0px">Method-of-Moments-Inference-for-GLMs</h1>
<h3 align="center" style="margin-bottom:0px; border-bottom:0px; padding-bottom:0px">Method-of-Moments Inference for GLMs and Doubly Robust Functionals under Proportional Asymptotics</h3>
<p align="center" style="margin-bottom:0px; border-bottom:0px; padding-bottom:0px">Xingyu Chen, Lin Liu, Rajarshi Mukherjee</p>

<p align="center">
    <a style="text-decoration:none !important;" href=" " alt="arXiv"><img src="https://img.shields.io/badge/paper-arXiv-red" /></a>
</p>

We introduce moments-based identification strategies for statistical functionals within high-dimensional generalized linear models (GLMs), where the dimension \( p \) is proportional to the sample size \( n \) and **the covariance matrix of covariates is known**. Key advantages of our methods include:

1. **Computational Efficiency**: Relying on only a few low-dimensional moments of the data, our methods are computationally efficient.

2. **Generalization to Observational Studies**: Our strategies extend to inferential techniques for average treatment effects (ATE) and mean estimands under missing data without requiring sample-splitting or cross-fitting.

3. **Root-n-consistent and Normal Estimators**: The estimators are Root-n-consistent and asymptotically normal (CAN) under Gaussian covariates and demonstrate universal applicability beyond Gaussian cases.


We implemented two simulation settings.

#### Estimating Logistic Model Coefficients and Quadratic Forms
The first involves estimating GLM (logistic model) coefficients and quadratic forms, where we compare our methods with that of **Bellec (2022)**[^1].

- Our method's code is located in the folder `code_simulation/code_glm_mom`. The code is consistent across all files, with variations only in the initial settings section.
- Bellec's method's code can be found in the folder `utils/function_of_bellec.R` and `code_simulation/code_glm_bellec`. The main function resides in `utils/function_of_bellec.R`, complete with detailed notes. In `code_simulation/code_glm_bellec`, the files are identical except for the initial settings section.

To generate Figures 1, 3, 7, 9, and 11, use `glm_clean_cluster_bellec.R`, `glm_clean_cluster_mom.R`, and `glm_plot_lambda.R` in the folder `code_clean_plot` after obtaining the data.

To generate Figures 2, 4, 6, 8, and 12, use `glm_clean_cluster_bellec.R`, `glm_clean_cluster_mom.R`, and `glm_plot_mom_hist_qqnorm.R` in the folder `code_clean_plot` after obtaining the data.

#### Estimating the Mean of a Response Under Missing Data
The second simulation setting involves estimating the mean of a response under random missingness, where the outcome model is a linear model and the missingness mechanism is a logistic model. We compare our methods with those of **Celentano and Wainwright (2023)**[^2].

- Our method's code is in the folder `code_simulation/code_mar_mom`. The code is the same in all files except for the initial setting section.
- Celentano and Wainwright's method's code is in the folder `code_simulation/code_mar_celentano`, with minor adjustments made to their original code available at [this GitHub repository](https://github.com/mcelentano/Debiasing_for_missing_data). The code is the same in all files except for the initial setting section.

To generate Figures 5, 10, 13, and 15, use `mar_clean_cluster_mom.R` and `mar_plot_lambda.R` in the folder `code_clean_plot` after obtaining the data.

To generate Figures 6, 12, 14, and 16, use `mar_clean_cluster_mom.R` and `mar_plot_mom_hist_qqnorm.R` in the folder `code_clean_plot` after obtaining the data.

Any questions about the code, feel free to contact me at xingyuchen0714@sjtu.edu.cn. By the way, you should really encapsulate functions that are used multiple times！ I'm already a mess anyway. I will fight this mountain of crappy code later. 



[^1]: **Bellec P C. Observable adjustments in single-index models for regularized M-estimators**[J]. arXiv preprint arXiv:2204.06990, 2022.

[^2]: **Celentano M, Wainwright M J. Challenges of the inconsistency regime: Novel debiasing methods for missing data models**[J]. arXiv preprint arXiv:2309.01362, 2023.