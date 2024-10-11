# Load function needed for glm_mom_unknown_sigma
source("demo_glm_MoM/function_of_glm_mom.R")
# After first run, You can use following load for fast start
# source("demo_glm_MoM/dependecy.R")
# load("demo_glm_MoM/all_functions_and_objects_mom.RData")

# Define parameters
n <- 2000
p <- 600

X <- matrix(rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p)
# Generate alpha vector, alpha %*% Sigma %*% alpha should be between 0 and 2 preferably
set.seed(1234)
alpha <- runif(p, 0, 1)
alpha <- alpha / sqrt(sum(alpha ^ 2))
z <- X %*% alpha
A <- rbinom(n, 1, 1 / (1 + exp(-z)))
A <- as.numeric(A)

###################################################
# par_initial is the initial estimator of alpha %*% mu_x, alpha %*% Sigma %*% alpha, we will fix this initial problem later
# B_bootstrap is the number of bootstrap iterations used to calculate the variance
###################################################
par_initial <- c(0.1, 1.1)
model_results <- glm_mom_unknown_sigma(X, A, B_bootstrap = 1000, par_initial = par_initial)
str(model_results)
###################################################
# alpha_em_S_4 is estimator of alpha, when unknown mu_x = 0, similarly, suffixes with 4 are all for unknown mu_x = 0
# alpha_em_S_1 is estimator of alpha, when known mu_x = 0, similarly, suffixes with 1 are all for known mu_x = 0
# sols_em_S_1 is the estimator of alpha %*% Sigma %*% alpha
# sols_em_S_4 is the estimator of alpha %*% mu_x, alpha %*% Sigma %*% alpha
# f_em_S_1, f_em_S_4 are estimates of 0th, 1st, and 2nd order derivatives of the logistic function
# moments_em_S_4 are estimates of 4 moments: mA, mX_2, mAX_1, mAX_2
# Variables starting with 'var' are variance estimates using the bootstrap method
# Following are the corresponding values and variances estimates
###################################################
plot(model_results$alpha_em_S_4)
plot(model_results$var_alpha_split_sample_4)

print(model_results$sols_em_S_4)
print(model_results$var_paras_split_sample_4)

print(model_results$f_em_S_4)
print(model_results$var_f_paras_split_sample_4)

print(model_results$moments_em_S_4)
print(model_results$var_moments_split_sample_4[c(1:4)])