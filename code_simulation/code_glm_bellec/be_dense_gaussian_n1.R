rm(list = ls())

library("parallel")
library('doParallel')
library('doRNG')
library('foreach')
library('glmnet')
library('MASS')
library("SMUT")

# Set parameters -----

Is_sparse =  FALSE
Is_sparse_only_1 =  FALSE
Is_Rad =  FALSE

nsmall = 1000
nbig = 5000
ns.length = 5
ns <- seq(nsmall, nbig, length.out = ns.length)

# Uncomment to use specific sections of 'ns'
# ns <- ns[c(1:4)]
ns <- ns[5]

N.replicates = 50

p.n.ratio = 1.2
lambda_value = exp(seq(log(10), log(.05), length.out = 12))
mu_x = 0

R.version <- version
glmnet.version <- packageVersion("glmnet")
system.version <- Sys.info()
sim.time <- format(Sys.time(), "%Y%m%d_%H%M%S")

filename.out = paste0(
  "data/bellec_n",
  ns[length(ns)],
  "_p_",
  p.n.ratio,
  "_sparse_",
  as.integer(Is_sparse),
  "_one_",
  as.integer(Is_sparse_only_1),
  "_Rad_",
  as.numeric(Is_Rad),
  "_",
  sim.time
)
filename.tracking = paste0(
  "track/bellec_n",
  ns[length(ns)],
  "_p_",
  p.n.ratio,
  "_sparse_",
  as.integer(Is_sparse),
  "_one_",
  as.integer(Is_sparse_only_1),
  "_Rad_",
  as.numeric(Is_Rad),
  "_",
  sim.time
)

# Execution section -----

set.seed(123)  # Fix the random seed to ensure consistency in alpha generation
# Define a function to generate alpha
generate_alpha <- function(p, omega_11, Is_sparse, Is_sparse_only_1) {
  if (Is_sparse & Is_sparse_only_1) {
    alpha <- rep(0, p)
    alpha[1] <- 1
  } else if (Is_sparse) {
    s <- round(sqrt(p))
    alpha <- rep(0, p)
    alpha[1:s] <- sqrt(omega_11 / s)
  } else {
    alpha <- runif(p, -sqrt(3 * omega_11 / p), sqrt(3 * omega_11 / p))
  }
  return(alpha)
}

# Generate alpha for each n
alpha_list <- lapply(ns, function(n) {
  p <- round(n * p.n.ratio)
  generate_alpha(p, omega_11 = 1, Is_sparse, Is_sparse_only_1)  # Assume omega_11 as 1, adjust as needed
})

# Generate a unique integer as seed using the current time
unique_seed <- as.integer(Sys.time())
set.seed(unique_seed)
# Generate random integers
random_replicates <- sample(1:1000000, N.replicates)
# Generate specs containing n, p, and alpha
specs <- list()
index <- 1
for (i in seq_along(ns)) {
  n <- ns[i]
  p <- round(n * p.n.ratio)
  for (replicate in random_replicates) {
    specs[[index]] <- list(n = n, p = p, alpha = alpha_list[[i]], seed = replicate)
    index <- index + 1
  }
}
# Create and register parallel cores
numCores <- detectCores()
cl <- makeCluster(numCores)

clusterEvalQ(cl, {
  source("utils/function_of_bellec.R")
  library('glmnet')
  library('MASS')
  library("SMUT")
})
# Export necessary objects and functions to the cluster
clusterExport(
  cl,
  c(
    "filename.tracking",
    "Is_Rad",
    "Is_sparse",
    "Is_sparse_only_1",
    "mu_x",
    "lambda_value",
    "p.n.ratio"
  )
)

experiment.data <- parLapply(cl, specs, function(spec) {
  run.replicate(
    n = spec$n,
    p = spec$p,
    alpha = spec$alpha,
    omega_11 = ifelse(1 == 1, 1, round(spec[2] * p.n.ratio)),
    Is_Rad = Is_Rad,
    Is_sparse = Is_sparse,
    Is_sparse_only_1 = Is_sparse_only_1,
    mu_x = mu_x,
    lambda_value = lambda_value,
    seed = spec$seed,
    filename.tracking = paste0(filename.tracking, ".txt")
  )
})
save(
  experiment.data,
  R.version,
  glmnet.version,
  system.version,
  sim.time,
  unique_seed,
  random_replicates,
  file = paste0(filename.out, ".Rda")
)

stopCluster(cl)
