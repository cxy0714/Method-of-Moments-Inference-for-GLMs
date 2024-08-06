rm(list = ls())


library('doParallel')
library('doRNG')
library('foreach')
library('glmnet')
library('MASS')
library(mvtnorm)
library(lsbclust)
library(rje)
library(nleqslv)
library(rootSolve)
library("SMUT")

N.replicates <- 100
n.p.ratio <- 1.2
n_list <-  seq(1000, 5000, length.out = 5)

Is_sparse =  TRUE
Is_sparse_only_1 =  TRUE
Is_Rad =  TRUE
mu_x = 0

indices_first_1 <- 1:10
indices_middle_2 <- 100:109
indices_selected <- c(indices_first_1, indices_middle_2)

# 使用当前时间生成一个唯一的整数作为种子
unique_seed <- as.integer(Sys.time())
set.seed(unique_seed)
# 生成随机整数
seeds <- sample(1:1000000, N.replicates)

numCores <- detectCores()

R.version <- version
glmnet.version <- packageVersion("glmnet")
system.version <- Sys.info()
sim.time <- format(Sys.time(), "%Y%m%d_%H%M%S")


# function ----
epsilon <- 0
integrand_forward_mapping <- function (x, mu, sigma_sq) {
  (epsilon + (1 - epsilon) * exp(x) / (1 + exp(x))) * dnorm(x, mu, sqrt(sigma_sq))
}
integrand_forward_mapping_quad <- function (x, mu, sigma_sq) {
  (1 - epsilon) * exp(x) / (1 + exp(x)) ^ 2 * dnorm(x, mu, sqrt(sigma_sq))
}

integrand_forward_mapping_hess <- function (x, mu, sigma_sq) {
  (1 - epsilon) * exp(x) / (1 + exp(x)) ^ 2 * (1 - 2 * exp(x) / (1 + exp(x))) * dnorm(x, mu, sqrt(sigma_sq))
}

# plot(seq(-10,10,0.01),sapply(seq(-10,10,0.01),function(i)(integrand_forward_mapping_quad(i,0,1))))
# plot(seq(-10,10,0.01),sapply(seq(-10,10,0.01),function(i)(integrand_forward_mapping_hess(i,0,1))))

forward_mapping_0 <- function (mu, sigma_sq) {
  integrate(
    integrand_forward_mapping,
    -100,
    100,
    mu = mu,
    sigma_sq = sigma_sq,
    stop.on.error = FALSE
  )$value
}
forward_mapping_1_mean <- function (mu, sigma_sq) {
  integrate(
    integrand_forward_mapping_quad,
    -100,
    100,
    mu = mu,
    sigma_sq = sigma_sq,
    stop.on.error = FALSE
  )$value * mu
}
forward_mapping_1_square <- function (mu, sigma_sq) {
  integrate(
    integrand_forward_mapping_quad,
    -100,
    100,
    mu = mu,
    sigma_sq = sigma_sq,
    stop.on.error = FALSE
  )$value ^ 2 * sigma_sq
}
forward_mapping_1 <- function (mu, sigma_sq) {
  integrate(
    integrand_forward_mapping_quad,
    -100,
    100,
    mu = mu,
    sigma_sq = sigma_sq,
    stop.on.error = FALSE
  )$value
}
forward_mapping_2 <- function (mu, sigma_sq) {
  integrate(
    integrand_forward_mapping_hess,
    -100,
    100,
    mu = mu,
    sigma_sq = sigma_sq,
    stop.on.error = FALSE
  )$value
}

## simu ----  
# 创建并行核并注册
cl <- makeCluster(numCores)
registerDoParallel(cl)

for (j in c(1:length(n_list))) {
  set.seed(random_replicates[j])
  
  n <- n_list[j]
  p <- round(n.p.ratio * n)
  omega_11 <- ifelse(1 == 1, 1, p)
  
  
  if (Is_sparse & Is_sparse_only_1) {
    alpha <- rep(0, p)
    alpha[1] <- 1
  } else  if (Is_sparse) {
    s <- round(sqrt(p))
    alpha <- rep(0, p)
    alpha[c(1:s)] <- sqrt(omega_11 / s)
  } else {
    alpha <- runif(p, -sqrt(3 * omega_11 / p), sqrt(3 * omega_11 / p))
  }
  s <- sum(abs(alpha) > 0)
  N.replicates <- N.replicates
  mu <- rep(mu_x, p)
  
  
  par1_truth <- sum(alpha * mu)
  par2_truth <- sum(alpha * alpha / omega_11)
  
  par_truth <- c(par1_truth, par2_truth)
  
  term0_truth <- forward_mapping_0(par1_truth, par2_truth)
  term1_truth <- forward_mapping_1(par1_truth, par2_truth)
  
  
  mA_truth <- term0_truth
  mX_2_truth <- as.numeric(sum(mu * mu) * omega_11)
  mXA_1_truth <- as.numeric(mA_truth * mX_2_truth + term1_truth * par1_truth)
  mXA_2_truth <- as.numeric(
    mA_truth ^ 2 * mX_2_truth + term1_truth ^ 2 * par2_truth + 2 * mA_truth * term1_truth * par1_truth
  )
  moments_truth <-  c(mA_truth, mX_2_truth, mXA_1_truth, mXA_2_truth)
  fn <- function (x) {
    term0 <- forward_mapping_0(x[1], x[2])
    term1 <- forward_mapping_1(x[1], x[2])
    
    a <- -mXA_1_truth + term1 * x[1] + mA_truth * mX_2_truth
    b <- -mXA_2_truth +  term1 ^ 2 * x[2] + 2 * mA_truth * term1 * x[1] + mA_truth ^
      2 * mX_2_truth
    return(c(a, b))
  }
  epsilon_par <- rnorm(2, 0, 1 / p)
  sols_oracle <- nleqslv(par_truth + epsilon_par, fn)$x
  
  
  sim_res <-  foreach(
    i = 1:N.replicates,
    .packages = c("glmnet", "MASS", "mvtnorm", "lsbclust", "nleqslv", "SMUT")
  ) %dorng% {
    result <- tryCatch({
      if (Is_Rad) {
        X_1 <- matrix(rbinom(n * p, size = 1, prob = 0.5),
                      nrow = n,
                      ncol = p)
        
        X_1 <- (2 * X_T - 1) / sqrt(omega_11)
        
      } else  {
        X_1 <- matrix(rnorm(n * p, mean = 0, sd = 1 / sqrt(omega_11)),
                      nrow = n,
                      ncol = p)
      }
      
      X_1 <- X_1 + mu_x
      
      z_1 <- X_1 %*% alpha
      A_1 <- rbinom(n, 1, 1 / (1 + exp(-z_1)))
      A_1 <- as.numeric(A_1)
      
      H <-  eigenMapMatMult(X_1, t(X_1)) * omega_11
      H <- H - diag(diag(H))
      vector_one <- rep(1, n)
      
      mA_em <- mean(A_1)
      mX_2_em <- as.numeric(sum(vector_one * eigenMapMatMult(H, vector_one)) / (n * (n - 1)))
      mXA_1_em <- as.numeric(sum(A_1 * eigenMapMatMult(H, vector_one)) / (n * (n - 1)))
      mXA_2_em <- as.numeric(sum(A_1 * eigenMapMatMult(H, A_1)) / (n * (n - 1)))
      
      
      moments_em <- c(mA_em, mX_2_em, mXA_1_em, mXA_2_em)
      fn_em <- function (x) {
        term0 <- forward_mapping_0(x[1], x[2])
        term1 <- forward_mapping_1(x[1], x[2])
        
        
        a <- -mXA_1_em + term1 * x[1] + mA_em * mX_2_em
        b <- -mXA_2_em +  term1 ^ 2 * x[2] + 2 * mA_em * term1 * x[1] + mA_em ^
          2 * mX_2_em
        return(c(a, b))
      }
      epsilon_par <- rnorm(2, 0, 1 / sqrt(p))
      sols_em <- nleqslv(par_truth + epsilon_par, fn_em)$x
      
      mu_alpha_est <- sols_em[1]
      alpha_L2_est_N <- sols_em[2]
      
      Ax_1_mean <- colMeans(X_1 * A_1)
      f_1_est <- forward_mapping_1(mu_alpha_est, alpha_L2_est_N)
      f_0_est <- mA_em
      mu_est <- colMeans(X_1)
      alpha_est_N <-  if (f_1_est > 0) {
        omega_11 * (Ax_1_mean - f_0_est * mu_est) / f_1_est
      } else {
        rep(NA, p)
      }
      
      
      alpha_est_N <- alpha_est_N[indices_selected]
      return(
        list(
          success = TRUE,
          alpha_est_N = alpha_est_N,
          mu_alpha_est = mu_alpha_est,
          alpha_L2_est_N = alpha_L2_est_N,
          sols_em = sols_em,
          moments_em = moments_em
        )
      )
    }, error = function(e) {
      # 如果出错，记录错误并返回失败的信号
      cat("Error: ", e$message, "\n")
      return(list(
        success = FALSE,
        message = e$message,
        seed = random_replicates[j],
        n = n,
        p = p
      ))
      
    })
    
    return(result)
  }
  
  ## save ----
  mu_alpha_est <- alpha_L2_est_N <- rep(NA, N.replicates)
  moments_em <-  matrix(NA, nrow = N.replicates, ncol = 4)
  sols_em <- matrix(NA, nrow = N.replicates, ncol = 2)
  alpha_est_N  <-  matrix(NA, nrow = N.replicates, ncol = length(indices_selected))
  for (i in 1:N.replicates) {
    
    if (sim_res[[i]]$success){
      alpha_L2_est_N[i] <- sim_res[[i]]$alpha_L2_est_N
      alpha_est_N[i, ] <- sim_res[[i]]$alpha_est_N
      mu_alpha_est[i] <- sim_res[[i]]$mu_alpha_est
      sols_em[i, ] <- sim_res[[i]]$sols_em
      moments_em[i, ] <- sim_res[[i]]$moments_em
    }
  }
  # options(scipen = 999)
  save(
    n,
    p,
    s,
    alpha,
    N.replicates,
    omega_11,
    
    
    moments_truth,
    moments_em,
    par_truth,
    sols_em,
    
    mu_alpha_est,
    alpha_L2_est_N,
    
    alpha_est_N,
    indices_selected,
    R.version,
    glmnet.version,
    system.version,
    sim.time,
    unique_seed,
    random_replicates,
    omega_11,
    file = paste0(
      "data/MoM_n_",
      n,
      "_p_",
      p,
      "_iter_",
      N.replicates,
      "_sparse_",
      as.numeric(s),
      "_omegma",
      as.numeric(omega_11),
      "_Ra_",
      as.numeric(Is_Rad),
      "_",
      sim.time,
      ".Rda"
    )
  )
}

# 关闭聚核
stopCluster(cl)
