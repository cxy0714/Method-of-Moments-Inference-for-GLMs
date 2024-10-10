library('glmnet')
library('MASS')
library(mvtnorm)
library(lsbclust)
library(rje)
library(nleqslv)
library(rootSolve)
library("SMUT")
# function -----
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
    integrand_forward_mapping,-100,
    100,
    mu = mu,
    sigma_sq = sigma_sq,
    stop.on.error = FALSE
  )$value
}
forward_mapping_1 <- function (mu, sigma_sq) {
  integrate(
    integrand_forward_mapping_quad,-100,
    100,
    mu = mu,
    sigma_sq = sigma_sq,
    stop.on.error = FALSE
  )$value
}
forward_mapping_2 <- function (mu, sigma_sq) {
  integrate(
    integrand_forward_mapping_hess,-100,
    100,
    mu = mu,
    sigma_sq = sigma_sq,
    stop.on.error = FALSE
  )$value
}
expectation_G <- function(sigma_sq) {
  g <- forward_mapping_1(0, sigma_sq)
  G <- g ^ 2 * sigma_sq
  return(G)
}

# 使用parSapply进行并行计算.
tol <- 1e-5
sigma_L2_values <- seq(0.01, 2, by = tol)
expectation_values <- sapply(sigma_L2_values, expectation_G)

# 创建反函数
inverse_expectation <- approxfun(expectation_values, sigma_L2_values, rule = 2)

### Jacobian----

library(nleqslv)
# 定义目标函数
moments_to_paras_4_reduced_func <- function(moments_reduced, par_initial) {
  
  fn_4 <- function(x) {
    term0 <- forward_mapping_0(x[1], x[2])
    term1 <- forward_mapping_1(x[1], x[2])
    
    a <- -moments_reduced[1] + term0
    b <- -moments_reduced[2] + term1 ^ 2 * x[2] 
    return(c(a, b))
  }
  
  # 尝试不同的初始值，直到找到一个有效的解
  max_attempts <- 10
  attempt <- 1
  while (attempt <= max_attempts) {
    tryCatch({
      sols <- nleqslv(par_initial, fn_4, control = list(maxit = 100))$x
      return(sols)
    }, error = function(e) {
      # 如果解不出来，对初始值进行微小扰动
      par_initial <- par_initial + rnorm(2, mean = 0, sd = 0.1)
      attempt <- attempt + 1
    })
  }
  
  # 如果尝试了多次仍然无法找到解，返回NA
  return(c(NA, NA))
}

moments_to_paras_func <- function(moments, par_initial, type) {
  if (type == "1") {
    mXA_2 <- moments
    
    sols <- inverse_expectation(mXA_2)
    
  } else if (type == "4") {
    mA <- moments[1]
    mX_2 <- ifelse(moments[2]>0,moments[2],1e-5)
    mXA_1 <- moments[3]
    mXA_2 <- ifelse(moments[4]>0,moments[4],1e-5)
    
    fn_4 <- function(x) {
      term0 <- forward_mapping_0(x[1], x[2])
      term1 <- forward_mapping_1(x[1], x[2])
      term2 <- forward_mapping_2(x[1], x[2])
      
      a <- -mXA_1 + term1 * x[1] + mA * mX_2
      b <- -mXA_2 + term1 ^ 2 * x[2] + 2 * mA * term1 * x[1] + mA ^ 2 * mX_2
      return(c(a, b))
    }
    
    sols <- nleqslv(par_initial, fn_4, control = list(maxit = 1000))$x
    
  } else if (type == "6") {
    mA <- moments[1]
    mX_2 <- ifelse(moments[2]>0,moments[2],1e-5)
    mXA_1 <- moments[3]
    mXA_2 <- ifelse(moments[4]>0,moments[4],1e-5)
    mAY <- moments[5]
    mXAY_1 <- moments[6]
    
    
    fn_6 <- function (x) {
      term0 <- forward_mapping_0(x[1], x[2])
      term1 <- forward_mapping_1(x[1], x[2])
      term2 <- forward_mapping_2(x[1], x[2])
      
      a <- -mXA_1 + term1 * x[1] + mA * mX_2
      b <- -mXA_2 +  term1 ^ 2 * x[2] + 2 * mA * term1 * x[1] + mA ^
        2 * mX_2
      c <- -mAY + mA * x[4] + term1 * x[3]
      d <- -mXAY_1 + (mA + mXA_1) * x[4] + (mX_2 * term1 + term2 * x[1]) * x[3]
      
      return(c(a, b, c, d))
    }
    
    sols <- nleqslv(par_initial, fn_6)$x
  } else if (type == "7") {
    mA <- moments[1]
    mX_2 <- ifelse(moments[2]>0,moments[2],1e-5)
    mXA_1 <- moments[3]
    mXA_2 <- ifelse(moments[4]>0,moments[4],1e-5)
    mAY <- moments[5]
    mXAY_1 <- moments[6]
    mXAYXA <- moments[7]
    
    
    fn_7 <- function (x) {
      term0 <- forward_mapping_0(x[1], x[2])
      term1 <- forward_mapping_1(x[1], x[2])
      term2 <- forward_mapping_2(x[1], x[2])
      
      a <- -mXA_1 + term1 * x[1] + mA * mX_2
      b <- -mXA_2 +  term1 ^ 2 * x[2] + 2 * mA * term1 * x[1] + mA ^
        2 * mX_2
      c <- -mAY + mA * x[4] + term1 * x[3]
      d <- -mXAY_1 + (mA + mXA_1) * x[4] + (mX_2 * term1 + term2 * x[1]) * x[3]
      e <- -mXAYXA + (mA ^ 2 * mX_2 + 2 * mA * term1 *
                        x[1] + mA ^ 2 + term1 ^ 2 * x[2]) * x[4] +
        (mA * (term1 * mX_2 + term2 * x[1]) +
           term1 ^ 2 * x[1] + term1 * term2 * x[2]) * x[3]
      return(c(a, b, c, d, e))
    }
    
    sols <- nleqslv(par_initial, fn_7)$x
  } else if (type == "6_new") {
    mA <- moments[1]
    mX_2 <- moments[2]
    mXA_1 <- moments[3]
    mXA_2 <- moments[4]
    mXAY_1 <- moments[5]
    mXAYXA <- moments[6]
    
    fn_6_new <- function (x) {
      term0 <- forward_mapping_0(x[1], x[2])
      term1 <- forward_mapping_1(x[1], x[2])
      term2 <- forward_mapping_2(x[1], x[2])
      
      a <- -mXA_1 + term1 * x[1] + mA * mX_2
      b <- -mXA_2 +  term1 ^ 2 * x[2] + 2 * mA * term1 * x[1] + mA ^
        2 * mX_2
      d <- -mXAY_1 + (mA + mXA_1) * x[4] + (mX_2 * term1 + term2 * x[1]) * x[3]
      e <- -mXAYXA + (mA ^ 2 * mX_2 + 2 * mA * term1 *
                        x[1] + mA ^ 2 + term1 ^ 2 * x[2]) * x[4] +
        (mA * (term1 * mX_2 + term2 * x[1]) +
           term1 ^ 2 * x[1] + term1 * term2 * x[2]) * x[3]
      return(c(a, b, d, e))
    }
    
    sols <- nleqslv(par_initial, fn_6_new)$x
    
  }
  
  
  return(sols)
}

moments_to_f_paras_func <- function(moments, par_initial, type) {
  
  sols <- moments_to_paras_func(moments,par_initial,type)
  if (type == "1") {
    
    f_2 <- forward_mapping_2(0, sols)
    f_1 <- forward_mapping_1(0, sols)
    f_0 <- forward_mapping_0(0, sols)
    
  } else {
    f_2 <- forward_mapping_2(sols[1], sols[2])
    f_1 <- forward_mapping_1(sols[1], sols[2])
    f_0 <- forward_mapping_0(sols[1], sols[2])
    
  }
  
  f <- c(f_0, f_1, f_2)
  
  return(f)
}

moments_to_alpha_func <- function(moments, par_initial, type, Omega_multi_mean_X, Omega_multi_mean_AX) {
  
  f <-  moments_to_f_paras_func(moments, par_initial, type)
  
  alpha <-  if (f[2] > 0) {
    (Omega_multi_mean_AX - f[1] * Omega_multi_mean_X) / f[2]
  } else {
    rep(NA, p)
  }
  
  return(alpha)
}
moments_to_paras_split_sample_func <- function(type, par_initial, moments_1,  moments_2) {
  
  
  paras_1 <- moments_to_paras_func(moments_1, par_initial, type) 
  
  paras_2 <- moments_to_paras_func(moments_2, par_initial, type) 
  
  paras <- (paras_1 + paras_2)/2
  return(paras)
}
moments_to_f_paras_split_sample_func <- function(type, par_initial, moments_1,  moments_2) {
  
  
  paras_1 <- moments_to_f_paras_func(moments_1, par_initial, type) 
  
  paras_2 <- moments_to_f_paras_func(moments_2, par_initial, type) 
  
  paras <- (paras_1 + paras_2)/2
  return(paras)
}
moments_to_alpha_split_sample_func <- function(type, par_initial, moments_1, Omega_multi_mean_X_1, Omega_multi_mean_AX_1, moments_2, Omega_multi_mean_X_2, Omega_multi_mean_AX_2) {
  
  
  alpha_split_1 <- moments_to_alpha_func(moments_1, par_initial, type, Omega_multi_mean_X_1, Omega_multi_mean_AX_1) 
  
  alpha_split_2 <- moments_to_alpha_func(moments_2, par_initial, type, Omega_multi_mean_X_2, Omega_multi_mean_AX_2) 
  
  alpha <- (alpha_split_1 + alpha_split_2)/2
  return(alpha)
}

call_f_function <- function(moments, par_initial, type, f_function_name, Omega_multi_mean_X = NA, Omega_multi_mean_AX = NA, moments_2 = NA, Omega_multi_mean_X_2 = NA, Omega_multi_mean_AX_2 = NA) {
  if (f_function_name == "moments_to_paras_func") {
    return(moments_to_paras_func(moments, par_initial, type))
  } else if (f_function_name == "moments_to_f_paras_func") {
    return(moments_to_f_paras_func(moments, par_initial, type))
  } else if (f_function_name == "moments_to_alpha_func") {
    return(moments_to_alpha_func(moments, par_initial, type, Omega_multi_mean_X, Omega_multi_mean_AX))
  } else if (f_function_name == "moments_to_alpha_split_sample_func") {
    return(moments_to_alpha_split_sample_func(type, par_initial, moments, Omega_multi_mean_X, Omega_multi_mean_AX, moments_2, Omega_multi_mean_X_2, Omega_multi_mean_AX_2))
  } else if (f_function_name == "moments_to_paras_split_sample_func"){
    return(moments_to_paras_split_sample_func(type, par_initial, moments, moments_2))
  } else if (f_function_name == "moments_to_f_paras_split_sample_func"){
    return(moments_to_f_paras_split_sample_func(type, par_initial, moments, moments_2))
  } else {
    stop("Unknown function name")
  }
}
jacobian_central_difference_main <- function(row_num, var_name, moments, par_initial, type, f_function_name, Omega_multi_mean_X = NA, Omega_multi_mean_AX = NA, moments_2 = NA, Omega_multi_mean_X_2 = NA, Omega_multi_mean_AX_2 = NA, epsilon = 1e-6) {
  var_used <- get(var_name)
  n <- length(var_used)
  jacobian <- matrix(0, nrow = row_num, ncol = n)
  for (i in 1:n ){
    
    var_used_plus <- var_used
    var_used_plus[i] <-var_used_plus[i] + epsilon
    
    var_used_minus <- var_used
    var_used_minus[i] <- var_used_minus[i] - epsilon
    
    if (var_name == "moments"){
      f_plus <- call_f_function( var_used_plus, par_initial, type, f_function_name, Omega_multi_mean_X = Omega_multi_mean_X, Omega_multi_mean_AX = Omega_multi_mean_AX, moments_2 = moments_2,  Omega_multi_mean_X_2 = Omega_multi_mean_X_2, Omega_multi_mean_AX_2 = Omega_multi_mean_AX_2)  
      f_minus <- call_f_function( var_used_minus, par_initial, type,f_function_name, Omega_multi_mean_X = Omega_multi_mean_X, Omega_multi_mean_AX = Omega_multi_mean_AX, moments_2 = moments_2, Omega_multi_mean_X_2 = Omega_multi_mean_X_2, Omega_multi_mean_AX_2 = Omega_multi_mean_AX_2)  
      
    } else if (var_name == "moments_2"){
      f_plus <- call_f_function( moments, par_initial, type, f_function_name,Omega_multi_mean_X, Omega_multi_mean_AX, var_used_plus,  Omega_multi_mean_X_2, Omega_multi_mean_AX_2)  
      f_minus <- call_f_function( moments, par_initial, type,f_function_name, Omega_multi_mean_X, Omega_multi_mean_AX, var_used_minus, Omega_multi_mean_X_2, Omega_multi_mean_AX_2)  
      
    } else if (var_name == "Omega_multi_mean_X") {
      f_plus <- call_f_function( moments, par_initial, type, f_function_name,var_used_plus, Omega_multi_mean_AX, moments_2, Omega_multi_mean_X_2, Omega_multi_mean_AX_2)  
      f_minus <- call_f_function( moments, par_initial, type,f_function_name, var_used_minus, Omega_multi_mean_AX, moments_2, Omega_multi_mean_X_2, Omega_multi_mean_AX_2)  
      
    } else if (var_name == "Omega_multi_mean_X_2") {
      f_plus <- call_f_function( moments, par_initial, type, f_function_name,Omega_multi_mean_X, Omega_multi_mean_AX, moments_2,  var_used_plus, Omega_multi_mean_AX_2)  
      f_minus <- call_f_function( moments, par_initial, type,f_function_name, Omega_multi_mean_X, Omega_multi_mean_AX, moments_2,  var_used_minus, Omega_multi_mean_AX_2)  
      
    } else if (var_name == "Omega_multi_mean_AX") {
      f_plus <- call_f_function( moments, par_initial, type, f_function_name,Omega_multi_mean_X,  var_used_plus, moments_2,  Omega_multi_mean_X_2, Omega_multi_mean_AX_2)  
      f_minus <- call_f_function( moments, par_initial, type,f_function_name, Omega_multi_mean_X,  var_used_minus, moments_2,  Omega_multi_mean_X_2, Omega_multi_mean_AX_2)  
      
    } else if (var_name == "Omega_multi_mean_AX_2"){
      f_plus <- call_f_function( moments, par_initial, type, f_function_name,Omega_multi_mean_X,  Omega_multi_mean_AX, moments_2,  Omega_multi_mean_X_2, var_used_plus)  
      f_minus <- call_f_function( moments, par_initial, type,f_function_name, Omega_multi_mean_X,  Omega_multi_mean_AX, moments_2,  Omega_multi_mean_X_2, var_used_minus)  
      
    }
    
    jacobian[, i] <- (f_plus - f_minus) / (2 * epsilon)
  }
  
  
  return(jacobian)
}

jacobian_central_difference <- function(moments, par_initial, type, f_function_name, Omega_multi_mean_X = NA, Omega_multi_mean_AX = NA, moments_2 = NA, Omega_multi_mean_X_2 = NA, Omega_multi_mean_AX_2 = NA, epsilon = 1e-6) {
  if (!is.character(f_function_name) || length(f_function_name) != 1) {
    stop("f_function_name must be 'moments_to_paras_func', 'moments_to_f_paras_func' , 'moments_to_alpha_func' or 'moments_to_alpha_split_sample_func' ")
  }
  # 修正：使用 any() 函数来处理向量中的 NA 检查
  if (f_function_name == "moments_to_alpha_func" & (any(is.na(Omega_multi_mean_X)) | any(is.na(Omega_multi_mean_AX)))) {
    stop("'Omega_multi_mean_X' and 'Omega_multi_mean_AX' must be provided (non-NA) when using 'moments_to_alpha_func'")
  }
  
  if (f_function_name == "moments_to_alpha_func"){
    
    m <- length(call_f_function(moments, par_initial, type, f_function_name, Omega_multi_mean_X, Omega_multi_mean_AX))  
    
    jacobian1 <- jacobian_central_difference_main(row_num = m, var_name = "moments", moments, par_initial, type, f_function_name, Omega_multi_mean_X, Omega_multi_mean_AX)
    jacobian2 <- jacobian_central_difference_main(row_num = m, var_name = "Omega_multi_mean_X", moments, par_initial, type, f_function_name, Omega_multi_mean_X, Omega_multi_mean_AX)
    jacobian3 <- jacobian_central_difference_main(row_num = m, var_name = "Omega_multi_mean_AX", moments, par_initial, type, f_function_name, Omega_multi_mean_X, Omega_multi_mean_AX)
    
    jacobian <- cbind(jacobian1,jacobian2,jacobian3)
    return(jacobian)
    
  } else if (f_function_name == "moments_to_paras_func" | f_function_name == "moments_to_f_paras_func"){
    m <- length(call_f_function(moments, par_initial, type, f_function_name, Omega_multi_mean_X, Omega_multi_mean_AX))  
    
    jacobian <- jacobian_central_difference_main(row_num = m, var_name = "moments", moments, par_initial, type, f_function_name, Omega_multi_mean_X, Omega_multi_mean_AX)
    return(jacobian)
    
  } else if (f_function_name == "moments_to_alpha_split_sample_func") {
    m <- length(call_f_function(moments, par_initial, type, f_function_name, Omega_multi_mean_X, Omega_multi_mean_AX, moments_2, Omega_multi_mean_X_2, Omega_multi_mean_AX_2))  
    
    jacobian1 <- jacobian_central_difference_main(row_num = m, var_name = "moments", moments, par_initial, type, f_function_name, Omega_multi_mean_X, Omega_multi_mean_AX, moments_2, Omega_multi_mean_X_2, Omega_multi_mean_AX_2)
    jacobian2 <- jacobian_central_difference_main(row_num = m, var_name = "Omega_multi_mean_X", moments, par_initial, type, f_function_name, Omega_multi_mean_X, Omega_multi_mean_AX, moments_2, Omega_multi_mean_X_2, Omega_multi_mean_AX_2)
    jacobian3 <- jacobian_central_difference_main(row_num = m, var_name = "Omega_multi_mean_AX", moments, par_initial, type, f_function_name, Omega_multi_mean_X, Omega_multi_mean_AX, moments_2, Omega_multi_mean_X_2, Omega_multi_mean_AX_2)
    jacobian4 <- jacobian_central_difference_main(row_num = m, var_name = "moments_2", moments, par_initial, type, f_function_name, Omega_multi_mean_X, Omega_multi_mean_AX, moments_2, Omega_multi_mean_X_2, Omega_multi_mean_AX_2)
    jacobian5 <- jacobian_central_difference_main(row_num = m, var_name = "Omega_multi_mean_X_2", moments, par_initial, type, f_function_name, Omega_multi_mean_X, Omega_multi_mean_AX, moments_2, Omega_multi_mean_X_2, Omega_multi_mean_AX_2)
    jacobian6 <- jacobian_central_difference_main(row_num = m, var_name = "Omega_multi_mean_AX_2", moments, par_initial, type, f_function_name, Omega_multi_mean_X, Omega_multi_mean_AX, moments_2, Omega_multi_mean_X_2, Omega_multi_mean_AX_2)
    
    
    jacobian <- cbind(jacobian1,jacobian2,jacobian3,jacobian4,jacobian5,jacobian6)
    
    return(jacobian)
    
  } else if (f_function_name == "moments_to_paras_split_sample_func" | f_function_name == "moments_to_f_paras_split_sample_func") {
    m <- length(call_f_function(moments, par_initial, type, f_function_name, Omega_multi_mean_X, Omega_multi_mean_AX, moments_2, Omega_multi_mean_X_2, Omega_multi_mean_AX_2))  
    
    jacobian1 <- jacobian_central_difference_main(row_num = m, var_name = "moments", moments, par_initial, type, f_function_name, Omega_multi_mean_X, Omega_multi_mean_AX, moments_2, Omega_multi_mean_X_2, Omega_multi_mean_AX_2)
    jacobian4 <- jacobian_central_difference_main(row_num = m, var_name = "moments_2", moments, par_initial, type, f_function_name, Omega_multi_mean_X, Omega_multi_mean_AX, moments_2, Omega_multi_mean_X_2, Omega_multi_mean_AX_2)
    
    jacobian <- cbind(jacobian1,jacobian4)
    
    return(jacobian)
    
  }
  
}
# 计算一阶泰勒展开
taylor_approximation_1st <- function(moments,  par_initial, type, f_function_name, delta_x, Omega_multi_mean_X = NA, Omega_multi_mean_AX = NA, moments_2 = NA, Omega_multi_mean_X_2 = NA, Omega_multi_mean_AX_2 = NA) {
  jacobian <- jacobian_central_difference(moments, par_initial, type, f_function_name, Omega_multi_mean_X, Omega_multi_mean_AX,moments_2, Omega_multi_mean_X_2, Omega_multi_mean_AX_2)
  
  f_x <- call_f_function(moments, par_initial, type, f_function_name, Omega_multi_mean_X, Omega_multi_mean_AX, moments_2, Omega_multi_mean_X_2, Omega_multi_mean_AX_2)  
  f_x_plus_delta_x <- f_x + as.numeric(jacobian %*% delta_x)
  
  return(f_x_plus_delta_x)
}
# moments_truth_4 <- moments_truth_7[c(1:4)]
# par_truth_4 <- par_truth_7[c(1:2)]
# par_initial_4 <- par_truth_4 + rep(1e-2,2)
# A_mat <- jacobian_central_difference(moments_truth_4, par_initial_4, type = "4", f_function_name = "moments_to_alpha_func", Omega_multi_mean_X = mu_X_em[1,], Omega_multi_mean_AX = Ax_T_mean[1,])
# 
# call_f_function(moments_truth_4, par_initial_4, type = "4", f_function_name = "moments_to_alpha_func", Omega_multi_mean_X = mu_X_em[1,], Omega_multi_mean_AX = Ax_T_mean[1,])
# 
# call_f_function(moments_truth_4, par_initial_4, type = "4", f_function_name = "moments_to_alpha_split_sample_func", Omega_multi_mean_X = mu_X_em[1,], Omega_multi_mean_AX = Ax_T_mean[1,], moments_2 = moments_truth_4 + rep(0.02,4),  Omega_multi_mean_X_2 = mu_X_em[2,]  + rep(1e-6,p), Omega_multi_mean_AX_2 = Ax_T_mean[2,] + rep(1e-6,p))
# jacobian_central_difference(moments_truth_4, par_initial_4, type = "4", f_function_name = "moments_to_alpha_split_sample_func", Omega_multi_mean_X = mu_X_em[1,], Omega_multi_mean_AX = Ax_T_mean[1,], moments_2 = moments_truth_4 + rep(0.02,4),  Omega_multi_mean_X_2 = mu_X_em[2,], Omega_multi_mean_AX_2 = Ax_T_mean[2,])
# # # 
# par_initial_4 <- par_truth_4 + 0.1
# delta_x_4 <- rep(0.01581139,4)
# approx_alpha <- taylor_approximation_1st(moments_truth_4, par_initial_4, type = "4", f_function_name = "moments_to_paras_func", delta_x = delta_x_4)
# plot(approx_alpha)
# actrual_alpha <- moments_to_paras_func(moments_truth_4 + delta_x_4[c(1:4)], par_initial_4, type = "4")
# plot(actrual_alpha)
# plot(approx_alpha - actrual_alpha)
# approx_alpha
# actrual_alpha
# approx_alpha - actrual_alpha
# summary(abs(approx_alpha - actrual_alpha))
# 
# par_initial_1 <- par_truth_1 + 0.1
# delta_x_1 <- rep(0.01,1)
# approx_alpha <- taylor_approximation_1st(moments_truth_1, par_initial_1, type = "1", f_function_name = "moments_to_paras_func", delta_x = delta_x_1, Omega_multi_mean_X = Omega_multi_mean_X_truth, Omega_multi_mean_AX = Omega_multi_mean_AX_truth)
# plot(approx_alpha)
# actrual_alpha <- moments_to_paras_func(moments_truth_1 + delta_x_1[1], par_initial_1, type = "1")
# plot(actrual_alpha)
# plot(approx_alpha - actrual_alpha)
# summary(abs(approx_alpha - actrual_alpha))
# 
# par_initial_1 <- par_truth_1 + 0.1
# delta_x_1 <- rep(0.1,1 + 2*p)
# approx_alpha <- taylor_approximation_1st(moments_truth_1, par_initial_1, type = "1", f_function_name = "moments_to_alpha_func", delta_x = delta_x_1, Omega_multi_mean_X = Omega_multi_mean_X_truth, Omega_multi_mean_AX = Omega_multi_mean_AX_truth)
# plot(approx_alpha)
# actrual_alpha_1 <- moments_to_alpha_func(moments_truth_1 + delta_x_1[1], par_initial_1, type = "1", Omega_multi_mean_X = Omega_multi_mean_X_truth + delta_x_1[c(1:p)], Omega_multi_mean_AX = Omega_multi_mean_AX_truth + delta_x_1[c(1:p)])
# plot(actrual_alpha_1)
# plot(approx_alpha - actrual_alpha_1)
# summary(abs(approx_alpha - actrual_alpha_1))
# # # 
# par_initial_4 <- par_truth_4 + 0.1
# delta_x_4 <- rep(0.1,4+2*p)
# approx_alpha <- taylor_approximation_1st(moments_truth_4, par_initial_4, type = "4", f_function_name = "moments_to_alpha_func", delta_x = delta_x_4, Omega_multi_mean_X = Omega_multi_mean_X_truth, Omega_multi_mean_AX = Omega_multi_mean_AX_truth)
# plot(approx_alpha)
# actrual_alpha <- moments_to_alpha_func(moments_truth_4 + delta_x_4[c(1:4)], par_initial_4, type = "4", Omega_multi_mean_X = Omega_multi_mean_X_truth + delta_x_4[c(1:p)], Omega_multi_mean_AX = Omega_multi_mean_AX_truth + delta_x_4[c(1:p)])
# plot(actrual_alpha)
# plot(approx_alpha - actrual_alpha)
# summary(abs(approx_alpha - actrual_alpha))
# # # 
# # par_initial_4 <- par_truth_4
# # delta_x_4 <- rep(0.1,4 + 2*p + 4 + 2*p)
# # approx_alpha <- taylor_approximation_1st(moments_truth_4, par_initial_4, type = "4", f_function_name = "moments_to_alpha_split_sample_func", delta_x = delta_x_4, Omega_multi_mean_X = Omega_multi_mean_X_truth, Omega_multi_mean_AX = Omega_multi_mean_AX_truth, moments_2 = moments_truth_4,  Omega_multi_mean_X_2 = Omega_multi_mean_X_truth, Omega_multi_mean_AX_2 = Omega_multi_mean_AX_truth)
# # plot(approx_alpha)
# # actrual_alpha <- call_f_function(moments_truth_4 + delta_x_4[c(1:4)], par_initial_4, type = "4", f_function_name = "moments_to_alpha_split_sample_func", Omega_multi_mean_X = Omega_multi_mean_X_truth + delta_x_4[c(1:p)], Omega_multi_mean_AX = Omega_multi_mean_AX_truth + delta_x_4[c(1:p)], moments_2 = moments_truth_4 + delta_x_4[c(1:4)],  Omega_multi_mean_X_2 = Omega_multi_mean_X_truth + delta_x_4[c(1:p)], Omega_multi_mean_AX_2 = Omega_multi_mean_AX_truth + delta_x_4[c(1:p)])
# # plot(actrual_alpha)
# # plot(approx_alpha - actrual_alpha)
# # summary(abs(approx_alpha - actrual_alpha))
# 
# par_initial_4 <- par_truth_4
# delta_x_4 <- rep(0.1,4 + 4)
# approx_alpha <- taylor_approximation_1st(moments_truth_4, par_initial_4, type = "4", f_function_name = "moments_to_paras_split_sample_func", delta_x = delta_x_4, Omega_multi_mean_X = Omega_multi_mean_X_truth, Omega_multi_mean_AX = Omega_multi_mean_AX_truth, moments_2 = moments_truth_4,  Omega_multi_mean_X_2 = Omega_multi_mean_X_truth, Omega_multi_mean_AX_2 = Omega_multi_mean_AX_truth)
# plot(approx_alpha)
# actrual_alpha <- call_f_function(moments_truth_4 + delta_x_4[c(1:4)], par_initial_4, type = "4", f_function_name = "moments_to_paras_split_sample_func", Omega_multi_mean_X = Omega_multi_mean_X_truth + delta_x_4[c(1:p)], Omega_multi_mean_AX = Omega_multi_mean_AX_truth + delta_x_4[c(1:p)], moments_2 = moments_truth_4 + delta_x_4[c(1:4)],  Omega_multi_mean_X_2 = Omega_multi_mean_X_truth + delta_x_4[c(1:p)], Omega_multi_mean_AX_2 = Omega_multi_mean_AX_truth + delta_x_4[c(1:p)])
# plot(actrual_alpha)
# plot(approx_alpha - actrual_alpha)
# summary(abs(approx_alpha - actrual_alpha))
# 
# par_initial_4 <- par_truth_4
# delta_x_4 <- rep(0.1,4 + 4)
# approx_alpha <- taylor_approximation_1st(moments_truth_4, par_initial_4, type = "4", f_function_name = "moments_to_f_paras_split_sample_func", delta_x = delta_x_4, Omega_multi_mean_X = Omega_multi_mean_X_truth, Omega_multi_mean_AX = Omega_multi_mean_AX_truth, moments_2 = moments_truth_4,  Omega_multi_mean_X_2 = Omega_multi_mean_X_truth, Omega_multi_mean_AX_2 = Omega_multi_mean_AX_truth)
# plot(approx_alpha)
# actrual_alpha <- call_f_function(moments_truth_4 + delta_x_4[c(1:4)], par_initial_4, type = "4", f_function_name = "moments_to_f_paras_split_sample_func", Omega_multi_mean_X = Omega_multi_mean_X_truth + delta_x_4[c(1:p)], Omega_multi_mean_AX = Omega_multi_mean_AX_truth + delta_x_4[c(1:p)], moments_2 = moments_truth_4 + delta_x_4[c(1:4)],  Omega_multi_mean_X_2 = Omega_multi_mean_X_truth + delta_x_4[c(1:p)], Omega_multi_mean_AX_2 = Omega_multi_mean_AX_truth + delta_x_4[c(1:p)])
# plot(actrual_alpha)
# plot(approx_alpha - actrual_alpha)
# summary(abs(approx_alpha - actrual_alpha))
# 
# par_initial_1 <- par_truth_1
# delta_x_1 <- rep(0.1,1 + 1)
# approx_alpha <- taylor_approximation_1st(moments_truth_1, par_initial_1, type = "1", f_function_name = "moments_to_f_paras_split_sample_func", delta_x = delta_x_1, Omega_multi_mean_X = Omega_multi_mean_X_truth, Omega_multi_mean_AX = Omega_multi_mean_AX_truth, moments_2 = moments_truth_1,  Omega_multi_mean_X_2 = Omega_multi_mean_X_truth, Omega_multi_mean_AX_2 = Omega_multi_mean_AX_truth)
# plot(approx_alpha)
# # delta_x_1 <- rep(0,1 + 1)
# actrual_alpha <- call_f_function(moments_truth_1 + delta_x_1[1], par_initial_1, type = "1", f_function_name = "moments_to_f_paras_split_sample_func", Omega_multi_mean_X = Omega_multi_mean_X_truth + delta_x_1[c(1:p)], Omega_multi_mean_AX = Omega_multi_mean_AX_truth + delta_x_1[c(1:p)], moments_2 = moments_truth_1 + delta_x_1[1],  Omega_multi_mean_X_2 = Omega_multi_mean_X_truth + delta_x_1[c(1:p)], Omega_multi_mean_AX_2 = Omega_multi_mean_AX_truth + delta_x_1[c(1:p)])
# plot(actrual_alpha)
# plot(approx_alpha - actrual_alpha)
# summary(abs(approx_alpha - actrual_alpha))

# ## feasible initial ----
# # 定义参数范围
# par1_range <- seq(-2, 2, length.out = 200)  # 第一个参数的取值范围
# par2_range <- seq(0.1, 4, length.out = 200)  # 第二个参数的取值范围，确保 > 0
# 
# # 创建存储结果的矩阵
# feasible_region <-  matrix(0, nrow = length(par1_range), ncol = length(par2_range))
# alpha_error <- matrix(NA, nrow = length(par1_range), ncol = length(par2_range))
# # 遍历所有参数组合
# actrual_alpha <- list()
# for (i in 1:length(par1_range)) {
#   for (j in 1:length(par2_range)) {
#     par_initial_4 <- c(par1_range[i], par2_range[j])
#     
#     # 调用 moments_to_alpha_func 函数计算解
#     actrual_alpha[[j*(i-1)+1]] <- tryCatch({
#       moments_to_paras_func(moments_truth_4 , 
#                             par_initial_4, 
#                             type = "4"
#                             )
#     }, error = function(e) NA)  # 如果出错，返回 NA
# 
#     # 检查是否为 NA
#     if (!any(is.na(actrual_alpha[[j*(i-1)+1]]))) {
#       alpha_error[i,j] <- max(actrual_alpha[[j*(i-1)+1]] - par_truth_4, na.rm = TRUE)
#       if (alpha_error[i,j] < 0.3) {
#         feasible_region[i, j] <- 1  # 符合条件的区域标记为 1
#       }
#     }
#   }
# }
# f_error <-  lapply(actrual_alpha, function(x) {
#   tryCatch({
#     fn_4(x)
#   }, error = function(e) {
#     NA
#   })
# })
# # 筛选满足条件的点
# feasible_points <- list()
# other_solutions <- list()
# infeasible_points <- list()
# error_bound <- 1e-5
# for (i in 1:length(actrual_alpha)) {
#   alpha <- actrual_alpha[[i]]
#   error <- f_error[[i]]
#   
#   if (is.na(alpha[1]) || is.na(error[1])) {
#     infeasible_points[[i]] <- c(NA, NA)
#   } else {
#     alpha_diff <- abs(alpha - c(0, 1))
#     error_diff <- abs(error)
#     
#     if (all(alpha_diff < error_bound) && all(error_diff < error_bound)) {
#       feasible_points[[i]] <- alpha
#     } else if (all(error_diff < error_bound)) {
#       other_solutions[[i]] <- alpha
#     } else {
#       infeasible_points[[i]] <- alpha
#     }
#   }
# }
# 
# # 转换为矩阵
# feasible_points_matrix <- matrix(0, nrow = length(par1_range), ncol = length(par2_range))
# other_solutions_matrix <- matrix(0, nrow = length(par1_range), ncol = length(par2_range))
# infeasible_points_matrix <- matrix(0, nrow = length(par1_range), ncol = length(par2_range))
# 
# for (i in 1:length(par1_range)) {
#   for (j in 1:length(par2_range)) {
#     index <- j*(i-1) + 1
#     if (index <= length(feasible_points) && !is.null(feasible_points[[index]])) {
#       feasible_points_matrix[i, j] <- 1
#     }
#     if (index <= length(other_solutions) && !is.null(other_solutions[[index]])) {
#       other_solutions_matrix[i, j] <- 1
#     }
#     if (index <= length(infeasible_points) && !is.null(infeasible_points[[index]])) {
#       infeasible_points_matrix[i, j] <- 1
#     }
#   }
# }
# 
# # 可视化结果
# library(ggplot2)
# library(reshape2)
# 
# # 转换矩阵为可绘制的数据框
# feasible_data <- melt(feasible_points_matrix)
# other_solutions_data <- melt(other_solutions_matrix)
# infeasible_data <- melt(infeasible_points_matrix)
# 
# colnames(feasible_data) <- c("par1", "par2", "feasible")
# colnames(other_solutions_data) <- c("par1", "par2", "other_solutions")
# colnames(infeasible_data) <- c("par1", "par2", "infeasible")
# 
# feasible_data$par1 <- par1_range[feasible_data$par1]
# feasible_data$par2 <- par2_range[feasible_data$par2]
# other_solutions_data$par1 <- par1_range[other_solutions_data$par1]
# other_solutions_data$par2 <- par2_range[other_solutions_data$par2]
# infeasible_data$par1 <- par1_range[infeasible_data$par1]
# infeasible_data$par2 <- par2_range[infeasible_data$par2]
# 
# # 合并数据
# combined_data <- merge(feasible_data, other_solutions_data, by = c("par1", "par2"))
# combined_data <- merge(combined_data, infeasible_data, by = c("par1", "par2"))
# 
# # 绘制图表
# ggplot(combined_data, aes(x = par1, y = par2)) +
#   geom_tile(aes(fill = as.factor(feasible))) +
#   geom_tile(aes(fill = as.factor(other_solutions)), alpha = 0.5) +
#   geom_tile(aes(fill = as.factor(infeasible)), alpha = 0.5) +
#   scale_fill_manual(values = c("0" = "green", "1" = "yellow"), 
#                     labels = c("Infeasible", "Feasible", "Other Solutions")) +
#   labs(title = "Feasible Region of Initial Values",
#        x = expression(alpha^T * mu),
#        y = expression(alpha^T * Sigma * alpha),
#        fill = "Feasibility") +
#   theme_minimal()
# 
# # Visualize the results
# library(ggplot2)
# library(reshape2)
# 
# 
# # Plot the feasible region
# ggplot(feasible_data, aes(x = par1, y = par2)) +
#   geom_tile(aes(fill = as.factor(feasible))) +
#   scale_fill_manual(values = c("0" = "red", "1" = "green"), 
#                     labels = c("Infeasible", "Feasible")) +
#   labs(title = "Feasible Region of Initial Values",
#        x = expression(alpha^T * mu),
#        y = expression(alpha^T * Sigma * alpha),
#        fill = "Feasibility") +
#   theme_minimal()
# 
# 
# # Convert the matrix to a plotable data frame
# alpha_error <- abs(alpha_error)
# 
# alpha_error_data <- melt(alpha_error)
# colnames(alpha_error_data) <- c("par1", "par2", "alpha_error")
# alpha_error_data$par1 <- par1_range[alpha_error_data$par1]
# alpha_error_data$par2 <- par2_range[alpha_error_data$par2]
# 
# # Fill NA values and categorize alpha_error
# alpha_error_data$alpha_error_filled <- ifelse(is.na(alpha_error_data$alpha_error), 2, alpha_error_data$alpha_error)
# alpha_error_data$alpha_error_filled <- ifelse(alpha_error_data$alpha_error_filled > 1, 1.5, alpha_error_data$alpha_error_filled)  # Mark values >1 as 1.1
# 
# # Plot the alpha_error heatmap
# ggplot(alpha_error_data, aes(x = par1, y = par2)) +
#   geom_tile(aes(fill = alpha_error_filled)) +
#   scale_fill_gradientn(
#     colors = c( "green", "yellow", "red"),
#     values = scales::rescale(c(0, 1, 1.5, 2)),
#     name = expression(alpha_error),
#     breaks = c(0, 1, 1.5, 2),
#     labels = c("0", "1", ">1","Infeasible" ),
#     na.value = "red"  # Set NA values to red
#   ) +
#   labs(title = "Alpha Error Heatmap  Region of Initial Values",
#        x = expression(alpha^T * mu),
#        y = expression(alpha^T * Sigma * alpha),
#        fill = "Alpha Error") +
#   theme_minimal()
# 
# # # 定义中心差分法计算海森矩阵
# # central_difference_hessian <- function(moments_truth_4, epsilon = error_bound) {
# #   n <- length(moments_truth_4)
# #   m <- length(moments_to_paras_func(moments_truth_4))
# #   hessian <- array(0, dim = c(m, n, n))
# #   
# #   for (i in 1:n) {
# #     for (j in 1:n) {
# #       # 中心差分
# #       moments_plus_i <- moments_truth_4
# #       moments_plus_i[i] <- moments_plus_i[i] + epsilon
# #       moments_plus_j <- moments_plus_i
# #       moments_plus_j[j] <- moments_plus_j[j] + epsilon
# #       f_plus_plus <- moments_to_paras_func(moments_plus_j)
# #       
# #       moments_minus_i <- moments_truth_4
# #       moments_minus_i[i] <- moments_minus_i[i] - epsilon
# #       moments_minus_j <- moments_minus_i
# #       moments_minus_j[j] <- moments_minus_j[j] - epsilon
# #       f_minus_minus <- moments_to_paras_func(moments_minus_j)
# #       
# #       moments_plus_minus <- moments_plus_i
# #       moments_plus_minus[j] <- moments_plus_minus[j] - epsilon
# #       f_plus_minus <- moments_to_paras_func(moments_plus_minus)
# #       
# #       moments_minus_plus <- moments_minus_i
# #       moments_minus_plus[j] <- moments_minus_plus[j] + epsilon
# #       f_minus_plus <- moments_to_paras_func(moments_minus_plus)
# #       
# #       hessian[, i, j] <- (f_plus_plus - f_plus_minus - f_minus_plus + f_minus_minus) / (4 * epsilon ^
# #                                                                                           2)
# #     }
# #   }
# #   
# #   return(hessian)
# # }
# # 
# # # 计算二阶泰勒展开
# # taylor_approximation_2nd_order <- function(moments_truth_4, delta_x) {
# #   # 计算目标函数在 moments_truth_4 处的值
# #   f_x <- moments_to_paras_func(moments_truth_4)
# #   
# #   # 计算雅可比矩阵
# #   jacobian <- central_difference_jacobian(moments_truth_4)
# #   
# #   # 计算海森矩阵
# #   hessian <- central_difference_hessian(moments_truth_4)
# #   
# #   # 计算一阶泰勒展开
# #   first_order_term <- f_x + jacobian %*% delta_x
# #   
# #   # 计算二阶泰勒展开
# #   second_order_term <- 0
# #   for (i in 1:length(delta_x)) {
# #     for (j in 1:length(delta_x)) {
# #       second_order_term <- second_order_term + 0.5 * delta_x[i] * delta_x[j] * hessian[, i, j]
# #     }
# #   }
# #   
# #   f_x_plus_delta_x <- first_order_term + second_order_term
# #   return(f_x_plus_delta_x)
# # }

## bootstrap for variance ----
# # Function to generate bootstrap vectors
# generate_bootstrap_vectors <- function(A, Y, X, H_S, inv_sigma, n_half, B_bootstrap, p) {
#   weights_matrix_half <- matrix(rmultinom(B_bootstrap * n_half, n_half, rep(1/n_half, n_half)), nrow = B_bootstrap, byrow = TRUE)
#   vector_one_half <- rep(1, n_half)
#   
#   mA_bootstrap <- sapply(1:B_bootstrap, function(i) mean(A * weights_matrix_half[i,]))
#   mAY_bootstrap <- sapply(1:B_bootstrap, function(i) mean(A * Y * weights_matrix_half[i,]))
#   
#   mX_2_bootstrap_1 <- sapply(1:B_bootstrap, function(i) sum((vector_one_half * weights_matrix_half[i,]) * eigenMapMatMult(H_S, vector_one_half * weights_matrix_half[i,])) / (n_half * (n_half - 1)))
#   mXA_1_bootstrap_1 <- sapply(1:B_bootstrap, function(i) sum((A * weights_matrix_half[i,]) * eigenMapMatMult(H_S, vector_one_half * weights_matrix_half[i,])) / (n_half * (n_half - 1)))
#   mXA_2_bootstrap_1 <- sapply(1:B_bootstrap, function(i) sum((A * weights_matrix_half[i,]) * eigenMapMatMult(H_S, A * weights_matrix_half[i,])) / (n_half * (n_half - 1)))
#   
#   mX_2_bootstrap_2 <- sapply(1:B_bootstrap, function(i) sum((vector_one_half * (weights_matrix_half[i,] - 1)) * eigenMapMatMult(H_S, vector_one_half * (weights_matrix_half[i,] - 1))) / (n_half * (n_half - 1)))
#   mXA_1_bootstrap_2 <- sapply(1:B_bootstrap, function(i) sum((A * (weights_matrix_half[i,] - 1)) * eigenMapMatMult(H_S, vector_one_half * (weights_matrix_half[i,] - 1))) / (n_half * (n_half - 1)))
#   mXA_2_bootstrap_2 <- sapply(1:B_bootstrap, function(i) sum((A * (weights_matrix_half[i,] - 1)) * eigenMapMatMult(H_S, A * (weights_matrix_half[i,] - 1))) / (n_half * (n_half - 1)))
#   
#   t_X <- t(X)
#   Omega_multi_mean_X <- t(sapply(1:B_bootstrap, function(i) as.numeric(eigenMapMatMult(inv_sigma, eigenMapMatMult(t_X, weights_matrix_half[i,])) / n_half)))
#   Omega_multi_mean_AX <- t(sapply(1:B_bootstrap, function(i) as.numeric(eigenMapMatMult(inv_sigma, eigenMapMatMult(t_X, A * weights_matrix_half[i,])) / n_half)))
#   
#   m1_U1 <- mA_bootstrap
#   m2_U2_1 <- cbind(mX_2_bootstrap_1, mXA_1_bootstrap_1, mXA_2_bootstrap_1)
#   m2_U2_2 <- cbind(mX_2_bootstrap_2, mXA_1_bootstrap_2, mXA_2_bootstrap_2)
#   m3_U1 <- cbind(Omega_multi_mean_X, Omega_multi_mean_AX)
#   
#   list(m1_U1 = m1_U1, m2_U2_1 = m2_U2_1, m2_U2_2 = m2_U2_2, m3_U1 = m3_U1)
# }
# 
# # Function to calculate covariance matrix from bootstrap vectors
# calculate_cov_matrix <- function(vectors_a, vectors_b = NULL) {
#   if (is.null(vectors_b)) vectors_b <- vectors_a
#   
#   B_bootstrap <- length(vectors_a$m1_U1)
#   p <- ncol(vectors_a$m3_U1)/2
#   
#   cov_matrix <- matrix(NA, 4 + 2*p, 4 + 2*p)
#   
#   cov_matrix[1,1] <- cov(vectors_a$m1_U1, vectors_b$m1_U1)
#   cov_matrix[1,2:4] <- cov(vectors_a$m1_U1, vectors_b$m2_U2_1)
#   cov_matrix[2:4,1] <- cov(vectors_a$m2_U2_1, vectors_b$m1_U1)
#   
#   cov_matrix[2:4,2:4] <- cov_bootstrap_U_stats(vectors_a$m2_U2_1, vectors_a$m2_U2_2, vectors_b$m2_U2_1, vectors_b$m2_U2_2)
#   
#   cov_matrix[1,(5:(4+2*p))] <- cov(vectors_a$m1_U1, vectors_b$m3_U1)
#   cov_matrix[(5:(4+2*p)),1] <- cov(vectors_a$m3_U1, vectors_b$m1_U1)
#   
#   cov_matrix[2:4,(5:(4+2*p))] <- cov(vectors_a$m2_U2_1, vectors_b$m3_U1)
#   cov_matrix[(5:(4+2*p)),2:4] <- cov(vectors_a$m3_U1, vectors_b$m2_U2_1)
#   
#   cov_matrix[(5:(4+2*p)), (5:(4+2*p))] <- eigenMapMatMult(t(vectors_a$m3_U1), vectors_b$m3_U1) / (B_bootstrap - 1) - 
#     eigenMapMatMult(colMeans(vectors_a$m3_U1), t(colMeans(vectors_b$m3_U1))) * B_bootstrap / (B_bootstrap - 1)
#   
#   cov_matrix
# }


cov_bootstrap_U_stats <- function(A_1, A_2, B_1, B_2 ){
  
  
  term_1 <- cov(A_1,B_1)
  term_2 <- cov(A_2,B_2)
  final <- term_1 - 2 * term_2
  
  return(final)
  
}

# summary ----

glm_summary <- function(coefficients, variances, data, formula) {
  
  # 提取因变量
  response_var <- all.vars(formula)[1]
  y <- data[[response_var]]
  
  # 提取所有变量
  all_vars <- all.vars(formula)
  
  # 创建模型矩阵
  # 如果公式中包含 `.`，则使用所有变量（不包括因变量）
  if (length(all_vars) == 2 && all_vars[2] == ".") {
    # 使用所有变量（不包括因变量）
    predictors <- setdiff(names(data), response_var)
    formula <- as.formula(paste(response_var, "~", paste(predictors, collapse = " + ")))
  }
  
  X <- model.matrix(formula, data)[, -1]
  
  n <- nrow(data)
  p <- length(coefficients)
  
  # 创建一个类似glm对象的列表
  glm_like <- list(
    coefficients = coefficients,
    vcov = diag(variances),
    df.residual = n - p,
    call = match.call(),
    family = binomial(link = "logit")
  )
  class(glm_like) <- "glm"
  
  # 计算必要的统计量
  linear.predictors <- as.vector(X %*% coefficients)
  fitted.values <- 1 / (1 + exp(-linear.predictors))
  residuals <- y - fitted.values
  
  deviance <- -2 * sum(y * log(fitted.values + 1e-10) + (1 - y) * log(1 - fitted.values + 1e-10))
  null.deviance <- -2 * sum(y * log(mean(y)) + (1 - y) * log(1 - mean(y)))
  
  # 创建summary对象
  summary_output <- list()
  summary_output$call <- glm_like$call
  summary_output$family <- glm_like$family
  summary_output$deviance <- deviance
  summary_output$aic <- deviance + 2 * p
  summary_output$contrasts <- NULL
  summary_output$df.residual <- n - p
  summary_output$null.deviance <- null.deviance
  summary_output$df.null <- n - 1
  summary_output$iter <- NA
  summary_output$deviance.resid <- sign(residuals) * sqrt(-2 * (y * log(fitted.values + 1e-10) + (1 - y) * log(1 - fitted.values + 1e-10)))
  
  # 计算系数表，不加入显著性标记，并确保结果是数值型
  z_values <- coefficients / sqrt(variances)
  p_values <- 2 * (1 - pnorm(abs(z_values)))
  
  coef_table <- cbind(
    Estimate = as.numeric(coefficients),
    `Std. Error` = as.numeric(sqrt(variances)),
    `z value` = as.numeric(z_values),
    `Pr(>|z|)` = as.numeric(p_values)
  )
  
  summary_output$coefficients <- coef_table
  
  # 计算aliased
  qr_decomp <- qr(X)
  aliased <- rep(FALSE, p)
  names(aliased) <- names(coefficients)
  aliased[qr_decomp$pivot[seq(p)]] <- qr_decomp$rank < p
  summary_output$aliased <- aliased
  
  summary_output$dispersion <- 1  # 对于binomial family，dispersion总是1
  summary_output$df <- c(p, n - p, p)
  
  # 计算协方差矩阵
  summary_output$cov.unscaled <- summary_output$cov.scaled <- glm_like$vcov
  
  # 设置terms属性
  summary_output$terms <- terms(formula)
  
  class(summary_output) <- c("custom_summary", "summary.glm")
  
  return(summary_output)
}



# 自定义打印方法
print.custom_summary <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients:\n")
  
  # 手动添加显著性标记
  signif_codes <- cut(x$coefficients[, 4], breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                      labels = c("***", "**", "*", ".", " "), right = FALSE)
  
  # 将因子类型转换为字符类型
  signif_codes <- as.character(signif_codes)
  
  # 生成带显著性标记的系数表
  coef_table <- x$coefficients
  coef_table <- cbind(coef_table, Signif. = signif_codes)  # 添加显著性标记列
  
  # 格式化系数表
  coef_table_formatted <- format(coef_table, digits = 3, justify = "right")
  
  # 手动调整每一列的宽度
  coef_table_formatted[, 1] <- format(coef_table_formatted[, 1], width = 8)
  coef_table_formatted[, 2] <- format(coef_table_formatted[, 2], width = 8)
  coef_table_formatted[, 3] <- format(coef_table_formatted[, 3], width = 8)
  coef_table_formatted[, 4] <- format(coef_table_formatted[, 4], width = 8)
  coef_table_formatted[, 5] <- format(coef_table_formatted[, 5], width = 4)
  
  # 打印系数表
  print(coef_table_formatted, quote = FALSE)
  
  cat("---\n")
  cat("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1\n")
  cat("\n(Dispersion parameter for binomial family taken to be", x$dispersion, ")\n\n")
  cat("    Null deviance:", format(x$null.deviance, digits = 5), "on", x$df.null, "degrees of freedom\n")
  cat("Residual deviance:", format(x$deviance, digits = 5), "on", x$df.residual, "degrees of freedom\n")
  cat("AIC:", format(x$aic, digits = 5), "\n\n")
  cat("Number of Fisher Scoring iterations:", x$iter, "\n")
}


# 
# coefficients <- c(0.5, -0.3, 0.2)
# variances <- c(0.01, 0.02, 0.015)
# names(coefficients) <- c("(Intercept)", "V1", "V2")
# 
# # 数据集
# data <- data.frame(A_T = c(0,1,0,1,1), V1 = c(1,2,3,4,5), V2 = c(5,4,3,2,1))
# 
# # 模型公式
# formula <- A_T ~ .
# 
# # 使用自定义函数
# result <- custom_glm_summary(coefficients, variances, data, formula)
# print(result)
# result$coefficients


# generate data -----
generate_data_XA <- function(n,p,omega_11,seed,alpha,mu_x,Is_Rad){
  set.seed(seed)
  if (Is_Rad) {
    X_T <- sqrt(2 / omega_11) * matrix(2 * rbinom(n * p, size = 1, prob = 0.5) - 1 ,
                                       nrow = n,
                                       ncol = p) + mu_x
  } else  {
    X_T <- matrix(rnorm(n * p, mean = mu_x, sd = 1 / omega_11),
                  nrow = n,
                  ncol = p)
  }
  
  z_T <- eigenMapMatMult(X_T , alpha)
  A_T <- rbinom(n, 1, 1 / (1 + exp(-z_T)))
  A_T <- as.numeric(A_T)
  
  return(list(X_T = X_T,
              A_T = A_T))
}
generate_data_XAY <- function(n,p,omega_11,seed,alpha,beta,mu_x,sd_mar_linear,Is_Rad){
  set.seed(seed)
  if (Is_Rad) {
    X_T <- sqrt(2 / omega_11) * matrix(2 * rbinom(n * p, size = 1, prob = 0.5) - 1 ,
                                       nrow = n,
                                       ncol = p) + mu_x
  } else  {
    X_T <- matrix(rnorm(n * p, mean = mu_x, sd = 1 / omega_11),
                  nrow = n,
                  ncol = p)
  }
  
  z_T <- eigenMapMatMult(X_T , alpha)
  A_T <- rbinom(n, 1, 1 / (1 + exp(-z_T)))
  A_T <- as.numeric(A_T)
  
  epsilon_1 <- rnorm(n, 0, sd_mar_linear)
  Y_T <- eigenMapMatMult(X_T , beta) + epsilon_1
  Y_T <- Y_T * A_T
  
  return(list(X_T = X_T,
              A_T = A_T,
              Y_T = Y_T))
}
# unknown sigma -----
glm_mom_unknown_sigma <- function(X_T,A_T,B_bootstrap,par_initial){
  n <- nrow(X_T)
  p <- ncol(X_T)
  

  
  n_half <- n / 2
  
  X_1 <- X_T[1:n_half, ]
  X_2 <- X_T[(n_half + 1):n, ]
  A_1 <- A_T[1:n_half]
  A_2 <- A_T[(n_half + 1):n ]
  
  par_initial_1 <- par_initial[1]
  par_initial_4 <- par_initial
  ## momoents ----
  sigma_split_a <- eigenMapMatMult(t(X_1), X_1) / n_half
  chol_sigma_split_a <- chol(sigma_split_a)
  inv_sigma_split_a <- chol2inv(chol_sigma_split_a) * (n_half -  p - 1) / n_half
  
  sigma_split_b <- eigenMapMatMult(t(X_2), X_2) / n_half
  chol_sigma_split_b <- chol(sigma_split_b)
  inv_sigma_split_b <- chol2inv(chol_sigma_split_b) * (n_half -  p - 1) / n_half
  
  
  H_S_a <- eigenMapMatMult(X_2, eigenMapMatMult(inv_sigma_split_a, t(X_2)))
  H_S_a <- H_S_a - diag(diag(H_S_a))
  
  H_S_b <- eigenMapMatMult(X_1, eigenMapMatMult(inv_sigma_split_b, t(X_1)))
  H_S_b <- H_S_b - diag(diag(H_S_b))
  
  vector_one_half <- rep(1, n_half)
  
  mA_em_a <- mean(A_2)
  mX_2_em_a <- as.numeric(sum(vector_one_half * eigenMapMatMult(H_S_a, vector_one_half)) / (n_half * (n_half - 1)))
  mXA_1_em_a <- as.numeric(sum(A_2 * eigenMapMatMult(H_S_a, vector_one_half)) / (n_half * (n_half - 1)))
  mXA_2_em_a <- as.numeric(sum(A_2 * eigenMapMatMult(H_S_a, A_2)) / (n_half * (n_half - 1)))
  
  moments_em_a_4 <-  c(mA_em_a, mX_2_em_a, mXA_1_em_a, mXA_2_em_a)
  
  mA_em_b <- mean(A_1)
  mX_2_em_b <- as.numeric(sum(vector_one_half * eigenMapMatMult(H_S_b, vector_one_half)) / (n_half * (n_half - 1)))
  mXA_1_em_b <- as.numeric(sum(A_1 * eigenMapMatMult(H_S_b, vector_one_half)) / (n_half * (n_half - 1)))
  mXA_2_em_b <- as.numeric(sum(A_1 * eigenMapMatMult(H_S_b, A_1)) / (n_half * (n_half - 1)))
  
  moments_em_b_4 <-  c(mA_em_b, mX_2_em_b, mXA_1_em_b, mXA_2_em_b)
  
  moments_em_S_1 <- (mXA_2_em_b +  mXA_2_em_a) / 2
  moments_em_S_4 <- (moments_em_b_4 +  moments_em_a_4) / 2
  
  Omega_multi_mean_AX_em_a <- colMeans(X_2 * A_2)
  Omega_multi_mean_AX_em_a <- eigenMapMatMult(inv_sigma_split_a, Omega_multi_mean_AX_em_a)
  
  
  Omega_multi_mean_AX_em_b <- colMeans(X_1 * A_1)
  Omega_multi_mean_AX_em_b <- eigenMapMatMult(inv_sigma_split_b, Omega_multi_mean_AX_em_b)
  
  Omega_multi_mean_X_em_a <- eigenMapMatMult(inv_sigma_split_a, colMeans(X_2))
  Omega_multi_mean_X_em_b <- eigenMapMatMult(inv_sigma_split_b, colMeans(X_1))
  
  ## paras ----
  sols_em_a_4 <- moments_to_paras_func(moments_em_a_4, par_initial_4, type = "4")
  sols_em_a_1 <- inverse_expectation(mXA_2_em_a)
  
  
  f_2_em_a_4 <- forward_mapping_2(sols_em_a_4[1], sols_em_a_4[2])
  f_1_em_a_4 <- forward_mapping_1(sols_em_a_4[1], sols_em_a_4[2])
  f_0_em_a_4 <- forward_mapping_0(sols_em_a_4[1], sols_em_a_4[2])
  
  f_2_em_a_1 <- forward_mapping_2(0, sols_em_a_1)
  f_1_em_a_1 <- forward_mapping_1(0, sols_em_a_1)
  f_0_em_a_1 <- forward_mapping_0(0, sols_em_a_1)
  
  f_em_a_4 <- c(f_0_em_a_4, f_1_em_a_4, f_2_em_a_4)
  f_em_a_1 <- c(f_0_em_a_1, f_1_em_a_1, f_2_em_a_1)
  
  sols_em_b_4 <- moments_to_paras_func(moments_em_b_4, par_initial_4, type = "4")
  sols_em_b_1 <- inverse_expectation(mXA_2_em_b)
  
  f_2_em_b_4 <- forward_mapping_2(sols_em_b_4[1], sols_em_b_4[2])
  f_1_em_b_4 <- forward_mapping_1(sols_em_b_4[1], sols_em_b_4[2])
  f_0_em_b_4 <- forward_mapping_0(sols_em_b_4[1], sols_em_b_4[2])
  
  f_2_em_b_1 <- forward_mapping_2(0, sols_em_b_1)
  f_1_em_b_1 <- forward_mapping_1(0, sols_em_b_1)
  f_0_em_b_1 <- forward_mapping_0(0, sols_em_b_1)
  
  f_em_b_4 <- c(f_0_em_b_4, f_1_em_b_4, f_2_em_b_4)
  f_em_b_1 <- c(f_0_em_b_1, f_1_em_b_1, f_2_em_b_1)
  ## alphas ----
  
  sols_em_S_1 <- (sols_em_b_1 + sols_em_a_1) / 2
  sols_em_S_4 <- (sols_em_b_4 + sols_em_a_4) / 2
  
  f_em_S_4 <- (f_em_b_4 + f_em_a_4) / 2
  f_em_S_1 <- (f_em_b_1 + f_em_a_1) / 2
  
  Omega_multi_mean_X_em_S <- (Omega_multi_mean_X_em_a + Omega_multi_mean_X_em_b) / 2
  Omega_multi_mean_AX_em_S <- (Omega_multi_mean_AX_em_a + Omega_multi_mean_AX_em_b) / 2
  
  alpha_em_S_1_a <-  if (f_em_a_1[2] > 0) {
    as.numeric(Omega_multi_mean_AX_em_a / f_em_a_1[2])
  } else {
    rep(NA_real_, p)
  }
  alpha_em_S_1_b <-  if (f_em_b_1[2] > 0) {
    as.numeric(Omega_multi_mean_AX_em_b / f_em_b_1[2])
  } else {
    rep(NA_real_, p)
  }
  
  alpha_em_S_1 <-  (alpha_em_S_1_a + alpha_em_S_1_b) / 2
  
  alpha_em_S_4_a <-  if (f_em_a_4[2] > 0) {
    as.numeric((
      Omega_multi_mean_AX_em_a - f_em_a_4[1] *  Omega_multi_mean_X_em_a
    ) / f_em_a_4[2])
  } else {
    rep(NA_real_, p)
  }
  alpha_em_S_4_b <-  if (f_em_b_4[2] > 0) {
    as.numeric((
      Omega_multi_mean_AX_em_b - f_em_b_4[1] *  Omega_multi_mean_X_em_b
    ) / f_em_b_4[2])
  } else {
    rep(NA_real_, p)
  }
  
  alpha_em_S_4 <-  (alpha_em_S_4_a + alpha_em_S_4_b) / 2
  
  
  ## jacobian ----
  jacobian_em_S_paras_1 <- jacobian_central_difference(moments_em_S_1,
                                                       par_initial_1,
                                                       type = "1",
                                                       f_function_name =  "moments_to_paras_func")
  jacobian_em_S_paras_4 <- jacobian_central_difference(moments_em_S_4,
                                                       par_initial_4,
                                                       type = "4",
                                                       f_function_name = "moments_to_paras_func")
  
  jacobian_em_S_f_paras_1 <- jacobian_central_difference(moments_em_S_1,
                                                         par_initial_1,
                                                         type = "1",
                                                         f_function_name = "moments_to_f_paras_func")
  jacobian_em_S_f_paras_4 <- jacobian_central_difference(moments_em_S_4,
                                                         par_initial_4,
                                                         type = "4",
                                                         f_function_name = "moments_to_f_paras_func")
  
  jacobian_em_S_alpha_1 <- jacobian_central_difference(
    moments_em_S_1,
    par_initial_1,
    type = "1",
    f_function_name = "moments_to_alpha_func",
    Omega_multi_mean_X = Omega_multi_mean_X_em_S,
    Omega_multi_mean_AX = Omega_multi_mean_AX_em_S
  )
  jacobian_em_S_alpha_4 <- jacobian_central_difference(
    moments_em_S_4,
    par_initial_4,
    type = "4",
    f_function_name = "moments_to_alpha_func",
    Omega_multi_mean_X = Omega_multi_mean_X_em_S,
    Omega_multi_mean_AX = Omega_multi_mean_AX_em_S
  )
  
  
  jacobian_em_S_alpha_split_sample_1 <- cbind(jacobian_em_S_alpha_1/2,jacobian_em_S_alpha_1/2)
  jacobian_em_S_alpha_split_sample_4 <- cbind(jacobian_em_S_alpha_4/2,jacobian_em_S_alpha_4/2)
  
  jacobian_em_S_paras_split_sample_1 <- cbind(jacobian_em_S_paras_1/2,jacobian_em_S_paras_1/2)
  jacobian_em_S_paras_split_sample_4 <- cbind(jacobian_em_S_paras_4/2,jacobian_em_S_paras_4/2)
  
  jacobian_em_S_f_paras_split_sample_1 <- cbind(jacobian_em_S_f_paras_1/2,jacobian_em_S_f_paras_1/2)
  jacobian_em_S_f_paras_split_sample_4 <- cbind(jacobian_em_S_f_paras_4/2,jacobian_em_S_f_paras_4/2)
  
  ## boostrap for split sample ------
  
  num_bootstrap <- n_half
  # 初始化权重矩阵
  weights_matrix_half <- matrix(0, nrow = B_bootstrap, ncol = num_bootstrap)
  # 生成多项式分布的样本
  for (i in 1:B_bootstrap) {
    weights_matrix_half[i, ] <- rmultinom(1, num_bootstrap, rep(1 / num_bootstrap, num_bootstrap))
  }
  cov_em_for_split_sample_4 <- matrix(NA_real_, 4 + 2 *
                                        p + 4 + 2 * p, 4 + 2 * p + 4 + 2 * p)
  cov_em_for_split_sample_4_a <- cov_em_for_split_sample_4_b <- cov_em_for_split_sample_4_ab <-  matrix(NA_real_, 4 + 2 *
                                                                                                          p , 4 + 2 * p)
  vector_one_half <- rep(1, n_half)
  
  ##### a-----
  mA_bootstrap_split_sample_a <- sapply(c(1:B_bootstrap), function(i)
    mean(A_2 * weights_matrix_half[i, ]))

  
  mX_2_bootstrap_split_sample_a_1 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (vector_one_half * weights_matrix_half[i, ]) * eigenMapMatMult(H_S_a, vector_one_half * weights_matrix_half[i, ])
    ) / (
      n_half * (n_half - 1)
    ))))
  mXA_1_bootstrap_split_sample_a_1 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_2 * weights_matrix_half[i, ]) * eigenMapMatMult(H_S_a, vector_one_half * weights_matrix_half[i, ])
    ) / (
      n_half * (n_half - 1)
    ))))
  mXA_2_bootstrap_split_sample_a_1 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_2 * weights_matrix_half[i, ]) * eigenMapMatMult(H_S_a, A_2 * weights_matrix_half[i, ])
    ) / (
      n_half * (n_half - 1)
    ))))

  
  mX_2_bootstrap_split_sample_a_2 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (vector_one_half * (weights_matrix_half[i, ] - 1)) * eigenMapMatMult(H_S_a, vector_one_half * (weights_matrix_half[i, ] -
                                                                                                       1))
    ) / (
      n_half * (n_half - 1)
    ))))
  mXA_1_bootstrap_split_sample_a_2 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_2 * (weights_matrix_half[i, ] - 1)) * eigenMapMatMult(H_S_a, vector_one_half * (weights_matrix_half[i, ] -
                                                                                           1))
    ) / (
      n_half * (n_half - 1)
    ))))
  mXA_2_bootstrap_split_sample_a_2 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_2 * (weights_matrix_half[i, ] - 1)) * eigenMapMatMult(H_S_a, A_2 * (weights_matrix_half[i, ] -
                                                                               1))
    ) / (
      n_half * (n_half - 1)
    ))))
 
  
  t_X_2 <- t(X_2)
  Omega_multi_mean_X_split_sample_a <- t(sapply(c(1:B_bootstrap), function(i)
    (
      as.numeric(
        eigenMapMatMult(
          inv_sigma_split_a,
          eigenMapMatMult(t_X_2, weights_matrix_half[i, ])
        ) / n_half
      )
    )))
  Omega_multi_mean_AX_split_sample_a <- t(sapply(c(1:B_bootstrap), function(i)
    (
      as.numeric(
        eigenMapMatMult(
          inv_sigma_split_a,
          eigenMapMatMult(t_X_2, A_2 * weights_matrix_half[i, ])
        ) / n_half
      )
    )))
  
  
  m1_U1_a <- mA_bootstrap_split_sample_a
  m2_U2_a_1 <- cbind(
    mX_2_bootstrap_split_sample_a_1,
    mXA_1_bootstrap_split_sample_a_1,
    mXA_2_bootstrap_split_sample_a_1
  )
  m2_U2_a_2 <- cbind(
    mX_2_bootstrap_split_sample_a_2,
    mXA_1_bootstrap_split_sample_a_2,
    mXA_2_bootstrap_split_sample_a_2
  )
  m3_U1_a <- cbind(Omega_multi_mean_X_split_sample_a,
                   Omega_multi_mean_AX_split_sample_a)
  
  cov_em_for_split_sample_4_a[1, 1] <- var(m1_U1_a)
  cov_em_for_split_sample_4_a[1, c(2:4)] <- cov(m1_U1_a, m2_U2_a_1)
  cov_em_for_split_sample_4_a[c(2:4), 1] <- cov(m2_U2_a_1, m1_U1_a)
  
  cov_em_for_split_sample_4_a[c(2:4), c(2:4)] <- cov_bootstrap_U_stats(m2_U2_a_1, m2_U2_a_2, m2_U2_a_1, m2_U2_a_2)
  
  cov_em_for_split_sample_4_a[1, c((4 + 1):(4 + 2 *
                                              p))] <- cov(m1_U1_a, m3_U1_a)
  cov_em_for_split_sample_4_a[c((4 + 1):(4 + 2 *
                                           p)), 1] <- cov(m3_U1_a, m1_U1_a)
  
  cov_em_for_split_sample_4_a[c(2:4), c((4 + 1):(4 + 2 *
                                                   p))] <- cov(m2_U2_a_1, m3_U1_a)
  cov_em_for_split_sample_4_a[c((4 + 1):(4 + 2 *
                                           p)), c(2:4)] <- cov(m3_U1_a, m2_U2_a_1)
  
  cov_em_for_split_sample_4_a[c((4 + 1):(4 + 2 *
                                           p)), c((4 + 1):(4 + 2 * p))] <- eigenMapMatMult(t(m3_U1_a), m3_U1_a) / (B_bootstrap - 1) - eigenMapMatMult(colMeans(m3_U1_a), t(colMeans(m3_U1_a))) *
    B_bootstrap / (B_bootstrap - 1)
  
  #### b ------
  mA_bootstrap_split_sample_b <- sapply(c(1:B_bootstrap), function(i)
    mean(A_1 * weights_matrix_half[i, ]))
 
  
  mX_2_bootstrap_split_sample_b_1 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (vector_one_half * weights_matrix_half[i, ]) * eigenMapMatMult(H_S_b, vector_one_half * weights_matrix_half[i, ])
    ) / (
      n_half * (n_half - 1)
    ))))
  mXA_1_bootstrap_split_sample_b_1 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_1 * weights_matrix_half[i, ]) * eigenMapMatMult(H_S_b, vector_one_half * weights_matrix_half[i, ])
    ) / (
      n_half * (n_half - 1)
    ))))
  mXA_1_bootstrap_split_sample_b_1 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_1 * weights_matrix_half[i, ]) * eigenMapMatMult(H_S_b, A_1 * weights_matrix_half[i, ])
    ) / (
      n_half * (n_half - 1)
    ))))
 
  mX_2_bootstrap_split_sample_b_2 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (vector_one_half * (weights_matrix_half[i, ] - 1)) * eigenMapMatMult(H_S_b, vector_one_half * (weights_matrix_half[i, ] -
                                                                                                       1))
    ) / (
      n_half * (n_half - 1)
    ))))
  mXA_1_bootstrap_split_sample_b_2 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_1 * (weights_matrix_half[i, ] - 1)) * eigenMapMatMult(H_S_b, vector_one_half * (weights_matrix_half[i, ] -
                                                                                           1))
    ) / (
      n_half * (n_half - 1)
    ))))
  mXA_1_bootstrap_split_sample_b_2 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_1 * (weights_matrix_half[i, ] - 1)) * eigenMapMatMult(H_S_b, A_1 * (weights_matrix_half[i, ] -
                                                                               1))
    ) / (
      n_half * (n_half - 1)
    ))))
 
  
  t_X_1 <- t(X_1)
  Omega_multi_mean_X_split_sample_b <- t(sapply(c(1:B_bootstrap), function(i)
    (
      as.numeric(
        eigenMapMatMult(
          inv_sigma_split_b,
          eigenMapMatMult(t_X_1, weights_matrix_half[i, ])
        ) / n_half
      )
    )))
  Omega_multi_mean_AX_split_sample_b <- t(sapply(c(1:B_bootstrap), function(i)
    (
      as.numeric(
        eigenMapMatMult(
          inv_sigma_split_b,
          eigenMapMatMult(t_X_1, A_1 * weights_matrix_half[i, ])
        ) / n_half
      )
    )))
  
  m1_U1_b <- mA_bootstrap_split_sample_b
  m2_U2_b_1 <- cbind(
    mX_2_bootstrap_split_sample_b_1,
    mXA_1_bootstrap_split_sample_b_1,
    mXA_1_bootstrap_split_sample_b_1
  )
  m2_U2_b_2 <- cbind(
    mX_2_bootstrap_split_sample_b_2,
    mXA_1_bootstrap_split_sample_b_2,
    mXA_1_bootstrap_split_sample_b_2
  )
  m3_U1_b <- cbind(Omega_multi_mean_X_split_sample_b,
                   Omega_multi_mean_AX_split_sample_b)
  
  cov_em_for_split_sample_4_b[1, 1] <- var(m1_U1_b)
  cov_em_for_split_sample_4_b[1, c(2:4)] <- cov(m1_U1_b, m2_U2_b_1)
  cov_em_for_split_sample_4_b[c(2:4), 1] <- cov(m2_U2_b_1, m1_U1_b)
  
  cov_em_for_split_sample_4_b[c(2:4), c(2:4)] <- cov_bootstrap_U_stats(m2_U2_b_1, m2_U2_b_2, m2_U2_b_1, m2_U2_b_2)
  
  cov_em_for_split_sample_4_b[1, c((4 + 1):(4 + 2 *
                                              p))] <- cov(m1_U1_b, m3_U1_b)
  cov_em_for_split_sample_4_b[c((4 + 1):(4 + 2 *
                                           p)), 1] <- cov(m3_U1_b, m1_U1_b)
  
  cov_em_for_split_sample_4_b[c(2:4), c((4 + 1):(4 + 2 *
                                                   p))] <- cov(m2_U2_b_1, m3_U1_b)
  cov_em_for_split_sample_4_b[c((4 + 1):(4 + 2 *
                                           p)), c(2:4)] <- cov(m3_U1_b, m2_U2_b_1)
  
  cov_em_for_split_sample_4_b[c((4 + 1):(4 + 2 *
                                           p)), c((4 + 1):(4 + 2 * p))] <- eigenMapMatMult(t(m3_U1_b), m3_U1_b) / (B_bootstrap - 1) - eigenMapMatMult(colMeans(m3_U1_b), t(colMeans(m3_U1_b))) *
    B_bootstrap / (B_bootstrap - 1)
  
  
  #### ab-----
  
  cov_em_for_split_sample_4_ab[1, 1] <- cov(m1_U1_a, m1_U1_b)
  cov_em_for_split_sample_4_ab[1, c(2:4)] <- cov(m1_U1_a, m2_U2_b_1)
  cov_em_for_split_sample_4_ab[c(2:4), 1] <- cov(m2_U2_b_1, m1_U1_a)
  
  cov_em_for_split_sample_4_ab[c(2:4), c(2:4)] <- cov_bootstrap_U_stats(m2_U2_a_1, m2_U2_a_2, m2_U2_b_1, m2_U2_b_2)
  
  cov_em_for_split_sample_4_ab[1, c((4 + 1):(4 + 2 *
                                               p))] <- cov(m1_U1_a, m3_U1_b)
  cov_em_for_split_sample_4_ab[c((4 + 1):(4 + 2 *
                                            p)), 1] <- cov(m3_U1_b, m1_U1_a)
  
  cov_em_for_split_sample_4_ab[c(2:4), c((4 + 1):(4 + 2 *
                                                    p))] <- cov(m2_U2_a_1, m3_U1_b)
  cov_em_for_split_sample_4_ab[c((4 + 1):(4 + 2 *
                                            p)), c(2:4)] <- cov(m3_U1_b, m2_U2_a_1)
  
  cov_em_for_split_sample_4_ab[c((4 + 1):(4 + 2 *
                                            p)), c((4 + 1):(4 + 2 * p))] <- eigenMapMatMult(t(m3_U1_a), m3_U1_b) / (B_bootstrap - 1) - eigenMapMatMult(colMeans(m3_U1_a), t(colMeans(m3_U1_b))) *
    B_bootstrap / (B_bootstrap - 1)
  
  
  #### summary -----
  cov_em_for_split_sample_4 <- rbind(
    cbind(
      cov_em_for_split_sample_4_a,
      cov_em_for_split_sample_4_ab
    ),
    cbind(
      t(cov_em_for_split_sample_4_ab),
      cov_em_for_split_sample_4_b
    )
  )
  cov_alpha_split_sample_4 <- eigenMapMatMult(
    eigenMapMatMult(
      jacobian_em_S_alpha_split_sample_4,
      cov_em_for_split_sample_4
    ),
    t(jacobian_em_S_alpha_split_sample_4)
  )
  
  
  var_moments_split_sample_4 <- diag(cov_em_for_split_sample_4)
  var_alpha_split_sample_4 <- diag(cov_alpha_split_sample_4)
  cov_em_for_paras_split_sample_4 <- rbind(
    cbind(
      cov_em_for_split_sample_4_a[c((1):(4)), c((1):(4))],
      cov_em_for_split_sample_4_ab[c((1):(4)), c((1):(4))]
    ),
    cbind(
      t(cov_em_for_split_sample_4_ab[c((1):(4)), c((1):(4))]),
      cov_em_for_split_sample_4_b[c((1):(4)), c((1):(4))]
    )
  )
  
  var_paras_split_sample_4 <- diag(eigenMapMatMult(
    eigenMapMatMult(
      jacobian_em_S_paras_split_sample_4,
      cov_em_for_paras_split_sample_4
    ),
    t(jacobian_em_S_paras_split_sample_4)
  ))
  var_f_paras_split_sample_4 <- diag(eigenMapMatMult(
    eigenMapMatMult(
      jacobian_em_S_f_paras_split_sample_4,
      cov_em_for_paras_split_sample_4
    ),
    t(jacobian_em_S_f_paras_split_sample_4)
  ))
  
  
  cov_em_for_split_sample_1 <- rbind(
    cbind(
      cov_em_for_split_sample_4_a[c((4):(4 + 2 * p)), c((4):(4 + 2 * p))],
      cov_em_for_split_sample_4_ab[c((4):(4 + 2 * p)), c((4):(4 + 2 * p))]
    ),
    cbind(
      t(cov_em_for_split_sample_4_ab[c((4):(4 + 2 * p)), c((4):(4 + 2 * p))]),
      cov_em_for_split_sample_4_b[c((4):(4 + 2 * p)), c((4):(4 + 2 * p))]
    )
  )
  
  var_alpha_split_sample_1 <- diag(eigenMapMatMult(
    eigenMapMatMult(
      jacobian_em_S_alpha_split_sample_1,
      cov_em_for_split_sample_1
    ),
    t(jacobian_em_S_alpha_split_sample_1)
  ))
  
  cov_em_for_paras_split_sample_1 <- rbind(
    cbind(
      cov_em_for_split_sample_4_a[4, 4],
      cov_em_for_split_sample_4_ab[4, 4]
    ),
    cbind(
      t(cov_em_for_split_sample_4_ab[4, 4]),
      cov_em_for_split_sample_4_b[4, 4]
    )
  )
  
  var_paras_split_sample_1 <- diag(eigenMapMatMult(
    eigenMapMatMult(
      jacobian_em_S_paras_split_sample_1,
      cov_em_for_paras_split_sample_1
    ),
    t(jacobian_em_S_paras_split_sample_1)
  ))
  var_f_paras_split_sample_1 <- diag(eigenMapMatMult(
    eigenMapMatMult(
      jacobian_em_S_f_paras_split_sample_1,
      cov_em_for_paras_split_sample_1
    ),
    t(jacobian_em_S_f_paras_split_sample_1)
  ))
  ### return ----
  return(list(
    sols_em_S_1 = sols_em_S_1,
    sols_em_S_4 = sols_em_S_4,
    
    f_em_S_1 = f_em_S_1,
    f_em_S_4 = f_em_S_4,
    
    alpha_em_S_1 = alpha_em_S_1,
    alpha_em_S_4 = alpha_em_S_4,
    
    moments_em_S_4 = moments_em_S_4,
    moments_em_S_1 = moments_em_S_1,
    
    var_moments_split_sample_4 = var_moments_split_sample_4,
    
    var_alpha_split_sample_4 = var_alpha_split_sample_4,
    var_paras_split_sample_4 = var_paras_split_sample_4,
    var_f_paras_split_sample_4 = var_f_paras_split_sample_4,
    
    var_alpha_split_sample_1 = var_alpha_split_sample_1,
    var_paras_split_sample_1 = var_paras_split_sample_1,
    var_f_paras_split_sample_1 = var_f_paras_split_sample_1
    
  ))
}
glm_mom_unknown_sigma_mar <- function(X_T,A_T,Y_T,B_bootstrap,par_initial_6){
  n <- nrow(X_T)
  p <- ncol(X_T)
  
  
  
  n_half <- n / 2
  
  X_1 <- X_T[1:n_half, ]
  X_2 <- X_T[(n_half + 1):n, ]
  A_1 <- A_T[1:n_half]
  A_2 <- A_T[(n_half + 1):n ]
  Y_1 <- Y_T[1:n_half]
  Y_2 <- Y_T[(n_half + 1):n ]
  
  par_initial_1 <- par_initial_6[1]
  par_initial_4 <- par_initial_6[1:2]
  
  ## momoents ----
  sigma_split_a <- eigenMapMatMult(t(X_1), X_1) / n_half
  chol_sigma_split_a <- chol(sigma_split_a)
  inv_sigma_split_a <- chol2inv(chol_sigma_split_a) * (n_half -  p - 1) / n_half
  
  sigma_split_b <- eigenMapMatMult(t(X_2), X_2) / n_half
  chol_sigma_split_b <- chol(sigma_split_b)
  inv_sigma_split_b <- chol2inv(chol_sigma_split_b) * (n_half -  p - 1) / n_half
  
  
  H_S_a <- eigenMapMatMult(X_2, eigenMapMatMult(inv_sigma_split_a, t(X_2)))
  H_S_a <- H_S_a - diag(diag(H_S_a))
  
  H_S_b <- eigenMapMatMult(X_1, eigenMapMatMult(inv_sigma_split_b, t(X_1)))
  H_S_b <- H_S_b - diag(diag(H_S_b))
  
  vector_one_half <- rep(1, n_half)
  
  mA_em_a <- mean(A_2)
  mX_2_em_a <- as.numeric(sum(vector_one_half * eigenMapMatMult(H_S_a, vector_one_half)) / (n_half * (n_half - 1)))
  mXA_1_em_a <- as.numeric(sum(A_2 * eigenMapMatMult(H_S_a, vector_one_half)) / (n_half * (n_half - 1)))
  mXA_2_em_a <- as.numeric(sum(A_2 * eigenMapMatMult(H_S_a, A_2)) / (n_half * (n_half - 1)))
  mAY_em_a <-  as.numeric(mean(A_2 * Y_2))
  mXAY_1_em_a <-  as.numeric(sum(A_2 * Y_2 * eigenMapMatMult(H_S_a, vector_one_half)) / (n_half * (n_half - 1)))
  
  moments_em_a_6 <-  c(mA_em_a,
                       mX_2_em_a,
                       mXA_1_em_a,
                       mXA_2_em_a,
                       mAY_em_a,
                       mXAY_1_em_a)
  moments_em_a_4 <-  c(mA_em_a, mX_2_em_a, mXA_1_em_a, mXA_2_em_a)
  
  mA_em_b <- mean(A_1)
  mX_2_em_b <- as.numeric(sum(vector_one_half * eigenMapMatMult(H_S_b, vector_one_half)) / (n_half * (n_half - 1)))
  mXA_1_em_b <- as.numeric(sum(A_1 * eigenMapMatMult(H_S_b, vector_one_half)) / (n_half * (n_half - 1)))
  mXA_2_em_b <- as.numeric(sum(A_1 * eigenMapMatMult(H_S_b, A_1)) / (n_half * (n_half - 1)))
  mAY_em_b <-  as.numeric(mean(A_1 * Y_1))
  mXAY_1_em_b <-  as.numeric(sum(A_1 * Y_1 * eigenMapMatMult(H_S_b, vector_one_half)) / (n_half * (n_half - 1)))
  
  moments_em_b_6 <-  c(mA_em_b,
                       mX_2_em_b,
                       mXA_1_em_b,
                       mXA_2_em_b,
                       mAY_em_b,
                       mXAY_1_em_b)
  moments_em_b_4 <-  c(mA_em_b, mX_2_em_b, mXA_1_em_b, mXA_2_em_b)
  
  moments_em_S_1 <- (mXA_2_em_b +  mXA_2_em_a) / 2
  moments_em_S_4 <- (moments_em_b_4 +  moments_em_a_4) / 2
  moments_em_S_6 <-  (moments_em_b_6 +  moments_em_a_6) / 2
  
  Omega_multi_mean_AX_em_a <- colMeans(X_2 * A_2)
  Omega_multi_mean_AX_em_a <- eigenMapMatMult(inv_sigma_split_a, Omega_multi_mean_AX_em_a)
  
  
  Omega_multi_mean_AX_em_b <- colMeans(X_1 * A_1)
  Omega_multi_mean_AX_em_b <- eigenMapMatMult(inv_sigma_split_b, Omega_multi_mean_AX_em_b)
  
  Omega_multi_mean_X_em_a <- eigenMapMatMult(inv_sigma_split_a, colMeans(X_2))
  Omega_multi_mean_X_em_b <- eigenMapMatMult(inv_sigma_split_b, colMeans(X_1))
  
  ## paras ----
  sols_em_a_6 <- moments_to_paras_func(moments_em_a_6, par_initial_6, type = "6")
  sols_em_a_4 <- moments_to_paras_func(moments_em_a_4, par_initial_4, type = "4")
  sols_em_a_1 <- inverse_expectation(mXA_2_em_a)
  
  f_2_em_a_6 <- forward_mapping_2(sols_em_a_6[1], sols_em_a_6[2])
  f_1_em_a_6 <- forward_mapping_1(sols_em_a_6[1], sols_em_a_6[2])
  f_0_em_a_6 <- forward_mapping_0(sols_em_a_6[1], sols_em_a_6[2])
  
  f_2_em_a_4 <- forward_mapping_2(sols_em_a_4[1], sols_em_a_4[2])
  f_1_em_a_4 <- forward_mapping_1(sols_em_a_4[1], sols_em_a_4[2])
  f_0_em_a_4 <- forward_mapping_0(sols_em_a_4[1], sols_em_a_4[2])
  
  f_2_em_a_1 <- forward_mapping_2(0, sols_em_a_1)
  f_1_em_a_1 <- forward_mapping_1(0, sols_em_a_1)
  f_0_em_a_1 <- forward_mapping_0(0, sols_em_a_1)
  
  f_em_a_4 <- c(f_0_em_a_4, f_1_em_a_4, f_2_em_a_4)
  f_em_a_6 <- c(f_0_em_a_6, f_1_em_a_6, f_2_em_a_6)
  f_em_a_1 <- c(f_0_em_a_1, f_1_em_a_1, f_2_em_a_1)
  
  sols_em_b_6 <- moments_to_paras_func(moments_em_b_6, par_initial_6, type = "6")
  sols_em_b_4 <- moments_to_paras_func(moments_em_b_4, par_initial_4, type = "4")
  sols_em_b_1 <- inverse_expectation(mXA_2_em_b)
  
  f_2_em_b_6 <- forward_mapping_2(sols_em_b_6[1], sols_em_b_6[2])
  f_1_em_b_6 <- forward_mapping_1(sols_em_b_6[1], sols_em_b_6[2])
  f_0_em_b_6 <- forward_mapping_0(sols_em_b_6[1], sols_em_b_6[2])
  
  f_2_em_b_4 <- forward_mapping_2(sols_em_b_4[1], sols_em_b_4[2])
  f_1_em_b_4 <- forward_mapping_1(sols_em_b_4[1], sols_em_b_4[2])
  f_0_em_b_4 <- forward_mapping_0(sols_em_b_4[1], sols_em_b_4[2])
  
  f_2_em_b_1 <- forward_mapping_2(0, sols_em_b_1)
  f_1_em_b_1 <- forward_mapping_1(0, sols_em_b_1)
  f_0_em_b_1 <- forward_mapping_0(0, sols_em_b_1)
  
  f_em_b_4 <- c(f_0_em_b_4, f_1_em_b_4, f_2_em_b_4)
  f_em_b_6 <- c(f_0_em_b_6, f_1_em_b_6, f_2_em_b_6)
  f_em_b_1 <- c(f_0_em_b_1, f_1_em_b_1, f_2_em_b_1)
  ## alphas ----
  
  sols_em_S_1 <- (sols_em_b_1 + sols_em_a_1) / 2
  sols_em_S_6 <- (sols_em_b_6 + sols_em_a_6) / 2
  sols_em_S_4 <- (sols_em_b_4 + sols_em_a_4) / 2
  
  f_em_S_4 <- (f_em_b_4 + f_em_a_4) / 2
  f_em_S_6 <- (f_em_b_6 + f_em_a_6) / 2
  f_em_S_1 <- (f_em_b_1 + f_em_a_1) / 2
  
  Omega_multi_mean_X_em_S <- (Omega_multi_mean_X_em_a + Omega_multi_mean_X_em_b) / 2
  Omega_multi_mean_AX_em_S <- (Omega_multi_mean_AX_em_a + Omega_multi_mean_AX_em_b) / 2
  
  alpha_em_S_1_a <-  if (f_em_a_1[2] > 0) {
    as.numeric(Omega_multi_mean_AX_em_a / f_em_a_1[2])
  } else {
    rep(NA_real_, p)
  }
  alpha_em_S_1_b <-  if (f_em_b_1[2] > 0) {
    as.numeric(Omega_multi_mean_AX_em_b / f_em_b_1[2])
  } else {
    rep(NA_real_, p)
  }
  
  alpha_em_S_1 <-  (alpha_em_S_1_a + alpha_em_S_1_b) / 2
  
  alpha_em_S_4_a <-  if (f_em_a_4[2] > 0) {
    as.numeric((
      Omega_multi_mean_AX_em_a - f_em_a_4[1] *  Omega_multi_mean_X_em_a
    ) / f_em_a_4[2])
  } else {
    rep(NA_real_, p)
  }
  alpha_em_S_4_b <-  if (f_em_b_4[2] > 0) {
    as.numeric((
      Omega_multi_mean_AX_em_b - f_em_b_4[1] *  Omega_multi_mean_X_em_b
    ) / f_em_b_4[2])
  } else {
    rep(NA_real_, p)
  }
  
  alpha_em_S_4 <-  (alpha_em_S_4_a + alpha_em_S_4_b) / 2
  
  
  alpha_em_S_6_a <-  if (f_em_a_6[2] > 0) {
    as.numeric((
      Omega_multi_mean_AX_em_a - f_em_a_6[1] *  Omega_multi_mean_X_em_a
    ) / f_em_a_6[2])
  } else {
    rep(NA_real_, p)
  }
  alpha_em_S_6_b <-  if (f_em_b_6[2] > 0) {
    as.numeric((
      Omega_multi_mean_AX_em_b - f_em_b_6[1] *  Omega_multi_mean_X_em_b
    ) / f_em_b_6[2])
  } else {
    rep(NA_real_, p)
  }
  
  alpha_em_S_6 <-  (alpha_em_S_6_a + alpha_em_S_6_b) / 2
  
  ## jacobian ----
  jacobian_em_S_paras_1 <- jacobian_central_difference(moments_em_S_1,
                                                       par_initial_1,
                                                       type = "1",
                                                       f_function_name =  "moments_to_paras_func")
  jacobian_em_S_paras_4 <- jacobian_central_difference(moments_em_S_4,
                                                       par_initial_4,
                                                       type = "4",
                                                       f_function_name = "moments_to_paras_func")
  jacobian_em_S_paras_6 <- jacobian_central_difference(moments_em_S_6,
                                                       par_initial_6,
                                                       type = "6",
                                                       f_function_name = "moments_to_paras_func")
  
  jacobian_em_S_f_paras_1 <- jacobian_central_difference(moments_em_S_1,
                                                         par_initial_1,
                                                         type = "1",
                                                         f_function_name = "moments_to_f_paras_func")
  jacobian_em_S_f_paras_4 <- jacobian_central_difference(moments_em_S_4,
                                                         par_initial_4,
                                                         type = "4",
                                                         f_function_name = "moments_to_f_paras_func")
  jacobian_em_S_f_paras_6 <- jacobian_central_difference(moments_em_S_6,
                                                         par_initial_6,
                                                         type = "6",
                                                         f_function_name = "moments_to_f_paras_func")
  
  jacobian_em_S_alpha_1 <- jacobian_central_difference(
    moments_em_S_1,
    par_initial_1,
    type = "1",
    f_function_name = "moments_to_alpha_func",
    Omega_multi_mean_X = Omega_multi_mean_X_em_S,
    Omega_multi_mean_AX = Omega_multi_mean_AX_em_S
  )
  jacobian_em_S_alpha_4 <- jacobian_central_difference(
    moments_em_S_4,
    par_initial_4,
    type = "4",
    f_function_name = "moments_to_alpha_func",
    Omega_multi_mean_X = Omega_multi_mean_X_em_S,
    Omega_multi_mean_AX = Omega_multi_mean_AX_em_S
  )
  jacobian_em_S_alpha_6 <- jacobian_central_difference(
    moments_em_S_6,
    par_initial_6,
    type = "6",
    f_function_name = "moments_to_alpha_func",
    Omega_multi_mean_X = Omega_multi_mean_X_em_S,
    Omega_multi_mean_AX = Omega_multi_mean_AX_em_S
  )
  
  jacobian_em_S_alpha_split_sample_1 <- cbind(jacobian_em_S_alpha_1/2,jacobian_em_S_alpha_1/2)
  jacobian_em_S_alpha_split_sample_4 <- cbind(jacobian_em_S_alpha_4/2,jacobian_em_S_alpha_4/2)
  jacobian_em_S_alpha_split_sample_6 <- cbind(jacobian_em_S_alpha_6/2,jacobian_em_S_alpha_6/2)
  
  jacobian_em_S_paras_split_sample_1 <- cbind(jacobian_em_S_paras_1/2,jacobian_em_S_paras_1/2)
  jacobian_em_S_paras_split_sample_4 <- cbind(jacobian_em_S_paras_4/2,jacobian_em_S_paras_4/2)
  jacobian_em_S_paras_split_sample_6 <- cbind(jacobian_em_S_paras_6/2,jacobian_em_S_paras_6/2)
  
  jacobian_em_S_f_paras_split_sample_1 <- cbind(jacobian_em_S_f_paras_1/2,jacobian_em_S_f_paras_1/2)
  jacobian_em_S_f_paras_split_sample_4 <- cbind(jacobian_em_S_f_paras_4/2,jacobian_em_S_f_paras_4/2)
  jacobian_em_S_f_paras_split_sample_6 <- cbind(jacobian_em_S_f_paras_6/2,jacobian_em_S_f_paras_6/2)
  
  ## boostrap for split sample ------
  
  num_bootstrap <- n_half
  # 初始化权重矩阵
  weights_matrix_half <- matrix(0, nrow = B_bootstrap, ncol = num_bootstrap)
  # 生成多项式分布的样本
  for (i in 1:B_bootstrap) {
    weights_matrix_half[i, ] <- rmultinom(1, num_bootstrap, rep(1 / num_bootstrap, num_bootstrap))
  }
  cov_em_for_split_sample_4 <- matrix(NA_real_, 4 + 2 *
                                        p + 4 + 2 * p, 4 + 2 * p + 4 + 2 * p)
  cov_em_for_split_sample_4_a <- cov_em_for_split_sample_4_b <- cov_em_for_split_sample_4_ab <-  matrix(NA_real_, 4 + 2 *
                                                                                                          p , 4 + 2 * p)
  vector_one_half <- rep(1, n_half)
  
  ##### a-----
  mA_bootstrap_split_sample_a <- sapply(c(1:B_bootstrap), function(i)
    mean(A_2 * weights_matrix_half[i, ]))
  mAY_bootstrap_split_sample_a <- sapply(c(1:B_bootstrap), function(i)
    mean(A_2 * Y_2 * weights_matrix_half[i, ]))
  
  mX_2_bootstrap_split_sample_a_1 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (vector_one_half * weights_matrix_half[i, ]) * eigenMapMatMult(H_S_a, vector_one_half * weights_matrix_half[i, ])
    ) / (
      n_half * (n_half - 1)
    ))))
  mXA_1_bootstrap_split_sample_a_1 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_2 * weights_matrix_half[i, ]) * eigenMapMatMult(H_S_a, vector_one_half * weights_matrix_half[i, ])
    ) / (
      n_half * (n_half - 1)
    ))))
  mXA_2_bootstrap_split_sample_a_1 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_2 * weights_matrix_half[i, ]) * eigenMapMatMult(H_S_a, A_2 * weights_matrix_half[i, ])
    ) / (
      n_half * (n_half - 1)
    ))))
  mXAY_1_bootstrap_split_sample_a_1 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_2 * Y_2 * weights_matrix_half[i, ]) * eigenMapMatMult(H_S_a, vector_one_half * weights_matrix_half[i, ])
    ) / (
      n_half * (n_half - 1)
    ))))
  mXAYXA_bootstrap_split_sample_a_1 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_2 * Y_2 * weights_matrix_half[i, ]) * eigenMapMatMult(H_S_a, A_2 * weights_matrix_half[i, ])
    ) / (
      n_half * (n_half - 1)
    ))))
  
  mX_2_bootstrap_split_sample_a_2 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (vector_one_half * (weights_matrix_half[i, ] - 1)) * eigenMapMatMult(H_S_a, vector_one_half * (weights_matrix_half[i, ] -
                                                                                                       1))
    ) / (
      n_half * (n_half - 1)
    ))))
  mXA_1_bootstrap_split_sample_a_2 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_2 * (weights_matrix_half[i, ] - 1)) * eigenMapMatMult(H_S_a, vector_one_half * (weights_matrix_half[i, ] -
                                                                                           1))
    ) / (
      n_half * (n_half - 1)
    ))))
  mXA_2_bootstrap_split_sample_a_2 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_2 * (weights_matrix_half[i, ] - 1)) * eigenMapMatMult(H_S_a, A_2 * (weights_matrix_half[i, ] -
                                                                               1))
    ) / (
      n_half * (n_half - 1)
    ))))
  mXAY_1_bootstrap_split_sample_a_2 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_2 * Y_2 * (weights_matrix_half[i, ] - 1)) * eigenMapMatMult(H_S_a, vector_one_half * (weights_matrix_half[i, ] -
                                                                                                 1))
    ) / (
      n_half * (n_half - 1)
    ))))
  mXAYXA_bootstrap_split_sample_a_2 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_2 * Y_2 * (weights_matrix_half[i, ] - 1)) * eigenMapMatMult(H_S_a, A_2 * (weights_matrix_half[i, ] -
                                                                                     1))
    ) / (
      n_half * (n_half - 1)
    ))))
  
  t_X_2 <- t(X_2)
  Omega_multi_mean_X_split_sample_a <- t(sapply(c(1:B_bootstrap), function(i)
    (
      as.numeric(
        eigenMapMatMult(
          inv_sigma_split_a,
          eigenMapMatMult(t_X_2, weights_matrix_half[i, ])
        ) / n_half
      )
    )))
  Omega_multi_mean_AX_split_sample_a <- t(sapply(c(1:B_bootstrap), function(i)
    (
      as.numeric(
        eigenMapMatMult(
          inv_sigma_split_a,
          eigenMapMatMult(t_X_2, A_2 * weights_matrix_half[i, ])
        ) / n_half
      )
    )))
  
  
  m1_U1_a <- mA_bootstrap_split_sample_a
  m2_U2_a_1 <- cbind(
    mX_2_bootstrap_split_sample_a_1,
    mXA_1_bootstrap_split_sample_a_1,
    mXA_2_bootstrap_split_sample_a_1
  )
  m2_U2_a_2 <- cbind(
    mX_2_bootstrap_split_sample_a_2,
    mXA_1_bootstrap_split_sample_a_2,
    mXA_2_bootstrap_split_sample_a_2
  )
  m3_U1_a <- cbind(Omega_multi_mean_X_split_sample_a,
                   Omega_multi_mean_AX_split_sample_a)
  
  cov_em_for_split_sample_4_a[1, 1] <- var(m1_U1_a)
  cov_em_for_split_sample_4_a[1, c(2:4)] <- cov(m1_U1_a, m2_U2_a_1)
  cov_em_for_split_sample_4_a[c(2:4), 1] <- cov(m2_U2_a_1, m1_U1_a)
  
  cov_em_for_split_sample_4_a[c(2:4), c(2:4)] <- cov_bootstrap_U_stats(m2_U2_a_1, m2_U2_a_2, m2_U2_a_1, m2_U2_a_2)
  
  cov_em_for_split_sample_4_a[1, c((4 + 1):(4 + 2 *
                                              p))] <- cov(m1_U1_a, m3_U1_a)
  cov_em_for_split_sample_4_a[c((4 + 1):(4 + 2 *
                                           p)), 1] <- cov(m3_U1_a, m1_U1_a)
  
  cov_em_for_split_sample_4_a[c(2:4), c((4 + 1):(4 + 2 *
                                                   p))] <- cov(m2_U2_a_1, m3_U1_a)
  cov_em_for_split_sample_4_a[c((4 + 1):(4 + 2 *
                                           p)), c(2:4)] <- cov(m3_U1_a, m2_U2_a_1)
  
  cov_em_for_split_sample_4_a[c((4 + 1):(4 + 2 *
                                           p)), c((4 + 1):(4 + 2 * p))] <- eigenMapMatMult(t(m3_U1_a), m3_U1_a) / (B_bootstrap - 1) - eigenMapMatMult(colMeans(m3_U1_a), t(colMeans(m3_U1_a))) *
    B_bootstrap / (B_bootstrap - 1)
  
  #### b ------
  mA_bootstrap_split_sample_b <- sapply(c(1:B_bootstrap), function(i)
    mean(A_1 * weights_matrix_half[i, ]))
  mAY_bootstrap_split_sample_b <- sapply(c(1:B_bootstrap), function(i)
    mean(A_1 * Y_1 * weights_matrix_half[i, ]))
  
  mX_2_bootstrap_split_sample_b_1 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (vector_one_half * weights_matrix_half[i, ]) * eigenMapMatMult(H_S_b, vector_one_half * weights_matrix_half[i, ])
    ) / (
      n_half * (n_half - 1)
    ))))
  mXA_1_bootstrap_split_sample_b_1 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_1 * weights_matrix_half[i, ]) * eigenMapMatMult(H_S_b, vector_one_half * weights_matrix_half[i, ])
    ) / (
      n_half * (n_half - 1)
    ))))
  mXA_1_bootstrap_split_sample_b_1 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_1 * weights_matrix_half[i, ]) * eigenMapMatMult(H_S_b, A_1 * weights_matrix_half[i, ])
    ) / (
      n_half * (n_half - 1)
    ))))
  mXAY_1_bootstrap_split_sample_b_1 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_1 * Y_1 * weights_matrix_half[i, ]) * eigenMapMatMult(H_S_b, vector_one_half * weights_matrix_half[i, ])
    ) / (
      n_half * (n_half - 1)
    ))))
  mXAYXA_bootstrap_split_sample_b_1 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_1 * Y_1 * weights_matrix_half[i, ]) * eigenMapMatMult(H_S_b, A_1 * weights_matrix_half[i, ])
    ) / (
      n_half * (n_half - 1)
    ))))
  
  mX_2_bootstrap_split_sample_b_2 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (vector_one_half * (weights_matrix_half[i, ] - 1)) * eigenMapMatMult(H_S_b, vector_one_half * (weights_matrix_half[i, ] -
                                                                                                       1))
    ) / (
      n_half * (n_half - 1)
    ))))
  mXA_1_bootstrap_split_sample_b_2 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_1 * (weights_matrix_half[i, ] - 1)) * eigenMapMatMult(H_S_b, vector_one_half * (weights_matrix_half[i, ] -
                                                                                           1))
    ) / (
      n_half * (n_half - 1)
    ))))
  mXA_1_bootstrap_split_sample_b_2 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_1 * (weights_matrix_half[i, ] - 1)) * eigenMapMatMult(H_S_b, A_1 * (weights_matrix_half[i, ] -
                                                                               1))
    ) / (
      n_half * (n_half - 1)
    ))))
  mXAY_1_bootstrap_split_sample_b_2 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_1 * Y_1 * (weights_matrix_half[i, ] - 1)) * eigenMapMatMult(H_S_b, vector_one_half * (weights_matrix_half[i, ] -
                                                                                                 1))
    ) / (
      n_half * (n_half - 1)
    ))))
  mXAYXA_bootstrap_split_sample_b_2 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_1 * Y_1 * (weights_matrix_half[i, ] - 1)) * eigenMapMatMult(H_S_b, A_1 * (weights_matrix_half[i, ] -
                                                                                     1))
    ) / (
      n_half * (n_half - 1)
    ))))
  
  t_X_1 <- t(X_1)
  Omega_multi_mean_X_split_sample_b <- t(sapply(c(1:B_bootstrap), function(i)
    (
      as.numeric(
        eigenMapMatMult(
          inv_sigma_split_b,
          eigenMapMatMult(t_X_1, weights_matrix_half[i, ])
        ) / n_half
      )
    )))
  Omega_multi_mean_AX_split_sample_b <- t(sapply(c(1:B_bootstrap), function(i)
    (
      as.numeric(
        eigenMapMatMult(
          inv_sigma_split_b,
          eigenMapMatMult(t_X_1, A_1 * weights_matrix_half[i, ])
        ) / n_half
      )
    )))
  
  m1_U1_b <- mA_bootstrap_split_sample_b
  m2_U2_b_1 <- cbind(
    mX_2_bootstrap_split_sample_b_1,
    mXA_1_bootstrap_split_sample_b_1,
    mXA_1_bootstrap_split_sample_b_1
  )
  m2_U2_b_2 <- cbind(
    mX_2_bootstrap_split_sample_b_2,
    mXA_1_bootstrap_split_sample_b_2,
    mXA_1_bootstrap_split_sample_b_2
  )
  m3_U1_b <- cbind(Omega_multi_mean_X_split_sample_b,
                   Omega_multi_mean_AX_split_sample_b)
  
  cov_em_for_split_sample_4_b[1, 1] <- var(m1_U1_b)
  cov_em_for_split_sample_4_b[1, c(2:4)] <- cov(m1_U1_b, m2_U2_b_1)
  cov_em_for_split_sample_4_b[c(2:4), 1] <- cov(m2_U2_b_1, m1_U1_b)
  
  cov_em_for_split_sample_4_b[c(2:4), c(2:4)] <- cov_bootstrap_U_stats(m2_U2_b_1, m2_U2_b_2, m2_U2_b_1, m2_U2_b_2)
  
  cov_em_for_split_sample_4_b[1, c((4 + 1):(4 + 2 *
                                              p))] <- cov(m1_U1_b, m3_U1_b)
  cov_em_for_split_sample_4_b[c((4 + 1):(4 + 2 *
                                           p)), 1] <- cov(m3_U1_b, m1_U1_b)
  
  cov_em_for_split_sample_4_b[c(2:4), c((4 + 1):(4 + 2 *
                                                   p))] <- cov(m2_U2_b_1, m3_U1_b)
  cov_em_for_split_sample_4_b[c((4 + 1):(4 + 2 *
                                           p)), c(2:4)] <- cov(m3_U1_b, m2_U2_b_1)
  
  cov_em_for_split_sample_4_b[c((4 + 1):(4 + 2 *
                                           p)), c((4 + 1):(4 + 2 * p))] <- eigenMapMatMult(t(m3_U1_b), m3_U1_b) / (B_bootstrap - 1) - eigenMapMatMult(colMeans(m3_U1_b), t(colMeans(m3_U1_b))) *
    B_bootstrap / (B_bootstrap - 1)
  
  
  #### ab-----
  
  cov_em_for_split_sample_4_ab[1, 1] <- cov(m1_U1_a, m1_U1_b)
  cov_em_for_split_sample_4_ab[1, c(2:4)] <- cov(m1_U1_a, m2_U2_b_1)
  cov_em_for_split_sample_4_ab[c(2:4), 1] <- cov(m2_U2_b_1, m1_U1_a)
  
  cov_em_for_split_sample_4_ab[c(2:4), c(2:4)] <- cov_bootstrap_U_stats(m2_U2_a_1, m2_U2_a_2, m2_U2_b_1, m2_U2_b_2)
  
  cov_em_for_split_sample_4_ab[1, c((4 + 1):(4 + 2 *
                                               p))] <- cov(m1_U1_a, m3_U1_b)
  cov_em_for_split_sample_4_ab[c((4 + 1):(4 + 2 *
                                            p)), 1] <- cov(m3_U1_b, m1_U1_a)
  
  cov_em_for_split_sample_4_ab[c(2:4), c((4 + 1):(4 + 2 *
                                                    p))] <- cov(m2_U2_a_1, m3_U1_b)
  cov_em_for_split_sample_4_ab[c((4 + 1):(4 + 2 *
                                            p)), c(2:4)] <- cov(m3_U1_b, m2_U2_a_1)
  
  cov_em_for_split_sample_4_ab[c((4 + 1):(4 + 2 *
                                            p)), c((4 + 1):(4 + 2 * p))] <- eigenMapMatMult(t(m3_U1_a), m3_U1_b) / (B_bootstrap - 1) - eigenMapMatMult(colMeans(m3_U1_a), t(colMeans(m3_U1_b))) *
    B_bootstrap / (B_bootstrap - 1)
  
  
  #### summary -----
  cov_em_for_split_sample_4 <- rbind(
    cbind(
      cov_em_for_split_sample_4_a,
      cov_em_for_split_sample_4_ab
    ),
    cbind(
      t(cov_em_for_split_sample_4_ab),
      cov_em_for_split_sample_4_b
    )
  )
  cov_alpha_split_sample_4 <- eigenMapMatMult(
    eigenMapMatMult(
      jacobian_em_S_alpha_split_sample_4,
      cov_em_for_split_sample_4
    ),
    t(jacobian_em_S_alpha_split_sample_4)
  )
  
  
  var_moments_split_sample_4 <- diag(cov_em_for_split_sample_4)
  var_alpha_split_sample_4 <- diag(cov_alpha_split_sample_4)
  cov_em_for_paras_split_sample_4 <- rbind(
    cbind(
      cov_em_for_split_sample_4_a[c((1):(4)), c((1):(4))],
      cov_em_for_split_sample_4_ab[c((1):(4)), c((1):(4))]
    ),
    cbind(
      t(cov_em_for_split_sample_4_ab[c((1):(4)), c((1):(4))]),
      cov_em_for_split_sample_4_b[c((1):(4)), c((1):(4))]
    )
  )
  
  var_paras_split_sample_4 <- diag(eigenMapMatMult(
    eigenMapMatMult(
      jacobian_em_S_paras_split_sample_4,
      cov_em_for_paras_split_sample_4
    ),
    t(jacobian_em_S_paras_split_sample_4)
  ))
  var_f_paras_split_sample_4 <- diag(eigenMapMatMult(
    eigenMapMatMult(
      jacobian_em_S_f_paras_split_sample_4,
      cov_em_for_paras_split_sample_4
    ),
    t(jacobian_em_S_f_paras_split_sample_4)
  ))
  
  
  cov_em_for_split_sample_1 <- rbind(
    cbind(
      cov_em_for_split_sample_4_a[c((4):(4 + 2 * p)), c((4):(4 + 2 * p))],
      cov_em_for_split_sample_4_ab[c((4):(4 + 2 * p)), c((4):(4 + 2 * p))]
    ),
    cbind(
      t(cov_em_for_split_sample_4_ab[c((4):(4 + 2 * p)), c((4):(4 + 2 * p))]),
      cov_em_for_split_sample_4_b[c((4):(4 + 2 * p)), c((4):(4 + 2 * p))]
    )
  )
  
  var_alpha_split_sample_1 <- diag(eigenMapMatMult(
    eigenMapMatMult(
      jacobian_em_S_alpha_split_sample_1,
      cov_em_for_split_sample_1
    ),
    t(jacobian_em_S_alpha_split_sample_1)
  ))
  
  cov_em_for_paras_split_sample_1 <- rbind(
    cbind(
      cov_em_for_split_sample_4_a[4, 4],
      cov_em_for_split_sample_4_ab[4, 4]
    ),
    cbind(
      t(cov_em_for_split_sample_4_ab[4, 4]),
      cov_em_for_split_sample_4_b[4, 4]
    )
  )
  
  var_paras_split_sample_1 <- diag(eigenMapMatMult(
    eigenMapMatMult(
      jacobian_em_S_paras_split_sample_1,
      cov_em_for_paras_split_sample_1
    ),
    t(jacobian_em_S_paras_split_sample_1)
  ))
  var_f_paras_split_sample_1 <- diag(eigenMapMatMult(
    eigenMapMatMult(
      jacobian_em_S_f_paras_split_sample_1,
      cov_em_for_paras_split_sample_1
    ),
    t(jacobian_em_S_f_paras_split_sample_1)
  ))
  ### return ----
  return(list(
    sols_em_S_1 = sols_em_S_1,
    sols_em_S_6 = sols_em_S_6,
    sols_em_S_4 = sols_em_S_4,
    
    f_em_S_1 = f_em_S_1,
    f_em_S_6 = f_em_S_6,
    f_em_S_4 = f_em_S_4,
    
    
    alpha_em_S_1 = alpha_em_S_1,
    alpha_em_S_4 = alpha_em_S_4,
    alpha_em_S_6 = alpha_em_S_6,
    
    moments_em_S_6 = moments_em_S_6,
    moments_em_S_4 = moments_em_S_4,
    moments_em_S_1 = moments_em_S_1,
    
    var_moments_split_sample_4 = var_moments_split_sample_4,
    
    var_alpha_split_sample_4 = var_alpha_split_sample_4,
    var_paras_split_sample_4 = var_paras_split_sample_4,
    var_f_paras_split_sample_4 = var_f_paras_split_sample_4,
    
    var_alpha_split_sample_1 = var_alpha_split_sample_1,
    var_paras_split_sample_1 = var_paras_split_sample_1,
    var_f_paras_split_sample_1 = var_f_paras_split_sample_1
    
  ))
}
glm_mom_known_sigma_mar <- function(X_T,A_T,Y_T,omega_11,B_bootstrap,par_initial_6){
  n <- nrow(X_T)
  p <- ncol(X_T)
  
  H <- eigenMapMatMult(X_T, t(X_T)) * omega_11
  H <- H - diag(diag(H))
  vector_one <- rep(1, n)
  
  mA_em <- mean(A_T)
  mX_2_em <- as.numeric(sum(vector_one * eigenMapMatMult(H, vector_one)) / (n * (n - 1)))
  mXA_1_em <- as.numeric(sum(A_T * eigenMapMatMult(H, vector_one)) / (n * (n - 1)))
  mXA_2_em <- as.numeric(sum(A_T * eigenMapMatMult(H, A_T)) / (n * (n - 1)))
  mAY_em <-  as.numeric(mean(A_T * Y_T))
  mXAY_1_em <-  as.numeric(sum(A_T * Y_T * eigenMapMatMult(H, vector_one)) / (n * (n - 1)))
  
  moments_em_6 <-  c(mA_em, mX_2_em, mXA_1_em, mXA_2_em, mAY_em, mXAY_1_em)
  moments_em_4 <-  c(mA_em, mX_2_em, mXA_1_em, mXA_2_em)
  moments_em_1 <- mXA_2_em
  Omega_multi_mean_X_em <- colMeans(X_T) * omega_11
  Omega_multi_mean_AX_em <- omega_11 * colMeans(X_T * A_T)
  
  par_initial_1 <- par_initial_6[1]
  par_initial_4 <- par_initial_6[1:2]
  
  jacobian_em_paras_1 <- jacobian_central_difference(moments_em_1,
                                                     par_initial_1,
                                                     type = "1",
                                                     f_function_name =  "moments_to_paras_func")
  jacobian_em_paras_4 <- jacobian_central_difference(moments_em_4,
                                                     par_initial_4,
                                                     type = "4",
                                                     f_function_name = "moments_to_paras_func")
  jacobian_em_paras_6 <- jacobian_central_difference(moments_em_6,
                                                     par_initial_6,
                                                     type = "6",
                                                     f_function_name = "moments_to_paras_func")
  
  jacobian_em_f_paras_1 <- jacobian_central_difference(moments_em_1,
                                                       par_initial_1,
                                                       type = "1",
                                                       f_function_name = "moments_to_f_paras_func")
  jacobian_em_f_paras_4 <- jacobian_central_difference(moments_em_4,
                                                       par_initial_4,
                                                       type = "4",
                                                       f_function_name = "moments_to_f_paras_func")
  jacobian_em_f_paras_6 <- jacobian_central_difference(moments_em_6,
                                                       par_initial_6,
                                                       type = "6",
                                                       f_function_name = "moments_to_f_paras_func")
  
  jacobian_em_alpha_1 <- jacobian_central_difference(
    moments_em_1,
    par_initial_1,
    type = "1",
    f_function_name = "moments_to_alpha_func",
    Omega_multi_mean_X = Omega_multi_mean_X_em,
    Omega_multi_mean_AX = Omega_multi_mean_AX_em
  )
  jacobian_em_alpha_4 <- jacobian_central_difference(
    moments_em_4,
    par_initial_4,
    type = "4",
    f_function_name = "moments_to_alpha_func",
    Omega_multi_mean_X = Omega_multi_mean_X_em,
    Omega_multi_mean_AX = Omega_multi_mean_AX_em
  )
  jacobian_em_alpha_6 <- jacobian_central_difference(
    moments_em_6,
    par_initial_6,
    type = "6",
    f_function_name = "moments_to_alpha_func",
    Omega_multi_mean_X = Omega_multi_mean_X_em,
    Omega_multi_mean_AX = Omega_multi_mean_AX_em
  )
  
  ## paras,alpha ----
  
  sols_em_6 <- moments_to_paras_func(moments_em_6, par_initial_6, type = "6")
  sols_em_4 <- moments_to_paras_func(moments_em_4, par_initial_4, type = "4")
  sols_em_1 <- inverse_expectation(mXA_2_em)
  
  f_2_em_6 <- forward_mapping_2(sols_em_6[1], sols_em_6[2])
  f_1_em_6 <- forward_mapping_1(sols_em_6[1], sols_em_6[2])
  f_0_em_6 <- forward_mapping_0(sols_em_6[1], sols_em_6[2])
  
  f_2_em_4 <- forward_mapping_2(sols_em_4[1], sols_em_4[2])
  f_1_em_4 <- forward_mapping_1(sols_em_4[1], sols_em_4[2])
  f_0_em_4 <- forward_mapping_0(sols_em_4[1], sols_em_4[2])
  
  f_2_em_1 <- forward_mapping_2(0, sols_em_1)
  f_1_em_1 <- forward_mapping_1(0, sols_em_1)
  f_0_em_1 <- forward_mapping_0(0, sols_em_1)
  
  f_em_1 <- c(f_0_em_1, f_1_em_1, f_2_em_1)
  f_em_4 <- c(f_0_em_4, f_1_em_4, f_2_em_4)
  f_em_6 <- c(f_0_em_6, f_1_em_6, f_2_em_6)
  
  
  alpha_em_6 <-  if (f_1_em_6 > 0) {
    (Omega_multi_mean_AX_em - f_0_em_6 * Omega_multi_mean_X_em) / f_1_em_6
  } else {
    rep(NA_real_, p)
  }
  
  alpha_em_4 <-  if (f_1_em_4 > 0) {
    (Omega_multi_mean_AX_em - f_0_em_4 * Omega_multi_mean_X_em) / f_1_em_4
  } else {
    rep(NA_real_, p)
  }
  
  alpha_em_1 <-  if (f_1_em_1 > 0) {
    Omega_multi_mean_AX_em / f_1_em_1
  } else {
    rep(NA_real_, p)
  }
  ## boostrap for variance ----
  num_bootstrap <- n
  weights_matrix <- matrix(0, nrow = B_bootstrap, ncol = num_bootstrap)
  for (i in 1:B_bootstrap) {
    weights_matrix[i, ] <- rmultinom(1, num_bootstrap, rep(1 / num_bootstrap, num_bootstrap))
  }
  
  
  mA_bootstrap <- sapply(c(1:B_bootstrap), function(i)
    mean(A_T * weights_matrix[i, ]))
  mAY_bootstrap <- sapply(c(1:B_bootstrap), function(i)
    mean(A_T * Y_T * weights_matrix[i, ]))
  
  mX_2_bootstrap_1 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (vector_one * weights_matrix[i, ]) * eigenMapMatMult(H, vector_one * weights_matrix[i, ])
    ) / (n * (
      n - 1
    )))))
  mXA_1_bootstrap_1 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_T * weights_matrix[i, ]) * eigenMapMatMult(H, vector_one * weights_matrix[i, ])
    ) / (n * (
      n - 1
    )))))
  mXA_2_bootstrap_1 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_T * weights_matrix[i, ]) * eigenMapMatMult(H, A_T * weights_matrix[i, ])
    ) / (n * (
      n - 1
    )))))
  mXAY_1_bootstrap_1 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_T * Y_T * weights_matrix[i, ]) * eigenMapMatMult(H, vector_one * weights_matrix[i, ])
    ) / (n * (
      n - 1
    )))))
  mXAYXA_bootstrap_1 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_T * Y_T * weights_matrix[i, ]) * eigenMapMatMult(H, A_T * weights_matrix[i, ])
    ) / (n * (
      n - 1
    )))))
  
  mX_2_bootstrap_2 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (vector_one * (weights_matrix[i, ] - 1)) * eigenMapMatMult(H, vector_one * (weights_matrix[i, ] -
                                                                                    1))
    ) / (n * (
      n - 1
    )))))
  mXA_1_bootstrap_2 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_T * (weights_matrix[i, ] - 1)) * eigenMapMatMult(H, vector_one * (weights_matrix[i, ] -
                                                                             1))
    ) / (n * (
      n - 1
    )))))
  mXA_2_bootstrap_2 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_T * (weights_matrix[i, ] - 1)) * eigenMapMatMult(H, A_T * (weights_matrix[i, ] -
                                                                      1))
    ) / (n * (
      n - 1
    )))))
  mXAY_1_bootstrap_2 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_T * Y_T * (weights_matrix[i, ] - 1)) * eigenMapMatMult(H, vector_one * (weights_matrix[i, ] -
                                                                                   1))
    ) / (n * (
      n - 1
    )))))
  mXAYXA_bootstrap_2 <- sapply(c(1:B_bootstrap), function(i)
    (as.numeric(sum(
      (A_T * Y_T * (weights_matrix[i, ] - 1)) * eigenMapMatMult(H, A_T * (weights_matrix[i, ] -
                                                                            1))
    ) / (n * (
      n - 1
    )))))
  
  Omega_multi_mean_X_bootstrap <- t(sapply(c(1:B_bootstrap), function(i)
    (as.numeric(
      (eigenMapMatMult(t(weights_matrix[i, ]), X_T)) * (omega_11  / n)
    ))))
  Omega_multi_mean_AX_bootstrap <- t(sapply(c(1:B_bootstrap), function(i)
    (as.numeric((eigenMapMatMult(t(
      A_T * weights_matrix[i, ]
    ), X_T)) * (omega_11  / n)
    ))))
  
  m1_U1 <- mA_bootstrap
  m2_U2_1 <- cbind(mX_2_bootstrap_1, mXA_1_bootstrap_1, mXA_2_bootstrap_1)
  m2_U2_2 <- cbind(mX_2_bootstrap_2, mXA_1_bootstrap_2, mXA_2_bootstrap_2)
  m3_U1 <- cbind(Omega_multi_mean_X_bootstrap,
                 Omega_multi_mean_AX_bootstrap)
  
  cov_em_4 <- matrix(NA_real_, 4 + 2 * p , 4 + 2 * p)
  cov_em_4[1, 1] <- var(m1_U1)
  cov_em_4[1, c(2:4)] <- cov(m1_U1, m2_U2_1)
  cov_em_4[c(2:4), 1] <- cov(m2_U2_1, m1_U1)
  
  cov_em_4[c(2:4), c(2:4)] <- cov_bootstrap_U_stats(m2_U2_1, m2_U2_2, m2_U2_1, m2_U2_2)
  
  cov_em_4[1, c((4 + 1):(4 + 2 * p))] <- cov(m1_U1, m3_U1)
  cov_em_4[c((4 + 1):(4 + 2 * p)), 1] <- cov(m3_U1, m1_U1)
  
  cov_em_4[c(2:4), c((4 + 1):(4 + 2 * p))] <- cov(m2_U2_1, m3_U1)
  cov_em_4[c((4 + 1):(4 + 2 * p)), c(2:4)] <- cov(m3_U1, m2_U2_1)
  
  cov_em_4[c((4 + 1):(4 + 2 * p)), c((4 + 1):(4 + 2 * p))] <- eigenMapMatMult(t(m3_U1), m3_U1) / (B_bootstrap - 1) - eigenMapMatMult(colMeans(m3_U1), t(colMeans(m3_U1))) *
    B_bootstrap / (B_bootstrap - 1)
  
  var_moments_em_4 <- diag(cov_em_4)
  
  
  cov_em_1 <- cov_em_4[4, 4]
  cov_em_1_alpha <- cov_em_4[c((4):(4 + 2 * p)), c((4):(4 + 2 * p))]
  
  var_paras_em_1 <- cov_em_1 * jacobian_em_paras_1[1, 1] ^ 2
  var_f_paras_em_1 <- diag(cov_em_1 * jacobian_em_f_paras_1 %*% t(jacobian_em_f_paras_1))
  var_alpha_em_1 <- diag(eigenMapMatMult(jacobian_em_alpha_1, eigenMapMatMult(cov_em_1_alpha, t(jacobian_em_alpha_1))))
  
  var_paras_em_4 <- diag(eigenMapMatMult(jacobian_em_paras_4, eigenMapMatMult(cov_em_4[c(1:4), c(1:4)], t(jacobian_em_paras_4))))
  var_f_paras_em_4 <- diag(eigenMapMatMult(jacobian_em_f_paras_4, eigenMapMatMult(cov_em_4[c(1:4), c(1:4)], t(jacobian_em_f_paras_4))))
  var_alpha_em_4 <- diag(eigenMapMatMult(jacobian_em_alpha_4, eigenMapMatMult(cov_em_4, t(jacobian_em_alpha_4))))
  
  return(list(
    
    var_moments_em_4 = var_moments_em_4,
    
    var_paras_em_1 = var_paras_em_1,
    var_f_paras_em_1 = var_f_paras_em_1,
    var_alpha_em_1 = var_alpha_em_1,
    
    var_paras_em_4 = var_paras_em_4,
    var_f_paras_em_4 = var_f_paras_em_4,
    var_alpha_em_4 = var_alpha_em_4,
    
    sols_em_1 = sols_em_1,
    sols_em_6 = sols_em_6,
    sols_em_4 = sols_em_4,
    
    f_em_1 = f_em_1,
    f_em_6 = f_em_6,
    f_em_4 = f_em_4,
    
    alpha_em_1 = alpha_em_1,
    alpha_em_4 = alpha_em_4,
    alpha_em_6 = alpha_em_6
  ))
}

# run.relicate -----
run_replicate_mom <- function(n,
                              p,
                              seed,
                              alpha,
                              beta,
                              Is_sparse,
                              Is_sparse_only_1,
                              Is_Rad,
                              mu_x,
                              sd_mar_linear,
                              omega_11,
                              epsilon_initial,
                              B_bootstrap,
                              filename.tracking) {
  cat(n, seed, "\n", file = filename.tracking, append = TRUE)
  result <- tryCatch({
    mu <- rep(mu_x, p)
    
    # paras_truth ----
    par1_truth <- sum(alpha * mu)
    par2_truth <- sum(alpha * alpha / omega_11)
    cov_truth <- sum(alpha * beta / omega_11)
    psi_truth <- sum(beta * mu)
    
    par_truth_6 <- c(par1_truth, par2_truth, cov_truth, psi_truth)
    par_truth_4 <- c(par1_truth, par2_truth)
    par_truth_1 <- par2_truth
    
    f_2_truth  <- forward_mapping_2(par_truth_6[1], par_truth_6[2])
    f_1_truth <- forward_mapping_1(par_truth_6[1], par_truth_6[2])
    f_0_truth  <- forward_mapping_0(par_truth_6[1], par_truth_6[2])
    f_truth <- c(f_0_truth, f_1_truth, f_2_truth)
    
    
    mA_truth <- f_0_truth
    mX_2_truth <- as.numeric(sum(mu * mu) * omega_11)
    mXA_1_truth <- as.numeric(mA_truth * mX_2_truth + f_1_truth * par1_truth)
    mXA_2_truth <- as.numeric(
      mA_truth ^ 2 * mX_2_truth + f_1_truth ^ 2 * par2_truth + 2 * mA_truth * f_1_truth * par1_truth
    )
    mAY_truth <-  as.numeric(mA_truth * psi_truth + f_1_truth * cov_truth)
    mXAY_1_truth <- as.numeric(
      (mA_truth + mXA_1_truth) * psi_truth + (mX_2_truth * f_1_truth + f_2_truth * par1_truth) * cov_truth
    )
    moments_truth_6 <-  c(mA_truth,
                          mX_2_truth,
                          mXA_1_truth,
                          mXA_2_truth,
                          mAY_truth,
                          mXAY_1_truth)
    
    moments_truth_4 <-  c(mA_truth, mX_2_truth, mXA_1_truth, mXA_2_truth)
    moments_truth_1 <- mXA_2_truth
    
    Omega_multi_mean_X_truth <- mu * omega_11
    Omega_multi_mean_AX_truth <- alpha * f_1_truth * omega_11
    
    par_initial_6 <- par_truth_6 + epsilon_initial
    par_initial_4 <- par_initial_6[c(1:2)]
    par_initial_1 <- par_initial_6[2]
    sols_oracle_1 <- moments_to_paras_func(moments_truth_1, par_initial_1, type = "1")
    sols_oracle_4 <-  moments_to_paras_func(moments_truth_4, par_initial_4, type = "4")
    sols_oracle_6 <-  moments_to_paras_func(moments_truth_6, par_initial_6, type = "6")
    
    # simulation ----
    ## generate data ----
    
    data_XAY <- generate_data_XAY(n, p, omega_11, seed, alpha, beta, mu_x = mu_x, sd_mar_linear, Is_Rad)
    X_T <- data_XAY$X_T
    A_T <- data_XAY$A_T
    Y_T <- data_XAY$Y_T
    
    known_sigma_results <- glm_mom_known_sigma_mar(
      X_T,
      A_T,
      Y_T,
      omega_11 = omega_11,
      B_bootstrap = B_bootstrap,
      par_initial_6 = par_initial_6
    )
    if (p < n) {
      unknown_sigma_results <- glm_mom_unknown_sigma_mar(X_T, A_T, Y_T, B_bootstrap, par_initial_6 = par_initial_6)
    }
    
    # return ----
    return(c(
      list(
        success = TRUE,
        n = n,
        p = p,
        seed = seed,
        Is_sparse = Is_sparse,
        Is_sparse_only_1 = Is_sparse_only_1,
        Is_Rad = Is_Rad,
        
        B_bootstrap = B_bootstrap,
        
        alpha = alpha,
        beta = beta,
        omega_11 = omega_11,
        mu_x = mu_x,
        
        sd_mar_linear = sd_mar_linear,
        epsilon_initial = epsilon_initial,
        
        moments_truth_6 = moments_truth_6,
        par_truth_6 = par_truth_6,
        sols_oracle_6 = sols_oracle_6,
        f_truth = f_truth,
        Omega_multi_mean_X_truth = Omega_multi_mean_X_truth,
        Omega_multi_mean_AX_truth = Omega_multi_mean_AX_truth
      ),
      known_sigma_results,
      unknown_sigma_results
    ))
  }, error = function(e) {
    # 如果出错，记录错误并返回失败的信号
    cat("Error: ", e$message, "\n")
    return(list(
      success = FALSE,
      message = e$message,
      seed = seed,
      n = n,
      p = p
    ))
  })
  
  
}
run_replicate_glm <- function(n,
                              p,
                              seed,
                              alpha,
                              beta,
                              Is_sparse,
                              Is_sparse_only_1,
                              Is_Rad,
                              mu_x,
                              sd_mar_linear,
                              omega_11,
                              epsilon_initial,
                              B_bootstrap,
                              filename.tracking) {
  # simulation ----
  ## generate data ----
  
  data_XAY <- generate_data_XAY(n, p, omega_11, seed, alpha, beta, mu_x = mu_x, sd_mar_linear, Is_Rad)
  X_T <- data_XAY$X_T
  A_T <- data_XAY$A_T
  Y_T <- data_XAY$Y_T
  
  df <- as.data.frame(data_XAY)
  glm_results <- glm(formula = A_T~X_T-1, df,family= binomial)
  # return ----
  return(
    list(
      n = n,
      p = p,
      seed = seed,
      Is_sparse = Is_sparse,
      Is_sparse_only_1 = Is_sparse_only_1,
      Is_Rad = Is_Rad,
      glm_converged = glm_results$converged,
      glm_boundary = glm_results$boundary,
      glm_coefficients = glm_results$coefficients
    ))
}



# 保存当前环境中的所有对象
save(list = ls(), file = "all_functions_and_objects_mom.RData")
