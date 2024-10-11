source("demo_glm_MoM/dependecy.R")
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
#
# plot(seq(-10,10,0.01),sapply(seq(-10,10,0.01),function(i)(integrand_forward_mapping_quad(i,0,1))))
# plot(seq(-10,10,0.01),sapply(seq(-10,10,0.01),function(i)(integrand_forward_mapping_hess(i,0,1))))

forward_mapping_0 <- function (mu, sigma_sq) {
  integrate(
    integrand_forward_mapping,
    mu - 10,
    mu + 10,
    mu = mu,
    sigma_sq = sigma_sq,
    stop.on.error = FALSE
  )$value
}
forward_mapping_1 <- function (mu, sigma_sq) {
  integrate(
    integrand_forward_mapping_quad,
    mu - 10,
    mu + 10,
    mu = mu,
    sigma_sq = sigma_sq,
    stop.on.error = FALSE
  )$value
}
forward_mapping_2 <- function (mu, sigma_sq) {
  integrate(
    integrand_forward_mapping_hess,
    mu - 10,
    mu + 10,
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
tol <- 5e-6
sigma_L2_values <- seq(0.01, 2, by = tol)
expectation_values <- sapply(sigma_L2_values, expectation_G)

# 创建反函数
inverse_expectation <- approxfun(expectation_values, sigma_L2_values, rule = 2)

### Jacobian----
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
    mX_2 <- ifelse(moments[2] > 0, moments[2], 1e-5)
    mXA_1 <- moments[3]
    mXA_2 <- ifelse(moments[4] > 0, moments[4], 1e-5)
    
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
    mX_2 <- ifelse(moments[2] > 0, moments[2], 1e-5)
    mXA_1 <- moments[3]
    mXA_2 <- ifelse(moments[4] > 0, moments[4], 1e-5)
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
    mX_2 <- ifelse(moments[2] > 0, moments[2], 1e-5)
    mXA_1 <- moments[3]
    mXA_2 <- ifelse(moments[4] > 0, moments[4], 1e-5)
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
  sols <- moments_to_paras_func(moments, par_initial, type)
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

moments_to_alpha_func <- function(moments,
                                  par_initial,
                                  type,
                                  Omega_multi_mean_X,
                                  Omega_multi_mean_AX) {
  f <-  moments_to_f_paras_func(moments, par_initial, type)
  
  alpha <-  if (f[2] > 0) {
    (Omega_multi_mean_AX - f[1] * Omega_multi_mean_X) / f[2]
  } else {
    rep(NA, p)
  }
  
  return(alpha)
}
moments_to_paras_split_sample_func <- function(type, par_initial, moments_1, moments_2) {
  paras_1 <- moments_to_paras_func(moments_1, par_initial, type)
  
  paras_2 <- moments_to_paras_func(moments_2, par_initial, type)
  
  paras <- (paras_1 + paras_2) / 2
  return(paras)
}
moments_to_f_paras_split_sample_func <- function(type, par_initial, moments_1, moments_2) {
  paras_1 <- moments_to_f_paras_func(moments_1, par_initial, type)
  
  paras_2 <- moments_to_f_paras_func(moments_2, par_initial, type)
  
  paras <- (paras_1 + paras_2) / 2
  return(paras)
}
moments_to_alpha_split_sample_func <- function(type,
                                               par_initial,
                                               moments_1,
                                               Omega_multi_mean_X_1,
                                               Omega_multi_mean_AX_1,
                                               moments_2,
                                               Omega_multi_mean_X_2,
                                               Omega_multi_mean_AX_2) {
  alpha_split_1 <- moments_to_alpha_func(moments_1,
                                         par_initial,
                                         type,
                                         Omega_multi_mean_X_1,
                                         Omega_multi_mean_AX_1)
  
  alpha_split_2 <- moments_to_alpha_func(moments_2,
                                         par_initial,
                                         type,
                                         Omega_multi_mean_X_2,
                                         Omega_multi_mean_AX_2)
  
  alpha <- (alpha_split_1 + alpha_split_2) / 2
  return(alpha)
}

call_f_function <- function(moments,
                            par_initial,
                            type,
                            f_function_name,
                            Omega_multi_mean_X = NA,
                            Omega_multi_mean_AX = NA,
                            moments_2 = NA,
                            Omega_multi_mean_X_2 = NA,
                            Omega_multi_mean_AX_2 = NA) {
  if (f_function_name == "moments_to_paras_func") {
    return(moments_to_paras_func(moments, par_initial, type))
  } else if (f_function_name == "moments_to_f_paras_func") {
    return(moments_to_f_paras_func(moments, par_initial, type))
  } else if (f_function_name == "moments_to_alpha_func") {
    return(
      moments_to_alpha_func(
        moments,
        par_initial,
        type,
        Omega_multi_mean_X,
        Omega_multi_mean_AX
      )
    )
  } else if (f_function_name == "moments_to_alpha_split_sample_func") {
    return(
      moments_to_alpha_split_sample_func(
        type,
        par_initial,
        moments,
        Omega_multi_mean_X,
        Omega_multi_mean_AX,
        moments_2,
        Omega_multi_mean_X_2,
        Omega_multi_mean_AX_2
      )
    )
  } else if (f_function_name == "moments_to_paras_split_sample_func") {
    return(moments_to_paras_split_sample_func(type, par_initial, moments, moments_2))
  } else if (f_function_name == "moments_to_f_paras_split_sample_func") {
    return(moments_to_f_paras_split_sample_func(type, par_initial, moments, moments_2))
  } else {
    stop("Unknown function name")
  }
}
#############################################
# Too small epsilon will make some approximate function can not separate the values
# like (inverse_expectation(0.09+1e-6) - inverse_expectation(0.09))/(1e-6)
#############################################
jacobian_central_difference_main <- function(row_num,
                                             var_name,
                                             moments,
                                             par_initial,
                                             type,
                                             f_function_name,
                                             Omega_multi_mean_X = NA,
                                             Omega_multi_mean_AX = NA,
                                             moments_2 = NA,
                                             Omega_multi_mean_X_2 = NA,
                                             Omega_multi_mean_AX_2 = NA,
                                             epsilon = 1e-5) {
  var_used <- get(var_name)
  n <- length(var_used)
  jacobian <- matrix(0, nrow = row_num, ncol = n)
  for (i in 1:n) {
    var_used_plus <- var_used
    var_used_plus[i] <- var_used_plus[i] + epsilon
    
    var_used_minus <- var_used
    var_used_minus[i] <- var_used_minus[i] - epsilon
    
    if (var_name == "moments") {
      f_plus <- call_f_function(
        var_used_plus,
        par_initial,
        type,
        f_function_name,
        Omega_multi_mean_X = Omega_multi_mean_X,
        Omega_multi_mean_AX = Omega_multi_mean_AX,
        moments_2 = moments_2,
        Omega_multi_mean_X_2 = Omega_multi_mean_X_2,
        Omega_multi_mean_AX_2 = Omega_multi_mean_AX_2
      )
      f_minus <- call_f_function(
        var_used_minus,
        par_initial,
        type,
        f_function_name,
        Omega_multi_mean_X = Omega_multi_mean_X,
        Omega_multi_mean_AX = Omega_multi_mean_AX,
        moments_2 = moments_2,
        Omega_multi_mean_X_2 = Omega_multi_mean_X_2,
        Omega_multi_mean_AX_2 = Omega_multi_mean_AX_2
      )
      
    } else if (var_name == "moments_2") {
      f_plus <- call_f_function(
        moments,
        par_initial,
        type,
        f_function_name,
        Omega_multi_mean_X,
        Omega_multi_mean_AX,
        var_used_plus,
        Omega_multi_mean_X_2,
        Omega_multi_mean_AX_2
      )
      f_minus <- call_f_function(
        moments,
        par_initial,
        type,
        f_function_name,
        Omega_multi_mean_X,
        Omega_multi_mean_AX,
        var_used_minus,
        Omega_multi_mean_X_2,
        Omega_multi_mean_AX_2
      )
      
    } else if (var_name == "Omega_multi_mean_X") {
      f_plus <- call_f_function(
        moments,
        par_initial,
        type,
        f_function_name,
        var_used_plus,
        Omega_multi_mean_AX,
        moments_2,
        Omega_multi_mean_X_2,
        Omega_multi_mean_AX_2
      )
      f_minus <- call_f_function(
        moments,
        par_initial,
        type,
        f_function_name,
        var_used_minus,
        Omega_multi_mean_AX,
        moments_2,
        Omega_multi_mean_X_2,
        Omega_multi_mean_AX_2
      )
      
    } else if (var_name == "Omega_multi_mean_X_2") {
      f_plus <- call_f_function(
        moments,
        par_initial,
        type,
        f_function_name,
        Omega_multi_mean_X,
        Omega_multi_mean_AX,
        moments_2,
        var_used_plus,
        Omega_multi_mean_AX_2
      )
      f_minus <- call_f_function(
        moments,
        par_initial,
        type,
        f_function_name,
        Omega_multi_mean_X,
        Omega_multi_mean_AX,
        moments_2,
        var_used_minus,
        Omega_multi_mean_AX_2
      )
      
    } else if (var_name == "Omega_multi_mean_AX") {
      f_plus <- call_f_function(
        moments,
        par_initial,
        type,
        f_function_name,
        Omega_multi_mean_X,
        var_used_plus,
        moments_2,
        Omega_multi_mean_X_2,
        Omega_multi_mean_AX_2
      )
      f_minus <- call_f_function(
        moments,
        par_initial,
        type,
        f_function_name,
        Omega_multi_mean_X,
        var_used_minus,
        moments_2,
        Omega_multi_mean_X_2,
        Omega_multi_mean_AX_2
      )
      
    } else if (var_name == "Omega_multi_mean_AX_2") {
      f_plus <- call_f_function(
        moments,
        par_initial,
        type,
        f_function_name,
        Omega_multi_mean_X,
        Omega_multi_mean_AX,
        moments_2,
        Omega_multi_mean_X_2,
        var_used_plus
      )
      f_minus <- call_f_function(
        moments,
        par_initial,
        type,
        f_function_name,
        Omega_multi_mean_X,
        Omega_multi_mean_AX,
        moments_2,
        Omega_multi_mean_X_2,
        var_used_minus
      )
      
    }
    
    jacobian[, i] <- (f_plus - f_minus) / (2 * epsilon)
  }
  
  
  return(jacobian)
}

jacobian_central_difference <- function(moments,
                                        par_initial,
                                        type,
                                        f_function_name,
                                        Omega_multi_mean_X = NA,
                                        Omega_multi_mean_AX = NA,
                                        moments_2 = NA,
                                        Omega_multi_mean_X_2 = NA,
                                        Omega_multi_mean_AX_2 = NA,
                                        epsilon = 1e-5) {
  if (!is.character(f_function_name) ||
      length(f_function_name) != 1) {
    stop(
      "f_function_name must be 'moments_to_paras_func', 'moments_to_f_paras_func' , 'moments_to_alpha_func' or 'moments_to_alpha_split_sample_func' "
    )
  }
  if (f_function_name == "moments_to_alpha_func" &
      (any(is.na(Omega_multi_mean_X)) |
       any(is.na(Omega_multi_mean_AX)))) {
    stop(
      "'Omega_multi_mean_X' and 'Omega_multi_mean_AX' must be provided (non-NA) when using 'moments_to_alpha_func'"
    )
  }
  
  if (f_function_name == "moments_to_alpha_func") {
    m <- length(
      call_f_function(
        moments,
        par_initial,
        type,
        f_function_name,
        Omega_multi_mean_X,
        Omega_multi_mean_AX
      )
    )
    
    jacobian1 <- jacobian_central_difference_main(
      row_num = m,
      var_name = "moments",
      moments,
      par_initial,
      type,
      f_function_name,
      Omega_multi_mean_X,
      Omega_multi_mean_AX
    )
    jacobian2 <- jacobian_central_difference_main(
      row_num = m,
      var_name = "Omega_multi_mean_X",
      moments,
      par_initial,
      type,
      f_function_name,
      Omega_multi_mean_X,
      Omega_multi_mean_AX
    )
    jacobian3 <- jacobian_central_difference_main(
      row_num = m,
      var_name = "Omega_multi_mean_AX",
      moments,
      par_initial,
      type,
      f_function_name,
      Omega_multi_mean_X,
      Omega_multi_mean_AX
    )
    
    jacobian <- cbind(jacobian1, jacobian2, jacobian3)
    return(jacobian)
    
  } else if (f_function_name == "moments_to_paras_func" |
             f_function_name == "moments_to_f_paras_func") {
    m <- length(
      call_f_function(
        moments,
        par_initial,
        type,
        f_function_name,
        Omega_multi_mean_X,
        Omega_multi_mean_AX
      )
    )
    
    jacobian <- jacobian_central_difference_main(
      row_num = m,
      var_name = "moments",
      moments,
      par_initial,
      type,
      f_function_name,
      Omega_multi_mean_X,
      Omega_multi_mean_AX
    )
    return(jacobian)
    
  } else if (f_function_name == "moments_to_alpha_split_sample_func") {
    m <- length(
      call_f_function(
        moments,
        par_initial,
        type,
        f_function_name,
        Omega_multi_mean_X,
        Omega_multi_mean_AX,
        moments_2,
        Omega_multi_mean_X_2,
        Omega_multi_mean_AX_2
      )
    )
    
    jacobian1 <- jacobian_central_difference_main(
      row_num = m,
      var_name = "moments",
      moments,
      par_initial,
      type,
      f_function_name,
      Omega_multi_mean_X,
      Omega_multi_mean_AX,
      moments_2,
      Omega_multi_mean_X_2,
      Omega_multi_mean_AX_2
    )
    jacobian2 <- jacobian_central_difference_main(
      row_num = m,
      var_name = "Omega_multi_mean_X",
      moments,
      par_initial,
      type,
      f_function_name,
      Omega_multi_mean_X,
      Omega_multi_mean_AX,
      moments_2,
      Omega_multi_mean_X_2,
      Omega_multi_mean_AX_2
    )
    jacobian3 <- jacobian_central_difference_main(
      row_num = m,
      var_name = "Omega_multi_mean_AX",
      moments,
      par_initial,
      type,
      f_function_name,
      Omega_multi_mean_X,
      Omega_multi_mean_AX,
      moments_2,
      Omega_multi_mean_X_2,
      Omega_multi_mean_AX_2
    )
    jacobian4 <- jacobian_central_difference_main(
      row_num = m,
      var_name = "moments_2",
      moments,
      par_initial,
      type,
      f_function_name,
      Omega_multi_mean_X,
      Omega_multi_mean_AX,
      moments_2,
      Omega_multi_mean_X_2,
      Omega_multi_mean_AX_2
    )
    jacobian5 <- jacobian_central_difference_main(
      row_num = m,
      var_name = "Omega_multi_mean_X_2",
      moments,
      par_initial,
      type,
      f_function_name,
      Omega_multi_mean_X,
      Omega_multi_mean_AX,
      moments_2,
      Omega_multi_mean_X_2,
      Omega_multi_mean_AX_2
    )
    jacobian6 <- jacobian_central_difference_main(
      row_num = m,
      var_name = "Omega_multi_mean_AX_2",
      moments,
      par_initial,
      type,
      f_function_name,
      Omega_multi_mean_X,
      Omega_multi_mean_AX,
      moments_2,
      Omega_multi_mean_X_2,
      Omega_multi_mean_AX_2
    )
    
    
    jacobian <- cbind(jacobian1,
                      jacobian2,
                      jacobian3,
                      jacobian4,
                      jacobian5,
                      jacobian6)
    
    return(jacobian)
    
  } else if (f_function_name == "moments_to_paras_split_sample_func" |
             f_function_name == "moments_to_f_paras_split_sample_func") {
    m <- length(
      call_f_function(
        moments,
        par_initial,
        type,
        f_function_name,
        Omega_multi_mean_X,
        Omega_multi_mean_AX,
        moments_2,
        Omega_multi_mean_X_2,
        Omega_multi_mean_AX_2
      )
    )
    
    jacobian1 <- jacobian_central_difference_main(
      row_num = m,
      var_name = "moments",
      moments,
      par_initial,
      type,
      f_function_name,
      Omega_multi_mean_X,
      Omega_multi_mean_AX,
      moments_2,
      Omega_multi_mean_X_2,
      Omega_multi_mean_AX_2
    )
    jacobian4 <- jacobian_central_difference_main(
      row_num = m,
      var_name = "moments_2",
      moments,
      par_initial,
      type,
      f_function_name,
      Omega_multi_mean_X,
      Omega_multi_mean_AX,
      moments_2,
      Omega_multi_mean_X_2,
      Omega_multi_mean_AX_2
    )
    
    jacobian <- cbind(jacobian1, jacobian4)
    
    return(jacobian)
    
  }
  
}
# 计算一阶泰勒展开
taylor_approximation_1st <- function(moments,
                                     par_initial,
                                     type,
                                     f_function_name,
                                     delta_x,
                                     Omega_multi_mean_X = NA,
                                     Omega_multi_mean_AX = NA,
                                     moments_2 = NA,
                                     Omega_multi_mean_X_2 = NA,
                                     Omega_multi_mean_AX_2 = NA) {
  jacobian <- jacobian_central_difference(
    moments,
    par_initial,
    type,
    f_function_name,
    Omega_multi_mean_X,
    Omega_multi_mean_AX,
    moments_2,
    Omega_multi_mean_X_2,
    Omega_multi_mean_AX_2
  )
  
  f_x <- call_f_function(
    moments,
    par_initial,
    type,
    f_function_name,
    Omega_multi_mean_X,
    Omega_multi_mean_AX,
    moments_2,
    Omega_multi_mean_X_2,
    Omega_multi_mean_AX_2
  )
  f_x_plus_delta_x <- f_x + as.numeric(jacobian %*% delta_x)
  
  return(f_x_plus_delta_x)
}

cov_bootstrap_U_stats <- function(A_1, A_2, B_1, B_2) {
  term_1 <- cov(A_1, B_1)
  term_2 <- cov(A_2, B_2)
  final <- term_1 - 2 * term_2
  
  return(final)
  
}
# unknown sigma -----
split_matrix <- function(mat) {
  half_mat <- mat / 2
  cbind(half_mat, half_mat)
}
calculate_covariance_matrix <- function(X,
                                        A,
                                        H_S,
                                        B_bootstrap,
                                        weights_matrix_half,
                                        inv_sigma_split) {
  n_half <- length(A)
  p <- ncol(X)
  cov_em_for_split_sample_4 <- matrix(NA_real_, 4 + 2 * p, 4 + 2 * p)
  vector_one_half <- rep(1, n_half)
  
  # Calculate bootstrapped means and statistics
  mA_bootstrap <- sapply(1:B_bootstrap, function(i)
    mean(A * weights_matrix_half[i, ]))
  
  calculate_m <- function(vec1, vec2, weights) {
    sum((vec1 * weights) * eigenMapMatMult(H_S, vec2 * weights)) / (n_half * (n_half - 1))
  }
  
  mX_2_1 <- sapply(1:B_bootstrap, function(i)
    calculate_m(vector_one_half, vector_one_half, weights_matrix_half[i, ]))
  mXA_1_1 <- sapply(1:B_bootstrap, function(i)
    calculate_m(A, vector_one_half, weights_matrix_half[i, ]))
  mXA_2_1 <- sapply(1:B_bootstrap, function(i)
    calculate_m(A, A, weights_matrix_half[i, ]))
  
  mX_2_2 <- sapply(1:B_bootstrap, function(i)
    calculate_m(vector_one_half, vector_one_half, weights_matrix_half[i, ] - 1))
  mXA_1_2 <- sapply(1:B_bootstrap, function(i)
    calculate_m(A, vector_one_half, weights_matrix_half[i, ] - 1))
  mXA_2_2 <- sapply(1:B_bootstrap, function(i)
    calculate_m(A, A, weights_matrix_half[i, ] - 1))
  
  t_X <- t(X)
  Omega_multi_mean_X <- t(sapply(1:B_bootstrap, function(i)
    as.numeric(eigenMapMatMult(
      inv_sigma_split,
      eigenMapMatMult(t_X, weights_matrix_half[i, ])
    ) / n_half)))
  Omega_multi_meanX <- t(sapply(1:B_bootstrap, function(i)
    as.numeric(eigenMapMatMult(
      inv_sigma_split,
      eigenMapMatMult(t_X, A * weights_matrix_half[i, ])
    ) / n_half)))
  
  m1_U1 <- mA_bootstrap
  m2_U2_1 <- cbind(mX_2_1, mXA_1_1, mXA_2_1)
  m2_U2_2 <- cbind(mX_2_2, mXA_1_2, mXA_2_2)
  m3_U1 <- cbind(Omega_multi_mean_X, Omega_multi_meanX)
  
  # Calculate covariance matrix
  cov_em_for_split_sample_4[1, 1] <- var(m1_U1)
  cov_em_for_split_sample_4[1, 2:4] <- cov(m1_U1, m2_U2_1)
  cov_em_for_split_sample_4[2:4, 1] <- cov(m2_U2_1, m1_U1)
  cov_em_for_split_sample_4[2:4, 2:4] <- cov_bootstrap_U_stats(m2_U2_1, m2_U2_2, m2_U2_1, m2_U2_2)
  cov_em_for_split_sample_4[1, (5):(4 + 2 * p)] <- cov(m1_U1, m3_U1)
  cov_em_for_split_sample_4[(5):(4 + 2 * p), 1] <- cov(m3_U1, m1_U1)
  cov_em_for_split_sample_4[2:4, (5):(4 + 2 * p)] <- cov(m2_U2_1, m3_U1)
  cov_em_for_split_sample_4[(5):(4 + 2 * p), 2:4] <- cov(m3_U1, m2_U2_1)
  cov_em_for_split_sample_4[(5):(4 + 2 * p), (5):(4 + 2 * p)] <-
    eigenMapMatMult(t(m3_U1), m3_U1) / (B_bootstrap - 1) -
    eigenMapMatMult(colMeans(m3_U1), t(colMeans(m3_U1))) * B_bootstrap / (B_bootstrap - 1)
  
  return(
    list(
      cov_em_for_split_sample_4 = cov_em_for_split_sample_4,
      m1_U1 = m1_U1,
      m2_U2_1 = m2_U2_1,
      m2_U2_2 = m2_U2_2,
      m3_U1 = m3_U1
    )
  )
}
calculate_estimate_cov <- function(jacobian_main,
                                   cov_em_moments_a,
                                   cov_em_moments_ab,
                                   cov_em_moments_b) {
  t_jacobian_main <- t(jacobian_main)
  term_1 <- eigenMapMatMult(eigenMapMatMult(jacobian_main, cov_em_moments_a),
                            t_jacobian_main)
  term_2 <- eigenMapMatMult(eigenMapMatMult(jacobian_main, cov_em_moments_ab),
                            t_jacobian_main)
  term_4 <- eigenMapMatMult(eigenMapMatMult(jacobian_main, cov_em_moments_b),
                            t_jacobian_main)
  estimate_cov <- term_1 + term_2 + t(term_2) + term_4
}
glm_mom_unknown_sigma <- function(X, A, B_bootstrap = 1000, par_initial) {
  n <- nrow(X)
  p <- ncol(X)
  
  n_half <- n / 2
  
  X_1 <- X[1:n_half, ]
  X_2 <- X[(n_half + 1):n, ]
  A_1 <- A[1:n_half]
  A_2 <- A[(n_half + 1):n]
  
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
  
  ## boostrap for split sample ------
  
  # Initialize weight matrix
  weights_matrix_half <- matrix(0, nrow = B_bootstrap, ncol = n_half)
  for (i in 1:B_bootstrap) {
    weights_matrix_half[i, ] <- rmultinom(1, n_half, rep(1 / n_half, n_half))
  }
  cov_em_list <- calculate_covariance_matrix(X_2,
                                             A_2,
                                             H_S_a,
                                             B_bootstrap,
                                             weights_matrix_half,
                                             inv_sigma_split_a)
  cov_em_for_split_sample_4_a <- cov_em_list$cov_em_for_split_sample_4
  m1_U1_a <-  cov_em_list$m1_U1
  m2_U2_a_1 <-  cov_em_list$m2_U2_1
  m2_U2_a_2 <-  cov_em_list$m2_U2_2
  m3_U1_a <-  cov_em_list$m3_U1
  
  cov_em_list <- calculate_covariance_matrix(X_1,
                                             A_1,
                                             H_S_b,
                                             B_bootstrap,
                                             weights_matrix_half,
                                             inv_sigma_split_b)
  cov_em_for_split_sample_4_b <- cov_em_list$cov_em_for_split_sample_4
  m1_U1_b <-  cov_em_list$m1_U1
  m2_U2_b_1 <-  cov_em_list$m2_U2_1
  m2_U2_b_2 <-  cov_em_list$m2_U2_2
  m3_U1_b <-  cov_em_list$m3_U1
  #### ab-----
  cov_em_for_split_sample_4_ab <- matrix(NA_real_, nrow = 4 + 2 * p, ncol = 4 +
                                           2 * p)
  
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
  var_moments_split_sample_4 <- c(diag(cov_em_for_split_sample_4_a),
                                  diag(cov_em_for_split_sample_4_b))
  
  cov_alpha_split_sample_4 <-  calculate_estimate_cov(
    jacobian_em_S_alpha_4 / 2,
    cov_em_for_split_sample_4_a,
    cov_em_for_split_sample_4_ab,
    cov_em_for_split_sample_4_b
  )
  var_alpha_split_sample_4 <- diag(cov_alpha_split_sample_4)
  
  cov_paras_split_sample_4 <- calculate_estimate_cov(
    jacobian_em_S_paras_4 / 2,
    cov_em_for_split_sample_4_a[c((1):(4)), c((1):(4))],
    cov_em_for_split_sample_4_ab[c((1):(4)), c((1):(4))],
    cov_em_for_split_sample_4_b[c((1):(4)), c((1):(4))]
  )
  var_paras_split_sample_4 <- diag(cov_paras_split_sample_4)
  
  cov_f_paras_split_sample_4 <- calculate_estimate_cov(
    jacobian_em_S_f_paras_4 / 2,
    cov_em_for_split_sample_4_a[c((1):(4)), c((1):(4))],
    cov_em_for_split_sample_4_ab[c((1):(4)), c((1):(4))],
    cov_em_for_split_sample_4_b[c((1):(4)), c((1):(4))]
  )
  var_f_paras_split_sample_4 <- diag(cov_f_paras_split_sample_4)
  
  cov_alpha_split_sample_1 <-  calculate_estimate_cov(
    jacobian_em_S_alpha_1 / 2,
    cov_em_for_split_sample_4_a[c((4):(4 + 2 * p)), c((4):(4 + 2 * p))],
    cov_em_for_split_sample_4_ab[c((4):(4 + 2 * p)), c((4):(4 + 2 * p))],
    cov_em_for_split_sample_4_b[c((4):(4 + 2 * p)), c((4):(4 + 2 * p))]
  )
  var_alpha_split_sample_1 <- diag(cov_alpha_split_sample_1)
  
  cov_paras_split_sample_1 <- calculate_estimate_cov(
    jacobian_em_S_paras_1 / 2,
    cov_em_for_split_sample_4_a[4, 4],
    cov_em_for_split_sample_4_ab[4, 4],
    cov_em_for_split_sample_4_b[4, 4]
  )
  var_paras_split_sample_1 <- diag(cov_paras_split_sample_1)
  
  cov_f_paras_split_sample_1 <- calculate_estimate_cov(
    jacobian_em_S_f_paras_1 / 2,
    cov_em_for_split_sample_4_a[4, 4],
    cov_em_for_split_sample_4_ab[4, 4],
    cov_em_for_split_sample_4_b[4, 4]
  )
  var_f_paras_split_sample_1 <- diag(cov_f_paras_split_sample_1)
  ### return ----
  return(
    list(
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
      
    )
  )
}

# 保存当前环境中的所有对象
save(list = ls(), file = "demo_glm_MoM/all_functions_and_objects_mom.RData")
