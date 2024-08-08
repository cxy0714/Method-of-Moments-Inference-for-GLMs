
# function ----

##########################################################
# socre and variance function for logistic model in bellec's method
##########################################################
logit_variance <- function(u) {
  exp_u <- exp(u)
  exp_u / (1 + exp_u) ^ 2
}
score_func <- function(u, y) {
  sigmoid <- exp(u) / (1 + exp(u))
  sigmoid - y
}

##########################################################
# Main Function for Implementing Bellec's Method
# This implementation assumes Sigma_X = (Identity / omega_11).
# The code is based on Sections 3, 4.1 (Theorem 4.1), and 6.3 of Bellec's paper.
# Parameter naming conventions follow Section 3: Derivatives and Observable Adjustments.
#
# When omega_11 == p:
#   This setup corresponds to Section 6.3 of Bellec's paper, where a simplified formula is provided.
#   Here, alpha_de_p2 is the unbiased estimate of alpha, computed by dividing the initial
#   Ridge regression alpha by the factor given in formula (6.8).
#
# When omega_11 == 1:
#   The implementation is based on formulas from Sections 3 and 4.1, with parameters suffixed by _Ge (Generalized).
#   - alpha_de_p2 represents the preliminary unbiased version as per formula (3.9) in Section 3.
#   - alpha_de_N_p2 is the final unbiased version, obtained by dividing alpha_de_p2
#     by the factor specified in Bellec's Theorem 4.1.
#
# Z_test is always the parameter to test Bellec's normal resutls formula (3.9) or formula (6.8)
#
# The approach uses two samples, split to calculate the quadratic form of alpha and compare it with our methods, which will be reported later. 
# Sample 1 is suffixed with _T, and Sample 2 with _H  (Tail and Head, well, a quirky nod—puns intended!)
#
# PS: Utilizing Cholesky decomposition to invert a positive definite matrix (here, A_hat) can be (nearly three times?) faster than using SVD.
#
##########################################################

run.replicate <- function(n,
                          p,
                          alpha = alpha,
                          omega_11,
                          Is_Rad,
                          Is_sparse,
                          Is_sparse_only_1,
                          mu_x,
                          lambda_value,
                          seed,
                          filename.tracking) {
  set.seed(seed)
  cat(n, seed, "\n", file = filename.tracking, append = TRUE)
  
  result <- tryCatch({
   
    s <- sum(abs(alpha) > 0)
    if (Is_Rad) {
      X_T <- matrix(rbinom(n * p, size = 1, prob = 0.5),
                    nrow = n,
                    ncol = p)
      
      X_T <- (2 * X_T - 1) / sqrt(omega_11)
      
      X_H <- matrix(rbinom(n * p, size = 1, prob = 0.5),
                    nrow = n,
                    ncol = p)
      
      X_H <- (2 * X_H - 1) / sqrt(omega_11)
    } else  {
      X_T <- matrix(rnorm(n * p, mean = 0, sd = 1 / sqrt(omega_11)),
                    nrow = n,
                    ncol = p)
      X_H <- matrix(rnorm(n * p, mean = 0, sd = 1 / sqrt(omega_11)),
                    nrow = n,
                    ncol = p)
    }
    X_T <- X_T + mu_x
    X_H <- X_H + mu_x
    
    z_T <- X_T %*% alpha
    A_T <- rbinom(n, 1, 1 / (1 + exp(-z_T)))
    
    Z_H <- X_H %*% alpha
    A_H <- rbinom(n, 1, 1 / (1 + exp(-Z_H)))
    
    
    indices_first_1 <- 1:10
    indices_middle_2 <- 100:109
    indices_selected <- c(indices_first_1, indices_middle_2)
    
    ## ridge penalized bellec------once------------
    fit_T <- glmnet(
      X_T,
      A_T,
      family = "binomial",
      alpha = 0,
      lambda = lambda_value,
      maxit = 10000,
      intercept = FALSE
    )
    X_T_trans <- t(X_T)
    alpha_p2_T <- fit_T$beta
    lambda_true <- fit_T$lambda
    if (omega_11 == p) {
      u_T <- X_T %*% alpha_p2_T
      D_T <- sapply(c(1:length(lambda_true)), function(i)
        (logit_variance(u_T[, i])))
      score_T <-  sapply(c(1:length(lambda_true)), function(i)
        (-score_func(u_T[, i], A_T)))
      r_T <-  colMeans(score_T ^ 2) ^ (1 / 2)
      X_score_T <-  t(X_T) %*% score_T
      lambda_true <- colMeans(X_score_T / alpha_p2_T) * omega_11 / n
      
      V_T <- rep(NA, length(lambda_true))
      
      for (i in 1:length(lambda_true)) {
        XDX_T <- eigenMapMatMult(X_T_trans, eigenMapMatMult(diag(D_T[, i]), X_T))
        
        L <- chol(XDX_T + diag(lambda_true[i] * n / p, p))
        
        A_hat_T  <- chol2inv(L)
        
        V_T[i] <-  sum(D_T[, i]) / n - sum(diag(eigenMapMatMult(
          A_hat_T,
          eigenMapMatMult(X_T_trans, eigenMapMatMult(diag(D_T[, i] ^ 2), X_T))
        ))) / n
        
      }
      sign_p2_T <- sapply(c(1:length(lambda_true)), function(i)
        (sign(crossprod(
          alpha, alpha_p2_T[, i]
        ))))
      t_square_T <- sapply(c(1:length(lambda_true)), function(i)
        ((V_T[i] + lambda_true[i]) ^ 2 * sum(alpha_p2_T[, i] ^ 2) / p - (p / n) * r_T[i] ^
           2
        ))
      Z_test_p2_T  <- sapply(1:length(lambda_true), function(i) {
        if (t_square_T[i] > 0) {
          (sqrt(n / p) / r_T[i]) * (as.vector((V_T[i] + lambda_true[i])) * alpha_p2_T[, i] - alpha * as.vector(sign_p2_T[i] * sqrt(t_square_T[i])))
        } else {
          rep(NA, length(alpha_p2_T[, i]))
        }
      })
      alpha_de_p2_T <- sapply(1:length(lambda_true), function(i) {
        if (t_square_T[i] > 0) {
          (V_T[i] + lambda_true[i]) / (sign_p2_T[i] * sqrt(t_square_T[i])) * alpha_p2_T[, i]
        } else {
          rep(NA, length(alpha_p2_T[, i]))
        }
      })
    } else {
      alpha_p2_T <- fit_T$beta
      lambda_true <- fit_T$lambda
      
      u_T <- X_T %*% alpha_p2_T
      D_T <- sapply(c(1:length(lambda_true)), function(i)
        (logit_variance(u_T[, i])))
      score_T <-  sapply(c(1:length(lambda_true)), function(i)
        (-score_func(u_T[, i], A_T)))
      r_T <-  colMeans(score_T ^ 2) ^ (1 / 2)
      X_score_T <-  eigenMapMatMult(X_T_trans, score_T)
      lambda_true <- colMeans(X_score_T / alpha_p2_T) / n
      
      
      V_T <- df_T <- gamma_T <-  rep(NA, length(lambda_true))
      for (i in 1:length(lambda_true)) {
        XDX_T <- eigenMapMatMult(X_T_trans, eigenMapMatMult(diag(D_T[, i]), X_T))
        
        L <- chol(XDX_T + diag(lambda_true[i] * n, p))
        
        A_hat_T  <- chol2inv(L)
        
        V_T[i] <-  sum(D_T[, i]) / n - sum(diag(eigenMapMatMult(
          A_hat_T,
          eigenMapMatMult(X_T_trans, eigenMapMatMult(diag(D_T[, i] ^ 2), X_T))
        ))) / n
        df_T[i] <- sum(diag(eigenMapMatMult(A_hat_T, XDX_T)))
        gamma_T[i] <- df_T[i] / (n * V_T[i])
        
      }
      t_square_T_Ge <- sqrt(omega_11) / n ^ 2 * colSums(X_score_T ^ 2)  +
        2 * V_T / n * colSums(alpha_p2_T * X_score_T) +
        V_T ^ 2 / n *  colSums((u_T - sweep(score_T, 2, gamma_T, "*")) ^
                                 2) -
        (p / n) * r_T ^ 2
      
      sign_p2_T_Ge <- sapply(c(1:length(lambda_true)), function(i)
        (sign(
          crossprod(alpha, V_T[i] * n / omega_11 * alpha_p2_T[, i] + X_score_T[, i])
        )))
      alpha_de_p2_T_Ge <- alpha_p2_T +  sweep(X_score_T, 2, (V_T * n) ^ -1 * omega_11, "*")
      Z_test_p2_T_Ge <-  sapply(1:length(lambda_true), function(i) {
        if (t_square_T_Ge[i] > 0) {
          (sqrt(n) / sqrt(omega_11)) * (
            V_T[i] / r_T[i] * alpha_de_p2_T_Ge[, i] -
              alpha * sign_p2_T_Ge[i] *  t_square_T_Ge[i] ^ (1 / 2) / r_T[i]
          )
        } else {
          rep(NA, p)
        }
      })
      
      de_p2_factor_T <- V_T / (sign_p2_T_Ge *  t_square_T_Ge ^ (1 / 2))
      alpha_de_N_p2_T_Ge <-   sapply(1:length(lambda_true), function(i) {
        if (t_square_T_Ge[i] > 0) {
          alpha_de_p2_T_Ge[, i] * de_p2_factor_T[i]
        } else {
          rep(NA, p)
        }
      })
      
    }
    
    ## ridge penalized bellec------twice------------
    fit_H <- glmnet(
      X_H,
      A_H,
      family = "binomial",
      alpha = 0,
      lambda = lambda_value,
      maxit = 10000,
      intercept = FALSE
    )
    
    X_H_trans <- t(X_H)
    
    alpha_p2_H <- fit_H$beta
    lambda_true <- fit_H$lambda
    if (omega_11 == p) {
      u_H <- X_H %*% alpha_p2_H
      D_H <- sapply(c(1:length(lambda_true)), function(i)
        (logit_variance(u_H[, i])))
      score_H <-  sapply(c(1:length(lambda_true)), function(i)
        (-score_func(u_H[, i], A_H)))
      r_H <-  colMeans(score_H ^ 2) ^ (1 / 2)
      X_score_H <-  t(X_H) %*% score_H
      lambda_true <- colMeans(X_score_H / alpha_p2_H) * omega_11 / n
      
      V_H <- rep(NA, length(lambda_true))
      for (i in 1:length(lambda_true)) {
        XDX_H <- eigenMapMatMult(X_H_trans, eigenMapMatMult(diag(D_H[, i]), X_H))
        
        L <- chol(XDX_H + diag(lambda_true[i] * n / p, p))
        
        A_hat_H  <- chol2inv(L)
        
        V_H[i] <-  sum(D_H[, i]) / n - sum(diag(eigenMapMatMult(
          A_hat_H,
          eigenMapMatMult(X_H_trans, eigenMapMatMult(diag(D_H[, i] ^ 2), X_H))
        ))) / n
        
      }
      sign_p2_H <- sapply(c(1:length(lambda_true)), function(i)
        (sign(crossprod(
          alpha, alpha_p2_H[, i]
        ))))
      t_square_H <- sapply(c(1:length(lambda_true)), function(i)
        ((V_H[i] + lambda_true[i]) ^ 2 * sum(alpha_p2_H[, i] ^ 2) / p - (p / n) * r_H[i] ^
           2
        ))
      Z_test_p2_H  <- sapply(1:length(lambda_true), function(i) {
        if (t_square_H[i] > 0) {
          (sqrt(n / p) / r_H[i]) * (as.vector((V_H[i] + lambda_true[i])) * alpha_p2_H[, i] - alpha * as.vector(sign_p2_H[i] * sqrt(t_square_H[i])))
        } else {
          rep(NA, length(alpha_p2_H[, i]))
        }
      })
      alpha_de_p2_H <- sapply(1:length(lambda_true), function(i) {
        if (t_square_H[i] > 0) {
          (V_H[i] + lambda_true[i]) / (sign_p2_H[i] * sqrt(t_square_H[i])) * alpha_p2_H[, i]
        } else {
          rep(NA, length(alpha_p2_H[, i]))
        }
      })
    } else {
      u_H <- X_H %*% alpha_p2_H
      D_H <- sapply(c(1:length(lambda_true)), function(i)
        (logit_variance(u_H[, i])))
      score_H <-  sapply(c(1:length(lambda_true)), function(i)
        (-score_func(u_H[, i], A_H)))
      r_H <-  colMeans(score_H ^ 2) ^ (1 / 2)
      X_score_H <-  eigenMapMatMult(X_H_trans, score_H)
      lambda_true <- colMeans(X_score_H / alpha_p2_H) / n
      
      V_H <- df_H <- gamma_H <-  rep(NA, length(lambda_true))
      for (i in 1:length(lambda_true)) {
        XDX_H <- eigenMapMatMult(X_H_trans, eigenMapMatMult(diag(D_H[, i]), X_H))
        
        L <- chol(XDX_H + diag(lambda_true[i] * n, p))
        
        A_hat_H  <- chol2inv(L)
        
        V_H[i] <-  sum(D_H[, i]) / n - sum(diag(eigenMapMatMult(
          A_hat_H,
          eigenMapMatMult(X_H_trans, eigenMapMatMult(diag(D_H[, i] ^ 2), X_H))
        ))) / n
        df_H[i] <- sum(diag(eigenMapMatMult(A_hat_H, XDX_H)))
        gamma_H[i] <- df_H[i] / (n * V_H[i])
        
      }
      t_square_H_Ge <- sqrt(omega_11) / n ^ 2 * colSums(X_score_H ^ 2)  +
        2 * V_H / n * colSums(alpha_p2_H * X_score_H) +
        V_H ^ 2 / n *  colSums((u_H - sweep(score_H, 2, gamma_H, "*")) ^
                                 2) -
        (p / n) * r_H ^ 2
      
      sign_p2_H_Ge <- sapply(c(1:length(lambda_true)), function(i)
        (sign(
          crossprod(alpha, V_H[i] * n / omega_11 * alpha_p2_H[, i] + X_score_H[, i])
        )))
      alpha_de_p2_H_Ge <- alpha_p2_H +  sweep(X_score_H, 2, (V_H * n) ^ -1 * omega_11, "*")
      Z_test_p2_H_Ge <-  sapply(1:length(lambda_true), function(i) {
        if (t_square_H_Ge[i] > 0) {
          (sqrt(n) / sqrt(omega_11)) * (
            V_H[i] / r_H[i] * alpha_de_p2_H_Ge[, i] -
              alpha * sign_p2_H_Ge[i] *  t_square_H_Ge[i] ^ (1 / 2) / r_H[i]
          )
        } else {
          rep(NA, p)
        }
      })
      
      de_p2_factor_H <- V_H / (sign_p2_H_Ge *  t_square_H_Ge ^ (1 / 2))
      alpha_de_N_p2_H_Ge <-   sapply(1:length(lambda_true), function(i) {
        if (t_square_H_Ge[i] > 0) {
          alpha_de_p2_H_Ge[, i] * de_p2_factor_H[i]
        } else {
          rep(NA, p)
        }
      })
      
    }
    ## combined----
    if (omega_11 == p) {
      alpha_L2_de_p2_TT <- sapply(c(1:length(lambda_true)), function(i)
        (sum(alpha_de_p2_T[, i] * alpha_de_p2_H[, i])))
      t_square_TT <- rbind(t_square_T, t_square_H)
      
      alpha_de_p2_TT <- list(alpha_de_p2_T[indices_selected, ], alpha_de_p2_H[indices_selected, ])
      Z_test_p2_TT <- list(Z_test_p2_T[indices_selected, ], Z_test_p2_T[indices_selected, ])
      
      return(
        list(
          success = TRUE,
          seed = seed,
          n = n,
          p = p,
          alpha = alpha,
          alpha_L2_de_p2_TT = alpha_L2_de_p2_TT,
          t_square_TT  =  t_square_TT,
          alpha_de_p2_TT = alpha_de_p2_TT,
          Z_test_p2_TT = Z_test_p2_TT
        )
      )
    } else {
      alpha_L2_de_N_Ge <- sapply(c(1:length(lambda_true)), function(i)
        (sum(
          alpha_de_N_p2_H_Ge[, i] * alpha_de_N_p2_T_Ge[, i]
        )))
      alpha_L2_de_Ge <- sapply(c(1:length(lambda_true)), function(i)
        (sum(
          alpha_de_p2_H_Ge[, i] * alpha_de_p2_T_Ge[, i]
        )))
      t_square_TT_Ge <- rbind(t_square_T_Ge, t_square_H_Ge)
      de_p2_factor_TT <- rbind(de_p2_factor_T, de_p2_factor_H)
      alpha_de_p2_TT_Ge <- list(alpha_de_p2_T_Ge[indices_selected, ], alpha_de_p2_H_Ge[indices_selected, ])
      alpha_de_N_p2_TT_Ge <- list(alpha_de_N_p2_T_Ge[indices_selected, ], alpha_de_N_p2_H_Ge[indices_selected, ])
      Z_test_p2_TT_Ge <- list(Z_test_p2_T_Ge[indices_selected, ], Z_test_p2_T_Ge[indices_selected, ])
      
      return(
        list(
          success = TRUE,
          seed = seed,
          n = n,
          p = p,
          alpha = alpha,
          alpha_L2_de_N_Ge = alpha_L2_de_N_Ge,
          alpha_L2_de_Ge = alpha_L2_de_Ge,
          
          t_square_TT_Ge  =  t_square_TT_Ge,
          de_p2_factor_TT = de_p2_factor_TT,
          alpha_de_p2_TT_Ge = alpha_de_p2_TT_Ge,
          alpha_de_N_p2_TT_Ge = alpha_de_N_p2_TT_Ge,
          Z_test_p2_TT_Ge = Z_test_p2_TT_Ge
        )
      )
    }
    
  }, error = function(e) {
    # if error, ruturn message and get off
    cat("Error: ", e$message, "\n")
    return(list(success = FALSE, message = e$message, seed = seed,
                n = n,
                p = p))
  })
  
  return(result)
}