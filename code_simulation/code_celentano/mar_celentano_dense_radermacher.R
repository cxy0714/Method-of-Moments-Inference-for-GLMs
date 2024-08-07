# run from the simulations directory

library(glmnet)
library(boot)
library(parallel)

# Define functions we need
d.inv.logit <- function(eta) {
  return(
    exp(eta)/(1+exp(eta))^2
  )
}
dd.inv.logit <- function(eta) {
  return(
    exp(eta)*(exp(eta)-1)/(exp(eta)+1)^3
  )
}
prop.score <- function(eta,epsilon) {
  return(
    epsilon + (1-epsilon) * inv.logit(eta)
  )
}
d.prop.score <- function(eta,epsilon) {
  return(
    (1-epsilon) * d.inv.logit(eta)
  )
}
dd.prop.score <- function(eta,epsilon) {
  return(
    (1-epsilon) * dd.inv.logit(eta)
  )
}
properr <- function(eta.d,eta.d.hat,epsilon) {
  return(
    1 - prop.score(eta.d,epsilon) / prop.score(eta.d.hat,epsilon)
  )
}
d.properr.d.etad <- function(eta.d,eta.d.hat,epsilon) {
  return(
    -d.prop.score(eta.d,epsilon) / prop.score(eta.d.hat,epsilon)
  )
}
d.properr.d.etadhat <- function(eta.d,eta.d.hat,epsilon) {
  return(
    prop.score(eta.d,epsilon) * d.prop.score(eta.d.hat,epsilon) / prop.score(eta.d.hat,epsilon)^2
  )
}
dd.properr.dd.etad <- function(eta.d,eta.d.hat,epsilon) {
  return(
    -dd.prop.score(eta.d,epsilon) / prop.score(eta.d.hat,epsilon)
  )
}
dd.properr.d.etad.d.etadhat <- function(eta.d,eta.d.hat,epsilon) {
  return(
    d.prop.score(eta.d,epsilon)*d.prop.score(eta.d.hat,epsilon) / prop.score(eta.d.hat,epsilon)^2
  )
}
dd.properr.dd.etadhat <- function(eta.d,eta.d.hat,epsilon) {
  return(
    prop.score(eta.d,epsilon) * dd.prop.score(eta.d.hat,epsilon) / prop.score(eta.d.hat,epsilon)^2 -
      2 * prop.score(eta.d,epsilon) * d.prop.score(eta.d.hat,epsilon)^2 / prop.score(eta.d.hat,epsilon)^3
  )
}
sqrtm <- function(K) {
  SVD <- svd(K)
  return(
    SVD$u %*% diag(sqrt(SVD$d)) %*% t(SVD$v)
  )
}
# Propensity score link and its derivatives
p.score <- function(eta) {
  epsilon + (1-epsilon) * inv.logit(eta)
}
d.p.score <- function(eta) {
  (1-epsilon)*d.inv.logit(eta)
}
dd.p.score <- function(eta) {
  (1-epsilon)*dd.inv.logit(eta)
}

# Expected value of a function under a guassian
gaussian.expectation <- function(
    f,
    mu,
    std.dev,
    width
) {
  if (std.dev < .01) {return(f(mu))}
  integrate(
    function(eta) {
      dnorm(eta,mean=mu,sd=std.dev)*f(eta)
    },
    lower = mu - width*std.dev , upper = mu + width * std.dev
  )[[1]]
}

# Binary search for root of an increasing function f
binary.solve <- function(f,x.min,x.max,tol) {
  total.iter = 0
  
  while (f(x.max) <= 0 & total.iter < 1e4) { # initialize max
    x.max = 2 * x.max
    total.iter = total.iter + 1
  }
  while ( x.max - x.min > tol & total.iter < 1e4) {
    x.temp = ( x.min + x.max ) / 2
    total.iter = total.iter + 1
    if ( f(x.temp) > 0 ) {
      x.max = x.temp
    } else {
      x.min = x.temp
    }
  }
  return((x.min + x.max)/2)
}

# Estimate mu.d and gamma.d
compute.prop.mu.gamma <- function(pi.hat,gamma.dstar.hat,
                                  gamma.reg,max.iter) {
  
  # initialize
  mu.d.hat = 0
  gamma.d.hat = gamma.dstar.hat
  
  for (i in 1:max.iter) {
    # Alternatinv coordinate optimization
    mu.d.hat.next <- binary.solve(
      function(mu.d) {
        gaussian.expectation(f=p.score,mu=mu.d,std.dev=gamma.d.hat,width=10) - pi.hat
      },
      x.min=-10,x.max=10,tol=1e-6
    )
    gamma.d.hat.next <- binary.solve(
      function(gamma.d) {
        gamma.d * gaussian.expectation(f=d.p.score,mu=mu.d.hat.next,std.dev=gamma.d,width=10) - pi.hat*gamma.dstar.hat + gamma.reg * gamma.d # Add some regularization to guarantee existence of solution
      },
      x.min=.01,x.max=1,tol=1e-6
    )
    if (abs(mu.d.hat.next - mu.d.hat)+abs(gamma.d.hat.next - gamma.d.hat) < 1e-6) {
      return(list(mu.d.hat = mu.d.hat.next,gamma.d.hat = gamma.d.hat.next))
    } else {
      mu.d.hat = mu.d.hat.next
      gamma.d.hat = gamma.d.hat.next
    }
  }
  return(list(mu.d.hat = mu.d.hat.next,gamma.d.hat = gamma.d.hat.next))
}

# Functions for computing zetas (specific to ridge)
zeta.theta.rootfnc <- function(zeta.theta,weights,lambda,p,n) {
  zeta.eta = (p/n) / (zeta.theta + lambda)
  return(
    zeta.theta - sum(weights/(zeta.eta*weights+1))/n
  )
}

run.replicate <- function(
    n,
    p,
    sigma,
    Is_Rad,
    theta.d.0,
    theta.y.0,
    theta.d,
    theta.y,
    mu.x,
    epsilon,
    lambda.y.seq,
    seed,
    filename.out,
    filename.tracking
) {
  
  set.seed(seed)
  cat(n,seed,"\n",file=filename.tracking,append=TRUE)
  
  mu.y = theta.y.0 + t(mu.x) %*% theta.y
  mu.d = theta.d.0 + t(mu.x) %*% theta.d
  gamma.y = sqrt(sum(theta.y^2))
  gamma.d = sqrt(sum(theta.d^2))
  pi.bar = gaussian.expectation(p.score,mu.d,gamma.d,width=10)
  alpha.1 = gaussian.expectation(d.p.score,mu.d,gamma.d,width=10) / pi.bar
  alpha.2 = gaussian.expectation(dd.p.score,mu.d,gamma.d,width=10) / pi.bar
  mu.x.cfd = mu.x + alpha.1 * theta.d
  
  # Generate data
  if (Is_Rad){
    X <- matrix(rbinom(n * p, size = 1, prob = 0.5),
                nrow = n,
                ncol = p)
    X <- (2 * X - 1)
  } else {
    X = matrix(rnorm(n*p),nrow=n,ncol=p)
  }
  z.y = .2*rnorm(n)
  y = theta.y.0 + X %*% theta.y + z.y  
  d = rbinom(n,1,p.score(theta.d.0 + X %*% theta.d))
  
  # Some random variables defined by the data
  n.1 = sum(d)
  eta.y = theta.y.0 + X %*% theta.y
  eta.d = theta.d.0 + X %*% theta.d
  w.ipw = 1/p.score(eta.d)
  
  # Some random variables defined by the data
  data.df <- data.frame()
  
  model.y = glmnet( # beware that glmnet scales by sum(d) rather than n
    X,y,
    weights=d,
    family = "gaussian",
    lambda = lambda.y.seq * n / n.1,
    alpha = 0, # ridge
    standardize = FALSE
  )
  model.y.ipw = glmnet( # beware that glmnet scales by sum(d) rather than n
    X,y,
    weights=d*w.ipw,
    lambda = lambda.y.seq * n / sum(d*w.ipw),
    family = "gaussian",
    alpha = 0, # ridge
    standardize = FALSE
  )
  
  for (lambda.y in lambda.y.seq) {
    
    # Fit outcome model, unweighted
    lambda.y.glmnet = lambda.y*n/n.1 # adjust because glmnet normalizes loss by n.1 rather than n
    theta.y.hat = as.matrix(predict(model.y,
                                    s=lambda.y.glmnet,
                                    weights=d,
                                    type="coefficients",
                                    exact=TRUE,
                                    x=X,y=y))
    theta.y.0.hat = theta.y.hat[1]
    theta.y.hat = theta.y.hat[2:(p+1)]
    eta.y.hat = theta.y.0.hat + X %*% theta.y.hat
    score.y = -d*(y-eta.y.hat)
    sub.grad.y = -t(X) %*% score.y / n
    lambda.y.emp = coef(lm(sub.grad.y ~ theta.y.hat - 1))[[1]] # sometimes glmnet returns solution for slightly different lambda
    theta.y.err <- theta.y.0.hat^2 + sum((theta.y.hat - theta.y)^2)
    
    # Fit outcome model, ipw weights
    lambda.y.glmnet.ipw = lambda.y*n/sum(d*w.ipw) # adjust because glmnet normalizes loss to have sum of weights = n
    theta.y.hat.ipw = as.matrix(predict(model.y.ipw,
                                        s=lambda.y.glmnet.ipw,
                                        weights=d*w.ipw,
                                        type="coefficients",
                                        exact=TRUE,
                                        x=X,y=y))
    theta.y.0.hat.ipw = theta.y.hat.ipw[1]
    theta.y.hat.ipw = theta.y.hat.ipw[2:(p+1)]
    eta.y.hat.ipw = theta.y.0.hat.ipw + X %*% theta.y.hat.ipw
    score.y.ipw = -d*w.ipw*(y-eta.y.hat.ipw)
    sub.grad.y.ipw = -t(X) %*% score.y.ipw /n
    lambda.y.emp.ipw = coef(lm(sub.grad.y.ipw ~ theta.y.hat.ipw - 1))[[1]] # sometimes glmnet returns solution for slightly different lambda
    theta.y.ipw.err <- theta.y.0.hat.ipw^2 + sum((theta.y.hat.ipw - theta.y)^2)
    
    # Estimates of feature means
    mu.x.hat = matrix(colMeans(X),nrow=p)
    mu.x.cfd.hat = matrix(colSums(d * X),nrow=p)/n.1
    
    # Compute summary statistics
    pi.hat = n.1/n
    zeta.y.theta.hat = binary.solve(
      function(zeta.theta){zeta.theta.rootfnc(zeta.theta,weights=d,lambda=lambda.y.emp,p,n)},
      x.min=0,x.max=1,tol=1e-10
    )
    zeta.y.eta.hat = (p/n)/(zeta.y.theta.hat + lambda.y.emp)
    zeta.y.theta.hat.ipw = binary.solve(
      function(zeta.theta) {zeta.theta.rootfnc(zeta.theta,weights=d*w.ipw,lambda=lambda.y.emp.ipw,p,n)},
      x.min=0,x.max=1,tol=1e-10
    )
    zeta.y.eta.hat.ipw = (p/n)/(zeta.y.theta.hat.ipw + lambda.y.emp.ipw)
    gamma.mu.hat = sqrt(max(sum(mu.x.hat^2) - p/n,0))
    gamma.dstar.hat = sqrt(max(sum((mu.x.cfd.hat - mu.x.hat)^2) - (p/n)*(1-pi.hat)/pi.hat,0))
    mu.d.hat = compute.prop.mu.gamma(pi.hat,gamma.dstar.hat,
                                     2/sqrt(n),max.iter=1000) # gamma.reg
    gamma.d.hat <- mu.d.hat$gamma.d.hat
    mu.d.hat <- mu.d.hat$mu.d.hat
    alpha.1.hat <- gaussian.expectation(d.p.score,mu.d.hat,gamma.d.hat,width=10) / pi.hat
    alpha.2.hat <- gaussian.expectation(dd.p.score,mu.d.hat,gamma.d.hat,width=10) / pi.hat
    c.Sigma.hat <- (alpha.2.hat-alpha.1.hat^2)/(1+(alpha.2.hat-alpha.1.hat^2)*gamma.d.hat^2)
    
    # Plug-in outcome modeling
    mu.y.hat = mean(eta.y.hat)
    mu.y.hat.ipw = mean(eta.y.hat.ipw)
    
    # Naive debiasing
    theta.y.0.hat.d.naive = theta.y.0.hat - sum( mu.x * sub.grad.y ) / zeta.y.theta.hat
    theta.y.hat.d.naive = theta.y.hat + sub.grad.y / zeta.y.theta.hat
    mu.y.hat.d.naive = theta.y.0.hat.d.naive + sum(mu.x.hat*theta.y.hat.d.naive)
    
    # Oracle debiasing via confounded distribution
    c.mu = alpha.1
    c.Sigma = (alpha.2 - alpha.1^2)/(1+(alpha.2-alpha.1^2)*gamma.d^2)
    theta.y.0.hat.d.cfd = theta.y.0.hat - sum(mu.x.cfd*sub.grad.y)/zeta.y.theta.hat + c.Sigma*sum(mu.x.cfd*theta.d)*sum(theta.d*sub.grad.y)/zeta.y.theta.hat
    theta.y.hat.d.cfd = theta.y.hat + sub.grad.y / zeta.y.theta.hat - c.Sigma*sum(theta.d*sub.grad.y)/zeta.y.theta.hat*theta.d
    mu.y.hat.d.cfd = theta.y.0.hat.d.cfd + sum(mu.x.hat*theta.y.hat.d.cfd)
    
    # Oracle debiasing with IPW weights
    theta.y.0.hat.d.ipw = theta.y.0.hat.ipw - sum( mu.x * sub.grad.y.ipw ) / zeta.y.theta.hat.ipw
    theta.y.hat.d.ipw = theta.y.hat.ipw + sub.grad.y.ipw / zeta.y.theta.hat.ipw
    mu.y.hat.d.ipw = theta.y.0.hat.d.ipw + sum(mu.x.hat*theta.y.hat.d.ipw)
    
    # Debiasing with estimated propensity model: estimate propensity model with moment method
    s.y.hat <- sqrt( sum( score.y^2 ) / (n*zeta.y.theta.hat^2) )
    s.d.moment.hat <- sqrt( (1-pi.hat)/(pi.hat*alpha.1.hat^2) )
    s.y.d.moment.hat <- 0
    s.x.d.moment.hat <- (1-pi.hat)/(pi.hat*alpha.1.hat)
    
    theta.y.hat.loo <- theta.y.hat + sub.grad.y / zeta.y.theta.hat
    theta.d.hat.loo.moment <- mu.x.cfd.hat - mu.x.hat
    
    theta.y.0.hat.d.emp.moment <- 
      theta.y.0.hat - 
      sum( mu.x.cfd.hat * sub.grad.y ) / zeta.y.theta.hat +
      c.Sigma.hat *
      ( sum(mu.x.cfd.hat*theta.d.hat.loo.moment)/alpha.1.hat-s.x.d.moment.hat*p/n) * 
      ( sum(theta.d.hat.loo.moment/alpha.1.hat*(theta.y.hat.loo-theta.y.hat) ) - s.y.d.moment.hat*(p-n*zeta.y.eta.hat*zeta.y.theta.hat)/n )
    theta.y.hat.d.emp.moment <- 
      theta.y.hat + sub.grad.y / zeta.y.theta.hat -
      (c.Sigma.hat*sum((theta.d.hat.loo.moment/alpha.1.hat)*sub.grad.y)/zeta.y.theta.hat) * theta.d.hat.loo.moment/alpha.1.hat +
      (c.Sigma.hat*s.y.d.moment.hat*(p-n*zeta.y.eta.hat*zeta.y.theta.hat)/n) * theta.d.hat.loo.moment/alpha.1.hat
    mu.y.hat.d.emp.moment <- theta.y.0.hat.d.emp.moment + sum(mu.x.hat*theta.y.hat.d.emp.moment)
    
    data.df <- rbind(
      data.df,
      data.frame(
        
        seed = seed,
        n = n,
        p = p,
        sigma = sigma,
        theta.d.0,
        theta.y.0,
        
        lambda.y=lambda.y,
        lambda.y.emp=lambda.y.emp,
        lambda.y.emp.ipw=lambda.y.emp.ipw,
        
        mu.y = mu.y,
        mu.d = mu.d,
        gamma.y = gamma.y,
        gamma.d = gamma.d,
        pi.bar = pi.bar,
        alpha.1 = alpha.1,
        alpha.2 = alpha.2,
        
        pi.hat = pi.hat,
        alpha.1.hat = alpha.1.hat,
        alpha.2.hat = alpha.2.hat,
        zeta.y.theta.hat = zeta.y.theta.hat,
        zeta.y.eta.hat = zeta.y.eta.hat,
        zeta.y.theta.hat.ipw = zeta.y.theta.hat.ipw,
        zeta.y.eta.hat.ipw =zeta.y.eta.hat.ipw,
        gamma.mu.hat = gamma.mu.hat,
        gamma.dstar.hat = gamma.dstar.hat,
        mu.d.hat = mu.d.hat,
        gamma.d.hat = gamma.d.hat,
        mu.d.hat = mu.d.hat,
        alpha.1.hat = alpha.1.hat,
        alpha.2.hat = alpha.2.hat,
        c.Sigma.hat = c.Sigma.hat,
        
        # Plug-in outcome modeling
        mu.y.hat = mu.y.hat,
        mu.y.hat.ipw = mu.y.hat.ipw,
        
        # Naive debiasing
        mu.y.hat.d.naive = mu.y.hat.d.naive,
        
        # Oracle debiasing via confounded distribution
        mu.y.hat.d.cfd = mu.y.hat.d.cfd,
        
        # Oracle debiasing with IPW weights
        mu.y.hat.d.ipw = mu.y.hat.d.ipw,
        
        # Debiasing with estimated propensity model: estimate propensity model with moment method
        mu.y.hat.d.emp.moment = mu.y.hat.d.emp.moment,
        
        # Error in theta.
        theta.y.err = theta.y.err,
        theta.y.ipw.err = theta.y.ipw.err
        
      )
    )
    
  }
  
  return(data.df)
  
}

# simulation meta data
R.version <- version
glmnet.version <- packageVersion("glmnet")
system.version <- Sys.info()
sim.time <- format(Sys.time(), "%Y%m%d_%H%M%S")
# Set parameters -----
Is_Rad <- TRUE
Is_sparse = FALSE
N.replicates = 500
delta = .8 
theta.d.0 = 0
theta.y.0 = 0
epsilon = .1 # overlap parameter
sigma = .2
lambda.y.seq = exp( seq( log(10) , log(.01) , length.out = 50 ) )
lambda.d=.25

nsmall = 100
nbig = 10000
ns.length = 10
ns <- round( exp( seq( log(nsmall) , log(nbig) , length.out = ns.length ) ) )

specs <- mapply(c,
                rep(1:N.replicates,times=length(ns)),
                ns,
                SIMPLIFY=FALSE)

# Output files
filename.out = paste0("data/celentano_lambda","_n_",ns[length(ns)], "_sparse_",as.integer(Is_sparse),"_Rad_",as.numeric(Is_Rad),"_",sim.time)
filename.tracking = paste0("track/celentano_lambda","_n_",ns[length(ns)], "_sparse_",as.integer(Is_sparse),"_Rad_",as.numeric(Is_Rad),"_",sim.time)


numCores <- min(detectCores(),50)
cl <- makeCluster(numCores)


# 预装包并导出函数和变量到集群节点
clusterEvalQ(cl, {
  library(glmnet)
  library(boot)
})

# 导出所有自定义函数和变量
clusterExport(cl, c("run.replicate", "d.inv.logit", "dd.inv.logit", "prop.score", "d.prop.score", 
                    "dd.prop.score", "properr", "d.properr.d.etad", "d.properr.d.etadhat", 
                    "dd.properr.dd.etad", "dd.properr.d.etad.d.etadhat", "dd.properr.dd.etadhat", 
                    "sqrtm", "p.score", "d.p.score", "dd.p.score", "gaussian.expectation", 
                    "binary.solve", "compute.prop.mu.gamma", "zeta.theta.rootfnc", 
                    "delta", "theta.d.0", "theta.y.0", "sigma", "epsilon", "lambda.y.seq", 
                    "Is_sparse", "Is_Rad", "filename.out", "filename.tracking"))


experiment.data <-  parLapply(cl, specs, function(spec) {
  run.replicate(
    n = spec[2], # n
    p = round(spec[2] / delta), # p
    sigma = sigma,
    Is_Rad = Is_Rad,
    theta.d.0 = theta.d.0,
    theta.y.0 = theta.y.0,
    theta.d = if (Is_sparse) {
      c(1, rep(0, round(spec[2] / delta) - 1))
    } else {
      runif(round(spec[2] / delta), -sqrt(3/(spec[2] / delta)), sqrt(3/(spec[2] / delta)))
    },
    theta.y = if (Is_sparse) {
      c(1, rep(0, round(spec[2] / delta) - 1))
    } else {
      runif(round(spec[2] / delta), -sqrt(3/(spec[2] / delta)), sqrt(3/(spec[2] / delta)))
    },
    mu.x = matrix(rep(0, round(spec[2] / delta)), nrow = round(spec[2] / delta)),
    epsilon = epsilon,
    lambda.y.seq = lambda.y.seq,
    seed = spec[1], # seed
    filename.out = paste(filename.out,".csv",sep=""),
    filename.tracking =  paste(filename.tracking,".txt",sep="")
  )
})

sim.date <- sim.time
experiment.df <- do.call(rbind,experiment.data)
save(experiment.df,R.version,glmnet.version,system.version,sim.time,
     file=paste(filename.out,".Rda",sep=""))
write.csv(
  experiment.df,
  paste0("data/celentano_lambda_sparse_",as.numeric(Is_sparse),"_Rad_",as.numeric(Is_Rad),"_",sim.time,".csv")
)
stopCluster(cl)



