best.sub.sel <- function(x, y, k_max, nruns=50, maxiter=1000, tol=1e-4) {
  # should be normalized
  x <- as.matrix(x)
  y <- as.numeric(y)
  
  # get largest eigenvalue
  xtx <- crossprod(x)
  L <- eigen(xtx)$values[1]
  
  n <- nrow(x)
  p <- ncol(x)
  
  ks <- 1:min(p, n, k_max)
  nk <- length(ks)
  
  # initialize beta
  beta0 <- rnorm(p, 0, 1)
  ids <- order(abs(beta0), decreasing=TRUE)
  beta0[-ids[1:ks[1]]] <- 0
  
  best.betas <- matrix(0,p,nk)
  
  # Run for each k
  for (i in 1:nk) {
    beta.k <- bs.one.k(x, y, ks[i], xtx, L, beta0, nruns, maxiter, tol, polish=polish)
    best.betas[,i] <- beta.k
    beta0 <- beta.k  # Use as a warm start for the next value of k
  }
  
  # TODO return best
  return(best.betas)
}


bs.one.k <- function(x, y, k, xtx, L, beta0, nruns=50, maxiter=1000, tol=1e-4, polish=TRUE) {
  n <- nrow(x)
  p <- ncol(x)
  # form <- ifelse(n < p, 2, 1)
  
  # Run the projected gradient method
  best.beta = bs.proj.grad(x, y, k, L, beta0, nruns, maxiter, tol=tol, polish=polish)
  
  # Set up and run the MIO solver from Cplex. The general form is
  #   min         x^T Q x + c^T x
  #   subject to  Ax <= b
  #               l <= x <= u
  #               some x_i's binary or integral

  # TODO Cplex solver
  return(best.beta)
}


bs.proj.grad <- function(x, y, k, L, beta0, nruns=50, maxiter=1000, tol=1e-4, polish=TRUE) {
  n <- nrow(x)
  p <- ncol(x)
  
  best.beta <- beta0
  best.crit <- Inf
  beta <- beta0
  
  for (r in 1:nruns) {
    for (i in 1:maxiter) {
      beta.old <- beta
      
      # Take gradient descent step
      grad <- -t(x) %*% (y - x %*% beta)
      beta <- beta - grad/L
      
      # Set to zero all but the top k
      ids <- order(abs(beta), decreasing=TRUE)
      beta[-ids[1:k]] <- 0
      
      # Perform least squares polishing, if we are asked to
      # if (polish) beta[ids[1:k]] = lsfit(x[,ids[1:k]],y,int=FALSE)$coef
      
      # Stop condition
      if (norm(beta - beta.old) / max(norm(beta),1) < tol) break
    }
    
    # Check if this run was better than the previous ones
    cur.crit = sum((y - x%*%beta)^2)
    if (cur.crit < best.crit) {
      best.crit <- cur.crit
      best.beta <- beta
    }
    
    beta <- rnorm(p, 0, 1)
    ids <- order(abs(beta), decreasing=TRUE)
    beta[-ids[1:k]] <- 0
  }
  
  return(best.beta)
}


beta_true <- as.vector(c(0.5, -0.3, 0.2))
x <- replicate(10, rnorm(100))
y <- x[, c(2, 4, 7)] %*% beta_true
results <- best.sub.sel(x, y, 10)
