library(gurobi)
library(Matrix)
library(ggplot2)


best.sub.sel <- function(x, y, k_max, nruns=50, maxiter=1000, tol=1e-4, polish=TRUE, mio=TRUE, time.limit=100) {
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
  beta0 <- init.sparse.beta(p, ks[1])
  best.betas <- matrix(0,p,nk)
  
  # Run for each k
  for (i in 1:nk) {
    beta.k <- bs.fixed.k(x, y, ks[i], xtx, L, beta0, nruns, maxiter, tol, polish=polish)
    best.betas[,i] <- beta.k
    beta0 <- beta.k  # Use as a warm start for the next value of k
  }
  
  return(best.betas)
}


bs.fixed.k <- function(x, y, k, xtx=NULL, L=NULL, beta0=NULL, nruns=50, maxiter=1000, tol=1e-4, polish=TRUE, mio=TRUE, time.limit=100) {
  n <- nrow(x)
  p <- ncol(x)
  if (is.null(xtx)) {
    xtx <- crossprod(x)
    L <- eigen(xtx)$values[1]
  }
  if (is.null(beta0)) beta0 <- init.sparse.beta(p, k)
  
  # Run the projected gradient method
  warm.beta <- bs.proj.grad(x, y, k, L, beta0, nruns, maxiter, tol, polish)
  if(mio){
    best.beta <- mio.solve(x, y, k, xtx, warm.beta, time.limit=time.limit)
    return(best.beta)
  }
  else{
    return(warm.beta)
  }
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
      
      # Gradient descent step
      grad <- -t(x) %*% (y - x %*% beta)
      beta <- beta - grad/L
      
      # Set to zero all but the top k
      ids <- order(abs(beta), decreasing=TRUE)
      beta[-ids[1:k]] <- 0
      
      # Stop condition
      if (norm(beta - beta.old) / max(norm(beta),1) < tol) break
    }
    # Least squares polishing
    if (polish) beta[ids[1:k]] = lsfit(x[,ids[1:k]],y,int=FALSE)$coef
    
    # Check if this run was better than the previous ones
    cur.crit = sum((y - x%*%beta)^2)
    if (cur.crit < best.crit) {
      best.crit <- cur.crit
      best.beta <- beta
    }
    beta <- init.sparse.beta(p, k) 
  }
  
  return(best.beta)
}


mio.solve <- function(x, y, k, xtx=NULL, warm.beta=NULL, time.limit=100) {
  n <- nrow(x)
  p <- ncol(x)
  I <- diag(1, p, p)
  
  x_norm <- rep(1, p)
  y_norm <- 1
  if (is.null(warm.beta)) {
    x_norm <- apply(x, 2, function(c){sqrt(sum(c^2))})
    x <- t(apply(x, 1, function(c){c/x_norm}))
    y_norm <- sqrt(sum(y^2))
    y <- y / y_norm
    xtx <- crossprod(x)
    mu <- max(abs(xtx[row(xtx) != col(xtx)]))
    nk <- max(1 - mu * (k-1), 1e-4)
    bigm.u <- min(sqrt(sum(order((t(x)%*%y)^2, decreasing=TRUE)[1:k]))/nk, sqrt((t(y)%*%y)/nk))
    # bigm.l <- sum((order(abs(t(x)%*%y), decreasing=TRUE)[1:k]))/(1-((k-1)*mu))
  }
  else {
    if (is.null(xtx)) xtx <- crossprod(x)
    bigm.u <- 2*max(abs(warm.beta))
    # bigm.l <- k * bigm.u
  }
  
  sum.z.vec <- c(rep(0,p), rep(1,p))
  model <- list()
  model$A <- rbind(cbind(I,-bigm.u*I), cbind(-I,-bigm.u*I), sum.z.vec) 
  model$sense <- rep("<=",2*p+1)
  model$rhs <- c(rep(0, 2*p), k)            # The vector b
  # rows 1-10:  betai - bigm * zi <= 0 which is betai <= bigm * zi
  # rows 11-21: -betai - bigm * zi <= 0 which is betai >= -bigm * zi
  # last row:   sum(z) <= k
  model$vtypes <- c(rep("C",p), rep("B",p)) # Variable types: p continous betas + p discrete z
  model$ub <- c(rep(bigm.u,p), rep(1,p))      # Upper bound Mu on betas and 1 on all z 
  model$lb <- c(rep(-bigm.u,p), rep(0,p))     # Lower bound -Mu on betas and 0 on all z
  model$obj <- c(-2*t(x)%*%y, rep(0,p))     # The vector c in the objective
  model$Q <- bdiag(xtx, matrix(0,p,p))
  # first p rows are multiplied by betas and are from objective function
  # last p rows are multiplied by z and are all 0 since z is not in objective
  
  if (!is.null(warm.beta)) {                # Warm start from proj gradient
    zvec <- as.numeric(warm.beta != 0)
    model$start <- c(warm.beta, zvec)
  }
  params <- list()
  params$TimeLimit <- time.limit
  gur.obj <- quiet(gurobi(model, params))
  
  return(gur.obj$x[1:p] * y_norm / x_norm)
}


quiet = function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}


init.sparse.beta <- function(p, k) {
  beta <- rnorm(p, 0, 1)
  ids <- order(abs(beta), decreasing=TRUE)
  beta[-ids[1:k]] <- 0
  return(beta)
}
  

custom.labeller <- function(variable, value){
  return(paste0(variable, '=', value))
}



plot.times.k <- function(n, p, k.true, k.test.list, reps, time.limit=200) {
  if (length(time.limit) == 1) {
    time.limit <- rep(time.limit, length(k.test.list))
  }
  stepi <- 1
  results <- list()
  progress.bar <- txtProgressBar(min = 1, max = length(k.test.list), initial = 1) 
  df <- data.frame(algorithm=character(), k=integer(), time=double())
  ratios <- data.frame(k=integer(), ratio=double())  
  
  for (j in 1:length(k.test.list)) {
    k <- k.test.list[j]
    tl <- time.limit[j]
    k.ratios <- c()
    for (i in 1:reps) {
      dat <- create.artif.data(n, p, k.true, noise.std = 0.2)
      mio.time <- system.time(mio.solve(dat$x, dat$y, k, time.limit=tl))["elapsed"]
      bs.time <- system.time(bs.fixed.k(dat$x, dat$y, k, time.limit=tl))["elapsed"]
      tmp.df <- data.frame(c('mio', 'bs'), c(k, k), c(mio.time, bs.time))
      names(tmp.df) <- c('algorithm', 'k', 'time')    
      df <- rbind(df, tmp.df)
      k.ratios <- c(k.ratios, bs.time / mio.time)
    }
    tmp.ratio.df <- data.frame(k, mean(k.ratios))
    names(tmp.ratio.df) <- c('k', 'ratio')    
    ratios <- rbind(ratios, tmp.ratio.df)
    
    stepi <- stepi + 1
    setTxtProgressBar(progress.bar, stepi)
  }
  results$p <- ggplot(df, aes(x=time, fill=algorithm)) +
    geom_density(alpha = 0.2) +
    facet_wrap(~k, scales='free', labeller=custom.labeller) +
    ggtitle(paste0("n=", n, ", p=", p, ", k_true=", k.true)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  results$df <- df
  results$ratios <- ratios
  return(results)
}


plot.times.n <- function(p, k.true, k, n.list, reps, time.limit=200) {
  if (length(time.limit) == 1) {
    time.limit <- rep(time.limit, length(n.list))
  }
  stepi <- 1
  results <- list()
  progress.bar <- txtProgressBar(min = 1, max = length(n.list), initial = 1) 
  df <- data.frame(algorithm=character(), n=integer(), time=double())
  ratios <- data.frame(n=integer(), ratio=double())  
  
  for (j in 1:length(n.list)) {
    n <- n.list[j]
    tl <- time.limit[j]
    n.ratios <- c()
    for (i in 1:reps) {
      dat <- create.artif.data(n, p, k.true, noise.std = 0.2)
      mio.time <- system.time(mio.solve(dat$x, dat$y, k, time.limit=tl))["elapsed"]
      bs.time <- system.time(bs.fixed.k(dat$x, dat$y, k, time.limit=tl))["elapsed"]
      tmp.df <- data.frame(c('mio', 'bs'), c(n, n), c(mio.time, bs.time))
      names(tmp.df) <- c('algorithm', 'n', 'time')    
      df <- rbind(df, tmp.df)
      n.ratios <- c(n.ratios, bs.time / mio.time)
    }
    tmp.ratio.df <- data.frame(n, mean(n.ratios))
    names(tmp.ratio.df) <- c('n', 'ratio')    
    ratios <- rbind(ratios, tmp.ratio.df)
    
    stepi <- stepi + 1
    setTxtProgressBar(progress.bar, stepi)
  }
  results$p <- ggplot(df, aes(x=time, fill=algorithm)) +
    geom_density(alpha = 0.2) +
    facet_wrap(~n, scales='free', labeller=custom.labeller) +
    ggtitle(paste0("p=", p, ", k=", k, ", k_true=", k.true)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  results$df <- df
  results$ratios <- ratios
  return(results)
}


plot.times.p <- function(n, k.true, k, p.list, reps, time.limit=200) {
  stepi <- 1
  results <- list()
  progress.bar <- txtProgressBar(min = 1, max = length(p.list), initial = 1) 
  df <- data.frame(algorithm=character(), p=integer(), time=double())
  ratios <- data.frame(p=integer(), ratio=double())  
  
  for (j in 1:length(p.list)) {
    p <- p.list[j]
    tl <- time.limit[j]
    p.ratios <- c()
    for (i in 1:reps) {
      dat <- create.artif.data(n, p, k.true, noise.std = 0.2)
      mio.time <- system.time(mio.solve(dat$x, dat$y, k, time.limit=tl))["elapsed"]
      bs.time <- system.time(bs.fixed.k(dat$x, dat$y, k, time.limit=tl))["elapsed"]
      tmp.df <- data.frame(c('mio', 'bs'), c(p, p), c(mio.time, bs.time))
      names(tmp.df) <- c('algorithm', 'p', 'time')    
      df <- rbind(df, tmp.df)
      p.ratios <- c(p.ratios, bs.time / mio.time)
    }
    tmp.ratio.df <- data.frame(p, mean(p.ratios))
    names(tmp.ratio.df) <- c('p', 'ratio')    
    ratios <- rbind(ratios, tmp.ratio.df)
    
    stepi <- stepi + 1
    setTxtProgressBar(progress.bar, stepi)
  }
  results$p <- ggplot(df, aes(x=time, fill=algorithm)) +
    geom_density(alpha = 0.2) +
    facet_wrap(~p, scales='free', labeller=custom.labeller) +
    ggtitle(paste0("n=", n, ", k=", k, ", k_true=", k.true)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  results$df <- df
  results$ratios <- ratios
  return(results)
}



compare.solvers <- function(dat, k) {
  mio.result <- system.time(mio.solve(dat$x, dat$y, k, time.limit = 3600))
  print(paste0('MIO: ', mio.result['elapsed']))
  
  bs.result <- system.time(bs.fixed.k(dat$x, dat$y, k, time.limit = 3600))
  print(paste0('BS: ', bs.result['elapsed']))
  results <- list()
  results$mio <- mio.result
  results$bs <- bs.result
  return(results)
}


create.artif.data <- function(n, p, k, noise.std=0.1) {
  dataset <- list()
  dataset$x <- as.matrix(replicate(p, rnorm(n)))
  dataset$active.set <- sort(sample(1:p, k, replace=F))
  dataset$beta <- rep(0, p) 
  dataset$beta[dataset$active.set] <- rnorm(k)
  dataset$y <- dataset$x %*% dataset$beta 
  dataset$y <- dataset$y + rnorm(nrow(dataset$x), 0, noise.std)
  return(dataset)
}


res.k <- plot.times.k(500, 200, 32, c(5, 20, 32, 64), reps=50, time.limit=200)
ggsave("times_k.png", plot = res.k$p)

res.n <- plot.times.n(200, 32, 20, c(250, 500, 2000, 10000), reps=50, time.limit=200)
ggsave("times_n.png", plot = res.n$p)

res.p <- plot.times.p(1000, 32, 20, c(50, 100, 200, 250), reps=20, time.limit=c(200, 200, 200, 500))
ggsave("times_.png", plot = res.p$p)



dat <- create.artif.data(n=1000, p=500, k=32, noise.std=0.2)
res <- compare.solvers(dat, 20)

# n=1000, p=500, k=32, k.test=20, noise.std=0.2
# [1] "MIO: 3600.43"
# [1] "BS: 704.529999999999"

# length(res$mio[which(res$mio != 0, arr.ind = TRUE)])
# res$bs[which(res$bs != 0, arr.ind = TRUE)]
# dat$beta[which(res$bs != 0, arr.ind = TRUE)]


# beta_true <- as.vector(c(0.5, -0.3, 0.2))
# x <- replicate(10, rnorm(100))
# y <- x[, c(2, 4, 7)] %*% beta_true + rnorm(nrow(x), 0, 0.2)


# results <- best.sub.sel(x, y, 6, polish=TRUE)
# mio.solve(x, y, 4)

