# Author: Leandro Loriato (llrt at impa.br)
#
# Auxiliary functions to simulate the path of the solution of a given SDE
#
# SDE.solve function based on sde.sim function from Stefano Iacus’s Simulation and Inference for Stochastic
# Differential Equations book (section 2.4, page 69)

# auxiliary function that simulates the solution using Euler-Maruyama discretization method
SDE.solve.euler = function(t0, T, N.t, X0, drift_, sigma_, Z){
  t = seq(t0, T, length=(N.t+1)) # spans the time range
  dt = (T-t0)/N.t # calculates time discretization step
  sqrt.dt = sqrt(dt) # for saving computation, pre-computes sqrt(dt)
  X = matrix(ncol=ncol(Z), nrow=(N.t+1)) # pre-allocates a matrix with N.t+1 rows for X
  X[1,] = X0 # first X value is X0
  for(i in 2:(N.t+1)){
    drift__ = drift_(t[i-1], X[i-1,])
    sigma__ = sigma_(t[i-1], X[i-1,])
    X[i,] = X[i-1,] + drift__*dt + sigma__*sqrt.dt*Z[i-1,]
  }
  invisible(X)
}
# for performance, pre-compile CPU intensive function
SDE.solve.euler = cmpfun(SDE.solve.euler)
# auxiliary function that simulates the solution using Milstein discretization method
SDE.solve.milstein = function(t0, T, N.t, X0, drift_, sigma_, sigma.x_, Z){
  t = seq(t0, T, length=(N.t+1)) # spans the time range

  dt = (T-t0)/N.t # calculates time discretization step
  sqrt.dt = sqrt(dt) # for saving computation, pre-computes sqrt(dt)
  X = matrix(ncol=ncol(Z), nrow=(N.t+1)) # pre-allocates a matrix with N.t+1 rows for X
  X[1,] = X0 # first X value is X0
  for(i in 2:(N.t+1)){
    drift__ = drift_(t[i-1], X[i-1,])
    sigma__ = sigma_(t[i-1], X[i-1,])
    sigma.x__ = sigma.x_(t[i-1], X[i-1,])
    42
    X[i,] = X[i-1,] + drift__*dt +
      sigma__*sqrt.dt*Z[i-1,] +
      (1/2)*sigma__*sigma.x__*(dt*Z[i-1,]^2 - dt)
  }
  invisible(X)
}
# for performance, pre-compile CPU intensive function
SDE.solve.milstein = cmpfun(SDE.solve.milstein)
# main function that, for an SDE of form dX(t) = drift(t, X) dt + sigma(t, X) dW(t), simulates its solution
SDE.solve = function(t0=0, T=1, X0=1, N.t=100, drift, sigma, method=c("euler", "milstein"), Z){
  if(missing(method)){ # if not provided, use Euler-Maruyama method
    method = "euler"
  } else{ # if provided, use chosen method
    method = match.arg(method)
  }
  # evaluates given drift, sigma and sigma.x expressions into functions
  drift_ = function(t, x){eval(drift)} # drift function
  sigma_ = function(t, x){eval(sigma)} # sigma function
  if(method=="milstein"){ # if Milstein method is chosen, sigma.x should have been provided
    sigma.x = D(sigma, "x") # derivative of sigma with respect to X
    sigma.x_ = function(t, x){eval(sigma.x)} # sigma.x function
  }
  # generates a sample path from SDE solution’s X with chosen method
  if(method == "euler"){
    X = SDE.solve.euler(t0, T, N.t, X0, drift_, sigma_, Z)
  } else if(method == "milstein"){
    X = SDE.solve.milstein(t0, T, N.t, X0, drift_, sigma_, sigma.x_, Z)
  }
  invisible(X)
}
# for performance, pre-compile CPU intensive function
SDE.solve = cmpfun(SDE.solve)