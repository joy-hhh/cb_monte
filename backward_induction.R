# Author: Leandro Loriato (llrt at impa.br)
#
# Auxiliary functions to calculate a convertible bond’s price using Monte Carlo method with
# backward induction techniques: Least-Squared Monte Carlo (LSMC) e Hedged Monte Carlo (HMC)

library(orthopolynom)

# given relevant parameters, calculates the convertible bond price via Black formula
price.black = function(r, q, t, T, S, sigma, bond.params){
  bond = bond.value.T(bond.params, r)
  
  d1 = function(t, s, bond, sigma){1/(sigma * sqrt(T-t)) * (log((bond.params$conversion$ratio((T-t), s) * s)/bond) + ((r - q) + (sigma^2)/2) * (T-t))}
  d2 = function(t, s, bond, sigma){d1(t, s, bond, sigma) - sigma*sqrt(T-t)}
  
  price = discounted.coupon(bond.params, t, r) +
    exp(-r*(T-t))*bond +
    bond.params$conversion$ratio((T-t), S) * S * exp(-q*(T-t)) * pnorm(d1(t, S, bond, sigma)) -
  bond * exp(-r*(T-t)) * pnorm(d2(t, S, bond, sigma))
  invisible(price)
}

# given relevant parameters, calculates the delta hedge from Black formula
delta.hedge.black = function(r, q, t, T, S, sigma, bond.params){
  bond = bond.value.T(bond.params, r)
  d1 = function(t, s, bond, sigma){1/(sigma * sqrt(T-t)) * (log((bond.params$conversion$ratio((T-t), s) * s)/
                                                                  bond) + ((r - q) + (sigma^2)/2) * (T-t))}
  price = bond.params$conversion$ratio((T-t), S) * exp(-q*(T-t)) * pnorm(d1(t, S, bond, sigma))
  invisible(price)
}
# function that evaluates desired basis functions and its derivatives
basis = function(
  type=c("laguerre", "chebyshev_1", "chebyshev_2", "chebyshev_3", "hermite", "legendre",
          "w_laguerre", "w_chebyshev_1", "w_chebyshev_2", "w_chebyshev_3", "w_hermite", "w_legendre",
          "black"),
  M=4, ...){
  if(missing(type)){ # if not provided, use Laguerre Polynomials
    type = "laguerre"
  } else{ # if provided, use chosen type
    type = match.arg(type)
  }
  # auxiliary function that, given a recurrences list for desired polynom,
  # evaluates the polynom functions and its derivative functions
  polynom.fun = function(recurrences){
    polynom = orthogonal.polynomials(recurrences)

    
    
    fun = polynomial.functions(polynom)
    polynom.x = polynomial.derivatives(polynom)
    fun.x = polynomial.functions(polynom.x)
    invisible(list(fun=fun, fun.x=fun.x))
  }
  # auxiliary function that, given a recurrences list for desired polynom and its
  # respective weight expression, evaluates the weighted polynom functions and its
  # derivative functions
  weighted.polynom.fun = function(recurrences, weight){
    weight.fun = function(x){eval(weight)}
    weight.x = D(weight, "x")
    weight.x.fun = function(x){eval(weight.x)}
    gen.fun = function(f){
      force(f)
      function(x){weight.fun(x)*f(x)}
    }
    gen.fun.x = function(f, f.x){
      force(f)
      force(f.x)
      function(x){weight.fun(x)*f.x(x) + weight.x.fun(x)*f(x)}
    }
    polynom = orthogonal.polynomials(recurrences)
    polynom.fun = polynomial.functions(polynom)
    fun = list()
    for(i in 1:length(polynom.fun)){
      fun[[i]] = gen.fun(f=polynom.fun[[i]])
    }
    polynom.x = polynomial.derivatives(polynom)
    polynom.x.fun = polynomial.functions(polynom.x)
    fun.x = list()
    for(i in 1:length(polynom.x.fun)){
      fun.x[[i]] = gen.fun.x(f=polynom.fun[[i]], f.x=polynom.x.fun[[i]])
    }
    invisible(list(fun=fun, fun.x=fun.x))
  }
  # for specified basis, evaluate its functions and derivative functions list
  switch(type,
         "laguerre"={ # Laguerre Polynomials
           recurrences = laguerre.recurrences(n=(M-1), normalized=TRUE)
           ret.fun = polynom.fun(recurrences)
           fun = ret.fun$fun
           fun.x = ret.fun$fun.x
         },
         "w_laguerre"={ # Weighted Laguerre Polynomials
           recurrences = laguerre.recurrences(n=(M-1), normalized=TRUE)
           weight = expression(exp(-x))
           
           ret.weighted = weighted.polynom.fun(recurrences, weight)
           fun = ret.weighted$fun
           fun.x = ret.weighted$fun.x
         },
         "chebyshev_1"={ # 1st Kind Chebyshev Polynomials
           recurrences = chebyshev.t.recurrences(n=(M-1), normalized=TRUE)
           ret.fun = polynom.fun(recurrences)
           fun = ret.fun$fun
           fun.x = ret.fun$fun.x
         },
         "w_chebyshev_1"={ # Weighted 1st Kind Chebyshev Polynomials
           recurrences = chebyshev.t.recurrences(n=(M-1), normalized=TRUE)
           weight = expression(1/sqrt(1-(x^2)/4))
           ret.weighted = weighted.polynom.fun(recurrences, weight)
           fun = ret.weighted$fun
           fun.x = ret.weighted$fun.x
         },
         "chebyshev_2"={ # 2nd Kind Chebyshev Polynomials
           recurrences = chebyshev.s.recurrences(n=(M-1), normalized=TRUE)
           ret.fun = polynom.fun(recurrences)
           fun = ret.fun$fun
           fun.x = ret.fun$fun.x
         },
         "w_chebyshev_2"={ # Weighted 2nd Kind Chebyshev Polynomials
           recurrences = chebyshev.s.recurrences(n=(M-1), normalized=TRUE)
           weight = expression(1/sqrt(1-(x^2)/4))
           ret.weighted = weighted.polynom.fun(recurrences, weight)
           fun = ret.weighted$fun
           fun.x = ret.weighted$fun.x
         },
         "chebyshev_3"={ # 3rd Kind Chebyshev Polynomials
           recurrences = chebyshev.c.recurrences(n=(M-1), normalized=TRUE)
           ret.fun = polynom.fun(recurrences)
           fun = ret.fun$fun
           fun.x = ret.fun$fun.x
         },
         "w_chebyshev_3"={ # Weighted 3rd Kind Chebyshev Polynomials
           recurrences = chebyshev.c.recurrences(n=(M-1), normalized=TRUE)
           weight = expression(1/sqrt(1-x^2))
           ret.weighted = weighted.polynom.fun(recurrences, weight)
           fun = ret.weighted$fun
           fun.x = ret.weighted$fun.x
         },
         "hermite"={ # Hermite Polynomials
           recurrences = hermite.h.recurrences(n=(M-1), normalized=TRUE)
           ret.fun = polynom.fun(recurrences)
           fun = ret.fun$fun
           fun.x = ret.fun$fun.x
         },
         "w_hermite"={ # Weighted Hermite Polynomials
           recurrences = hermite.h.recurrences(n=(M-1), normalized=TRUE)
           
           weight = expression(exp(-x^2))
           ret.weighted = weighted.polynom.fun(recurrences, weight)
           fun = ret.weighted$fun
           fun.x = ret.weighted$fun.x
         },
         "legendre"={ # Legendre Polynomials
           recurrences = legendre.recurrences(n=(M-1), normalized=TRUE)
           ret.fun = polynom.fun(recurrences)
           fun = ret.fun$fun
           fun.x = ret.fun$fun.x
         },
         "w_legendre"={ # Weighted Legendre Polynomials
           recurrences = legendre.recurrences(n=(M-1), normalized=TRUE)
           weight = expression(1)
           ret.weighted = weighted.polynom.fun(recurrences, weight)
           fun = ret.weighted$fun
           fun.x = ret.weighted$fun.x
         },
         "black"={ # Black Basis
           bond.T = bond.value.T(bond.params, r)
           fun = list()
           fun[[1]] = function(x){1 + x * 0}
           fun[[2]] = function(x){x - exp(-r*T)*bond.T/bond.params$conversion$ratio(T, x)}
           fun[[3]] = function(x){price.black(x, r=r, q=q, t=0, T=T, sigma=.sigma, bond.params=bond.params)}
           fun.x = list()
           fun.x[[1]] = function(x){0 + x * 0}
           fun.x[[2]] = function(x){1 + x * 0}
           fun.x[[3]] = function(x){delta.hedge.black(x, r=r, q=q, t=0, T=T, sigma=.sigma, bond.params)}
         }
  )
  ret = list(type=type, M=M, fun=fun, fun.x=fun.x)
  invisible(ret)
}
# generate a function that, given an evaluted basis function, calculates its value given an input vector X
funval = function(X){
  invisible(function(fun.eval){fun.eval(X)})
}
# sanitize a exercise matrix, leaving for each path only the entry that effected the
# exercise
clean.exercise = function(exercise, N.t){
  cleaned.exercise = matrix(exercise, nrow=nrow(exercise), ncol=ncol(exercise))
  stopping.time = rep(1, ncol(exercise))
  for(i in 1:ncol(exercise)){
    for(j in 1:N.t){
      if(cleaned.exercise[[j,i]] > 0){ # some exercise action was taken at current time
        stopping.time[[i]] = j
        cleaned.exercise[(j+1):(N.t+1),i] = 0 # as instrument already suffered an exercise action,
        
        # every other future exercise action does not take place
        break
      }
    }
  }
  invisible(list(cleaned.exercise=cleaned.exercise, stopping.time=stopping.time))
}
# for performance, pre-compile CPU intensive function
clean.exercise = cmpfun(clean.exercise)
# utilitary function for packing results from a Monte Carlo simulation into an appropriate object
as.montecarlo = function(C, S, S0, method, alpha=0.05){
  .S = S[1,] # initial stock price values
  indexes = which(.S==S0) # which column indexes correspond to stock price starting at S0 at time t=t0
  filtered.C = C[,indexes] # matrix of paths of price starting at S0 at time t=t0
  value = mean(filtered.C[1,]) # fair price at time t=t0
  sigma = sd(filtered.C[1,]) # standard deviation of price at time t=t0
  M = length(filtered.C[1,])
  Z.alpha = qnorm(alpha/2, lower.tail=F)
  confidence.interval = c(value-Z.alpha*sigma/sqrt(M), value+Z.alpha*sigma/sqrt(M))
  MC = (list(C=filtered.C, method=method, value=value, sigma=sigma, alpha=alpha, confidence.interval=
               confidence.interval))
  class(MC) = "montecarlo"
  return(MC)
}
# for performance, pre-compile CPU intensive function
as.montecarlo = cmpfun(as.montecarlo)
# utilitary functions for summarizing relevant data from a montecarlo object
summary.montecarlo = function(MC){
  s = list(value=MC$value, sigma=MC$sigma, confidence.interval=MC$confidence.interval)
  class(s) = "summary.montecarlo"
  return(s)
}
setGeneric("summary.montecarlo")
print.summary.montecarlo = function(s){
  print(paste("value: ", s$value))
  print(paste("sigma: ", s$sigma))
  print("confidence interval: ")
  print(s$confidence.interval)
}
setGeneric("print.summary.montecarlo")
# auxiliary function that evaluates regression args for LSMC method
regression.args.LSMC = function(S, C, j, rho, Ca.val){
  # evaluate regression args for LSMC method
  y = rho * C[j+1,]
  x = Ca.val
  
  ret = list(x=x, y=y)
  invisible(ret)
}
# for performance, pre-compile CPU intensive function
regression.args.LSMC = cmpfun(regression.args.LSMC)
# auxiliary function that evaluates regression args for HMC method
regression.args.HMC = function(S, C, j, rho, Ca.val, Fa.val){
  # evaluate regression args for HMC method
  y = rho * C[j+1,]
  x = Ca.val + Fa.val * (rho*S[j+1,] - S[j,])
  ret = list(x=x, y=y)
  invisible(ret)
}
# for performance, pre-compile CPU intensive function
regression.args.HMC = cmpfun(regression.args.HMC)
# auxiliary function that solves regression using regular method,
# then calculating projected expected value according to regression
regression.regular = function(x, y, Ca.val){
  # in regular method, least squares fit is solved using QR on
  # whole sample
  alpha = solve(qr(x, LAPACK=TRUE), y)
  # calculate expected value
  exp.val = tcrossprod(Ca.val, t(alpha)) # Ca.val %*% alpha
  invisible(exp.val)
}
# for performance, pre-compile CPU intensive function
regression.regular = cmpfun(regression.regular)
# auxiliary function that solves regression using Bouchard’s method,
# then calculating projected expected value according to regression
regression.bouchard = function(x, y, S, regression.intervals, Ca.val){
  # in Bouchard’s method, least squares fit is solved using several QR,
  # each on a partition of the sample
  exp.val = vector(mode="numeric", length=length(y))
  .df.regression = data.frame(S=S, indexes=1:length(y))
  .df.regression = .df.regression[order(.df.regression$S),]
  interval.length = length(y)/regression.intervals
  for(k in 1:regression.intervals){ # for each regression interval
    # evaluate indexes corresponding to current interval
    interval.indexes = ((k-1)*interval.length+1):(min(k*interval.length, interval.length*regression.intervals))
    # evaluates corresponding original indexes
    original.indexes = .df.regression$indexes[interval.indexes]
    # extract partitioned x and y based on interval indexes
    .y = y[original.indexes]
    .x = x[original.indexes,]
    
    # solve regression for partitioned x and y
    alpha = solve(qr(.x, LAPACK=TRUE), .y)
    # calculate expected value
    exp.val[original.indexes] = tcrossprod(Ca.val[original.indexes,], t(alpha)) # Ca.val %*% alpha
  }
  invisible(exp.val)
}
# for performance, pre-compile CPU intensive function
regression.bouchard = cmpfun(regression.bouchard)
# main function that, given relevant parameters, calculates the american contingent claim price via Monte
Carlo
price.MC = function(S, S0, r, q, t0, T, N.t, N.MC,
                    backward.induction.method=c("LSMC", "HMC"), price.fun, bond.params,
                    regression.method=c("regular", "bouchard"), regression.intervals){
  # sets backward induction method
  if(missing(backward.induction.method)){ # if not provided, use LSMC method
    backward.induction.method = "LSMC"
  } else{ # if provided, use chosen method
    backward.induction.method = match.arg(backward.induction.method)
  }
  # sets regression method
  if(missing(regression.method)){ # if not provided, use regular method
    regression.method = "regular"
  } else{ # if provided, use chosen method
    regression.method = match.arg(regression.method)
    368
    # if Bouchard’s regression method is chosen and number of regression intervals is not informed,
    # sets it in such a way that each interval contains 200 values
    if(regression.method=="bouchard" && missing(regression.intervals)){
      regression.intervals = N.MC/200
    }
  }
  dt = (T-t0)/N.t # calculates time discretization step
  t = seq(t0, T, by=dt) # spans the time range
  # pre-allocates a N.t+1 x N.MC matrix for C (contingent claim’s price at each time) and exercise
  C = matrix(numeric(0), nrow=(N.t+1), ncol=N.MC)
  exercise = matrix(numeric(0), nrow=(N.t+1), ncol=N.MC)
  # bond value at time T
  bond.T = bond.value.T(bond.params, r)
  exercise.time.bitmap.ret = exercise.time.bitmap(T, bond.params) # computes exercise time restriction
  bitmap
  for(i in 1:N.MC){ # for each path, updates C and exercise values at time with index j
    S.vector = S[1:(N.t+1), i]
    payoff.ret = payoff(T, S.vector, r, bond.T, exercise.time.bitmap.ret, bond.params) # intrinsic value info at time T
    C[[N.t+1,i]] = payoff.ret$value # initializes C values at time T, as payoff at time T is known
    exercise[[1+N.t,i]] = payoff.ret$action # evaluates exercise decision
    
  }
  # sets basis functions
  Ca.fun = price.fun$fun # basis functions for contingent claim’s price
  if(backward.induction.method=="HMC"){
    Fa.fun = price.fun$fun.x # basis functions for hedge function for contingent claim’s price
  }
  rho = exp(-r*dt) # pre-computing exp(-r*dt) to save computation
  # apply the chosen backward induction algorithm
  for(j in N.t:1){
    cat("*") # print "*" character to give sense of computation progress
    exercise.time.bitmap.ret = exercise.time.bitmap(t[[j]], bond.params) # computes exercise time restriction
    bitmap
    for(i in 1:N.MC){ # for each path, updates C and exercise values at time with index j
      S.vector = S[1:j, i]
      payoff.ret = payoff(t[[j]], S.vector, r, rho*C[[j+1,i]], exercise.time.bitmap.ret, bond.params)
      C[[j,i]] = payoff.ret$value # initialize C value with intrinsic value
      exercise[[j,i]] = payoff.ret$action # evaluates exercise decision
    }
    # evaluate regression args according to chosen backward induction method
    if(backward.induction.method=="LSMC"){
      # apply Ca functions over S values at time with index j
      Ca.val = sapply(Ca.fun, funval(S[j,]))
      # evaluate regression args
      ret.regression.args = regression.args.LSMC(S, C, j, rho, Ca.val)
    } else if(backward.induction.method=="HMC"){
      # apply Ca and Fa functions over S values at time with index j
      Ca.val = sapply(Ca.fun, funval(S[j,]))
      Fa.val = sapply(Fa.fun, funval(S[j,]))
      # evaluate regression args
      ret.regression.args = regression.args.HMC(S, C, j, rho, Ca.val, Fa.val)
    }
    y = ret.regression.args$y
    x = ret.regression.args$x
    # solve least-squares fit using chosen regression method and
    # calculate projected expected value at time with index j
    if(regression.method=="regular"){ # regular regression method
      exp.val = regression.regular(x, y, Ca.val)
    } else if (regression.method=="bouchard"){ # Bouchard’s regression method
      exp.val = regression.bouchard(x, y, S[j,], regression.intervals, Ca.val)
    }
    # for each component over time with index j:
    # -> C if the value of C is greater than or equal the expected value calculated
    # -> expected value if the value of C is less than the expected value calculated
    for(i in 1:N.MC){
      if(C[[j,i]] >= exp.val[[i]]){ # exercise takes place
        #C[[j,i]] = C[[j,i]] # update value to current intrinsic value
      } else{ # return next time’s discounted value
        C[[j,i]] = rho*C[[j+1,i]] # update value to discounted value of time j+1
        exercise[[j,i]] = exercise.legend[["continuation"]] # update current exercise action to continuation
        
      }
      # if current time is a coupon payment date, also add coupon to current value
      if(t[[j]] %in% bond.params$coupon$dates){
        C[[j,i]] = C[[j,i]] + bond.params$coupon$rate * bond.params$principal
      }
    }
  }
  # clean exercise matrix to retain only the first exercise action at each path, calculating its stopping time
  ret.clean = clean.exercise(exercise, N.t)
  stopping.time = ret.clean$stopping.time
  exercise = ret.clean$cleaned.exercise
  ret = list(MC=as.montecarlo(C, S, S0, method=backward.induction.method, alpha=0.05), exercise=exercise,
             stopping.time=stopping.time)
  invisible(ret)
}
price.MC = cmpfun(price.MC)