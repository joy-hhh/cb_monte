# Author: Leandro Loriato (llrt at impa.br)
#
# Calculates the price of a Convertible Bond in setting 3, in which we have:
# − Complex product (callable, puttable, american−style conversion)
# − Simplified model (only asset is stochastic, do not consider credit risk issues)
#

library(compiler)

# import auxiliary functions
source('bond.R')
source('payoff.r')
source('sde.r')
source('backward_induction.r')
source('boundary.r')

# general parameters
t0 = 0 # initial time
T = 2 # final time
N.t = 100 #252*2 # number of discretization steps in time
dt = (T-t0)/N.t # calculates time discretization step

# bond model
r = 0.05 # assuming deterministic risk−free interest rate
bond.principal = 100 # bond’s principal
bond.redemption.ratio = 1 # ratio at which principal is redempted at maturity
bond.coupon.frequency = 0.5 # frequency on which coupons are paid (assuming equally spaced coupons)
bond.coupon.dates = seq(t0 + bond.coupon.frequency, T, by=bond.coupon.frequency) # coupon dates, starting from payment period
# immediately after t0 to, and including, T
bond.coupon.rate = 0 # coupon rate to be paid on bond’s principal (assuming equally valued coupons and pre−fixed rate)
# combine coupon parameters into one data structure
bond.coupon = list(frequency=bond.coupon.frequency, dates=bond.coupon.dates, rate=bond.coupon.rate)
# number of stocks a convertible bond’s may be exchanged for
bond.conversion.ratio = function(t, S.vector){
  # if mean of stock prices of the last 20 last days (including current one) is greater
  # than 130% of initial stock price, then conversion ratio is reset to 0.8
  n = length(S.vector)
  mean.S = mean(S.vector[max(1, n-19):n])
  if(mean.S > 1.3*S0){
    conversion.ratio = 0.8 # sets conversion ratio to 0.8
  } else{
    conversion.ratio = 1 # sets conversion ratio to 1
  }
  invisible(conversion.ratio)
}
bond.conversion.dates = seq(t0, T, by=dt) # assuming american−style (conversion possible at each time)
bond.conversion.restriction = function(t, S.vector){ # embodies restrictions, beside those time−related,
  # that must be met for conversion action
  # to be allowed (can be modified to encompass other types of
  # restrictions, e.g. of Contingent Convertible Bonds)
  invisible(TRUE) # no conversion restriction
}
# combine coupon parameters into one data structure
bond.conversion = list(ratio=bond.conversion.ratio, dates=bond.conversion.dates, restriction=bond.conversion.restriction)
# call optionality:
# issuer may choose to end contract before maturity, forcing the investor the option to exchange the
convertible
# bond for given strike price or to convert it (if current time is one of conversion dates)
bond.call.present = TRUE # whether a call optionality is present or not
bond.call.strike = function(t, S.vector){ # call’s strike price
  call.strike = 110 # assuming constant over time
  invisible(call.strike)
}
bond.call.dates = seq(t0, T, by=dt) # dates where call optionality may take place
bond.call.restriction = function(t, S.vector){ # embodies restrictions, beside those time−related,
  # that must be met for call action to be allowed
  # (can be modified to encompass other types of restrictions)
  # if mean of stock prices of the last 20 days (including current one) is greater
  # than 110% of initial stock price, then call exercise is allowed
  n = length(S.vector)
  mean.S = mean(S.vector[max(1, n-19):n])

    if(mean.S > 1.1*S0){
    ret = TRUE # call exercise is allowed
  } else{
    ret = FALSE # call exercise is not allowed
  }
  invisible(ret)
}
# combine call parameters into one data structure
bond.call = list(present=bond.call.present, strike=bond.call.strike, dates=bond.call.dates, restriction=bond.call.restriction)

# put optionality:
# investor may choose to end contract before maturity, forcing the issuer to buy the convertible bond for
given
# strike price
bond.put.present = TRUE # whether a put optionality is present or not
bond.put.strike = function(t, S.vector){ # put’s strike price
  put.strike = 98 # assuming constant over time
  invisible(put.strike)
}
bond.put.dates = seq(t0, T, by=dt) # dates where put optionality may take place
bond.put.restriction = function(t, S.vector){ # embodies restrictions, beside those time−related,
  # that must be met for put action to be allowed
  # (can be modified to encompass other types of restrictions, e.g. of
  # Contingent Convertible Bonds)
  invisible(TRUE)
}
# combine put parameters into one data structure
bond.put = list(present=bond.put.present, strike=bond.put.strike, dates=bond.put.dates, restriction=bond.put.restriction)
# creates a data structure to contain all relevant bond parameters
bond.params = list(principal=bond.principal, redemption.ratio=bond.redemption.ratio,
                   coupon=bond.coupon, conversion=bond.conversion,
                   call=bond.call, put=bond.put)

# config about which prices, given by each method, to calculate
calculate.LSMC.euler = TRUE # whether to calculate price given by LSMC/Euler or not
calculate.LSMC.milstein = TRUE # whether to calculate price given by LSMC/Milstein or not
calculate.LSMC = list(euler=calculate.LSMC.euler,
                      milstein=calculate.LSMC.milstein)
calculate.HMC.euler = TRUE # whether to calculate price given by HMC/Euler or not
calculate.HMC.milstein = TRUE # whether to calculate price given by HMC/Milstein or not
calculate.HMC = list(euler=calculate.HMC.euler,
                     milstein=calculate.HMC.milstein)
# creates a data structure to contain all relevant graph config about which prices to calculate
calculate = list(LSMC=calculate.LSMC,
                 HMC=calculate.HMC)

# stock model
S0 = 100 # initial value for S
q = 0.1 # continuous dividend yield
# SDE of form dX(t) = drift(t, X) dt + sigma(t, X) dW(t)
.drift = (r-q) # stock drift equals the risk−free interest rate in risk−neutral measure
drift = expression(.drift*x) # expression for the drift
.sigma = 0.4 # stock volatility
sigma = expression(.sigma*x) # expression for the sigma

# Monte Carlo parameters
S0.delta = 0.8*S0 # distance of original S0 to maximum/minimum value in the array
S0.max = S0 + S0.delta # maximum S0 value in initial stock price’s array
S0.min = S0 - S0.delta # minimum S0 value in initial stock price’s array
N.S0.band = 100 # number of values in each band of the array
N.S0 = 2*N.S0.band # number of values in the whole array (excluding the central value, i.e. original one)
d.S0 = (S0.delta)/N.S0.band # stock price step in the array
S0.array = seq(S0.min, S0.max, by=d.S0) # initial stock price’s array
S0.array.without.central = S0.array[-(N.S0.band+1)] # initial stock price’s array without central value

N.MC.central = 10000 # number of Monte Carlo simulations for original initial stock price value
N.MC.band = 10 # number of Monte Carlo simulations for rest of initial stock price’s array

price.fun.LSMC = basis(type="hermite", M=5) # LSMC’s basis functions for instrument price
price.fun.HMC = basis(type="hermite", M=5) # HMC’s basis functions for instrument price

# generating paths
Z = matrix(rnorm(N.t*(N.MC.central + N.MC.band*N.S0)), ncol=(N.MC.central + N.MC.band*N.S0), nrow=N.t)
# generates N(0,1) samples

if(any(calculate$LSMC$euler, calculate$HMC$euler)){ # simulate Euler−Maruyama paths only if its
  # corresponding prices are to be calculated
  # simulates the solution for the given SDE using Euler−Maruyama method
  
  S.euler = matrix(ncol=(N.MC.central + N.MC.band*N.S0), nrow=(N.t+1)) # pre−allocates matrix for S
  
  for(.S0 in S0.array.without.central){ # for each value in the initial stock price’s array, except the original one
    # generates N.MC.band stock price paths
    j = which(S0.array.without.central==.S0) # obtain index of current stock price in array without central
    indexes = ((j-1)*N.MC.band+1):(min(j*N.MC.band, N.MC.band*N.S0))
    S.euler[,indexes] = SDE.solve(t0=t0, T=T, X0=.S0, N.t=N.t, drift=drift, sigma=sigma, method="euler", Z=Z[,indexes])
  }

  # finally, for the original initial stock price, generates N.MC.central stock price paths
  indexes = (N.MC.band*N.S0+1):(N.MC.central + N.MC.band*N.S0)
  S.euler[,indexes] = SDE.solve(t0=t0, T=T, X0=S0, N.t=N.t, drift=drift, sigma=sigma, method="euler", Z=Z[,indexes])
}

if(any(calculate$LSMC$milstein, calculate$HMC$milstein)){ # simulate Milstein paths only if its
  # corresponding prices are to be calculated
  # simulates the solution for the given SDE using Milstein method
  S.milstein = matrix(ncol=(N.MC.central + N.MC.band*N.S0), nrow=(N.t+1)) # pre−allocates matrix for S
  for(.S0 in S0.array.without.central){ # for each value in the initial stock price’s array, except the original one

    # generates N.MC.band stock price paths
    j = which(S0.array.without.central==.S0) # obtain index of current stock price in array without central
    indexes = ((j-1)*N.MC.band+1):(min(j*N.MC.band, N.MC.band*N.S0))
    S.milstein[,indexes] = SDE.solve(t0=t0, T=T, X0=.S0, N.t=N.t, drift=drift, sigma=sigma, method="milstein",
                                     Z=Z[,indexes])
  }
  # finally, for the original initial stock price, generates N.MC.central stock price paths
  indexes = (N.MC.band*N.S0+1):(N.MC.central + N.MC.band*N.S0)
  S.milstein[,indexes] = SDE.solve(t0=t0, T=T, X0=S0, N.t=N.t, drift=drift, sigma=sigma, method="milstein", Z
                                   =Z[,indexes])
}
# bond price at time T
bond.T = bond.value.T(bond.params, r)
# pricing
# pricing with LSMC algorithm
if(calculate$LSMC$euler){ # calculate LSMC/Euler price only if it is to be calculated
  ret.LSMC.euler = price.MC(S=S.euler, S0=S0, r=r, q=q, t0=t0, T=T, N.t=N.t, N.MC=(N.MC.central + N.MC.band*N.S0),
                            backward.induction.method="LSMC", price.fun=price.fun.LSMC, bond.params=bond.params,
                            regression.method="bouchard", regression.intervals=50)
  price.LSMC.euler = ret.LSMC.euler$MC
  exercise.LSMC.euler = ret.LSMC.euler$exercise
}
if(calculate$LSMC$milstein){ # calculate LSMC/Milstein price only if it is to be calculated
  ret.LSMC.milstein = price.MC(S=S.milstein, S0=S0, r=r, q=q, t0=t0, T=T, N.t=N.t, N.MC=(N.MC.central + N.MC.band*N.S0),
                               backward.induction.method="LSMC", price.fun=price.fun.LSMC, bond.params=bond.params,
                               regression.method="bouchard", regression.intervals=50)
  price.LSMC.milstein = ret.LSMC.milstein$MC
  exercise.LSMC.milstein = ret.LSMC.milstein$exercise
}
# pricing with HMC algorithm
if(calculate$HMC$euler){ # calculate HMC/Euler price only if it is to be calculated
  ret.HMC.euler = price.MC(S=S.euler, S0=S0, r=r, q=q, t0=t0, T=T, N.t=N.t, N.MC=(N.MC.central + N.MC.band*N.S0),
                           backward.induction.method="HMC", price.fun=price.fun.HMC, bond.params=bond.params,
                           regression.method="bouchard", regression.intervals=50)
  price.HMC.euler = ret.HMC.euler$MC
  exercise.HMC.euler = ret.HMC.euler$exercise
}
if(calculate$HMC$milstein){ # calculate HMC/Milstein price only if it is to be calculated
  ret.HMC.milstein = price.MC(S=S.milstein, S0=S0, r=r, q=q, t0=t0, T=T, N.t=N.t, N.MC=(N.MC.central + N.MC.band*N.S0),
                              backward.induction.method="HMC", price.fun=price.fun.HMC, bond.params=bond.params,
                              regression.method="bouchard", regression.intervals=50)
  price.HMC.milstein = ret.HMC.milstein$MC
  exercise.HMC.milstein = ret.HMC.milstein$exercise
}
  
# Black’s formula pricing
price.black_ = price.black(r=r, q=q, t=0, T=T, S=S0, sigma=.sigma, bond.params=bond.params)

# combine obtained prices into one data structure
price.LSMC_ = list()
if(calculate$LSMC$euler){ # adds LSMC/Euler price to structure if it was to be computed
  price.LSMC_$euler = price.LSMC.euler$value
}
if(calculate$LSMC$milstein){ # adds LSMC/Milstein price to structure if it was to be computed
  price.LSMC_$milstein = price.LSMC.milstein$value
}

price.HMC_ = list()
if(calculate$HMC$euler){ # adds HMC/Euler price to structure if it was to be computed
  price.HMC_$euler = price.HMC.euler$value
}
if(calculate$HMC$milstein){ # adds HMC/Milstein price to structure if it was to be computed
  price.HMC_$milstein = price.HMC.milstein$value
}

price = list(black=price.black_, LSMC=price.LSMC_, HMC=price.HMC_)
  
  
  

