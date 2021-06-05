# Author: Leandro Loriato (llrt at impa.br)
#
# Auxiliary functions to calculate the value of a bond at a given time
discounted.coupon = function(bond.params, t, r){
  payments = c() # initialize vector of made payments as empty
  for (t_ in bond.params$coupon$dates){
    # adds discounted value of accrued coupons to payments vector
    payments = c(payments, exp(-r*(t_-t)) * bond.params$coupon$rate * bond.params$principal)
  }
  value = sum(payments)
  invisible(value)
}
# calculates the redemption value of a bond, given relevant parameters
redemption = function(bond.params, r){
  # redemeed value is the principal times the redemption ratio
  value = bond.params$principal * bond.params$redemption.ratio
  invisible(value)
}
# for performance, pre−compile CPU intensive function
redemption = cmpfun(redemption)
# calculates the value of a bond at time T, given relevant parameters
bond.value.T = function(bond.params, r){
  # value of the bond at time T is the redemption value
  value = redemption(bond.params, r)
  invisible(value)
}