# Author: Leandro Loriato (llrt at impa.br)
#
# Auxiliary functions to calculate the payoff of a Convertible Bond given the relevant parameters

# legend of exercise possibilities
exercise.legend = c("continuation"=0,
                     "voluntary_conversion"=1,
                     "put"=2,
                     "call"=3,
                     "forced_conversion"=4,
                     "redemption"=5)

# given a time step, computes its exercise time restriction bitmap,
# evaluating which exercise actions' time restrictions are met
exercise.time.bitmap = function(t, bond.params){
  ret = list()
  ret$conversion = FALSE
  ret$put = FALSE
  ret$call = FALSE
  ret$redemption = FALSE

  # conversion
  if(t %in% bond.params$conversion$dates){
    ret$conversion = TRUE
  }
  # put
  if(bond.params$put$present && t %in% bond.params$put$dates){
    ret$put = TRUE
  }
  # call
  if(bond.params$call$present && t %in% bond.params$call$dates){
    ret$call = TRUE
  }
  # redemption
  if(t==T){
    ret$redemption = TRUE
  }

  invisible(ret)
}
# for performance, pre−compile CPU intensive function
exercise.time.bitmap = cmpfun(exercise.time.bitmap)

# calculates a payoff value for the convertible bond at a time, given a stock price,
# a reference value and other relevant parameters
payoff = function(t, S.vector, r, continuation.value, exercise.time.bitmap, bond.params){

  # initialize auxiliary variables
  S = tail(S.vector, n=1)
  action = 'continuation'
  conversion.value = bond.params$conversion$ratio(t, S.vector) * S
  call.value = bond.params$call$strike(t, S.vector)
  put.value = bond.params$put$strike(t, S.vector)
  redemption.value = redemption(bond.params, r)
  
  if(bond.params$conversion$restriction(t, S.vector) && exercise.time.bitmap$conversion && conversion.value > continuation.value){
    if(bond.params$put$restriction(t, S.vector) && exercise.time.bitmap$put && put.value > conversion.value)
    {
      action = 'put'
    } else{
      action = 'voluntary_conversion'
    }
  } else if(bond.params$put$restriction(t, S.vector) && exercise.time.bitmap$put && put.value > continuation.value){
    action = 'put'
  } else if(bond.params$call$restriction(t, S.vector) && exercise.time.bitmap$call && continuation.value > call.value){
    if(exercise.time.bitmap$conversion && conversion.value > call.value){
      action = 'forced_conversion'
    } else{
      action = 'call'
    }
  } else if(exercise.time.bitmap$redemption && redemption.value > conversion.value){
    action = 'redemption'
  }

  # based on which action takes place, evaluates corresponding payoff
  switch(action,
         'voluntary_conversion' = { # voluntary conversion takes place
           value = conversion.value
         },
         'put'={ # put exercise takes place
           value = put.value
         },
         'call'={ # call exercise takes place
           value = call.value
         },
         'forced_conversion'={ # forced conversion takes place
           value = conversion.value
         },
         'redemption'={ # redemption takes place
           value = redemption.value
         },
         'continuation'={ # none of other actions take place, investor holds convertible bond for one more period
           value = continuation.value
         }
  )

  action = exercise.legend[[action]]

  ret = list(value=value, action=action)
  invisible(ret)
}
# for performance, pre−compile CPU intensive function
payoff = cmpfun(payoff)