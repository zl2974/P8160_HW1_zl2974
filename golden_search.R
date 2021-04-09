library(tidyverse)

goose_egg = 
  function(
    fun,
    left = NULL,
    right = NULL,
    range = NULL,
    ratio = 0.618,
    tol = 10e-4,
    ...
  ){
    if (!any(left,right)){
      left = range[1]
      right = range[2]
    }
    
    mid_1 = left + ratio*(right - left)
    
    f_mid_1 = fun(mid_1)
    
    mid_2 = mid_1 + ratio*(right-mid_1)
    
    f_mid_2 = fun(mid_2)
    
    f_left = fun(left)
    
    f_right = fun(right)
    
    i = 1
    
    while (abs(f_left - f_right)>tol && i<1000){
      i = i + 1
      if (f_mid_1 < f_mid_2) {
        f_left = f_mid_1
        left = mid_1
      } else {
        f_right = f_mid_2
        right = mid_2
      }
      mid_1 = left + ratio * (right - left)
      f_mid_1 = fun(mid_1)
      mid_2 = mid_1 + ratio * (right - mid_1)
      f_mid_2 = fun(mid_2)
    }
    return(mean(mid_1,mid_2))
  }