
## soft threshold update
soft_threshold =
  function(beta, lambda) {
    if (abs(beta)>lambda){
      if (beta >0) return(beta - lambda)
      return(beta +lambda)
    }
    return(0)
  }


## Family Function

### Binomial Family
bfun =
  function(y, x, theta_vec) {
    if (is.factor(y))
      y = y %>% as.numeric() - 1
    
    p_prev = exp(x %*% theta_vec) / (1 + exp(x %*% theta_vec)) #  row matrix
    
    p_prev[p_prev<1e-4] = 1e-4
    
    p_prev[p_prev > 1-1e-4] = 0.999
    
    p_prev[p_prev == 0] = 1e-4
    
    w_prev = p_prev * (1 - p_prev) # row matrix
    
    w_prev[w_prev<1e-5] = 1e-5
    
    z_prev = x %*% theta_vec + (y - p_prev) / w_prev # row * 1 matrix
    
    loglink = 
      -sum(y * log(p_prev/(1-p_prev))+log(1-p_prev))
      #sum(w_prev * (z_prev - x %*% theta_vec) ^ 2) / (2 * length(y))
    
    return(list(
      loglink = loglink,
      p = p_prev,
      w = w_prev,
      z = z_prev
    ))
  }


### Gaussian Family

gfun = function(y,x,theta_vec){
  z = y
  w = rep(1/length(y),length(y))
  loglink = sum(w *(y - x %*% theta_vec) ^ 2) / 2
  return(
    list(z=z,w=w,loglink=loglink)
  )
}


## Lasso Function ###########
#############################
Qlasso =
  function(y,
           x,
           family_fun,
           lambda = 0,
           theta_vec = NA,
           maxiter = 500,
           tol = 1e-6,
           intercept = T,
           ...) {
    # data preparation
    if (is.factor(y))
      y = y %>% as.numeric() - 1
    y = y
    x = scale(x)
    xsd = attr(x, "scaled:scale") %>% as.vector()
    if (intercept) {
      x = cbind(rep(1, length(y)), x) # add alpha
      xsd = append(1, xsd)
    }
    
    # Manually choosing starting value
    if (!is.na(theta_vec)) {
      if (length(theta_vec) != ncol(x))
        theta_vec = append(rep(0, ncol(x) - length(theta_vec)), theta_vec)
    } else
      theta_vec = rep(0, ncol(x))
    
    #update part
    iter = 0
    
    cur_result = family_fun(y, x, theta_vec)$loglink
    
    if (abs(cur_result) == Inf|is.na(cur_result))
      stop("Diverge at starting value")
    
    prev_result = -Inf
    
    while (abs(cur_result - prev_result) > tol
           && iter < maxiter) {
      iter = iter + 1
      
      prev_result = cur_result
      
      for (i in 1:length(theta_vec)) {
        #weight
        fun_prev = family_fun(y, x, theta_vec)
        z_prev = fun_prev$z # working response
        w_prev = fun_prev$w # weight
        
        #update
        cur_theta =sum((z_prev - x[, -i] %*% theta_vec[-i]) * x[, i] * w_prev)
        
        #soft-threshold
        if (i > 1) {
          cur_theta =
            soft_threshold(cur_theta,
                           lambda) / sum(w_prev * (x[, i] ^ 2))
        } else {
          cur_theta = cur_theta / sum(w_prev * (x[, i] ^ 2))
        }
        
        #update theta
        theta_vec[[i]] = cur_theta
      }
      cur_result = family_fun(y, x, theta_vec)$loglink + lambda * sum(abs(theta_vec))
    }
    return(list(coefficient = theta_vec / xsd,
                loglikelihood = cur_result))
  }

## Cross validation

cv_Qlasso = 
  function(y,x,family_fun,lambda=NA,number = 5,intercept = T){
    
    result = tibble()
    
    lambda = as.vector(lambda)
    
    if (all(is.na(lambda))) {
      Y = y
      if (is.factor(y)) Y = y %>% as.numeric()-1
      Y = as.vector(Y)
      max_lambda = length(Y)*max(colSums(diag(Y) %*% scale(x)))
      lambda = exp(seq(log(max_lambda),-10,length=20))
    }
    
    block_length = length(y)/number
    
    test_index = sample(rep(1:number,len = length(y)))
    
    for (lam in lambda) {
      for (i in 1:number) {
        #create training partition
        trainX = x[-which(test_index==i),]
        testX = x[which(test_index==i),]
        trainY = y[-which(test_index==i)]
        testY = y[which(test_index==i)]
        
        #train
        trainF = Qlasso(trainY,trainX,family_fun,lambda = lam,intercept = intercept)
        if (intercept) {prev_coef = trainF$coefficient[-1]
        }else {prev_coef = trainF$coefficient}
        
        #test
        testF = family_fun(testY,testX,prev_coef)
        result = 
          result %>% 
          rbind(expand.grid(coefficient = trainF$coefficient,
                            lambda = lam,
                            loglikelihood = testF$logli+ lam*sum(abs(trainF$coefficient))) %>% 
                  cbind(tibble(term = (1:length(trainF$coefficient))-1) %>% 
                          mutate_all(as.character))
          )
      }
    }
    
    # take the mean of everything
    result = result %>% 
      group_by(lambda,term) %>% 
      summarise(across(c(coefficient,loglikelihood),mean))%>% 
      arrange((loglikelihood)) %>% 
      nest(coef = c(coefficient,term))
    
    coef = result[1,] %>% 
      unnest() %>% 
      pull(coefficient)
    
    lambda = result[1,"lambda"] %>% as.numeric()
    
    loglikelihood = result[1,"loglikelihood"] %>% as.numeric()
    
    return(list(
      coefficient = coef,
      loglikelihood = loglikelihood,
      lambda = lambda,
      cvtable = result
    ))
  }