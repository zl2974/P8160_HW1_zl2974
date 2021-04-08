
## soft threshold update
soft_threshold =
  function(beta, lambda) {
    if (abs(beta)>lambda){
      if (beta >0) return(beta - lambda)
      return(beta +lambda)
    }
    return(0)
  }

classification = 
  function(theta_vec,new_x,new_y = NA,...){
    
    x = new_x
    
    p_prev = exp(x %*% theta_vec) / (1 + exp(x %*% theta_vec)) #  row matrix
    
    p_prev[p_prev<1e-4] = 1e-4
    
    p_prev[p_prev > 1-1e-4] = 0.999
    
    p_prev[p_prev == 0] = 1e-4
    
    cl = ifelse(p_prev>0.5,1,0)
    
    if (!all(is.na(new_y))) {
      if (is.factor(new_y))
        new_y = new_y %>% as.numeric() - 1
      err = mean(cl != new_y)
    } else {err = NA}
    
    return(list(prob = p_prev,
                class = cl,
                error_rate = err))
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
    
    objective = 
      -sum(y * log(p_prev/(1-p_prev))+log(1-p_prev))
      #sum(w_prev * (z_prev - x %*% theta_vec) ^ 2) / (2 * length(y))
    
    return(list(
      objective = objective,
      p = p_prev,
      w = w_prev,
      z = z_prev
    ))
  }


### Gaussian Family

gfun = function(y,x,theta_vec){
  z = y
  w = rep(1/length(y),length(y))
  objective = sum(w *(y - x %*% theta_vec) ^ 2) / 2
  return(
    list(z=z,w=w,objective=objective)
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
           return_scale = T,
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
    
    cur_result = family_fun(y, x, theta_vec)$objective
    
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
        if (intercept & i == 1) {
          cur_theta =
            cur_theta / sum(w_prev * (x[, i] ^ 2))
        } else {
          cur_theta = soft_threshold(cur_theta,
                                     lambda) / sum(w_prev * (x[, i] ^ 2))
        }
        
        #update theta
        theta_vec[[i]] = cur_theta
      }
      cur_result = 
        (family_fun(y, x, theta_vec)$objective)/length(y) + lambda * sum(abs(theta_vec))
    }
    
    if (!return_scale) theta_vec = theta_vec / xsd
    return(list(coefficient =theta_vec,
                objective = cur_result))
  }

## Cross validation

cv_Qlasso = 
  function(y,
           x,
           family_fun,
           lambda = NA,
           number = 5,
           intercept = T,
           return_scale = T,
           obj = "test_error") {
    
    result = tibble()
    
    lambda = as.vector(lambda)
    
    if (all(is.na(lambda))) {
      Y = y
      if (is.factor(y)) Y = y %>% as.numeric()-1
      Y = as.vector(Y)
      max_lambda = length(Y)*max(colSums(diag(Y) %*% scale(x)))
      lambda = c(0,exp(seq(log(max_lambda),-5,length=20)),1e+3)
    }
    
    block_length = length(y)/number
    
    test_index = sample(rep(1:number,len = length(y)))
    
    for (lam in lambda) {
      for (i in 1:number) {
        #create training partition
        trainX = x[-which(test_index==i),]
        
        trainY = y[-which(test_index==i)]
        
        #train
        trainF = Qlasso(trainY,trainX,family_fun,lambda = lam,intercept = intercept,return_scale = return_scale)
        prev_coef = trainF$coefficient
        
        #test
        testX = x[which(test_index==i),]
        
        testY = y[which(test_index==i)]
      
        if (return_scale) testX = scale(testX)
        
        if (intercept) {testX = cbind(rep(1,length(testY)),testX)}
        
        testF = family_fun(testY,testX,prev_coef)
        
        penalty = sum(abs(prev_coef))
        
        if (obj == "test_error"){
          objective = classification(prev_coef,testX,testY)$err
        } else {
          if(obj == "deviance"){
            objective = testF$obj
          }
        }
        
        result = 
          result %>% 
          rbind(expand.grid(coefficient = trainF$coefficient,
                            lambda = lam,
                            objective = 
                              objective
                              ) %>% 
                  cbind(tibble(term = (1:length(trainF$coefficient))-1) %>% 
                          mutate_all(as.character))
          )
      }
    }
    
    # take the mean of everything
    result = result %>% 
      group_by(lambda,term) %>% 
      summarise(sd = sd(objective),
                across(c(coefficient,objective),mean))%>% 
      arrange((objective)) %>% 
      relocate(objective,lambda) %>% 
      nest(coef = c(coefficient,term))
    
    coef = result[1,] %>% 
      unnest(c(coef)) %>% 
      pull(coefficient)
    
    lambda = result[1,"lambda"] %>% as.numeric()
    
    objective = result[1,"objective"] %>% as.numeric()
    
    onese = min(result[,"objective"] %>% unlist())+result[1,"sd"] %>% as.numeric()

    result =result %>% 
      mutate(`1se` = objective <= onese)
    
    return(list(
      coefficient = coef,
      objective = objective,
      lambda = lambda,
      cvtable = result
    ))
  }
