


## 1. DATA GENERATING FUNCTION - returns dataframe of treatment groups, 1 covariate and outcomes
## generates data from two overlapping uniform distributions

create_data_unif = function(
                       n=50,
                       ratio = 2, # ration of controls to treated
                       min_T=0.1, # min value treated group
                       max_T=1, # max value treated group
                       min_C=0, 
                       max_C=0.9, 
                       err_sd = 0.5,
                       beta1 = 2, #slope of treated
                       beta2 = 2, #slope of control
                       effect = 1){
  # create covariate
  t<- runif(n,min_T,max_T)
  c<- runif(ratio*n,min_C, max_C )
  x1 <- c(t,c) 
  
  #create outcome variable
  y<- as.numeric(n+(ratio*n))
  err <- rnorm(length(y),0,err_sd)
  y[1:n] <- effect+ beta1*t  # intercept (effect) + linear coefficient + error
  y[(n+1):(n+ratio*n)] <- beta2*c  # 0 intercept + linear coefficient + error 
  y <- y + err # add error to all outcomes
  
  
  #create treatment classes
  treatment <- c(rep("T",n),rep("C",ratio*n)) # create treatment groups
  
  #combine as dataframe
  df <- as.data.frame(cbind(treatment,x1,y))
  df$x1 <- as.numeric(df$x1)
  df$y = as.numeric(df$y)
  
  return(df)
}



#### 2. GET EFFECT FUNCTION: function to estimate the effect of a treatment using the matched data using a linear model
#### Estimates the ATT (usually, when no treatments are discarded, otherwise ATM, effect of treatment on matched data)


get_effect <- function(n=50, #parameters are for the create data function
                       ratio = 2, 
                       min_T=0.1, 
                       max_T=1,
                       min_C=0, 
                       max_C=0.9, 
                       err_sd = 0.5,
                       beta1 = 2, 
                       beta2 = 2, 
                       effect = 1,
                       replace = FALSE,  
                       caliper = NULL,
                       distance = "mahalanobis",
                       method = "nearest"){ 
  
  df <- create_data_unif(n=n,
                         ratio = ratio,
                         min_T = min_T,
                         max_T = max_T,
                         min_C = min_C,
                         max_C = max_C,
                         err_sd = err_sd,
                         beta1 = beta1,
                         beta2 = beta2,
                         effect = effect)
  
  df$treatment <- as.factor(df$treatment) # change to factor to use caliper (assumes logit)
  
  m.out <- matchit(treatment~x1,data = df,
                   method = method, 
                   distance = distance, 
                   replace = replace, 
                   caliper = caliper) 
  
  # add new parts here to return number unique matches (i.e. num controls used) 
  # and num discarded from treatment (t not used)
  controls_used <- length(unique(m.out$match.matrix))
  discards <- sum(is.na(m.out$match.matrix))
  
  fit1 <- lm(y ~ treatment, data = match.data(m.out), weights = weights) #estimating the effect
  
  output <-c(fit1$coefficients[2], controls_used, discards)
  return(output)
}



### 3. MONTE CARLO FUNTCION: uses the data generating and get effect functions

g.x <- function(m=100, # number of simulations
                n=50, #other parameters are for the create data function
                ratio = 2, 
                min_T=0.1, 
                max_T=1,
                min_C=0, 
                max_C=0.9, 
                err_sd = 0.5,
                beta1 = 2, 
                beta2 = 2, 
                effect = 1,
                replace = FALSE,
                caliper = NULL,
                method = "nearest",
                distance = "mahalanobis"){
  
  reps <- replicate(m, get_effect(
    n=n,
    ratio = ratio,
    min_T = min_T,
    max_T = max_T,
    min_C = min_C,
    max_C = max_C,
    err_sd = err_sd,
    beta1 = beta1,
    beta2 = beta2,
    effect = effect,
    replace = replace,
    caliper = caliper,
    method = method,
    distance = distance))
  return(reps)
}



### 4. CREATE NORMAL DATA FUNCTION
create_norm_data = function(
  n=50,
  ratio = 2, # ration of controls to treated
  mean_C = 0,
  mean_T = 0.1,
  sd_T = sqrt(1/12), # same as for uniform(0,1)
  sd_C =sqrt(1/12),
  err_sd = 0.5,
  beta1 = 2, #slope of treated
  beta2 = 2, #slope of control
  effect = 1)
  {
  # create covariate
  t<- rnorm(n,mean = mean_T, sd= sd_T)
  c<- rnorm(n*ratio,mean = mean_C, sd= sd_C)
  x1 <- c(t,c) 
  
  #create outcome variable
  y<- numeric(n+ratio*n)
  err <- rnorm(length(y),0,err_sd)
  y[1:n] <- effect+ beta1*t  # intercept (effect) + linear coefficient + error
  y[(n+1):(n+ratio*n)] <- beta2*c  # 0 intercept + linear coefficient + error 
  y <- y + err # add error to all outcomes
  
  #create treatment classes
  treatment <- c(rep("T",n),rep("C",ratio*n)) # create treatment groups
  
  #combine as dataframe
  df <- as.data.frame(cbind(treatment,x1,y))
  df$x1 <- as.numeric(df$x1)
  df$y = as.numeric(df$y)
  
  return(df)
}


#### 5 NORMAL GET EFFECT FUNCTION

get_effect_norm <- function(n=50, #parameters are for the create data function
                       ratio = 2, 
                       mean_C = 0,
                       mean_T = 0.1,
                       sd_T = sqrt(1/12),
                       sd_C = sqrt(1/12),
                       err_sd = 0.5,
                       beta1 = 2, 
                       beta2 = 2, 
                       effect = 1,
                       replace = FALSE,  
                       caliper = NULL,
                       distance = "mahalanobis",
                       method = "nearest"){ 
  
  df <- create_norm_data(n=n,
                         ratio = ratio,
                         mean_C = mean_C,
                         mean_T = mean_T,
                         sd_T = sd_T,
                         sd_C = sd_C,
                         err_sd = err_sd,
                         beta1 = beta1,
                         beta2 = beta2,
                         effect = effect)
  
  df$treatment <- as.factor(df$treatment) # change to factor to use caliper (assumes logit)
  
  m.out <- matchit(treatment~x1,data = df,
                   method = method, 
                   distance = distance, 
                   replace = replace, 
                   caliper = caliper) 
  
  fit1 <- lm(y ~ treatment, data = match.data(m.out), weights = weights) #estimating the effect
  
  # add new parts here to return num unique matches (i.e. num controls used) 
  # and num discarded from treatment (t not used)
  controls_used <- length(unique(m.out$match.matrix))
  discards <- sum(is.na(m.out$match.matrix))
  
  fit1 <- lm(y ~ treatment, data = match.data(m.out), weights = weights) #estimating the effect
  
  output <-c(fit1$coefficients[2], controls_used, discards)
  return(output)
  #return(fit1$coefficients[2])
}


### 6. MONTE CARLO FUNTCION NORMAL: uses the data generating and get effect functions

g.x_norm <- function(m=100, # number of simulations
                n=50, #other parameters are for the create data function
                ratio = 2, 
                mean_C = 0,
                mean_T = 0.1,
                sd_T = sqrt(1/12),
                sd_C = sqrt(1/12),
                err_sd = 0.5,
                beta1 = 2, 
                beta2 = 2, 
                effect = 1,
                replace = FALSE,
                caliper = NULL,
                method = "nearest",
                distance = "mahalanobis"){
  
  reps <- replicate(m, get_effect_norm(
    n=n,
    ratio = ratio,
    mean_C = mean_C,
    mean_T = mean_T,
    sd_T = sd_T,
    sd_C =sd_C,
    err_sd = err_sd,
    beta1 = beta1,
    beta2 = beta2,
    effect = effect,
    replace = replace,
    caliper = caliper,
    method = method,
    distance = distance))
  return(reps)
}



### CREATE DATA BETA - create two skewed distributions, between 0 and 1, similar to propensity scores???

create_data_beta = function(
  n=10,
  ratio = 2, # ration of controls to treated
  shape1 = 5, 
  shape2 = 2, 
  shape3 = 2,
  shape4 = 5,
  err_sd = 0.5,
  beta1 = 2, #slope of treated
  beta2 = 2, #slope of control
  effect = 1){
  # create covariate
  t<- rbeta(n=n,shape1 = shape1, shape2 = shape2)
  c<- rbeta(ratio*n, shape1 = shape3, shape2 = shape4)
  x1 <- c(t,c) 
  
  #create outcome variable
  y<- numeric(n+ratio*n)
  err <- rnorm(length(y),0,err_sd)
  y[1:n] <- effect+ beta1*t  # intercept (effect) + linear coefficient + error
  y[(n+1):(n+ratio*n)] <- beta2*c  # 0 intercept + linear coefficient + error 
  y <- y + err # add error to all outcomes
  
  
  #create treatment classes
  treatment <- c(rep("T",n),rep("C",ratio*n)) # create treatment groups
  
  #combine as dataframe
  df <- as.data.frame(cbind(treatment,x1,y))
  df$x1 <- as.numeric(df$x1)
  df$y = as.numeric(df$y)
  
  return(df)
}

