rm(list=ls())
library(shiny)

options(digits = 3)



## Power Functions for Sample Size Approaches

##Testing Approach
## rho.0 <- value under null
## rho.A <- value under alternate

TestingSS.F<-function(rho.0, rho.A, k, n, alpha=0.05){
  t.0 <- 1+k*rho.0/(1-rho.0)
  t.1 <- 1+k*rho.A/(1-rho.A)
  t.R <- t.0/t.1
  est.power <- 1-pf(t.R*qf(1-alpha, n-1, n*(k-1)), n-1, n*(k-1))
  return(est.power)
}


TestingSS.ZF<-function(rho.0, rho.A, k, n, alpha=0.05){
  muZ.0 <- log(1+k*rho.0/(1-rho.0))/2
  muZ.1 <- log(1+k*rho.A/(1-rho.A))/2
  sdZ <- sqrt(k/(2*(k-1)*(n-1)))
  
  est.power <- 1-pnorm(qnorm(1-alpha)-(muZ.1-muZ.0)/sdZ, mean=0, sd=1)
  return(est.power)
}


TestingSS.Z <- function(rho.0, rho.A, k, n, alpha=0.05){
  muZ.0 <- log((1+rho.0)/(1-rho.0))/2
  muZ.1 <- log((1+rho.A)/(1-rho.A))/2
  
  sd.Swiger <- function(n,k,rho){sqrt(2*(n*k-1)*(1+(k-1)*rho)^2/(k^2*(1+rho)^2*(n*k-n)*(n-1)))}
  sd.Fisher <-function(n,k,rho){sqrt(2*(1+(k-1)*rho)^2/(n*k*(1+rho)^2*(k-1)))}
  sd.Zerbe <- function(n,k,rho){sqrt(2*(n*k-n)^2*(n*k-3)*(1+(k-1)*rho)^2/(k^2*(1+rho)^2*(n*k-n-2)^2*(n*k-n-4)*(n-1)))}
  
  est.power <- function(fun = sd.Swiger){
    1-pnorm(qnorm(1-alpha)-(muZ.1/fun(n,k,rho.A)-muZ.0/fun(n,k,rho.A)), mean=0, sd=1)}
  
  return(c(est.power(fun=sd.Swiger), est.power(fun=sd.Fisher), est.power(fun=sd.Zerbe)))
  
}



TestingSS.W <- function(rho.0, rho.A, k, n, alpha=0.05){
  sd.Swiger <- function(n,k,rho){sqrt(2*(n*k-1)*(1+(k-1)*rho)^2*(1-rho)^2/(k^2*(n*k-n)*(n-1)))}
  sd.Fisher <-function(n,k,rho){sqrt(2*(1+(k-1)*rho)^2*(1-rho)^2/(n*k*(k-1)))}
  sd.Zerbe <- function(n,k,rho){sqrt(2*(n*k-n)^2*(n*k-3)*(1+(k-1)*rho)^2*(1-rho)^2/(k^2*(n*k-n-2)^2*(n*k-n-4)*(n-1)))}
  
  est.power <- function(fun = sd.Swiger){
    1-pnorm((qnorm(1-alpha)-(rho.A-rho.0)/fun(n,k,rho.A)), mean=0, sd=1)}
  
  return(c(est.power(fun=sd.Swiger), est.power(fun=sd.Fisher), est.power(fun=sd.Zerbe)))
  
}











##Assurance Probability Approach
## rho <- target ICC
## width <- value of width around rho




AssuranceSS.F<-function(rho, width, k, n, alpha=0.05){
  d1 <- n-1
  d2<- n*(k-1)
  
  
  F.t<-function(rho) {1+k*rho/(1-rho)} 
  Fu <- qf(1-alpha/2, n-1, n*(k-1))
  Fl <- qf(alpha/2, n-1, n*(k-1))
  F.rho <- F.t(rho)
  F.L <- Fl * (k-width+k*width)
  F.U <- Fu * (k+width-k*width)
  
  
  sq.pt <-max(0,(Fu-Fl)*(F.U^2/Fu-F.L^2/Fl))
  
  f.r.1<- (F.U-F.L+sqrt(sq.pt))/(2*width)
  f.r.2<- (F.U-F.L-sqrt(sq.pt))/(2*width)
  est_power <-1-pf(f.r.1/F.t(rho),d1, d2) + pf(f.r.2/F.t(rho),d1, d2)
  return(est_power)
}


AssuranceSS.ZF<-function(rho, width, k, n, alpha=0.05){
  F.t<-function(rho) {0.5*log(1+k*rho/(1-rho))}
  var.z <- 0.5*(1/(n-1)+1/(n*(k-1)))
  Fu <- qf(1-alpha/2, n-1, n*(k-1))
  Fl <- qf(alpha/2, n-1, n*(k-1))
  F.L <- Fl * (k-width+k*width)
  F.U <- Fu * (k+width-k*width)
  
  sq.pt <-max(0,(Fu-Fl)*(F.U^2/Fu-F.L^2/Fl))
  
  f.r.1<- (F.U-F.L+sqrt(sq.pt))/(2*width)
  f.r.2<- (F.U-F.L-sqrt(sq.pt))/(2*width)
  
  est_power <- pnorm((F.t(rho)-0.5*log(f.r.1))/sqrt(var.z)) + 1- pnorm((F.t(rho)-0.5*log(f.r.2))/sqrt(var.z))
  
  return(est_power)
}


AssuranceSS.W <- function(rho, width, k, n, alpha=0.05){
  width = width/(2*qnorm(1-alpha/2))
  f.rho <- (1-rho)*(1+(k-1)*rho)
  f.d.rho <- abs(k-2-2*(k-1)*rho)
  f.Swiger <- function(n){(k^2*(n*k-n)*(n-1))/(2*(k*n-1))}
  f.Fisher <-function(n){n*k*(k-1)/2}
  f.Zerbe <- function(n){(k^2*(n*k-n-2)^2*(n*k-n-4)*(n-1))/(2*(n*k-n)^2*(n*k-3))}
  
  est.power <- function(fun){
    1-pnorm((f.rho*sqrt(fun(n=n))-width*fun(n=n))/(f.rho*f.d.rho), mean=0, sd=1)}
  return(c(est.power(fun=f.Swiger), est.power(fun=f.Fisher), est.power(fun=f.Zerbe)))
}




AssuranceSS.Z<-function(rho, width, k, n, alpha=0.05){
  F.t<-0.5*log((1+rho)/(1-rho))
  cosh.inv<-function(x){
    x<-abs(x)
    return(log(x+sqrt(abs(x^2-1))))
  }
  sd.Swiger <- function(n){sqrt(2*(n*k-1)*(1+(k-1)*rho)^2/(k^2*(1+rho)^2*(n*k-n)*(n-1)))}
  sd.Fisher <- function(n){sqrt(2*(1+(k-1)*rho)^2/(n*k*(1+rho)^2*(k-1)))}
  sd.Zerbe <- function(n){sqrt(2*(n*k-n)^2*(n*k-3)*(1+(k-1)*rho)^2/(k^2*(1+rho)^2*(n*k-n-2)^2*(n*k-n-4)*(n-1)))}
  
  
  est.power <- function(fun){
    del <- 2*qnorm(1-alpha/2)*fun(n=n)
    cosh_val <- max(1,2*sinh(del)/width-cosh(del))
    f.r<- (F.t-0.5*cosh.inv(cosh_val))/fun(n=n)
    return(pnorm(f.r))
  }
  return(c(est.power(fun=sd.Swiger), est.power(fun=sd.Fisher), est.power(fun=sd.Zerbe)))
}


### Sample sizes are obtained at which the estimated power reaches observed power

SampleSize<-function(criterion_specification, 
                     ## For Assurance Probability, 
                     # named list of (rho, target_width and gamma_value) is to be provided
                     # For Testing method,
                     # named list of(rho.A, rho.0 and beta_value) is to be provided,
                     k,  
                     alpha = 0.05, 
                     step=1,
                     nmax=5e4){
  procedure = NULL
  if ("target_width" %in% names(criterion_specification)){
    rho = criterion_specification[["rho"]]
    width = criterion_specification[["target_width"]]
    pow = criterion_specification[["gamma_value"]]
    procedure = "Assurance"
  } else if("rho.0" %in% names(criterion_specification)){
    rho.0 = criterion_specification[["rho.0"]]
    rho.A = criterion_specification[["rho.A"]]
    pow = criterion_specification[["beta_value"]]
    procedure = "Testing"
    if(rho.A<=rho.0){
      cat("Error in Parameter Specifications")
      break
    }
  } else{
    "Error in Parameter Specifications"
  }
  
  if(procedure=="Testing"){
    power<- function(n){
      c(TestingSS.W(rho.0=rho.0, rho.A=rho.A, k=k, n=n, alpha=alpha),
        TestingSS.F(rho.0=rho.0, rho.A=rho.A, k=k, n=n, alpha=alpha),
        TestingSS.Z(rho.0=rho.0, rho.A=rho.A, k=k, n=n, alpha=alpha),
        TestingSS.ZF(rho.0=rho.0, rho.A=rho.A, k=k, n=n, alpha=alpha))
    }
  }else{
    power<- function(n){
      c(AssuranceSS.W(rho=rho, width=width, k=k, n=n, alpha=alpha),
        AssuranceSS.F(rho=rho, width=width, k=k, n=n, alpha=alpha),
        AssuranceSS.Z(rho=rho, width=width, k=k, n=n, alpha=alpha),
        AssuranceSS.ZF(rho=rho, width=width, k=k, n=n, alpha=alpha))
    }
  }
  
  
  power_cont<- matrix(0, nrow=0, ncol=8)
  for(n in seq(4,nmax,step)){
    power_cont <- rbind(power_cont,power(n=n))
    if(all(tail(power_cont,1)>=(1-pow), na.rm=TRUE)){
      max_n <- n
      break
    }
  }
  
  samplesize<- apply(power_cont,2,function(x) seq(4,nmax,step)[which(x>(1-pow)&!is.na(x))][1])
  names(samplesize)<- c('Wald.S', 'Wald.F', 'Wald.Ze', 'F',
                        'Z.S', 'Z.F', 'Z.Ze', 'ZF')
  return(samplesize)
}










## Width of Confidence interval method

#Calculate ICC from variance
ICC_theory <- function(sigma_s,# Subject measurements standard deviation 
                       sigma_e# Error standard deviation
){
  return(sigma_s^2/(sigma_s^2+sigma_e^2))
}


#ICC Estimate from data
ICC_estimate<-function(data,  MLE=FALSE){
  data[is.na(data)]<-0
  k <- ncol(data)
  n <- nrow(data)
  
  MS_between <- 1/(n-1) * sum(k*(rowMeans(data)-mean(data))^2)
  MS_within <- 1/(n*k-n)* sum((data-rowMeans(data))^2)
  
  rho = (MS_between-MS_within)/(MS_between + (k-1)*MS_within)
  
  rho_MLE = (MS_between * (n-1)/n - MS_within)/
    (MS_between * (n-1)/n + (k-1)* MS_within)
  
  if (MLE==FALSE){return(rho)
  }else{return(rho_MLE)}
  
}

### Fisher approximation ###
var_rho.Fisher <- function(rho,
                           k,#number of raters
                           n #number of subjects
){
  return(2*(1-rho)^2 * (1+(k-1)*rho)^2/(k*n*(k-1)))
}


### Swiger approximation ###
var_rho.Swiger <- function(rho,
                           k,#number of raters
                           n #number of subjects
){
  return(2*(n*k-1)*(1-rho)^2 * (1+(k-1)*rho)^2/(k^2*(n*k-n)*(n-1)))
}

### Zerbe approximation ###
var_rho.Zerbe <- function(rho,
                          k,#number of raters
                          n #number of subjects
){
  return(2*(n*k-n)^2*(n*k-3)*(1-rho)^2 * 
           (1+(k-1)*rho)^2/(k^2*(n*k-n-2)^2*(n*k-n-4)*(n-1)))
}


### Wald Method (Wald.F, Wald.S, Wald.Ze) ###
Wald <- function(rho, ## Value of ICC estimate 
                 k, ## number of raters
                 n, ## number of subjects
                 alpha, ## level of uncertainty
                 variance_function = var_rho.Fisher
){
  zalpha_by_2 <- qnorm(1-alpha/2)
  
  Lower = rho - zalpha_by_2*sqrt(variance_function(rho,k,n))
  Upper = rho + zalpha_by_2*sqrt(variance_function(rho,k,n))
  
  
  return(c("Lower"=Lower, 
           "Upper"=Upper))
}


### F-distr Method (F) ###
F_method<- function(rho, ## Value of ICC estimate 
                    k, ## number of raters
                    n, ## number of subjects
                    alpha ## level of uncertainty 
){
  
  F_val <- (1+(k-1)*rho)/(1-rho)
  Fl <- qf(alpha/2,df1=n-1,df2=n*(k-1))
  Fu <- qf(1-alpha/2,df1=n-1,df2=n*(k-1))
  
  Lower = (F_val/Fu-1)/(F_val/Fu+k-1)
  Upper = (F_val/Fl-1)/(F_val/Fl+k-1)
  
  return(c("Lower"=Lower, 
           "Upper"=Upper))
}



### Method with Fisher Transformation on rho (Z.S,Z.F,Z.Ze) ###
Z_method <- function(rho, ## Value of ICC estimate 
                     k, ## number of raters
                     n, ## number of subjects
                     alpha, ## level of uncertainty
                     variance_function = var_rho.Fisher
                     
){
  
  zalpha_by_2 <- qnorm(1-alpha/2)
  Z_trans <- (1+rho)/(1-rho)
  VarZ = variance_function(rho,k,n)/((1+rho)^2*(1-rho)^2)
  
  a = 0.5* log(Z_trans)- zalpha_by_2*sqrt(VarZ)
  b = 0.5* log(Z_trans)+ zalpha_by_2*sqrt(VarZ)
  
  Lower = (exp(2*a)-1)/(exp(2*a)+1)
  Upper = (exp(2*b)-1)/(exp(2*b)+1)
  
  
  return(c("Lower"=Lower, 
           "Upper"=Upper))
}


### F-distr Method with Fisher Transformation (ZF.rho) ###
ZF <- function(rho, ## Value of ICC estimate 
               k, ## number of raters
               n, ## number of subjects
               alpha ## level of uncertainty
){
  
  zalpha_by_2 <- qnorm(1-alpha/2)
  F_val <- (1+(k-1)*rho)/(1-rho)
  VarF = 0.5*(1/(n-1) + 1/(n*k-n))
  
  a = 0.5* log(F_val)- zalpha_by_2*sqrt(VarF)
  b = 0.5* log(F_val)+ zalpha_by_2*sqrt(VarF)
  
  Lower = (exp(2*a)-1)/(exp(2*a)+k-1)
  Upper = (exp(2*b)-1)/(exp(2*b)+k-1)
  
  
  return(c("Lower"=Lower, 
           "Upper"=Upper))
}

CI_methods<-function(y, ## data matrix
                     alpha ## level of uncertainty
){
  rho = ICC_estimate(y)
  n <- dim(y)[1]
  k <- dim(y)[2]
  return(list("Parameters"=list("rho"=rho, "n"=n, "k"=k),
              "Confidence_intervals"=list("Wald.S" = Wald(rho = rho,k=k,n=n,alpha=alpha,variance_function = var_rho.Swiger),
                                          "Wald.F" = Wald(rho = rho,k=k,n=n,alpha=alpha,variance_function = var_rho.Fisher),
                                          "Wald.Ze" = Wald(rho = rho,k=k,n=n,alpha=alpha,variance_function = var_rho.Zerbe),
                                          "F" =  F_method(rho = rho,k=k,n=n,alpha=alpha),
                                          "Z.S" = Z_method(rho = rho,k=k,n=n,alpha=alpha,variance_function = var_rho.Swiger),
                                          "Z.F" = Z_method(rho = rho,k=k,n=n,alpha=alpha,variance_function = var_rho.Fisher),
                                          "Z.Ze" = Z_method(rho = rho,k=k,n=n,alpha=alpha,variance_function = var_rho.Zerbe),
                                          "ZF" = ZF(rho = rho,k=k,n=n,alpha=alpha))))
}

get_var_theory_distr<-function(param_dist,type="norm"){
  ## This function will help get the variance from the parameters 
  ##defined in functional call. These will be used later
  if(type == "norm"){
    variance = (param_dist[names(param_dist) == "sd"])^2
  }else if(type == "pois"){
    variance = param_dist[names(param_dist) == "lambda"]
  }else if(type == "chisq"){
    variance = 2*param_dist[names(param_dist) == "df"]
  }else if(type == "beta"){
    a<-param_dist[names(param_dist) == "shape1"]
    b<-param_dist[names(param_dist) == "shape2"]
    variance = unname(a*b/((a+b)^2 * (a+b+1)))
  }else if(type == "Beta.4P"){
    mI<-param_dist[names(param_dist) == "l"]
    mA<-param_dist[names(param_dist) == "u"]
    a <- param_dist[names(param_dist) == "alpha"]
    b <- param_dist[names(param_dist) == "beta"]
    variance = unname(a*b/((a+b)^2 * (a+b+1))*(mA-mI)^2)
  }else if(type == "lnorm"){
    sdl = ifelse("sdlog" %in% param_dist,
                 param_dist[names(param_dist) == "sdlog"],1)
    ml = ifelse("meanlog" %in% param_dist,
                param_dist[names(param_dist) == "meanlog"],0)
    variance = unname(exp(sdl^2)-1)*exp(2*ml+sdl^2)
  }else if(type == "gamma"){
    sh <- param_dist[names(param_dist) == "shape"]
    if("rate" %in% param_dist){
      variance = sh/(param_dist[names(param_dist) == "rate"]^2)
    }else if("scale" %in% param_dist){
      variance = sh*param_dist[names(param_dist) == "scale"]^2
    }else{
      variance = sh
    }
  }
  
  return(variance)
}

skew_kurt<-function(param_dist,type="norm"){
  if(type == "norm"){
    skew = 0
    kurt = 0
  }else if(type == "Beta.4P"){
    a <- param_dist[names(param_dist) == "alpha"]
    b <- param_dist[names(param_dist) == "beta"]
    skew = unname(2*(b-a)/(a+b+2) * sqrt((a+b+1)/(a*b)))
    kurt = unname(6*((a-b)^2*(a+b+1)-a*b*(a+b+2))/(a*b*(a+b+2)*(a+b+3)))
    #unname(3*(a+b+1)*(2*(a+b)^2+a*b*(a+b-6))/(a*b*(a+b+2)*(a+b+3)))
  }else if(type == "lnorm"){
    sdl = param_dist[names(param_dist) == "sdlog"]
    ml = ifelse("meanlog" %in% param_dist,
                param_dist[names(param_dist) == "meanlog"],0)
    skew = unname((exp(sdl^2)+2)*sqrt(exp(sdl^2)-1))
    kurt = unname(exp(4*sdl^2)+2*exp(3*sdl^2) + 3*exp(2*sdl^2)- 6)
  }else if(type == "gamma"){
    sh <- param_dist[names(param_dist) == "shape"]
    skew = unname(2/sqrt(sh))
    kurt = unname(6/sh)
  }  
  return(list("Skewness"=skew, "Kurtosis"=kurt))
}



generate_data <- function(param_dist, ## named vector of error sd/param matching 
                          #distribution parameters as in R functional calls
                          target_ICC, ## target ICC to generate data
                          form ="norm", ## distribution for subjects
                          const = 0, ## intercept term for the model 
                          alpha = 0.05, ##Corresponding to 
                          #(1-alpha)x100% confidence 
                          verbose = TRUE,
                          nk=NULL, ##If eg. nk=c(10,3) for n and k, optimization
                          # will be skipped
                          imbalance = 0, ## for unbalancedness
                          seed = 0
                          
){
  
  set.seed(seed)
  n=nk[1]
  k=nk[2]
  
  N <- n*k
  
  
  ## We use the theoritical value of the variance for calculating the error 
  #variance based on target ICC
  sigma_s <- sqrt(get_var_theory_distr(param_dist = param_dist, type=form))
  sigma_e <- sigma_s*sqrt((1 - target_ICC)/target_ICC)
  
  ## Subject values are drawn from defined distribution parameters 
  distr_sub<-paste0("r",form,"(n=",n,",", paste0(names(param_dist),"=", param_dist,collapse=","), ")")
  subject_measurements <- eval(parse(text=distr_sub))
  
  ## Errors are generated from a standard normal distrubution
  #distr_err <- paste0("r",form,"(n=",N,",",names(param_dist),"=", sigma_e,")")
  eij_data <- rnorm(n=N, sd=sigma_e)#eval(parse(text=distr_err))
  
  eij <- matrix(eij_data,nrow=n, ncol=k)
  
  y <- const + subject_measurements + eij
  
  if (imbalance>0){
    delete_resp = as.integer(imbalance*n/100)
    y[sample(nrow(y)*ncol(y), delete_resp)] <- NA
    n <- mean(colSums(!is.na(y)))
  }
  
  
  rho_est = ICC_estimate(y)
  bias_ANOVA <- function(rho){
    rt <- -2*(1-rho)*(rho+(1-rho)/k)*(rho+(1-rho)/(k*n))/(n-1)
    return(rt)
  }
  rho_bias = bias_ANOVA(rho=target_ICC)
  rho_ML =ICC_estimate(y,MLE = TRUE)
  Confidence_intervals <- list("Wald.S" = Wald(rho_est, k=k, n=n, alpha=alpha,
                                               variance_function = var_rho.Swiger),
                               "Wald.F" = Wald(rho_est, k=k,n=n, alpha=alpha,
                                               variance_function = var_rho.Fisher),
                               "Wald.Ze" = Wald(rho_est, k=k,n=n, alpha=alpha,
                                                variance_function = var_rho.Zerbe),
                               "F" =  F_method(rho_est, k=k,n=n, alpha=alpha),
                               "Z.S" = Z_method(rho_est, k=k,n=n, alpha=alpha,
                                                variance_function = var_rho.Swiger),
                               "Z.F" = Z_method(rho_est, k=k,n=n, alpha=alpha,
                                                variance_function = var_rho.Fisher),
                               "Z.Ze" = Z_method(rho_est, k=k,n=n, alpha=alpha,
                                                 variance_function = var_rho.Zerbe),
                               "ZF" = ZF(rho_est, k=k,n=n, alpha=alpha)
  )
  
  return(list("rho_est"=rho_est,
              "rho_est_MLE"=rho_ML,
              "rho_bias.ANOVA" = rho_bias,
              "Confidence_intervals"=Confidence_intervals,
              "Y"=y,
              "X"=subject_measurements,
              "error"=eij,
              "seed"=seed))
}

CIwidth <- function(rho, k, target_width, alpha=0.05, n_max = 1e4, param_dist=c("sd"=1), 
                    form="norm"){
  n<-4
  widths <- list()
  while (n<=n_max) {
    widths[[n-3]] =  c("Wald.S" = unname(diff(Wald(rho = rho,k=k,n=n,alpha=alpha,variance_function = var_rho.Swiger))),
                       "Wald.F" = unname(diff(Wald(rho = rho,k=k,n=n,alpha=alpha,variance_function = var_rho.Fisher))),
                       "Wald.Ze" = unname(diff(Wald(rho = rho,k=k,n=n,alpha=alpha,variance_function = var_rho.Zerbe))),
                       "F" =  unname(diff(F_method(rho = rho,k=k,n=n,alpha=alpha))),
                       "Z.S" = unname(diff(Z_method(rho = rho,k=k,n=n,alpha=alpha,variance_function = var_rho.Swiger))),
                       "Z.F" = unname(diff(Z_method(rho = rho,k=k,n=n,alpha=alpha,variance_function = var_rho.Fisher))),
                       "Z.Ze" = unname(diff(Z_method(rho = rho,k=k,n=n,alpha=alpha,variance_function = var_rho.Zerbe))),
                       "ZF" = unname(diff(ZF(rho = rho,k=k,n=n,alpha=alpha)))
    )
    if (all(widths[[n-3]]<=target_width, na.rm = TRUE)){
      break
    }else{
      n <- n+1
    }
  }
  width <- do.call(rbind.data.frame, widths)
  colnames(width) = names(widths[[n-3]])
  samplesize <- apply(width, 2, function(x) which(x<=target_width)[1]+3)
  names(samplesize)<- c('Wald.S', 'Wald.F', 'Wald.Ze', 'F',
                        'Z.S', 'Z.F', 'Z.Ze', 'ZF')
  return(samplesize)
}



observed_power<-function(samplesize, rho.0, rho.A, k, alpha=0.05, nsims= 2.5e4){
  ci_met <-  c('Wald.S', 'Wald.F', 'Wald.Ze', 'F',
               'Z.S', 'Z.F', 'Z.Ze', 'ZF')
  
  bin_code <- function(x, rho.0){
    return(unname(ifelse(x>rho.0,1,0)))
  }
  Gen <- sapply(1:8, function(j)
    mean(sapply(1:nsims, function(sed)
      bin_code(generate_data(param_dist = c("sd"=1), 
                             form="norm", 
                             alpha=2*alpha,
                             nk = c(samplesize[j],k), 
                             seed=sed*10, 
                             target_ICC=rho.A)$Confidence_intervals[[ci_met[j]]][1],
               rho.0 = rho.0))))
  names(Gen)<- c('Wald.S', 'Wald.F', 'Wald.Ze', 'F',
                 'Z.S', 'Z.F', 'Z.Ze', 'ZF')
  return(Gen)
}

observed_coverage<- function(samplesize, rho, width, k, alpha=0.05, nsims= 2.5e4){
  ci_met <-  c('Wald.S', 'Wald.F', 'Wald.Ze', 'F',
               'Z.S', 'Z.F', 'Z.Ze', 'ZF')
  
  bin_code <- function(ci, rho){
    return(unname(ifelse(rho>=ci[1]&rho<=ci[2],1,0)))
  }
  Gen <- sapply(1:8, function(j)
    mean(sapply(1:nsims, function(sed)
      bin_code(generate_data(param_dist = c("sd"=1), 
                             form="norm", 
                             alpha=alpha,
                             nk = c(samplesize[j],k), 
                             seed=sed*10, 
                             target_ICC=rho)$Confidence_intervals[[ci_met[j]]],
               rho = rho))))
  names(Gen)<- c('Wald.S', 'Wald.F', 'Wald.Ze', 'F',
                 'Z.S', 'Z.F', 'Z.Ze', 'ZF')
  return(Gen)
}


observed_assurance<- function(samplesize, rho, width, k, alpha=0.05, nsims= 2.5e4){
  ci_met <-  c('Wald.S', 'Wald.F', 'Wald.Ze', 'F',
               'Z.S', 'Z.F', 'Z.Ze', 'ZF')
  
  bin_code <- function(ci, width){
    return(unname(ifelse((ci[2]-ci[1])<=width,1,0)))
  }
  Gen <- sapply(1:8, function(j)
    mean(sapply(1:nsims, function(sed)
      bin_code(generate_data(param_dist = c("sd"=1), 
                             form="norm", 
                             alpha=alpha,
                             nk = c(samplesize[j],k), 
                             seed=sed*10, 
                             target_ICC=rho)$Confidence_intervals[[ci_met[j]]],
               width = width))))
  names(Gen)<- c('Wald.S', 'Wald.F', 'Wald.Ze', 'F',
                 'Z.S', 'Z.F', 'Z.Ze', 'ZF')
  return(Gen)
}
choices <- c('Wald.S', 'Wald.F', 'Wald.Ze', 'F',
             'Z.S', 'Z.F', 'Z.Ze', 'ZF')
names(choices) <- c("Wald method with Swiger variance",
                    "Wald method with Fisher variance",
                    "Wald method with Zerbe variance",
                    "Searle (F) method",
                    "Normalization with Swiger variance",
                    "Normalization with Fisher variance",
                    "Normalization with Zerbe variance",
                    "Normalization of Searle method")
color_plate<- c('blue', 'violet', 'red', 'black', 
                      'lightblue', 'pink', 'orange', 'grey')
                      names(color_plate)<-choices
                      
                      shinyApp(
                        ui <- navbarPage(h4(em("Sample Size App")),
                                         tabPanel(strong("Choice of method"),
                                                  fluidRow(column(4, offset = 3,
                                                                  h3("Sample Size Approach: "),
                                                                  radioButtons("SSizeProc", label=" ",
                                                                               choices = c("Width of Confidence Interval Approach",
                                                                                           "Assurance Probability Approach",
                                                                                           "Testing Approach"),
                                                                               selected = "Assurance Probability Approach"),
                                                                  checkboxInput("obs", label="Output observed coverage/assurance/power (simulation based)", value=FALSE),
                                                                  checkboxInput("cst", label="Find the optimal sample size minimizing a linear cost function", value=TRUE),
                                                  )),
                                                  fluidRow(column(5, offset = 2,
                                                                  h3("Description: "),
                                                                  mainPanel(conditionalPanel(condition = "input.SSizeProc == 'Width of Confidence Interval Approach'",
                                                                                             p("The", strong("Width of the Confidence Interval Approach"), "aims to find the minimum number of participants", em("n"),
                                                                                               "given the width of the confidence interval around ICC"),
                                                                                             conditionalPanel(condition="input.obs==1",
                                                                                                              p("Observed coverage (based on simulations) will be displayed."),
                                                                                                              p("Simulations will be based on generating standard normal replicates for the specified", em("k"),
                                                                                                                "and obtained sample sizes (i.e.,", em("n"),") for each confidence interval method.")
                                                                                             )
                                                                  ),
                                                                  conditionalPanel(condition = "input.SSizeProc == 'Assurance Probability Approach'",
                                                                                   p("The", strong("Assurance Probability Approach"), "aims to find the minimum number of participants", em("n"),
                                                                                     "such that the width of the confidence interval around ICC is less than the specified", em("target width"),
                                                                                     "with a certain level of", em("assurance")
                                                                                   ),
                                                                                   conditionalPanel(condition="input.obs==1",
                                                                                                    p("Observed assurance probability (based on simulations) will be displayed."),
                                                                                                    p("Simulations will be based on generating standard normal replicates for the specified", em("k"),
                                                                                                      "and obtained sample sizes (i.e.,", em("n"),") for each confidence interval method.")
                                                                                   ),
                                                                                   conditionalPanel(condition="input.cst==1",
                                                                                                    p("The optimal number of participants, (",em("n"),"),and rater/repititions(",em("k"),
                                                                                                      ") will be provided based on minimizing a linear cost function. Given the total cost comprising of the",
                                                                                                      "cost of hiring a rater or making a repitition,","c",tags$sub('1'),"cost of hiring a participant", 
                                                                                                      "c",tags$sub('2'), "and cost of making one observation", "c",tags$sub('3'),"the optimal choice of",
                                                                                                      "(n.k) to satify the given criterion for the approach will be provided."
                                                                                                    )
                                                                                   )
                                                                  ),
                                                                  conditionalPanel(condition = "input.SSizeProc == 'Testing Approach'",
                                                                                   p("The", strong("Testing Approach"), "aims to find the minimum number of participants", em("n"),
                                                                                     "such that the hypothesis test that the ICC is less than or equal to the ", em("ICC value at null hypothesis"), 
                                                                                     "(null hypothesis) against the alternative that the ICC is greater than the ICC value at null hypothesis, achieves a certain", em("power"),
                                                                                     " The value of ICC under alternative hypothesis is", em("ICC value at alternative hypothesis"),"."),
                                                                                   conditionalPanel(condition="input.obs==1",
                                                                                                    p("Observed power for the hypothesis test (based on simulations) will be displayed."),
                                                                                                    p("Simulations will be based on generating standard normal replicates for the specified", em("k"),
                                                                                                      "and obtained sample sizes (i.e.,", em("n"),") for each confidence interval method.")
                                                                                   )),
                                                                  )))
                                         ),
                                         tabPanel(strong("Parameters"),
                                                  fixedPage(
                                                    sidebarPanel(width = 4,
                                                                 h4("General Parameters"),
                                                                 numericInput("k", label = "Number of Raters:",
                                                                              min = 2, max=100, step= 1, value = 10),
                                                                 sliderInput("alpha", label="Confidence level", 
                                                                             min = 0.8, max = 0.99, step = 0.01, value=0.95),
                                                                 conditionalPanel(condition="input.obs==1",
                                                                                  numericInput("nsims_obs", 
                                                                                               label = "Number of simulations:",
                                                                                               min=1000, max = 2.5e5, step=500, value= 2000)
                                                                 ),
                                                                 conditionalPanel(condition="input.cst==1",
                                                                                  h4("Total cost (T) = c1xk + c2xn + c3xnk"),
                                                                                  list(numericInput("c1",
                                                                                                    label="Cost of recruiting a rater/making a repitition(c1)",
                                                                                                    min=0, max=2.5e5, step=0.05, value = 0),
                                                                                       numericInput("c2",
                                                                                                    label="Cost of recruiting a participant(c2)",
                                                                                                    min=0, max=2.5e5, step=0.05, value = 0),
                                                                                       numericInput("c3",
                                                                                                    label="Cost of obtaining one observation(c3)",
                                                                                                    min=0, max=2.5e5, step=0.05, value = 0))
                                                                 )),
                                                    sidebarPanel(width=4,offset=5, h4("Sample Size Approach Specific Parameters"),
                                                                 conditionalPanel(
                                                                   condition = "input.SSizeProc == 'Width of Confidence Interval Approach'",
                                                                   list(sliderInput("targetwidth", label="Target Width", min = 0.05, max = 1, step = 0.05, value=0.2),
                                                                        sliderInput("targetICC", label="Value of ICC", min = 0.0, max = 0.999, step = 0.025, value=0.8),
                                                                        numericInput("nsims", label = "Maximum value of n:",min=1000, max = 2.5e5, step=500, value= 25000)
                                                                   )
                                                                 ),
                                                                 conditionalPanel(
                                                                   condition = "input.SSizeProc == 'Assurance Probability Approach'",
                                                                   list(sliderInput("assurance", label="Assurance Probability", min = 0.5, max = 1, step = 0.05, value=0.5),
                                                                        sliderInput("targetwidthA", label="Target Width", min = 0.05, max = 1.0, step = 0.025, value=0.2),
                                                                        sliderInput("targetICCA", label="Value of ICC", min = 0.0, max = 1.0, step = 0.025, value=0.9),
                                                                        numericInput("nmaxA", label = "Maximum value of n:",min=1000, max = 2.5e5, step=1000, value= 5000)
                                                                   )
                                                                 ),
                                                                 conditionalPanel(
                                                                   condition = "input.SSizeProc == 'Testing Approach'",
                                                                   list(sliderInput("power", label="Power of Hyothesis Test", min = 0.5, max = 1, step = 0.05, value=0.9),
                                                                        sliderInput("ICCNull", label="ICC value at null hypothesis", min = 0.0, max = 1.0, step = 0.025, value=0.8),
                                                                        sliderInput("ICCAlt", label="ICC value at alternative hypothesis", min = 0.0, max = 1.0, step = 0.025, value=0.9),
                                                                        numericInput("nmaxT", label = "Maximum value of n:",min=1000, max = 2.5e5, step=1000, value= 5000)
                                                                   )
                                                                 ),
                                                    ))),
                                         
                                         tabPanel(strong("Results"),
                                                  fluidRow(
                                                    column(4, offset = 3,
                                                           h3("Number of Participants Required:"),
                                                           tableOutput('tableSS'))),
                                                  conditionalPanel(
                                                    condition = "input.obs==1",
                                                    fluidRow(
                                                      column(4, offset = 3,
                                                             conditionalPanel(condition="input.SSizeProc == 'Width of Confidence Interval Approach'",
                                                                              h3("Observed Coverage:")),
                                                             conditionalPanel(condition="input.SSizeProc == 'Assurance Probability Approach'",
                                                                              h3("Observed Assurance Probability:")),
                                                             conditionalPanel(condition="input.SSizeProc == 'Testing Approach'",
                                                                              h3("Observed Power:")),
                                                             tableOutput('tableobs')
                                                      ))),
                                                  fluidRow(
                                                    column(4, offset = 3,
                                                           conditionalPanel(condition="input.cst == 1",
                                                                            h3("Optimal Combination by minimizing total cost:"),
                                                                            tableOutput('tableCost')))),
                                                  fluidRow(
                                                    sidebarPanel(width = 3,
                                                                 h3("Choose parameters for plotting:"),
                                                                 checkboxGroupInput("select", "Select methods", 
                                                                                    choices=choices)),
                                                    mainPanel(width = 4,
                                                              h3("Contour Plot:"),
                                                              plotOutput('plt')))
                                         ),
                                         navbarMenu(strong("Explanation"),
                                                    tabPanel("Confidence Interval Methods",
                                                             fluidRow(
                                                               column(6, offset = 2,
                                                                      h3("Variance of ICC"),
                                                                      p("Three different approximations of the variance of ICC exists:",
                                                                        tags$ul(tags$li("Swiger variance:",withMathJax("$$var(\\hat{\\rho})_S=\\frac{2(nk-1)(1-\\rho)^2\\{1+(k - 1)\\rho\\}^2}{k^2n(k-1)(n-1)}$$")),
                                                                                tags$li("Fisher variance:",withMathJax("$$var(\\hat{\\rho})_F = \\frac{2(1-\\rho)^2\\{1+(k - 1)\\rho\\}^2}{nk(k-1)}$$")),
                                                                                tags$li("Zerbe variance:",withMathJax("$$var(\\hat{\\rho})_{Ze} =\\frac{2(1-\\rho)^2\\{1+(k - 1)\\rho\\}^2n^2(k-1)^2(nk-3)}{k^2(n-1)(nk-n-2)^2(nk-n-4)}$$"))
                                                                        )
                                                                      ),
                                                                      h3("Confidence Interval of ICC"),
                                                                      p("In the literature there are currently four methods to compute the lower(L) and upper(U) bounds of the confidence interval of","\u03c1",  
                                                                        "namely Wald method, Searle method, and their normalized version, conforming eight different confidence interval methods :",
                                                                        tags$ul(tags$li("Wald",tags$sub("S"), "  : Described as ", strong("Wald confidence interval method"), "with", strong("Swiger variance")),
                                                                                tags$li("Wald",tags$sub("F"), "  : Described as ", strong("Wald confidence interval method"), "with", strong("Fisher variance")),
                                                                                tags$li("Wald",tags$sub("Ze"), "  : Described as ", strong("Wald confidence interval method"), "with", strong("Zerbe variance")),
                                                                                tags$li("F", "    : Described as ", strong("Searle method"), "or", em("Exact"), "method"),
                                                                                tags$li("Z",tags$sub("S"), "   :  Described as ", strong("normalized ICC method"), "with", strong("Swiger variance")),
                                                                                tags$li("Z",tags$sub("F"), "   :  Described as ", strong("normalized ICC method"), "with", strong("Fisher variance")),
                                                                                tags$li("Z",tags$sub("Ze"), "   :  Described as ", strong("normalized ICC method"), "with", strong("Zerbe variance")),
                                                                                tags$li("ZF",tags$sub("\u03c1"), "   : Described as ", strong("normalized Searle method")),
                                                                        )),
                                                                      h4("Wald Confidence Interval method:"),
                                                                      p("Based on the central limit theorem, the lower (L) and upper (U) bounds of the confidence interval of the ICC can be written as,",
                                                                        withMathJax("$$L,U = \\hat{\\rho}\\mp z_{1-\\alpha/2}\\sqrt{var(\\hat{\\rho})}$$"),
                                                                        "where z",tags$sub("1-\u03b1/2"), "is the (1-\u03b1/2)100 percentile of the standard normal distribution."),
                                                                      h4("Searle method:"),
                                                                      p("Under the assumption of normality of outcomes of the ANOVA model, the F-statistic follows a distribution ",
                                                                        "\u03c4 X F",tags$sub("\u03bd1, \u03bd2"), "where \u03c4=(1+(k-1)\u03c1)/(1-\u03c1), and F",tags$sub("\u03bd1, \u03bd2"),
                                                                        " is an F distribution with \u03bd1 = n-1, \u03bd2=n(k-1). The ratio of the mean squares of the ANOVA is represented as,",
                                                                        withMathJax("$$F(\\hat{\\rho}) = \\frac{BMS}{WMS} = \\frac{1+(k-1)\\hat{\\rho}}{1-\\hat{\\rho}}$$"), 
                                                                        "Then the lower (L) and upper (U) bounds of the confidence interval of the ICC can be written as, ",
                                                                        withMathJax("$$L,U =  \\frac{F(\\hat{\\rho})/F_u-1}{F(\\hat{\\rho})/F_u+k-1},\\frac{F(\\hat{\\rho})/F_L-1}{F(\\hat{\\rho})/F_L+k-1}$$")
                                                                      ),
                                                                      h4("Normalized ICC method:"),
                                                                      p("Fisher transformation can be applied to the ICC such that the transformed ICC approximately follows a normal distribution.Applying this transformation leads to",
                                                                        withMathJax("$$Z(\\hat{\\rho}) = \\frac{1}{2}ln\\frac{1+(k-1)\\hat{\\rho}}{1-\\hat{\\rho}} \\sim \\mathcal{N}\\bigg\\{\\frac{1}{2}ln\\frac{1+\\rho}{1-{\\rho}}, var(Z(\\hat{\\rho}))\\bigg\\}$$"), 
                                                                        "The variance can be obtained using the Delta method to one of the variances defined above, as, ",
                                                                        withMathJax("$$var(Z(\\hat{\\rho}))=\\frac{var(\\hat{\\rho})}{(1-\\rho^2)^2}$$"),
                                                                        "We compute the confidence interval around normalized \u03c1 as,",
                                                                        withMathJax("$$L_{Z},U_{Z} = Z(\\hat{\\rho}) \\mp z_{1-\\alpha/2}{\\sqrt{var(Z(\\hat{\\rho}))}}$$"),
                                                                        "Finally, he confidence interval of \u03c1 is obtained by back-transformation,",
                                                                        withMathJax("$$ L,U = \\frac{exp(2L_Z)-1}{exp(2L_Z)+1}, \\frac{exp(2U_Z)-1}{exp(2U_Z)+1}$$")
                                                                      ),
                                                                      h4("Normalized Searle method:"),
                                                                      p("The F-statistic can be normalized as,",
                                                                        withMathJax("$$Z(F(\\hat{\\rho})) = \\frac{1}{2}ln\\frac{1+(k-1)\\hat{\\rho}}{1-\\hat{\\rho}} \\sim \\mathcal{N}\\bigg\\{\\frac{1}{2}ln\\frac{1+(k-1)\\rho}{1-\\rho}, \\frac{1}{2}\\bigg(\\frac{1}{n-1}+\\frac{1}{n(k-1)}\\bigg)\\bigg\\}$$"),
                                                                        "The confidence interval on this Fisher transformed scale is",
                                                                        withMathJax("$$L_{ZF},U_{ZF} = Z(F(\\hat{\\rho})) \\mp z_{1-\\alpha/2}{\\sqrt{var(Z(F(\\hat{\\rho})))}}$$"),
                                                                        "The confidence limits for \u03c1 can be obtained directly by back-transforming as,",
                                                                        withMathJax("$$L,U= \\frac{exp(2L_{ZF})-1}{exp(2L_{ZF})+k-1}, \\frac{exp(2U_{ZF})-1}{exp(2U_{ZF})+k-1}$$")
                                                                      )
                                                               )))#,
                                                    #    tabPanel("Width of Confidence Interval Approach",
                                                    #             fluidRow(
                                                    #               column(4, offset = 2,
                                                    #             p("The", strong("Width of Confidence Interval Approach"), "consists in finding the minimum number of participants", em("n"),
                                                    #             "for a given width, \u03c9, of the confidence interval around a planned value of \u03c1 and a given number of raters,"
                                                    #             ,em("k"),".")
                                                    #             )))#,
                                                    #  tabPanel("Assurance Probability Approach"),
                                                    #  tabPanel("Testing Approach")
                                         )
                        ),
                        server <- function(input, output){
                          dataSSize <- reactive({
                            alpha <- 1-as.numeric(input$alpha)
                            rho <- NULL
                            txt <- NULL
                            
                            df <- NULL
                            obs<-NULL
                            c_1<-as.numeric(input$c1)
                            c_2<-as.numeric(input$c2)
                            c_3<-as.numeric(input$c3)
                            
                            observeEvent(input$ICCAlt,{
                              updateSliderInput(inputId = "ICCNull", 
                                                min =0,
                                                max=as.numeric(input$ICCAlt))
                            })
                            observeEvent(input$ICCNull,{
                              updateSliderInput(inputId = "ICCAlt", 
                                                min =as.numeric(input$ICCNull),
                                                max=1)
                            })
                            
                            
                            if(input$SSizeProc =='Assurance Probability Approach'){
                              withProgress(message = 'Computing', style = 'notification', value = 0,{
                                df<-t(sapply(2:as.numeric(input$k), function(x)
                                  SampleSize(criterion_specification = c('rho' = as.numeric(input$targetICCA),
                                                                         'target_width'=as.numeric(input$targetwidthA),
                                                                         'gamma_value'= as.numeric(1-input$assurance)),
                                             k = x, 
                                             alpha = alpha,
                                             nmax=as.numeric(input$nmaxA))))
                                rho = as.numeric(input$targetICCA)
                                
                                if(input$obs){
                                  incProgress(0.5)
                                  obs <- observed_assurance(samplesize = df[nrow(df),], 
                                                            rho=rho,
                                                            width=as.numeric(input$targetwidthA),
                                                            k = as.numeric(input$k),
                                                            alpha=alpha,
                                                            nsims=as.numeric(input$nsims_obs))
                                  txt = "Assurance Probability"
                                }
                              })
                            }else if(input$SSizeProc =='Testing Approach'){
                              withProgress(message = 'Computing', style = 'notification', value = 0,{
                                df<-t(sapply(2:as.numeric(input$k), function(x)
                                  SampleSize(criterion_specification = c('rho.0' = as.numeric(input$ICCNull),
                                                                         'rho.A' = as.numeric(input$ICCAlt),
                                                                         'beta_value'= as.numeric(1-input$power)),
                                             k = x, 
                                             alpha = alpha,
                                             nmax=as.numeric(input$nmaxT))))
                                rho = as.numeric(input$ICCAlt)
                                if(input$obs){
                                  incProgress(0.5)
                                  obs <- observed_power(samplesize = df[nrow(df),], 
                                                        rho.0=as.numeric(input$ICCNull),
                                                        rho.A=as.numeric(input$ICCAlt),
                                                        k = as.numeric(input$k),
                                                        alpha=alpha,
                                                        nsims=as.numeric(input$nsims_obs))
                                  txt = "Power"
                                }
                              })
                            }else if(input$SSizeProc =='Width of Confidence Interval Approach'){
                              withProgress(message = 'Computing', style = 'notification', value = 0,{
                                df<-t(sapply(2:as.numeric(input$k), function(x)
                                  CIwidth(rho = as.numeric(input$targetICC), 
                                          k = x, 
                                          target_width = as.numeric(input$targetwidth),
                                          alpha = alpha,
                                          param_dist = param_dist,
                                          form = form)))
                                rho = as.numeric(input$targetICC)
                                if(input$obs){
                                  incProgress(0.5)
                                  obs <- observed_coverage(samplesize = df[nrow(df),], 
                                                           rho=rho,
                                                           width=as.numeric(input$targetwidthA),
                                                           k = as.numeric(input$k),
                                                           alpha=alpha,
                                                           nsims=as.numeric(input$nsims_obs))
                                  txt = "Coverage"
                                }
                              })
                            }
                            Costs=t(sapply(2:as.numeric(input$k), 
                                           function(i) 
                                             sapply(1:ncol(df), 
                                                    function(j) 
                                                      c_1*i+c_2*df[i-1,j]+c_3*i*df[i-1,j])))
                            Optimal.Comb<- sapply(1:ncol(Costs), function(x)
                              paste0("(",
                                     'k'= c(2:as.numeric(input$k))[which.min(Costs[,x])], ",", 
                                     'n'=df[,x][which.min(Costs[,x])], ")")
                            )
                            list('SampleSize'=df, 
                                 'Cost'= Optimal.Comb,
                                 'Obs'=obs, 'txt'=txt)
                            
                            
                          })
                          
                          output$tableSS <-renderTable({
                            t(data.frame("No."=as.integer(dataSSize()$SampleSize[(as.numeric(input$k)-1),]),
                                         row.names = c(HTML(paste0("Wald",tags$sub("S"))),
                                                       HTML(paste0("Wald",tags$sub("F"))),
                                                       HTML(paste0("Wald",tags$sub("Ze"))),
                                                       "F",
                                                       HTML(paste0("Z",tags$sub("S"))),
                                                       HTML(paste0("Z",tags$sub("F"))),
                                                       HTML(paste0("Z",tags$sub("Ze"))),
                                                       HTML(paste0("ZF",tags$sub(HTML("&rho;")))))))
                          },
                          sanitize.text.function = function(x) x,
                          rownames = TRUE)
                          
                          output$tableCost <-renderTable({
                            t(data.frame("Combination(k,n)"=dataSSize()$Cost,
                                         row.names = c(HTML(paste0("Wald",tags$sub("S"))),
                                                       HTML(paste0("Wald",tags$sub("F"))),
                                                       HTML(paste0("Wald",tags$sub("Ze"))),
                                                       "F",
                                                       HTML(paste0("Z",tags$sub("S"))),
                                                       HTML(paste0("Z",tags$sub("F"))),
                                                       HTML(paste0("Z",tags$sub("Ze"))),
                                                       HTML(paste0("ZF",tags$sub(HTML("&rho;")))))))
                          },
                          sanitize.text.function = function(x) x,
                          rownames = TRUE)
                          
                          output$tableobs <-renderTable({
                            t(data.frame("Observed"=format(round(dataSSize()$Obs, 3), nsmall = 3),
                                         row.names = c(HTML(paste0("Wald",tags$sub("S"))),
                                                       HTML(paste0("Wald",tags$sub("F"))),
                                                       HTML(paste0("Wald",tags$sub("Ze"))),
                                                       "F",
                                                       HTML(paste0("Z",tags$sub("S"))),
                                                       HTML(paste0("Z",tags$sub("F"))),
                                                       HTML(paste0("Z",tags$sub("Ze"))),
                                                       HTML(paste0("ZF",tags$sub(HTML("&rho;")))))))
                          },sanitize.text.function = function(x) x,
                          rownames = TRUE)
                          
                          output$plt <- renderPlot({
                            data <- dataSSize()$SampleSize[, colnames(dataSSize()$SampleSize) 
                                                           %in% input$select]
                            color<- color_plate[names(color_plate) 
                                                %in% input$select]
                            
                            matplot(x=2:as.numeric(input$k), col = color,
                                    y=data, type='b', pch=19, lty=1,
                                    xlab = "Number of Raters/Repititions (k)",
                                    ylab="Number of Participants (n)", cex=1.25)
                            if (length(input$select)>1){
                              legend('topright', c(expression('Wald'['S']),
                                                   expression('Wald'['F']),
                                                   expression('Wald'['Ze']),
                                                   "F",
                                                   expression('Z'['S']),
                                                   expression('Z'['F']),
                                                   expression('Z'['Ze']),
                                                   expression('ZF'[rho]))[
                                                     which(choices %in% input$select)],
                                     pch=19, cex=1.5, col=color,
                                     ncol=ifelse(length(input$select)>5,2,1))}
                          })
                        },
                        options = list(port = 6555,launch.browser = TRUE)
                        
                        
                      )
