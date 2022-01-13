#--------------------------------------------------------------------
#  Computing responsiveness measures using GLM models for the mediator 
#---------------------------------------------------------------------
M.mod <- function(mediator="gaussian", link="identity"){
  if(mediator=="gaussian"){
    mod.glm <- glm(HBV ~ log.HCV + gender + agegp2 + agegp3 + agegp4 + smoke + alcohol, 
                   family = gaussian, data=dados)
    Cov <- summary(mod.glm)$cov.scaled
    scl.parm <- mod.glm$deviance/mod.glm$df.residual
    scl.parm.std = NA
  } else if(mediator == "gamma"){
    mod.glm <- glm(HBV ~ log.HCV +  gender + agegp2 + agegp3 + agegp4 + smoke + alcohol, 
                   family = Gamma(link=link), data=dados)
    Cov <- summary(mod.glm)$cov.scaled
    scl.parm <- as.numeric(gamma.shape(mod.glm, verbose=FALSE)[1])
    scl.parm.std <- as.numeric(gamma.shape(mod.glm, verbose=FALSE)[2])
  }
  return(list(zeta.0=mod.glm$coeff[1], zeta.1=mod.glm$coeff[2], zeta.2=mod.glm$coeff[-c(1,2)], 
              Cov=Cov, scl.parm=scl.parm, scl.parm.std=scl.parm.std))
}


#------------------------------------------------------------
#  Setting mediator models
#  ----------------------------------------------------------

gamma.theta <- function(link="identity", a, x, zeta.0, zeta.1, zeta.2){
  if(link == "identity"){
    return(-1/(zeta.0 + zeta.1 * a + t(zeta.2) %*% x))
  } else if(link == "log"){
    return(-exp(-(zeta.0 + zeta.1 * a + t(zeta.2) %*% x)))
  } else if(link == "inverse"){
    return(-(zeta.0 + zeta.1 * a + t(zeta.2) %*% x))
  }
}

d.gamma.theta <- function(link="identity", a, x, zeta.0, zeta.1, zeta.2){
  if(link == "identity"){
    return(zeta.1 * gamma.theta(link="identity", a, x, zeta.0, zeta.1, zeta.2)^2)
  } else if(link == "log"){
    return(-zeta.1 * gamma.theta(link="log", a, x, zeta.0, zeta.1, zeta.2))
  } else if(link == "inverse"){
    return(-zeta.1)
  }
}

#------------------------------------------------------------
#  Fitting survival models (Cox and AFT) for the outcome 
#  ----------------------------------------------------------

cox.aft.sen <- function(method="cox", mediator="gaussian", link="identity", a=0, x=rep(0,p)){
  if(method == "cox"){
    mod <- coxph(Surv(hcc.time/365.25, status) ~ log.HCV + HBV + gender + agegp2 + 
                   agegp3 + agegp4 + smoke + alcohol, data = dados)
    beta.1 <- as.numeric(mod$coefficients[1])
    beta.2 <- as.numeric(mod$coefficients[2])
    beta.4 <- as.numeric(mod$coefficients[-c(1,2)])
    T.cov <- mod$var[1:2,1:2]
  } else {
    mod <- survreg(Surv(hcc.time/365.25, status) ~ log.HCV + HBV + gender + agegp2 + 
                     agegp3 + agegp4 + smoke + alcohol, dist='weibull', data = dados)
    beta.0 <- as.numeric(mod$coefficients[1])
    beta.1 <- as.numeric(mod$coefficients[2])
    beta.2 <- as.numeric(mod$coefficients[3])
    beta.4 <- as.numeric(mod$coefficients[-c(1,2,3)])
    T.cov <- mod$var[2:3,2:3]
    zeta.scl <- mod$scale
    eps.mgf <- gamma(1+zeta.scl)
  }
  
  nds <- nis <- rep(NA, length(a))
  ## NDE without interaction:
  if(mediator == "gaussian"){  
    M <- M.mod()
    zeta.0 <- M$zeta.0
    zeta.1 <- M$zeta.1
    zeta.2 <- M$zeta.2
    sigma2 <- M$scl.parm
    cont <- 1
    for(a0 in a){
      nds[cont] <- eps.mgf * beta.1 * 
        exp(beta.2*zeta.0 + beta.0 + 0.5*sigma2*beta.2^2 + (beta.2*zeta.1+beta.1)*a[cont] + t(beta.2*zeta.2 + beta.4)%*%x)
      cont <- cont + 1
    }
    nis <- nds * zeta.1 * beta.2 / beta.1
  } else if(mediator == "gamma"){
    M <- M.mod(mediator, link)
    zeta.0 <- M$zeta.0
    zeta.1 <- M$zeta.1
    cat("zeta1 = ", zeta.1, "\n")
    zeta.2 <- M$zeta.2
    nu <- M$scl.parm
    cont <- 1
    for(a0 in a){
      theta <- gamma.theta(link, a[cont], x, zeta.0, zeta.1, zeta.2)
      d.theta <- d.gamma.theta(link, a[cont], x, zeta.0, zeta.1, zeta.2)
      nds[cont] <- eps.mgf * beta.1 * 
        exp(beta.0 + beta.1*a[cont] + t(beta.4)%*%x + nu*log(theta/(theta+beta.2/nu)))
      nis[cont] <- nds[cont] * nu * d.theta * ((1/theta) - (1/(theta+beta.2/nu))) / beta.1
      cont <- cont + 1
    }
  }
  ts <- nds + nis
  
  return(list(NIS=as.numeric(nis), NDS=as.numeric(nds), TS=as.numeric(ts)))
}

