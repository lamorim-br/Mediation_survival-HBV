#--------------------------------------------------------------------
#  Computing responsiveness measures using GLM models for the mediator
#  and using the Aalen model for the outcome
#---------------------------------------------------------------------

#------------------------------------------------------------
#  Setting mediator models
#----------------------------------------------------------
M.mod <- function(mediator="gaussian", link="identity"){
  if(mediator=="gaussian"){
    mod.glm <- glm(HBV ~ log.HCV + gender + agegp2 + agegp3 + agegp4 + smoke + alcohol                   , 
                   family = gaussian, data=dados)
    Cov <- summary(mod.glm)$cov.scaled
    scl.parm <- mod.glm$deviance/mod.glm$df.residual
    scl.parm.std = NA
  } else if(mediator == "gamma"){
    mod.glm <- glm(HBV ~ log.HCV + gender + agegp2 + agegp3 + agegp4 + smoke + alcohol, 
                   family = Gamma(link=link), data=dados)
    Cov <- summary(mod.glm)$cov.scaled
    scl.parm <- as.numeric(gamma.shape(mod.glm, verbose=FALSE)[1])
    scl.parm.std <- as.numeric(gamma.shape(mod.glm, verbose=FALSE)[2])
  }
  return(list(zeta.0=mod.glm$coeff[1], zeta.1=mod.glm$coeff[2], zeta.2=mod.glm$coeff[-c(1,2)], 
              Cov=Cov, scl.parm=scl.parm, scl.parm.std=scl.parm.std))
}

#------------------------------------------------------------
#  Auxiliary functions for the mediator model
#----------------------------------------------------------
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
#  Auxiliary functions for the Aalen model
#----------------------------------------------------------
gamma.aalen.Q <- function(tm=1, link="identity", a, zeta.0, zeta.1, zeta.2, 
                          beta.2, beta.4, nu, approximation=TRUE){
  N <- nrow(X)
  B2 <- beta.2 * tm
  B4 <- beta.4 * tm
  Q <- rep(NA, N)
  dnm <- 0
  for(i in 1:N){
    x.aux <- as.numeric(X[i,1:p])
    theta.aux <- gamma.theta(link, a, x.aux, zeta.0, zeta.1, zeta.2)
    if(approximation){
      Q[i] <- X$Prob[i]
    } else {
      Q[i] <- (((theta.aux/(theta.aux + B2/nu))^nu) * exp(t(B4) %*% x.aux)) * X$Prob[i]
    }
    dnm <- dnm + Q[i]
  }
  return(Q / dnm)
}

gamma.aalen.NIS.psi1 <- function(tm=1, x, link="identity", a=1, 
                                 zeta.0, zeta.1, zeta.2, beta.2, 
                                 nu=NULL, Q=NULL){
  N <- nrow(X)
  E <- 0
  
  for(i in 1:N){
    x.aux <- as.numeric(X[i,1:p])
    theta <- gamma.theta(link, a, x.aux, zeta.0, zeta.1, zeta.2)
    d.theta <- d.gamma.theta(link, a, x.aux, zeta.0, zeta.1, zeta.2)
    E <- E + d.theta*((1/theta)-(1/(theta+beta.2*tm/nu))) * Q[i]
  }
  
  theta <- gamma.theta(link, a, x, zeta.0, zeta.1, zeta.2)
  d.theta <- d.gamma.theta(link, a, x, zeta.0, zeta.1, zeta.2)
  
  return(d.theta*((1/theta)-(1/(theta+beta.2*tm/nu)))-E)
}

gamma.aalen.NIS <- function(tm=1, link="identity", a=1, 
                            zeta.0, zeta.1, zeta.2, 
                            beta.2, beta.4=NULL, 
                            nu=NULL, Q=NULL, approximation=TRUE){
  N <- nrow(X)
  nis <- 0
  for(i in 1:N){
    x.aux <- as.numeric(X[i,1:p])
    theta <- gamma.theta(link, a, x.aux, zeta.0, zeta.1, zeta.2)
    d.theta <- d.gamma.theta(link, a, x.aux, zeta.0, zeta.1, zeta.2)
    if(approximation){
      nis <- nis + (beta.2 * d.theta / theta^2) * X$Prob[i]
    } else {
      psi <- gamma.aalen.NIS.psi1(tm, x.aux, link, a, 
                                  zeta.0, zeta.1, zeta.2, beta.2, nu, Q)
      nis <- nis + (t(beta.4) %*% x.aux - beta.2 / (theta + beta.2 * tm / nu)) * psi * Q[i]
      nis <- nis + beta.2 * d.theta * Q[i] / (theta + beta.2 * tm / nu)^2
    }
  }
  nis <- nu * nis
  return(nis)
}

#------------------------------------------------------------
#  Setting the Aalen model
#----------------------------------------------------------

aalen.sen <- function(mediator="gaussian", link="identity", a=1, x=rep(0,p), 
                      grid.lo=0, grid.up=15, by=0.05, approximation=TRUE){
  mod <- aalen(Surv(hcc.time/365.25, status) ~ const(log.HCV) + const(HBV) + const(gender) + 
                 const(agegp2) + const(agegp3) + const(agegp4) + const(smoke) + 
                 const(alcohol), data = dados)
  
  beta.1 <- as.numeric(mod$gamma[1])
  beta.2 <- as.numeric(mod$gamma[2])
  beta.4 <- as.numeric(mod$gamma[-c(1,2)])
  
  tm.grid <- seq(from=grid.lo, to=grid.up, by=by)
  tm.length <- length(tm.grid)
  
  if(approximation){
    nis <- matrix(NA, ncol=1, nrow=length(a))
  } else {
    nis <- matrix(NA, ncol=tm.length, nrow=length(a))
  }
  
  nds <- beta.1
  
  if(mediator == "gaussian"){
    M <- M.mod()
    nis <- beta.2 * M$zeta.1
  } else if(mediator == "gamma"){
    M <- M.mod(mediator, link)
    nu <- M$scl.parm
    scl.parm.std <- M$scl.parm.std
    cont <- 1
    for(a0 in a){
      if(approximation){
        nis[cont,1] <- gamma.aalen.NIS(1, link, a0, M$zeta.0, M$zeta.1, M$zeta.2, beta.2, nu=nu)
      } else {
        for(i in 1:tm.length){
          Q <- gamma.aalen.Q(tm.grid[i], link, a0, M$zeta.0, M$zeta.1, M$zeta.2, beta.2, beta.4, nu)
          nis[cont,i] <- gamma.aalen.NIS(tm.grid[i], link, a0, M$zeta.0, M$zeta.1, M$zeta.2, 
                                         beta.2, beta.4, nu, Q, approximation=FALSE)
        }
        cat("Countdown: ", length(a)-cont, "\n")
      }
      cont <- cont + 1
    }
  }
  ts <- nds + nis
  
  return(list(NIS=nis, NDS=as.numeric(nds), TS=as.numeric(ts)))
}