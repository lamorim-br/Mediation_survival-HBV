#----------------------------------------
#  Fitting GLM models for the mediator 
#  --------------------------------------
M.mod <- function(mediator="gaussian", link="identity"){
  if(mediator=="gaussian"){
    mod.glm <- glm(HBV ~ HCV.cat + agegp2 + agegp3 + agegp4 + gender + smoke +alcohol, 
                   family = gaussian, data=dados)
    Cov <- summary(mod.glm)$cov.scaled
    scl.parm <- NA
    scl.parm.std = NA
  } else if(mediator == "gamma"){
    mod.glm <- glm(HBV ~ HCV.cat + agegp2 + agegp3 + agegp4 + gender + smoke +alcohol, 
                   family = Gamma(link=link), data=dados)
    Cov <- summary(mod.glm)$cov.scaled
    scl.parm <- as.numeric(gamma.shape(mod.glm, verbose=FALSE)[1])
    scl.parm.std <- as.numeric(gamma.shape(mod.glm, verbose=FALSE)[2])
  }
  return(list(zeta.0=mod.glm$coeff[1], zeta.1=mod.glm$coeff[2], zeta.2=mod.glm$coeff[-c(1,2)], 
              Cov=Cov, scl.parm=scl.parm, scl.parm.std=scl.parm.std))
}

gamma.theta <- function(link="identity", a, x, zeta.0, zeta.1, zeta.2){
  if(link == "identity"){
    return(-1/(zeta.0 + zeta.1 * a + t(zeta.2) %*% x))
  } else if(link == "log"){
    return(-exp(-(zeta.0 + zeta.1 * a + t(zeta.2) %*% x)))
  } else if(link == "inverse"){
    return(-(zeta.0 + zeta.1 * a + t(zeta.2) %*% x))
  }
}

gamma.aalen.Q <- function(tm, link="identity", a.star, zeta.0, zeta.1, zeta.2, 
                          beta.2, beta.4, nu){
  N <- nrow(X)
  B2 <- beta.2 * tm
  B4 <- beta.4 * tm
  Q <- rep(NA, N)
  dnm <- 0
  for(i in 1:N){
    x.aux <- as.numeric(X[i,1:p])
    theta.aux <- gamma.theta(link, a.star, x.aux, zeta.0, zeta.1, zeta.2)
    Q[i] <- (((theta.aux/(theta.aux+B2/nu))^nu)*exp(t(beta.4) %*% x.aux)) * X$Prob[i]
    dnm <- dnm + Q[i]
  }
  return(Q / dnm)
}

gamma.aalen.NIE <- function(tm=1, link="identity", a=1, a.star=0, zeta.0, zeta.1, zeta.2, beta.2, nu=NULL, Q=NULL, approximation=TRUE){
  N <- nrow(X)
  nie <- 0
  for(i in 1:N){
    x.aux <- as.numeric(X[i,1:p])
    theta.star <- gamma.theta(link, a.star, x.aux, zeta.0, zeta.1, zeta.2)
    theta <- gamma.theta(link, a, x.aux, zeta.0, zeta.1, zeta.2)
    if(approximation){
      nie <- nie + beta.2 * ((1/theta.star) - (1/theta)) * X$Prob[i]
    } else {
      B2 <- beta.2 * tm
      nie <- nie + beta.2 * ((1/(theta.star + B2/nu)) - (1/(theta + B2/nu))) * Q[i]
    }
  }
  return(nie)
}

#------------------------------------------------------------
#  Fitting survival models (Cox and AFT) for the outcome 
#  ----------------------------------------------------------
cox.aft.eff <- function(method="cox", mediator="gaussian", link="identity", a.star=0, a=1, x=rep(0,p), G=10^4){
  if(method == "cox"){
    mod <- coxph(Surv(hcc.time/365.25, status) ~ HCV.cat + HBV + gender + agegp2 + 
                   agegp3 + agegp4 + smoke + alcohol, data = dados)
    beta.1 <- as.numeric(mod$coefficients[1])
    beta.2 <- as.numeric(mod$coefficients[2])
    beta.4 <- as.numeric(mod$coefficients[-c(1,2)])
    T.cov <- mod$var[1:2,1:2]
  } else {
    mod <- survreg(Surv(hcc.time/365.25, status) ~ HCV.cat + HBV + gender + agegp2 + 
                     agegp3 + agegp4 + smoke + alcohol, dist='weibull', data = dados)
    beta.1 <- as.numeric(mod$coefficients[2])
    beta.2 <- as.numeric(mod$coefficients[3])
    beta.4 <- as.numeric(mod$coefficients[-c(1,2,3)])
    T.cov <- mod$var[2:3,2:3]
  }
  
  #----------------------------------------------------------------------------------
  #  Computing NDE and NIE (models without interaction between exposure and mediator) 
  #-----------------------------------------------------------------------------------
  ## NDE without interaction:
  nde <- beta.1 * (a - a.star)
  if(mediator == "gaussian"){  
    M <- M.mod()
    zeta.0 <- M$zeta.0
    zeta.1 <- M$zeta.1
    zeta.2 <- M$zeta.2
    nie <- beta.2 * zeta.1 * (a - a.star)
    # Confidence intervals:
    M.cov <- M$Cov
    ci <- cox.aft.eff.ci(beta.1, beta.2, T.cov, zeta.0, zeta.1, zeta.2, -1, M.cov, -1, G, a, a.star, x, mediator, link)
  } else if(mediator == "gamma"){
    M <- M.mod(mediator, link)
    zeta.0 <- M$zeta.0
    zeta.1 <- M$zeta.1
    zeta.2 <- M$zeta.2
    theta.star <- gamma.theta(link, a.star, x, zeta.0, zeta.1, zeta.2)
    theta <- gamma.theta(link, a, x, zeta.0, zeta.1, zeta.2)
    nu <- M$scl.parm
    nie <- nu * log(((beta.2/nu + theta.star)/theta.star)/((beta.2/nu + theta)/theta))
    # Confidence intervals:
    M.cov <- M$Cov
    scl.parm.std <- M$scl.parm.std
    ci <- cox.aft.eff.ci(beta.1, beta.2, T.cov, zeta.0, zeta.1, zeta.2, nu, M.cov, scl.parm.std, G, a, a.star, x, mediator, link)
  }
  te <- nde + nie
  
  return(list(NIE=as.numeric(nie), NDE=as.numeric(nde), TE=as.numeric(te),  
              IE.lo=ci$IE.lo, IE.up=ci$IE.up, DE.lo=ci$DE.lo, DE.up=ci$DE.up, TE.lo=ci$TE.lo, TE.up=ci$TE.up))
}


#----------------------------------------
#  Fitting Aalen model for the outcome 
#  ---------------------------------------
aalen.eff <- function(mediator="gaussian", link="identity", a.star=0, a=1, x=rep(0,p), G=10^4, grid.lo=0, grid.up=15, by=0.05, approximation=TRUE){
  mod <- aalen(Surv(hcc.time/365.25, status) ~ const(HCV.cat) + const(HBV) + const(gender) + 
                 const(agegp2) + const(agegp3) + const(agegp4) + const(smoke) + const(alcohol), data = dados)
  beta.1 <- as.numeric(mod$gamma[1])
  beta.2 <- as.numeric(mod$gamma[2])
  beta.4 <- as.numeric(mod$gamma[-c(1,2)])
  T.cov <- mod$robvar.gamma[1:2,1:2]
  
  tm.grid <- seq(from=grid.lo, to=grid.up, by=by)
  tm.length <- length(tm.grid)
  
  nde <- beta.1
  if(mediator == "gaussian"){
    M <- M.mod()
    nie <- beta.2 * M$zeta.1
    M.cov <- M$Cov
    ci <- aalen.ci(beta.1, beta.2, T.cov, M$zeta.0, M$zeta.1, M$zeta.2, -1, M.cov, -1, G, a, a.star, x, mediator, link)
  } else if(mediator == "gamma"){
    M <- M.mod(mediator, link)
    M.cov <- M$Cov
    nu <- M$scl.parm
    scl.parm.std <- M$scl.parm.std
    if(approximation){
      nie <- gamma.aalen.NIE(1, link, a, a.star, M$zeta.0, M$zeta.1, M$zeta.2, beta.2)
      ci <- aalen.ci(beta.1, beta.2, T.cov, M$zeta.0, M$zeta.1, M$zeta.2, nu, M.cov, scl.parm.std, G, a, a.star, x, mediator, link)
    } else {
      nie <- rep(0, tm.length)
      for(i in 1:tm.length){
        Q <- gamma.aalen.Q(tm.grid[i], link, a.star, M$zeta.0, M$zeta.1, M$zeta.2, beta.2, beta.4, nu)
        nie[i] <- gamma.aalen.NIE(tm.grid[i], link, a, a.star, M$zeta.0, M$zeta.1, M$zeta.2, beta.2, nu, Q, approximation=FALSE)
      }
   
      ci <- aalen.ci(beta.1, beta.2, T.cov, M$zeta.0, M$zeta.1, M$zeta.2, nu, M.cov, scl.parm.std, G, a, a.star, x, mediator, link)
    }
  }
  te <- nde + nie
  
  return(list(NIE=as.numeric(nie), NDE=as.numeric(nde), TE=as.numeric(te), 
              IE.lo=ci$IE.lo, IE.up=ci$IE.up, DE.lo=ci$DE.lo, DE.up=ci$DE.up, TE.lo=ci$TE.lo, TE.up=ci$TE.up,
              zeta.0=M$zeta.0, zeta.1=M$zeta.1, zeta.2=M$zeta.2, scl.parm=M$scl.parm))
}


#----------------------------------------------------------------------------------
#  Computing 95%CIs for NDE e NIE varying models for mediator and outcome
#-----------------------------------------------------------------------------------
cox.aft.eff.ci <- function(beta.1, beta.2, T.cov, zeta.0, zeta.1, zeta.2, scl.parm=-1, M.cov, scl.parm.std=-1, G, a, a.star, x, mediator, link){
  require(mvtnorm)
  IE <- DE <- TE <- rep(0,G)
  set.seed(137)
  
  beta.sim <- rmvnorm(G, mean=c(beta.1, beta.2), sigma=T.cov) 
  zeta.sim <- rmvnorm(G, mean=c(zeta.0, zeta.1, zeta.2), sigma=M.cov)
  if(scl.parm>0 & scl.parm.std>0){
    scl.parm.sim <- rnorm(G, scl.parm, scl.parm.std)
  }
  
  DE <- beta.sim[,1] * (a-a.star)
  
  if(mediator=="gaussian"){
    for(g in 1:G){
      beta.2.sim <- beta.sim[g,2]
      zeta.1.sim <- zeta.sim[g,2]
      
      IE[g] <- beta.2.sim * zeta.1.sim * (a - a.star)
    }
  } else if(mediator=="gamma"){
    for(g in 1:G){
      beta.2.sim <- beta.sim[g,2]
      zeta.0.sim <- zeta.sim[g,1]
      zeta.1.sim <- zeta.sim[g,2]
      zeta.2.sim <- zeta.sim[g,-c(1,2)]
      nu.sim <- scl.parm.sim[g]
      
      theta.star <- gamma.theta(link, a.star, x, zeta.0.sim, zeta.1.sim, zeta.2.sim)
      theta <- gamma.theta(link, a, x, zeta.0.sim, zeta.1.sim, zeta.2.sim)
      
      IE[g] <- nu.sim*log(((beta.2.sim/nu.sim+theta.star)/theta.star)/((beta.2.sim/nu.sim+theta)/theta))
    }
  }
  
  TE <- IE + DE
  
  IE.ci <- quantile(IE, c(0.025, 0.975))  
  DE.ci <- quantile(DE, c(0.025, 0.975))  
  TE.ci <- quantile(TE, c(0.025, 0.975))
  
  return(list(IE.lo = IE.ci[1], IE.up = IE.ci[2], DE.lo = DE.ci[1], DE.up = DE.ci[2], TE.lo = TE.ci[1], TE.up = TE.ci[2]))
}

aalen.ci <- function(beta.1, beta.2, T.cov, zeta.0, zeta.1, zeta.2, scl.parm=-1, M.cov, scl.parm.std=-1, G, a, a.star, x=rep(0,p), mediator, link){
  require(mvtnorm)
  IE <- DE <- TE <- rep(0,G)
  set.seed(137)
  
  beta.sim <- rmvnorm(G, mean=c(beta.1, beta.2), sigma=T.cov) 
  zeta.sim <- rmvnorm(G, mean=c(zeta.0, zeta.1, zeta.2), sigma=M.cov)
  if(scl.parm>0 & scl.parm.std>0){
    scl.parm.sim <- rnorm(G, scl.parm, scl.parm.std)
  }
  
  DE <- beta.sim[,1] * (a-a.star)
  
  if(mediator=="gaussian"){
    for(g in 1:G){
      beta.2.sim <- beta.sim[g,2]
      zeta.1.sim <- zeta.sim[g,2]
      
      IE[g] <- beta.2.sim * zeta.1.sim * (a - a.star)
    }
  } else if(mediator=="gamma"){
    for(g in 1:G){
      beta.2.sim <- beta.sim[g,2]
      zeta.0.sim <- zeta.sim[g,1]
      zeta.1.sim <- zeta.sim[g,2]
      zeta.2.sim <- zeta.sim[g,-c(1,2)]
      nu.sim <- scl.parm.sim[g]
      
      theta.star <- gamma.theta(link, a.star, x, zeta.0.sim, zeta.1.sim, zeta.2.sim)
      theta <- gamma.theta(link, a, x, zeta.0.sim, zeta.1.sim, zeta.2.sim)
      
      IE[g] <- gamma.aalen.NIE(1, link, a, a.star, zeta.0.sim, zeta.1.sim, zeta.2.sim, beta.2.sim)
    }
  }
  
  TE <- IE + DE
  
  IE.ci <- quantile(IE, c(0.025, 0.975))  
  DE.ci <- quantile(DE, c(0.025, 0.975))  
  TE.ci <- quantile(TE, c(0.025, 0.975))
  
  return(list(IE.lo = IE.ci[1], IE.up = IE.ci[2], DE.lo = DE.ci[1], DE.up = DE.ci[2], TE.lo = TE.ci[1], TE.up = TE.ci[2]))
}
