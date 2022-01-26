# Mediation_survival-HBV
Data related to the manuscript entitled "Causal mediation for survival data: a unifying approach via GLM" (Taddeo &amp;  Amorim, 2022) published in the Revista Colombiana de Estadística - Applied Statistics, January 2022, volume 45, issue 1, pp. 161-191 (http://doi.org/10.15446/rce.v45n1.94553).  The R files replicate all data analyses presented in the manuscript. The data are drawn from Huang & Yang (2017), so please cite the original work if you ever use any of these data sets for research purposes. 


Information about dataset can be found at README.txt.

R codes for the analyses presented in Section 5.1 of the paper can be found under the topic "Causal-mediation-in-survival analysis", while R codes for the analyses as in Section 5.2 can be found under the topic "Responsiveness-causal-measure".

------------------------------------------------------
------------------------------------------------------
# Topic: Causal-mediation-in-survival analysis

------------------------------------------------------
**Causal mediation-survival models (Tables 4-5)-CIs-main.R**
-------------------------------------------------------

    _R file for reading data set and estimating NID and NIE (Tables 4 & 5) in Taddeo & Amorim (2022)._

############################################################################

DATA*: hcvhbv_liver.txt

CODE:  causal_effects_mediation_CI.R

GOAL:  Estimate NID and NIE using mediation with survival models

***   Point and intervalar estimates ***

#############################################################################

 * Link for download of complete data set provided by Huang & Yang (2017)

 _R code for the main analysis in Section 5.1 (Tables 4 & 5) (Taddeo & Amorim, 2022)_

--------------------------------------------------------
  **Required packages**
--------------------------------------------------------

rm(list=ls())

library(survival)

library(timereg)

library(data.table)

library(MASS)

-----------------------------------
  **Reading data set**
-----------------------------------
dados <- read.table("hcvhbv_liver.txt", header = TRUE)

head(dados)

attach(dados)

-----------------------------------
**Specifying covariate set**
-----------------------------------
X <- dados[,c("agegp2", "agegp3", "agegp4", "gender", "smoke", "alcohol")]

p <- ncol(X)

X <- setDT(X)[, .N, by = c(names(X))]

Prob <- X$N / sum(X$N)

X <- cbind(X, Prob)

-----------------------------------
**R functions to estimates NDE e NIE (and their 95%CIs) varying models for mediator and outcome**
-----------------------------------
source("Causal-NIE-NDE.R")


#### _Running the R function to print NDE e NIE varying models for mediator and outcome_
-----------------------------------


----------------------------------------------------------------------------------------------------------
### Case 1: Normal GLM for mediator and varying survival model for outcome (AFT, Cox, Aalen)
---------------------------------------------------------------------------------------------------------
#### AFT
nat.effects <- cox.aft.eff(mediator="gaussian", method="aft", link="identity")

c(nat.effects$NDE,nat.effects$DE.lo,nat.effects$DE.up)

c(nat.effects$NIE,nat.effects$IE.lo,nat.effects$IE.up)

####  Cox
nat.effects <- cox.aft.eff(mediator="gaussian", method="cox", link="identity")

c(exp(nat.effects$NDE),exp(nat.effects$DE.lo),exp(nat.effects$DE.up))

c(exp(nat.effects$NIE),exp(nat.effects$IE.lo),exp(nat.effects$IE.up))

####  Aalen
nat.effects <- aalen.eff(mediator="gaussian", link="identity")

c(nat.effects$NDE,nat.effects$DE.lo,nat.effects$DE.up)

c(nat.effects$NIE,nat.effects$IE.lo,nat.effects$IE.up)

------------------------------------------------------------------------------------------------------------------------------
### Case 2: Gamma GLM (log link function) for mediator and varying survival model for outcome (AFT, Cox, Aalen)
-------------------------------------------------------------------------------------------------------------------------------
####  AFT
nat.effects <- cox.aft.eff(mediator="gamma", method="aft", link="log")

c(nat.effects$NDE,nat.effects$DE.lo,nat.effects$DE.up)

c(nat.effects$NIE,nat.effects$IE.lo,nat.effects$IE.up)

####  Cox
nat.effects <- cox.aft.eff(mediator="gamma", method="cox", link="log")

c(exp(nat.effects$NDE),exp(nat.effects$DE.lo),exp(nat.effects$DE.up))

c(exp(nat.effects$NIE),exp(nat.effects$IE.lo),exp(nat.effects$IE.up))

####  Aalen (it can take some minutes to run)
nat.effects <- aalen.eff(mediator="gamma", link="log")

c(nat.effects$NDE,nat.effects$DE.lo,nat.effects$DE.up)

c(nat.effects$NIE,nat.effects$IE.lo,nat.effects$IE.up)

------------------------------------------------------------------------------------------------------------------------
### Case 3: Gamma GLM (inverse link function) for mediator and varying survival model for outcome (AFT, Cox, Aalen)
------------------------------------------------------------------------------------------------------------------------
####  AFT
nat.effects <- cox.aft.eff(mediator="gamma", method="aft", link="inverse")

c(nat.effects$NDE,nat.effects$DE.lo,nat.effects$DE.up)

c(nat.effects$NIE,nat.effects$IE.lo,nat.effects$IE.up)

####  Cox
nat.effects <- cox.aft.eff(mediator="gamma", method="cox", link="inverse")

c(exp(nat.effects$NDE),exp(nat.effects$DE.lo),exp(nat.effects$DE.up))

c(exp(nat.effects$NIE),exp(nat.effects$IE.lo),exp(nat.effects$IE.up))

####  Aalen (it can take some minutes to run)
nat.effects <- aalen.eff(mediator="gamma", link="inverse")

c(nat.effects$NDE,nat.effects$DE.lo,nat.effects$DE.up)

c(nat.effects$NIE,nat.effects$IE.lo,nat.effects$IE.up)

---------------------------------------------------------------------------------------------------------------------
### Case 4: Gamma GLM (log link function) for mediator and varying survival model for outcome (AFT, Cox, Aalen)
### Changing values for X
---------------------------------------------------------------------------------------------------------------------
####  AFT
nat.effects <- cox.aft.eff(mediator="gamma", method="aft", link="log",x=c(0,0,1,1,1,1))

c(nat.effects$NDE,nat.effects$DE.lo,nat.effects$DE.up)

c(nat.effects$NIE,nat.effects$IE.lo,nat.effects$IE.up)


------------------------
### Comments
------------------------
#(1) AFT, Cox and Aalen estimates using Gamma GLM for mediator are conditional on X
    (here X = 0. This setup can be changed in line 158 and 250)

#(2) NIE in Aalen model with Gamma GLM for mediator is time-dependent. An approximated value is presented.



------------------------------------------------------
------------------------------------------------------
# Topic: Responsiveness-causal-measure

------------------------------------------------------
**HBV_responsiveness-survival.R**
-------------------------------------------------------

    _R file for generating figures with responsiveness measures (Figures 4 & 5) in Taddeo & Amorim (2022)._


###########################################################################

DATA*: hcvhbv_liver.txt

CODE:  Responsiveness_AFT_Figure4.R

GOAL:  Estimate NDR and NIR using AFT models

***   Varying mediation model ***

###########################################################################

 * Link for download of complete data set provided by Huang & Yang (2017)

 _R code for the main analyses in Section 5.2 (Figure 4) (Taddeo & Amorim, 2022)_
 
 --------------------------------------------------------
  **Required packages**
--------------------------------------------------------

rm(list=ls())

library(survival)

library(timereg)

library(data.table)

library(MASS)

-----------------------------------
  **Reading data set**
-----------------------------------
dados <- read.table("hcvhbv_liver.txt", header = TRUE)

head(dados)

attach(dados)

-----------------------------------
**Specifying covariate set**
-----------------------------------
qq <- quantile(dados$logc[dados$logc>0])

dd <- qq["50%"] - qq["25%"]

dados$logc <- dados$logc / dd

X <- dados[,c("agegp2", "agegp3", "agegp4", "gender", "smoke", "alcohol")]

p <- ncol(X)

X <- setDT(X)[, .N, by = c(names(X))]

Prob <- X$N / sum(X$N)

X <- cbind(X, Prob)


----------------------------------------------------------------------------------
**R functions to estimate NDR e NIR for varying mediator models for AFT**
-----------------------------------------------------------------------------------
source("Responsiveness-AFT.R")


----------------------------------------------------------------------------------------------------
### NDR and NIR for varying mediator model using AFT for outcome
----------------------------------------------------------------------------------------------------
a <- seq(0,20,length.out=1000)

nat.sen.gama.log <- cox.aft.sen(mediator="gamma", method="aft", link="log", a=a)

nat.sen.gama.inv <- cox.aft.sen(mediator="gamma", method="aft", link="inverse", a=a)

nat.sen.gaus <- cox.aft.sen(mediator="gaussian", method="aft", link="identity", a=a)


----------------------------------------------------------------------------------------------------------
### Figure 4 A: NDR
---------------------------------------------------------------------------------------------------------
plot(a, nat.sen.gama.log$NDS, type="l", 
     xlab="log of HCV viral load", ylab="NDS",
     lwd=2, cex.axis=1.3, cex.lab=1.3)

lines(a, nat.sen.gama.inv$NDS, lty="dashed", lwd=2)

lines(a, nat.sen.gaus$NDS, lty="dotted", lwd=3)

----------------------------------------------------------------------------------------------------------
### Figure 4 B: NIR
---------------------------------------------------------------------------------------------------------
plot(a, nat.sen.gama.log$NIS, type="l", 
     xlab="log of HCV viral load", ylab="NIS",
     ylim=c(0,8), lwd=2, cex.axis=1.3, cex.lab=1.3)

lines(a, nat.sen.gama.inv$NIS, lty="dashed", lwd=2)

lines(a, nat.sen.gaus$NIS, lty="dotted", lwd=3)


###########################################################################
###########################################################################

DATA*: hcvhbv_liver.txt

CODE:  Responsiveness_Aalen_Figure5.R

GOAL:  Estimate NIR and NDR using Aalen models

***   Varying mediation model ***

###########################################################################

 * Link for download of complete data set provided by Huang & Yang (2017)

 _R code for the main analyses in Section 5.2 (Figure 5) (Taddeo & Amorim, 2022)_
 
 Same data as before.
 
----------------------------------------------------------------------------------
**R functions to estimate NDR e NIR for varying mediator models for AFT**
-----------------------------------------------------------------------------------
source("Responsiveness-Aalen.R")


----------------------------------------------------------------------------------------------------
### NDR and NIR for varying mediator model using AFT for outcome
----------------------------------------------------------------------------------------------------
a <- seq(0,20,length.out=20)

----------------------------------------------------------------------------------------------------------
### Figure 5 A: NIR using gamma mediator model, with log link function
---------------------------------------------------------------------------------------------------------
nat.sen <- aalen.sen(mediator="gamma", link="log", a=a, x=rep(0,p), 
                     grid.lo=0, grid.up=20, by=2, approximation=FALSE)

persp(a, seq(0,20,by=2), 1000*nat.sen$NIS, 
      theta=25, phi=10, r=15, d=.3, ticktype="detailed", zlab="\n\nNIS (x 1,000)", 
      ylab="\nYear", xlab="\nlog of HCV viral load", 
      zlim=c(-2.3,-0.5), cex.axis=.75, shade=.5, nticks=4)

----------------------------------------------------------------------------------------------------------
### Figure 5 B: NIR using gamma mediator model, with inverse link function
---------------------------------------------------------------------------------------------------------
nat.sen <- aalen.sen(mediator="gamma", link="inverse", a=a, x=rep(0,p),
          grid.lo=0, grid.up=20, by=2, approximation=FALSE)


persp(a, seq(0,20,by=2), 1000*nat.sen$NIS, 
      theta=25, phi=10, r=15, d=.3, ticktype="detailed", zlab="\n\nNIS (x 1,000)", 
      ylab="\nYear", xlab="\nlog of HCV viral load", 
      zlim=c(-2.3,-0.5), cex.axis=.75, shade=.5, nticks=4)
