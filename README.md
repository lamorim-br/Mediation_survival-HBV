# Mediation_survival-HBV
Data related to the manuscript entitled "Causal mediation for survival data: a unifying approach via GLM" (Taddeo &amp;  Amorim, 2021) published in the Revista Colombiana de Estadistica.  The R files replicate all data analyses presented in the manuscript. The data are drawn from Huang & Yang (2017), so please cite the original work if you ever use any of these data sets for research purposes. 

Information about dataset can be found at README.txt.

R codes for the analyses presented in Section 5.1 of the paper can be found under the topic "Causal-mediation-in-survival analysis", while R codes for the analyses as in Section 5.2 can be found under the topic "Responsiveness-causal-measure".


Topic: Causal-mediation-in-survival analysis

------------------------------------------------------
**Causal mediation-survival models (Tables 4-5)-CIs-main.R**
-------------------------------------------------------

    _R file for reading data set and estimating NID and NIE (Tables 4 & 5) in Taddeo & Amorim (2022)._

#########################################################
DATA*: hcvhbv_liver.txt
CODE:  causal_effects_mediation_CI.R
GOAL:  Estimate NID and NIE using mediation with survival models
       ***   Point and intervalar estimates ***
########################################################

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
