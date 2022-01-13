-----------------------------------
Software & Data
-----------------------------------

Software: R

Operating system: Microsoft Windows 10

Replication Data for Analyses in the manuscript "Causal mediation for survival data: a unifying approach via GLM" (2021-10-29)

-----------------------------------
File list
-----------------------------------

Causal mediation-survival models (Tables 4-5)-CIs-main.R        R file for reading data set and estimating NIE and NDE (Tables 4 & 5)
Causal-NIE-NDE.R                                                R functions used for estimation of NIE and NDE 
HCV_responsiveness_survival.R                                   R file for generating figures with responsiveness measures (Figures 4 & 5)
hcvhbv_liver.txt                                                HCV (parcial) data set from Huang and Yang* (2017)    

* Link for download of complete data set provided by Huang & Yang (2017)

-----------------------------------
Data dictionary
-----------------------------------

hcc_time:           time to liver cancer diagnosis (in days)
status:             indicator of event (1=case, 0=censored)
agegr2:             dummy for age group2 (40-49 years, reference: 30-39 years)
agegr3:             dummy for age group3 (50-59 years, reference: 30-39 years)
agegr4:             dummy for age group4 (60-69 years, reference: 30-39 years)
gender:             gender indicator (0= female, 1=male) 
smoke:              cigarette smoking indicator (1=yes, 0=no)
alcohol:            alcohol comsuption indicator (1=yes, 0=no)
HCV.cat             HCV viral load indicator (1=detected, 0=non detected)
HBV:                natural logarithm of HBV viral load

----------------------------------
Sources of proprietary data
-----------------------------------

Source: Huang, Y-T and Yang, H-I. (2017) Causal mediation analysis of survival outcomes with multiple mediators. Epidemiology, 28(3), 370–378.**

Huang, Y-T, Jen, C-Lm Yang H-I, Lee M-H, Lu S-N, Iloeje U.H, Chen C-J. (2011) Lifetime risk and sex difference of hepatocellular carcinoma among patients with
chronic hepatitis B and C. Journal of Clinical Oncology, 29(27), 3643–3650.


** The supplementary appendix provides details on the design of the hepatitis study and the R codes for mediation causal analysis using Cox and Aalen models. 

