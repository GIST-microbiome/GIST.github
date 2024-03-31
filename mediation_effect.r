library(psych)
library(lavaan)
library(ggplot2)
#library(readxl)
library(semPlot)
library(knitr)
library(dplyr)
library(MBESS)
library(tidyverse)

#AI_med<-read.csv("Risk-AI-LM-Mediation.csv",header=TRUE)
#AI_med<-"Risk-AI-LM-Mediation.csv" %>% read_csv()
#micro_AI_med<-read.csv("Risk-shannon-LM-Mediation.csv",header=TRUE)
AI_med %>% 
  headTail() %>% 
  kable()
#AI_med %>% 
#  select(Risk,Adjuvant_imatinib,Liver_metastasis) %>% 
#  pairs.panels()

correlations <- cor(AI_med)
print(correlations)

mediation_model <- '
  # a path
         Adjuvant_imatinib ~ a * Risk

         # b path
         Liver_metastasis ~ b * Adjuvant_imatinib

         # c prime path 
         Liver_metastasis ~ cp * Risk

         # indirect and total effects
         ab := a * b
         total := cp + ab
'
# Estimate the mediation model
mediation_results <- sem(mediation_model, data = AI_med)

# Summarize the results
summary(mediation_results, standardized = TRUE, fit.measures = TRUE)
parameterestimates(mediation_results, boot.ci.type = "bca.simple", standardized = TRUE) %>% 
  kable()

with(AI_med, mediation(x = Risk, mediator = Adjuvant_imatinib, dv = Liver_metastasis, bootstrap = TRUE, which.boot = "BCa", B = 1000))

