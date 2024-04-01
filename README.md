# GIST Intratumoral Microbiome Analysis
Code used in 'Intratumoral Shewanella algae promotes liver metastasis and dampens adjuvant imatinib treatment in gastrointestinal stromal tumor'
 
 Contains feature selection and feature ranking models and agrithom code and Illumina 16s rRNA amplicon sequencing original microbiome read count matrix

 ## The feature selection script
 
 ### The feature selection script related to the feature selection model is placed in the "modeling.R script", which is dirven by the SIAMCAT (https://siamcat.embl.de/) package in R.
 
  Including batch-correction protocol, feature selection model training and validation. 
  
  The training data used in the paper is placed in the Risk Straitification Model Data and Liver Metastasis Model Data folders. 
  
  You can also use your own data origin from tumor microbes to test the feature selection model.

  You can run "Rscript modeling.R" or activate R environment to run this script 'source("modeling.R").
  Of course you need to install all R packages dependencies.

  The feature selection LASSO Machine learning model related Figures in the paper will be reproduced.

## The feature engineering method

### The feature engineering method including Priority Score caculation and Time-dependent explanations of survival model

 The script "confident_interval_survshap.r" is the Time-dependent explanations of survival (SurvSHAP,https://github.com/MI2DataLab/survshap, https://doi.org/10.1016/j.knosys.2022.110234) model for average importance score caculation and High Priority feature genera (14) ranking.

 The script "randomforest_feature.py" is the random forest agrithom for feature genera's AI (Adjuvant imatinib) Priority Score caculation.

 ## Mediation analysis between Risk stratifications and liver metastasis

 The script "mediation_effect.r" is the original code for the mediation analysis.
 
 We hypothesize that the only intermediate variable between risk and liver metastasis is the adjuvant imatinib.
 
 For 95% CI caculation, we utilized the Bias Corrected and Accelerated (BCa) Bootstrap resampling method.


 
 ## Over all study design in this project
 
![image](https://github.com/GIST-microbiome/GIST.github/assets/143196047/43130b82-62bf-414a-b295-90ee1e3ce1d8)
