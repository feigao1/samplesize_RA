# samplesize_RA
This is the R code for sample size calculation for active-arm trial with counterfactual incidence based on recency assay.

samplesize_cal.R is the main code for replicating results in the arxiv paper: https://arxiv.org/abs/2011.00725.

setting_code.R includes parameters for two settings for MSM and women trials.

recencyassay_code.R includes main computation codes.
    samplesize_CF: a function that calculates sample size with given input parameters para, hypothesis R0 vs R1, type I error alpha, and power beta.
    power_CF: a function that calculates power with given input parameters para, hypothesis R0 vs R1, type I error alpha, and sample size N.
