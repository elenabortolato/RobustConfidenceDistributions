# RobustConfidenceDistributions
code associated to the paper "On approximate robust confidence distributions" (Bortolato and  Ventura 2023)

### Folder "simulation study section 5"
- RDS files containing the results of the simulation study
- results_ks_w.R: R code for producing tables 7 8 of the paper using RDS files in the folder
- simulations_ks_w.R: R code for running the simulation study


### Folder "simulation study section 4"
- R code for reproducing the results in each scenario
- R code for reproduce tables 3, 4, 5, 6 using RData files "results_explorer.R"
- R code for results based on Bootstrap methods

- RData files with the results of each simulation study, in particular

- files cont_40, cont_80, nocont_40, nocont_80 contain results for when the true parameter value is 2.6
- files cont_40_1, cont_80_1, nocont_40_1, nocont_80_1  contain results for when the true parameter value is 1

- "nocont": non contaminated scenario, "cont": contaminated scenario
- "40" and "80" are used to indicate the total sample size

### Tutorial based on the a real data example in section 4.2 
(soon)
