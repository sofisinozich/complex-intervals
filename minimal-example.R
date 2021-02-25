library(survey)

library(dplyr)

source("VarKish.R") # Used to get ingredients for Franco et al.'s proposed calculation of a Kish-type effective sample size
source("CIarrFcn.R") # Used to calculate confidence intervals, given weighted counts and effective smaple size

# Create an example survey design object ----

  ess_rf <- haven::read_dta("chapter_exercises_ess6rf.dta")
  ess_design <- svydesign(id=~psu,strata=~stratify,weights=~pspwght,data=ess_rf)

# Calculate point estimate ----
  prop_est <- svymean(~ voted_lastelection, ess_design, na.rm = TRUE)
  
# Confidence interval using Franco et al's method ----
  
  # Step 1: Calculate Kish-type variance ----
  
    Yfr <- ess_design$variables %>%
      group_by(stratify, psu) %>%
      summarize(
        Cluster_Total = sum(pspwght * voted_lastelection),
        Cluster_Size = sum(pspwght) # !! Check to see whether to just use unweighted count
      ) %>%
      ungroup() %>%
      mutate(Stratum_ID = as.integer(factor(stratify))) %>%
      select(psu, Cluster_Total, Stratum_ID, Cluster_Size)
  
    H <- n_distinct(Yfr[['Stratum_ID']])
    
    nCvec <- tapply(Yfr[['psu']],
                    Yfr[['Stratum_ID']],
                    n_distinct)
    
    knows_what_Kvec_should_be <- FALSE
    if (knows_what_Kvec_should_be) {
      Kvec <- ??? # Maybe just `inf`?
        
        var_kish <- VarKish(
          Yfr, # One row for each sampled cluster
          # Column (1): Cluster totals Y_{hk+},
          # Column (2): Index for clusters h in 1:H
          # Column (3): Cluster size M_{kh}
          H, # Known number of strata [perhaps not all sampled]
          nCvec, # nCvec[h] = count of sampled clusters in strat h
          Kvec # Kvec[h] = known number of pop clusters in strat h
        )
    }
  
  # Step 2: Calculate effective sample size
    
    var_prop_srs <- as.numeric(
      (1/n) * ( prop_est * (1 - prop_est) )
    )
    
    deff_var_kish <- var_kish/var_prop_srs # Not quite right since var_kish is variance of a total
    n_eff_var_kish <- n / deff_var_kish
    
  # Step 3: Calculate intervals
    
    weighted_count <- n_eff_var_kish * prop_est
    
    CIarrFcn(
      kstar = weighted_count,
      m = n_eff_var_kish,
      alpha = 0.05
    )
  
    