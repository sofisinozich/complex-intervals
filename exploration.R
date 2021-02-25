library(survey)
library(haven)

source("CIarrFcn.R")

# Integrating survey package output in with the existing function
surveyci <- function(formula,design,deff) {
  if(length(all.vars(formula))>1) {
    message("Nope.")
    return()
  }
  p <- svymean(formula,design,na.rm=TRUE,deff=deff)
  n_eff <- dim(design$variables[all.vars(formula)])[1]/deff(p)
  
  CIarrFcn(p*n_eff,n_eff,0.05)
}

# Example
ess_rf <- read_dta("chapter_exercises_ess6rf.dta")
ess_design <- svydesign(id=~psu,strata=~stratify,weights=~pspwght,data=ess_rf)
surveyci(~voted_lastelection,ess_design,deff="replace")

voted_mean <- svymean(~voted_lastelection,ess_design,na.rm=TRUE,deff="replace")
# Allow infinite df
# Wald
confint(voted_mean,df=Inf)
# Logit
confint(svyciprop(~voted_lastelection,ess_design,method="logit",df=Inf))
# Similar but not identical results so far...