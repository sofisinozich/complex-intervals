# Description from Supplementary Materials section of Franco et al. paper
# (Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6690503/)
#
# "The R function VarKish calculates the design-based and Kish-type variances for a survey-weighted total of a binary attribute Yhki
#  in the setting of a stratified, single-stage, cluster-sample (in which all units are taken from each sampled cluster)"

VarKish <- function (Yfr, H, nCvec, Kvec) {
  # input Yfr has 3 cols, with 1 row for each sampled cluster
  # H is known number of strata [perhaps not all sampled]
  # cluster tots Y_{hk+}, clusters h in 1:H, and sizes M_{kh}
  # nCvec[h] = count of sampled clusters in strat h
  nsc = sum(nCvec)
  
  nstr = sum(nCvec>0)
  # Kvec[h] = known number of pop clusters in strat h
  csiz1 = (sum(Yfr[[3]]-1)==0) ## T if all cluster-sizes 1 
  
  miss2zer = function(vec) replace(vec, is.na(vec),0)
  #  function to zero out NA elements of a numeric vector
  var0  = function(x) if(length(x)>1) var(x) else 0
  sd0   = function(x) if(length(x)>1) sd(x) else 0
  
  nt.st = miss2zer( tapply(Yfr[[3]],Yfr[[2]],sum) )
  Msq.est = miss2zer( tapply(Yfr[[3]]^2,Yfr[[2]],sum) )
  Mvr.est = miss2zer( tapply(Yfr[[3]],Yfr[[2]],var0) )
  # NB var's in size-1 clusters defined as 0.
  Ymn.st = miss2zer( tapply(Yfr[[1]],Yfr[[2]],mean) )
  Yvr.st = miss2zer( tapply(Yfr[[1]],Yfr[[2]],var0) )
  nt = sum(nt.st)     
  Yest = sum(Kvec*Ymn.st)
  tauhat = Ymn.st*nCvec/nt.st
  sighsq = ifelse(nt.st >= 2, 
                  nt.st*tauhat*(1-tauhat)/(nt.st-1), sgsqhat)
  
  ### rhohat corrected as of 5/14/18 in unequal-cluster-size case
  rhohat = if(csiz1) 0 else
    max(0, 1 - (1-H/nt)*sum(Yfr[[1]]*(Yfr[[3]]-Yfr[[1]]))/
          sum((tauhat*(1-tauhat))[Yfr[[2]]]*Yfr[[3]]*(Yfr[[3]]-1)) )
  Vdsgn = sum(Kvec*(Kvec/nCvec-1)*Yvr.st)
  Vdsgn = ifelse(Vdsgn!=0, Vdsgn,  
                 (Yest+1/2)*(N-Yest+1/2)*(N*(N-nt)/(N+1)^2) )
  VestM = sum(tauhat^2*Kvec*(Kvec/nCvec-1)*Mvr.est) + 
    sum(sighsq*(Kvec/nCvec-1)*(Kvec/nCvec)*
          ((1-rhohat)*nt.st + rhohat*Msq.est))
  list(Vdsgn = Vdsgn, VestM = VestM, Yest = Yest, 
       tauhat= tauhat, sighsq = sighsq, rhohat= rhohat) 
}
