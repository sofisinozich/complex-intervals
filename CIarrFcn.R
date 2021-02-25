CIarrFcn<-function (kstar, m, alpha) {     
  #  Inputs: kstar = vector of estimated (non-integer) counts 
  #    m = vector of (estimated effective) sample sizes
  klgth = length(kstar)
  m[m==0] = 1.e-5
  za = qnorm(1 - alpha/2)
  CIarr = array(0, c(klgth, 5, 8), dimnames = list(1:klgth, 
                                                   c("Lowr", "Upr", "Width", "Flag", "LoInd"), 
                                                   c("Wald", "JeffPr", "UnifPr", "ClPe", "Wils", "AgCo", 
                                                     "Assqr", "Logit")))
  # Flag column indicates case of CI-piece outside [0,1]
  # or (in Bayesian intervals) of p-hat outside CI.
  # LoInd column indicates CI-piece < 0 or p-hat < CI.
  IntDef = function(intA, outA, ptA=0, qtA=1) {
    tmpmat = cbind(ifelse(outA[, 1], ptA, intA[, 1]), 
                   ifelse(outA[, 2], qtA, intA[, 2]))
    cbind(tmpmat, tmpmat[,2]-tmpmat[,1], outA[,1] |
            outA[,2], outA[,1] )  }
  # IntDef defines output matrix: for fixed  kstar[i],m[i] & 
  # CI-type, creates: Lowr, Upr, Width, Flag, LoInd 
  
  tmp0 = za * sqrt(kstar * (m - kstar))/sqrt(m^3)
  int0 = cbind(kstar/m - tmp0, kstar/m + tmp0)
  out0 = cbind(int0[, 1] < 0, int0[, 2] > 1)
  CIarr[,, 1] = IntDef(int0, out0)
  
  ptJ = (kstar + 0.5)/(m + 1)
  int0 = cbind(qbeta(alpha/2, kstar + 0.5, m - kstar + 0.5), 
               qbeta(1 - alpha/2, kstar + 0.5, m - kstar + 0.5))
  out0 = cbind(ptJ < int0[, 1], ptJ > int0[, 2])
  CIarr[,, 2] = IntDef(int0, out0, ptJ, ptJ)
  
  ptU = (kstar + 1)/(m + 2)
  int0 = cbind(qbeta(alpha/2, kstar + 1, m - kstar + 1), 
               qbeta(1 - alpha/2, kstar + 1, m - kstar + 1))
  out0 = cbind(ptU < int0[, 1], ptU > int0[, 2])
  CIarr[,, 3] = IntDef(int0, out0, ptU, ptU)
  
  v = 2 * cbind(kstar, m - kstar + 1, kstar + 1, m - kstar)
  CIarr[, -3, 4] = cbind(ifelse(kstar == 0, 0, {
    tmp = v[, 1] * qf(alpha/2, v[, 1], v[, 2])
    tmp/(v[, 2] + tmp)}), ifelse(kstar == m, 1, {
      tmp = v[, 3] * qf(1 - alpha/2, v[, 3], v[, 4])
      tmp/(v[, 4] + tmp)}), kstar == 0 | kstar == m, kstar == 0)
  CIarr[, 3, 4] = CIarr[, 2, 4] - CIarr[, 1, 4]
  
  vt = kstar * (m - kstar)/m^2 + za^2/(4 * m)
  at = (kstar + za^2/2)/(m + za^2)
  bt = (za * (m * vt)^0.5)/(m + za^2)
  int0 = cbind(at - bt, at + bt)
  out0 = cbind(int0[, 1] < 0, int0[, 2] > 1)
  CIarr[,, 5] = IntDef(int0,out0)
  
  mt = m + za^2
  ptA = (kstar + za^2/2)/mt
  st = za * sqrt(ptA * (1 - ptA)/mt)
  int0 = cbind(ptA - st, ptA + st)
  out0 = cbind(int0[, 1] < 0, int0[, 2] > 1)
  CIarr[,, 6] = IntDef(int0,out0)
  
  ast = asin(sqrt(ptJ))
  st = za/sqrt(4 * m)
  int0 = cbind(ast - st, ast + st)
  out0 = cbind(int0[, 1] < 0, int0[, 2] > pi/2)
  CIarr[,, 7] = IntDef(sin(int0)^2,out0)
  
  lgt = ifelse(kstar==0, -Inf, ifelse(kstar==m, Inf, log(kstar/(m-kstar)))) 
  tmp1 = ifelse(kstar==0 | kstar==m, 0, za/sqrt(kstar*(1-kstar/m)))
  int0 = cbind(lgt - tmp1, lgt + tmp1)
  CIarr[,, 8] = IntDef( plogis(int0), array(0, c(klgth,2)) )
  CIarr   }
