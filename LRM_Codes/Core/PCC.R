source('Core/GOR.R')
source('Core/RandomGraphOverlapRate.R')
PCC = function(A1, A2) {
  n = nrow(A1);
  A1[lower.tri(A1,diag=TRUE)] = 0;
  A2[lower.tri(A2,diag=TRUE)] = 0;
  m1 = sum(sum(A1));
  m2 = sum(sum(A2));
  PCC = ((m1 + m2)/(2*sqrt(m1*m2*(1 - 2*m1/n^2)*(1 - 2*m2/n^2))))*(GOR(A1, A2) - RandomGraphOverlapRate(A1, A2));
  return(PCC)
}