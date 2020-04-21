RandomGraphOverlapRate = function(A1, A2) {
  n = nrow(A1);
  A1[lower.tri(A1,diag=TRUE)] = 0;
  A2[lower.tri(A2,diag=TRUE)] = 0;
  m1 = sum(sum(A1));
  m2 = sum(sum(A2));
  RG = (1 / n^2) * (4 * m1 * m2) / (m1 + m2);
  return(RG)
}