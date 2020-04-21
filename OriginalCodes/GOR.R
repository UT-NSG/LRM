GOR = function(A1, A2) {
  A1[lower.tri(A1,diag=TRUE)] = 0;
  A2[lower.tri(A2,diag=TRUE)] = 0;
  overlap_mat = A1 * A2;
  overlaps = sum(sum(overlap_mat));
  m1 = sum(sum(A1));
  m2 = sum(sum(A2));
  GOR = (2 * overlaps) / (m1 + m2);
  return(GOR)
}