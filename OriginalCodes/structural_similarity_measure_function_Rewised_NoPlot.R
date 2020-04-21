structural_similarity_measure_function_Rewised_NoPlot <- 
  function(SS, PercentOfLeading = 0.1){

KK = nrow(SS)

K_ToShow = KK*PercentOfLeading
Acumulative_Measure = 0
for (i in 1:(K_ToShow)) {
  index_Max = (which(SS == max(SS)))
  index_MaxRC = (which(SS == max(SS), arr.ind = TRUE))
  index_Max_C = index_Max[1]
  # print(index_Max_C)
  Acumulative_Measure = Acumulative_Measure + SS[index_Max_C]
  if(i > KK-2){
    break;
  }
  SS = SS[-index_MaxRC[1,1],]
  SS = SS[,-index_MaxRC[1,2]]
}


return ((Acumulative_Measure)/floor(K_ToShow))

}


