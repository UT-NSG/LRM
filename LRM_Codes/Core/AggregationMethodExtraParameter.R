source('Core/RandomGraphOverlapRate.R')
source('Core/PCC.R')
source('Core/GOR.R')
AggregationMethodExtraParameter = 
  function (aggregation_method, Train_Matrix, AUX_Matrix, EPType = extra_processing) {
    if(EPType == 'None') {
      if(aggregation_method == 'SimpleAddition') {
        AggregationMethodExtraParameter = 0;
      }
      else if(aggregation_method == 'PCC'){
        AggregationMethodExtraParameter = PCC(Train_Matrix, AUX_Matrix);
      }
      else if(aggregation_method == 'GOR'){
        AggregationMethodExtraParameter = GOR(Train_Matrix, AUX_Matrix);
      }
      else {
        # To add other types of extra processing later
        print(aggregation_method)
        print(' aggregation method is not support yet.')
        stop('Error: Error, undefined aggregation method parameter')
      }
    }
    else {
      # To add other types of extra processing later
      print(EPType)
      print(' extra processing is not support yet.')
      stop('Error: Error, undefined extra processing parameter')
    }
    return(AggregationMethodExtraParameter)
}