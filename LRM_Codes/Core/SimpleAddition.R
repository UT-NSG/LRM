SimpleAddition = 
  function (ScoreSoFar, NewPrediction, EPType = 'None', AggregationMethodExtraParameter=0) {
    if(EPType == 'None') {
      NewScore = ScoreSoFar + AggregationMethodExtraParameter * NewPrediction
    }
    else if(EPType == 'Laplacian') {
      NewScore = ScoreSoFar + (diag(diag(NewPrediction)) - NewPrediction)
    }
    else if(EPType == 'Laplacian_A') {
      NewScore = ScoreSoFar + NewPrediction
    }
    else if(EPType == 'Normalized_Laplacian') {
      NewScore = ScoreSoFar + (diag(diag(NewPrediction)) - NewPrediction)
    }
    else if(EPType == 'Normalized_Laplacian_A') {
      NewScore = ScoreSoFar + NewPrediction
    }
    else {
      # To add other types of extra processing later
      print(EPType)
      stop('Error: Error, undefined extra processing parameter')
    }
    return(NewScore)
}