#' @title Validate prediction.
#' @description This function checks colwise equality of a two assignments z.
#' E.g. a predicted and the true or changes between different predictions.
#' @param true a matrix containing the correct assignment.
#' @param pred a matrix containing a predicted assignment.
#' @author Sarah Scharfenberg
#' @return A list storing the percentage of correct assigned and false
#' assigned compounds. 
#' @export
validate<-function(true, pred){
  correct<-0
  for(m in 1:ncol(true)){
    # compare columnwise
    if(sum(pred[,m]!=true[,m])==0){correct <- correct + 1}
  }
  all <- ncol(true)
  false <- all - correct
  
  results <- list()
  results[["correct in percent"]]<-correct/all*100.0
  results[["false in percent"]]<-false/all*100.0
  return(results)
}