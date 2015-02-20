#' @title Get true assignment.
#' @description For testing only. 
#' @param p a row- and colnamed matrix.
#' @author Sarah Scharfenberg
#' @return A matrix of the same size as p containing 1 if colname equals
#' rowname and 0 else.
#' @export
get_true<-function(p){
  z<-p
  cols <- colnames(p)
  rows <- rownames(p)
  for(m in 1:ncol(p)){
    for(c in 1:nrow(p)){
      
      if(p[c,m]!=0){
        if(rows[c] == cols[m]){
          z[c,m]<-1
        }else{
          z[c,m]<-0
        }
      }
      
    }
  }
  return(z)
}