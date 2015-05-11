#' @title Validate prediction.
#' @description This function checks colwise equality of two assignments z.
#' E.g. a predicted and the true or changes between different predictions.
#' @param true a matrix containing the correct assignment.
#' @param pred a matrix containing a predicted assignment.
#' @author Sarah Scharfenberg
#' @return A list storing the percentage of correct assigned and false
#' assigned compounds. 
# #' @export
validate_matrix<-function(true, pred){
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


#' @title Validate prediction.
#' @description This function checks equality of a two assignments z 
#' reduced to the index of each assignment.
#' E.g. a predicted and the true or changes between different predictions.
#' @param true a vector containing the correct assignment.
#' @param pred a vector containing a predicted assignment.
#' @author Sarah Scharfenberg
#' @return A list storing the percentage of correct assigned and false
#' assigned compounds. 
# #' @export
validate_vector<-function(true, pred){
  
  det <-which(true!=0 & pred!=0)
  
  if(length(det) > 0 ){
    true_det <- true[det]
    pred_det <- pred[det]
    
    all <- length(true_det)
    
    correct<-sum(true_det == pred_det)
    false <- all - correct
  }else{
    correct<-0
    false <- 0
    all<-0
  }
  results <- list()
  if(all>0){
    results[["correct in percent"]]<-correct/all*100.0
    results[["false in percent"]]<-false/all*100.0
  }else{
    results[["correct in percent"]]<-0
    results[["false in percent"]]<-0
  }
  results[["true id lost in percent"]]<-(length(true)-all) / length(true) *100.0
  return(results)
}



#' @title Worst case ranking.
#' @description This function ranks a numeric (scoring) matrix colwise
#' assigning the worst case rank.
#' @param p a matrix containing colwise scoring lists to be ranked.
#' @param comp_cores number of cores to use for parallel processing
#' @details Decreasing order - high scores appear at low rank.
#' If e.g. the first two candidates of a list have the same score 
#' the rank is the worst possible, so in this case 2 for both. 
#' @author Sarah Scharfenberg
#' @return A matrix similar to input but storing the ranks. 
#' @import parallel
#' @import doParallel
#' @export
worst_case_ranking <- function(p, comp_cores=NULL){
  turnon<-FALSE
  if(!is.null(comp_cores)){
    if(comp_cores > detectCores()){comp_cores <- detectCores()}
    cl <- makeCluster(comp_cores)
    registerDoParallel(cl)
    turnon <- TRUE
  }
#########################################################
  ranking<- foreach::foreach( m=1:ncol(p), 
                              .errorhandling = "stop",
                              .init=NULL,
                              .combine = "cbind") %dopar%{ 
######################################################### 

#  ranking <- sapply(1:ncol(p),function(m){    
    p_col <- p[,m]
    names(p_col)<-rownames(p)
    p_col <- sort(p_col, decreasing = TRUE)
    order <- sapply(1:length(rownames(p)), function(x){which(names(p_col)==rownames(p)[x])})
    r_col <- vector("numeric",length(p_col))
    set <- vector("logical",length(p_col))
    names(r_col)<-names(p_col)
    r <- 0
    #calculate ranking
    for(c in 1:nrow(p)){
        if(!set[c]){
          entry <- p_col[c]
          red <- length(which(p_col==entry))
          r <- r+red
          r_col[which(p_col==entry)]=r
          set[which(p_col==entry)]=TRUE
        }
    }
    return(r_col[order])
  }
  if(turnon){stopCluster(cl)}
  rownames(ranking)<-rownames(p)
  colnames(ranking)<-colnames(p)
  return(ranking)
}


#' @title Average case ranking.
#' @description This function ranks a numeric (scoring) matrix colwise
#' assigning the average case rank.
#' @param p a matrix containing colwise scoring lists to be ranked.
#' @param comp_cores number of cores to use for parallel processing
#' @details Decreasing order - high scores appear at low rank. 
#' If e.g. the first two candidates of a list have the same score 
#' the rank is the mean of all involved scores, so in this case (1+2)/2 = 1.5.
#' @author Sarah Scharfenberg
#' @return A matrix similar to input but storing the ranks. 
#' @import parallel
#' @import doParallel
#' @export
average_case_ranking <- function(p, comp_cores=1){
  
  if(comp_cores > detectCores()){comp_cores <- detectCores()}
  cl <- makeCluster(comp_cores)
  registerDoParallel(cl)
#########################################################
  ranking<- foreach::foreach( m=1:ncol(p), 
                              .errorhandling = "stop",
                              .init=NULL,
                              .combine = "cbind") %dopar%{ 
######################################################### 
#  for(m in 1:ncol(p)){
    
    p_col <- p[,m]
    names(p_col)<-rownames(p)
    p_col <- sort(p_col, decreasing = TRUE)
    order <- sapply(1:length(rownames(p)), function(x){which(names(p_col)==rownames(p)[x])})
    r_col <- vector("numeric",length(p_col))
    set <- vector("logical",length(p_col))
    names(r_col)<-names(p_col)
    r <- 0
    #calculate ranking
    for(c in 1:nrow(p)){
        if(!set[c]){
          entry <- p_col[c]        
          no <- length(which(p_col==entry))
          red <- 1/no * sum((r+1):(r+no))
          r_col[which(p_col==entry)]=red
          set[which(p_col==entry)]=TRUE
          r <- r+no
        }
    }
    return(r_col[order])
  }
  stopCluster(cl)
  rownames(ranking)<-rownames(p)
  colnames(ranking)<-colnames(p)
  return(ranking)
}

#' @title Best case ranking.
#' @description This function ranks a numeric (scoring) matrix colwise
#' assigning the best case rank.
#' @param p a matrix containing colwise scoring lists to be ranked.
#' @param comp_cores number of cores to use for parallel processing
#' @details Decreasing order - high scores appear at low rank. 
#' If e.g. the first two candidates of a list have the same score 
#' the rank is the best possible, so in this case 1 for both.
#' @author Sarah Scharfenberg
#' @return A matrix similar to input but storing the ranks. 
#' @export
best_case_ranking <- function(p, comp_core){
  
  if(comp_cores > detectCores()){comp_cores <- detectCores()}
  cl <- makeCluster(comp_cores)
  registerDoParallel(cl)
#########################################################
  ranking<- foreach::foreach( m=1:ncol(p), 
                              .errorhandling = "stop",
                              .init=NULL,
                              .combine = "cbind") %dopar%{ 
######################################################### 
  #for(m in 1:ncol(p)){
    
    p_col <- p[,m]
    names(p_col)<-rownames(p)
    p_col <- sort(p_col, decreasing = TRUE)
    order <- sapply(1:length(rownames(p)), function(x){which(names(p_col)==rownames(p)[x])})
    r_col <- vector("numeric",length(p_col))
    set <- vector("logical",length(p_col))
    names(r_col)<-names(p_col)
    r <- 1
    #calculate ranking
    for(c in 1:nrow(p)){
        if(!set[c]){
          entry <- p_col[c]
          r_col[which(p_col==entry)]=r
          set[which(p_col==entry)]=TRUE
          no <- length(which(p_col==entry))
          r <- r+no
        }
    }
    return(r_col[order])
  }
  return(ranking)
}

#' @title Validation of two rankings.
#' @description This function validates a rank given a vector 
#' containing the true compund ids.
#' @param z_true a named numeric vector containg forach query (named) the
#' index of the correct compound regarding the matrix p
#' @param p a named matrix containing the intial ranking
#' @param p_new a named matrix containing a predicted ranking
#' @author Sarah Scharfenberg
#' @return A list storing
#' $'better rank' - number and percentage of queries having a better rank
#' $'worse rank' - number and percentage of queries having a worse rank
#' $'same rank' - number and percentage of queries having the same rank
#' $'true id lost' - number and percent of ids lost because true id unknown
#' @export
validate_rank<-function(z_true,p,p_new){
  
  if(ncol(p)!=ncol(p_new) || nrow(p)!=nrow(p_new)){
    print("p and p_new are not same size")
    return(NULL)
  }
  take <- z_true[z_true!=0]
  report <- diag(p[take,names(take)])-diag(p_new[take,names(take)])
  
  res <- vector("numeric",4)
  res[1] <- sum(report > 0)
  res[2] <- sum(report < 0)
  res[3] <- sum(report==0)
  res[4] <- length(z_true) - length(take)
#   for(m in 1:ncol(p)){
#     if(z_true[m]==0){
#       res[4]<-res[4]+1
#     }
#     else{
#       if(p[z_true[m],m]!=p_new[z_true[m],m]){
#         if(p[z_true[m],m]>p_new[z_true[m],m]){
#           res[1]<-res[1]+1
#         }
#         else{
#           res[2]<-res[2]+1
#         }
#       }
#       else{
#         res[3]<-res[3]+1
#       }
#     }
#   }
  all <- ncol(p)
  results <- list()
  results[["better rank"]]<- c(res[1],res[1]/all*100.0)
  results[["worse rank"]]<-c(res[2], res[2]/all*100.0)
  results[["same rank"]]<-c(res[3], res[3]/all*100.0)
  results[["true id lost"]]<-c(res[4],res[4]/all*100.0)
  return(results)
}

#' @title Calculate Scores from sample table.
#' @description This function calculates the mean posterior estimator
#' given a set of samples.
#' @param p named Likelihood matrix
#' @param class_table a matrix containing as much rows as p has columns
#' and a column for each sample storing a complete assignemtn of all
#' unknown variables
#' @param samples the number of sampels to consider for calculation
#' @param number of burn in samples to discard 
#' @param mode one of "mul" or "add" or nothing if the resulting score should
#' be multiplied with or added to p 
#' @param comp_cores number of cores to use for parallel processing
#' @author Sarah Scharfenberg
#' @return A matrix storing the new scores.
#' @import doParallel
#' @import parallel
#' @export
calculate_score <- function(p,class_table,samples,
                            burnin_samples=NULL,mode=NULL,
                            comp_cores=NULL){
  
  if(ncol(p)!=nrow(class_table)){
    stop(print("wrong dimension"))
  }
  if(is.null(burnin_samples)){
      start <- 1
      end <- samples
  }else{
      start <- 1+burnin_samples
      end <- burnin_samples + samples
  }
  
 # result <- matrix(0, ncol=ncol(p), nrow=nrow(p))
  turnon<-FALSE
  if(!is.null(comp_cores)){
    if(comp_cores > detectCores()){comp_cores <- detectCores()}
    cl <- makeCluster(comp_cores)
    registerDoParallel(cl)
    turnon <- TRUE
  }
#########################################################
  result<- foreach::foreach( m=1:ncol(p), 
                              .errorhandling = "stop",
                              .init=NULL,
                              .combine = "cbind") %dopar%{ 
######################################################### 
#  for(m in 1:ncol(p)){

#     counts <- as.data.frame(table(class_table[m,start:end]))
#     for(i in 1:nrow(counts)){
#       # entry in counts corresponds to the position in p
#       # as the rows are the same we can directly take this
#       res[as.integer(as.character(counts[i,"Var1"]))] <-as.integer(counts[i,"Freq"])
#     }
    tmp_ct <- as.integer(class_table[m,start:end])
    ct_uniq <- unique(tmp_ct)
    count_uniq <- sapply(1:length(ct_uniq),function(i){
      return(sum(tmp_ct == ct_uniq[i]))
    })
    res <- vector("numeric",nrow(p))
    res[ct_uniq]<-count_uniq
    # res is a column vector containing the occurence of ...
    # in range of class_table[m,start:end]
    return(res)
  }
  if(turnon){  stopCluster(cl) }
  rownames(result)<-rownames(p)
  colnames(result)<-colnames(p)
  result <- result/ncol(class_table)
  if(mode=="mul"){return(result*p)}
  if(mode=="add"){return(result+p)}
  return(result)
}


#' @title Get prediction table.
#' @description A function to get a human readable output of a prediction.
#' @param p named Likelihood matrix
#' @param ranking a character string specifying the rank to use.
#' until now only "worst case" possible.
#' @param max_rank the maximal rank that should be diplayed
#' @author Sarah Scharfenberg
#' @return A data.frame storing 'query', 'rank', 'score', 'candidate', 'cname'
#' foreach compound with rank <=max_rank.
#' @export
get_prediction_table <- function(p, ranking="worst case", max_rank=10){
  rank_table <- NULL
  if(ranking=="worst case"){rank_table <- worst_case_ranking(p)}
  pred <- data.frame()
  for(x in 1:ncol(p)){
    tmp <- rownames(rank_table)
    tmp_table <-data.frame()
    for(y in 1:max_rank){
      ids <- which(rank_table[,x]==y)
      if(length(ids>0)){
        cnames <- get.cid(tmp[ids])$IUPACName
        if(length(cnames)==0){cnames=NA}
        tmp_table <- rbind(tmp_table,data.frame(query=colnames(p)[x],
                                                rank=y, score=p[ids,x] ,candidate=tmp[ids],
                                                cname=cnames))
      }
    }
    pred<-rbind(pred,tmp_table)
  }
  rownames(pred)<-NULL
  return(pred)
}

