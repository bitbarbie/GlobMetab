#' @title Get relationmatrix 0.
#' @description This function creates a relationmatrix without
#' any connection. It was written to validate the influence of 
#' different connections.
#' @param p a rownamed matrix containing the likelihoods (dimension CxM)
#' @author Sarah Scharfenberg
#' @return A matrix CxC containing 0.
#' @export
get_relationmatrix_zero<-function(p)
{
  all_in_p <- rownames(p)
  w <- matrix(data = 0, nrow = length(all_in_p), ncol = length(all_in_p))
  rownames(w)<-all_in_p
  colnames(w)<-all_in_p
  return(w)
}


#' @title Get relationmatrix from reaction pairs.
#' @description This function creates a relationmatrix from a list of 
#' reaction pairs.
#' @param p a rownamed matrix containing the likelihoods (dimension CxM)
#' @param rp a data.table containing reaction pairs in colum V1 and V2,
#' respectively.
#' @author Sarah Scharfenberg
#' @details A list of reaction pairs can be retrieved for example from
#' the kegg reactions database.
#' The ids of rp should match the ids of p. If not try to translate them
#' with ...
#' @return A sparse matrix CxC containing 1 if there is a reaction pair containing 
#' these compounds, 0 else. Connections are not directed and do not exist
#' for a compound with itself. 
#' @import data.table
#' @import Matrix
#' @export
get_relationmatrix_rp<-function(p, rp){
  all_in_p <- rownames(p)
  C<-length(all_in_p)  
  wRows <- vector(mode = "numeric", length = 0)
  wCols <- vector(mode = "numeric", length = 0)
  
  for(i in 1:nrow(rp)){
    c1 <- rp[i,V1]
    c2 <- rp[i,V2]
    if(c1 != c2){
      if(is.element(c1,all_in_p) && is.element(c2,all_in_p)){
        c1_n <- which(all_in_p==c1)
        c2_n <- which(all_in_p==c2)
        if(c1_n > c2_n){
          wRows <- c(wRows, c2_n)
          wCols <- c(wCols, c1_n)
        }else{
          wRows <- c(wRows, c1_n)
          wCols <- c(wCols, c2_n)
        }
      }
    }
  }
  wVals <- rep(1, times=length(wRows))
  w <- sparseMatrix(i = wRows, j = wCols, x = wVals, dims = c(C,C),
                    symmetric = TRUE, index1 = TRUE)
  return(w)
}


#'@title Calculate similarity of 2 smiles.
#'@description This function calculates the chosen distance between two
#'rcdk smiles.
#'@param fp1 fingerprint of compound 1
#'@param fp2 fingerprint of compound 2
#'@param dist the distance to calculate as parameter method defined 
#'in rcdk::distance
#'@author Sarah Scharfenberg
#'@return The distance of the two fringerprints.
#'@import rcdk
#'@export
similarity <- function(fp1, fp2, dist="tanimoto") {

  sim <- NULL
  if(!is.null(fp1) && !is.null(fp2)){
    sim <- rcdk::distance(fp1, fp2, method = dist)
  }
  return(sim)
}


#'@title Get relationmatrix based on chemical similarity.
#'@description This function creates a relationmatrix based on the chemical
#'similarity. 
#'@param p a matrix rownamed by CID
#'@param dist a character string specifying the distance that should be
#'calculated. Default is tanimoto.
#'@param mode a character string specifying the matrix value. If dist was
#'chosen the matrix contains the distances foreach pair. If relation was 
#'chosen the matrix contains 1 if the distance was greater than the given
#'threshold and 0 else.
#'@param threshold a double value specifying the threshold for a chemical
#'similarity representing a relation. (see mode)
#'@param folder_PubChemID a character string specifying the location of a
#'prepared local pubchem mirror
#'@param verbose switch detailed output on
#'@return A (sparse) matrix CxC containing the specified output. (see mode)
#' Similarities are not directed and do not exist for a compound with itself.
#'@import rpubchem
#'@import data.table
#'@import rcdk
#'@import Matrix
#'@export
get_relationmatrix_distance <- function(p, dist="tanimoto", 
                                        mode="relation",
                                        threshold=0.7,
                                        folder_PubChemID,verbose=FALSE){
  
  
  if(!is.element(mode, c("dist","relation"))){
    stop("Mode must be one of c(\"dist\",\"relation\")")
  }
  
  all_in_p <- rownames(p)
  C<-length(all_in_p)  
  wRows <- vector(mode = "numeric", length = 0)
  wCols <- vector(mode = "numeric", length = 0)
  wVals <- vector(mode = "numeric", length = 0)
  # omit scientific notation in filename
  options(scipen=999)
  
  pubs <- data.frame(sort(as.integer(all_in_p)))
  colnames(pubs)<-"cid"
  
  fingerprints <- list()
  for(id in pubs[,1]){
    if(verbose){ print(paste0("next id: ",id))}
    if(id%%100==0){
      end <- id 
    }else{
      end <- (id%/%100)*100+100
    }
    start <- end-100+1
    file <- paste0("part_",start,"_",end,".csv")
    if(verbose){ print(paste0("reading ",file))}
    tmp <- data.table::fread(paste0(folder_PubChemID,file))
    if(verbose){ print("succeded")}
    s <- tmp[which(tmp[,cid]==id),openeye_can_smiles]
    if(length(s)!=0){
      m <- parse.smiles(s)
      fp <- get.fingerprint(m[[1]])
      fingerprints[[as.character(id)]]<-fp   
    }
  }

  for(c1 in 1:(C-1)){
    fp1 <- fingerprints[[all_in_p[c1]]]
    print(paste0(c1,"/",C))

    for(c2 in (c1+1):C){
      fp2 <- fingerprints[[all_in_p[c2]]]

      if(!is.null(fp1) && !is.null(fp2)){
        sim <- similarity(fp1,fp2,dist)
        if(mode=="relation" && sim>threshold){
          if(c1 > c2){
            wRows <- c(wRows, c2)
            wCols <- c(wCols, c1)
            wVals <- c(wVals,1)
          }else{
            wRows <- c(wRows, c1)
            wCols <- c(wCols, c2)
            wVals <- c(wVals,1)
          }
          
        }
        if(mode=="dist"){
          
          if(c1 > c2){
            wRows <- c(wRows, c2)
            wCols <- c(wCols, c1)
            wVals <- c(wVals,sim)
          }else{
            wRows <- c(wRows, c1)
            wCols <- c(wCols, c2)
            wVals <- c(wVals,sim)
          }
        }
      }
    }
  }
  w <- sparseMatrix(i = wRows, j = wCols, x = wVals, dims = c(C,C),
                    symmetric = TRUE, index1 = TRUE)
  return(w)
}


#'@title Get relationmatrix based on chemical similarity.
#'@description This function creates a relationmatrix based on the chemical
#'similarity. 
#'@param p a matrix rownamed by CID
#'@param dist a character string specifying the distance that should be
#'calculated. Default is tanimoto.
#'@param mode a character string specifying the matrix value. If dist was
#'chosen the matrix contains the distances foreach pair. If relation was 
#'chosen the matrix contains 1 if the distance was greater than the given
#'threshold and 0 else.
#'@param threshold a double value specifying the threshold for a chemical
#'similarity representing a relation. (see mode)
#'@param folder_PubChemID a character string specifying the location of a
#'prepared local pubchem mirror
#'@param verbose switch detailed output on
#'@param comp_cores number of cores to use for parallel processing
#'@return A matrix CxC containing the specified output. (see mode)
#' Similarities are not directed and do not exist for a compound with itself.
#' Careful with big C and similarity output this
#' is no sparse matrix!
#'@import rpubchem
#'@import data.table
#'@import rcdk
#'@import parallel
#'@import doParallel
#'@import foreach
#'@export
get_relationmatrix_distance_par <- function(p, dist="tanimoto", 
                                        mode="relation",
                                        threshold=0.7,
                                        folder_PubChemID,verbose=FALSE,
                                        comp_cores=1){
 
  if(comp_cores > parallel::detectCores())
    stop("You assigned more cores to the comp_cores argument than are availible on your machine.")
  
  
  if(!is.element(mode, c("dist","relation"))){
    stop("Mode must be one of c(\"dist\",\"relation\")")
  }
  
  all_in_p <- rownames(p)
  C<-length(all_in_p)  

  # omit scientific notation in filename
  options(scipen=999)
  
  pubs <- data.frame(sort(as.integer(all_in_p)))
  colnames(pubs)<-"cid"
  
  fingerprints <- list()
  for(id in pubs[,1]){
    if(verbose){ print(paste0("next id: ",id))}
    if(id%%100==0){
      end <- id 
    }else{
      end <- (id%/%100)*100+100
    }
    start <- end-100+1
    file <- paste0("part_",start,"_",end,".csv")
    if(verbose){ print(paste0("reading ",file))}
    tmp <- data.table::fread(paste0(folder_PubChemID,file))
    if(verbose){ print("succeded")}
    s <- tmp[which(tmp[,cid]==id),openeye_can_smiles]
    if(length(s)!=0){
      m <- parse.smiles(s)
      fp <- get.fingerprint(m[[1]])
      fingerprints[[as.character(id)]]<-fp   
    }
  }
  
  cl <- makeCluster(comp_cores)
  registerDoParallel(cl)
  print("Cores registered")
  
  #########################################################
  par_res<- foreach::foreach( c1=1:(C-1), .errorhandling = "stop",
                    .verbose = verbose,
                    .init=NULL,
                    .combine = "rbind",
                    .packages = "rcdk") %dopar%{ 
  ##########################################################                    
# for(c1 in 1:(C-1)){

    wRows <- vector(mode = "numeric", length = 0)
    wCols <- vector(mode = "numeric", length = 0)
    wVals <- vector(mode = "numeric", length = 0)

    fp1 <- fingerprints[[all_in_p[c1]]]
    print(paste0(c1,"/",C))
    
    for(c2 in (c1+1):C){
      fp2 <- fingerprints[[all_in_p[c2]]]
      
      if(!is.null(fp1) && !is.null(fp2)){
        sim <- GlobMetab::similarity(fp1,fp2,dist)
        if(mode=="relation" && sim>threshold){
          if(c1 > c2){
            wRows <- c(wRows, c2)
            wCols <- c(wCols, c1)
            wVals <- c(wVals,1)
          }else{
            wRows <- c(wRows, c1)
            wCols <- c(wCols, c2)
            wVals <- c(wVals,1)
          }
          
        }
        if(mode=="dist"){
          
          if(c1 > c2){
            wRows <- c(wRows, c2)
            wCols <- c(wCols, c1)
            wVals <- c(wVals,sim)
          }else{
            wRows <- c(wRows, c1)
            wCols <- c(wCols, c2)
            wVals <- c(wVals,sim)
          }
        }
      }
    }
#  }
    return(cbind(wRows,wCols,wVals))
  } #end dopar

  stopCluster(cl)
  wRows <- par_res[,1]
  wCols <- par_res[,2]
  wVals <- par_res[,3]
  w <- sparseMatrix(i = wRows, j = wCols, x = wVals, dims = c(C,C),
                    symmetric = TRUE, index1 = TRUE)
  return(w)
}
