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
#' @return A matrix CxC containing 1 if there is a reaction pair containing 
#' these compounds, 0 else. Connections are not directed and do not exist
#' for a compound with itself.
#' @import data.table
#' @export
get_relationmatrix_rp<-function(p, rp){
  all_in_p <- rownames(p)
  w <- matrix(data = 0, nrow = length(all_in_p), ncol = length(all_in_p))
  rownames(w)<-all_in_p
  colnames(w)<-all_in_p
  
  for(i in 1:nrow(rp)){
    c1 <- rp[i,V1]
    c2 <- rp[i,V2]
    if(c1 != c2){
      if(is.element(c1,all_in_p) && is.element(c2,all_in_p)){
        w[which(rownames(w)==c1),which(rownames(w)==c2)]<-1
        w[which(rownames(w)==c2),which(rownames(w)==c1)]<-1
      }
    }
  }
  return(w)
}


#'@title Calculate similarity of 2 smiles.
#'@description This function calculates the chosen distance between two
#'smiles strings.
#'@param s1 smiles string of compound 1
#'@param s2 smiles string of compound 2
#'@param dist the distance to calculate as parameter method defined 
#'in rcdk::distance
#'@author Sarah Scharfenberg
#'@return The distance of the two smiles.
#'@import rcdk
#'@export
similarity <- function(s1, s2, dist="tanimoto") {
  m1 <- parse.smiles(s1)
  m2 <- parse.smiles(s2)
  
  fp1 <- get.fingerprint(m1[[1]])
  fp2 <- get.fingerprint(m2[[1]])
  
  return(distance(fp1, fp2, method = dist))
}


#'@title Get relationmatrix based on chemical similarity.
#'@description This function creates a relationmatrix based on the chemical
#'similarity. 
#'@param p a matrix rownamed by CID
#'@param distance a character string specifying the distance that should be
#'calculated. Default is tanimoto.
#'@param mode a character string specifying the matrix value. If dist was
#'chosen the matrix contains the distances foreach pair. If relation was 
#'chosen the matrix contains 1 if the distance was greater than the given
#'threshold and 0 else.
#'@param threshold a double value specifying the threshold for a chemical
#'similarity representing a relation. (see mode)
#' @return A matrix CxC containing the specified output. (see mode)
#' Similarities are not directed and do not exist for a compound with itself.
#'@import rpubchem
#'@export
get_relationmatrix_distance <- function(p, distance="tanimoto", 
                                        mode="relation",
                                        threshold=0.7){
  
  
  if(!is.element(mode, c("dist","relation"))){
    stop("Mode must be one of c(\"dist\",\"relation\")")
  }
  
  all_in_p <- rownames(p)
  C<-length(all_in_p)
  w <- matrix(data = 0, nrow = C, ncol = C)
  rownames(w)<-all_in_p
  colnames(w)<-all_in_p
  
  for(c1 in 1:C){
    s1 <- ""
    hits <- get.cid(all_in_p[c1])
    # get smile of unqiue matches, we assume that if its matches it is unique
    if(nrow(hits)==1){
          s1 <- as.character(hits["CanonicalSmiles"])
    }
    for(c2 in (c1+1):C){
      s2 <- ""
      hits <- get.cid(all_in_p[c2])
      if(nrow(hits)==1){
        s2 <- as.character(hits["CanonicalSmiles"])
      }
      if(s1 != "" && s2 != ""){
        sim <- similarity(s1,s2,distance)
        if(mode=="relation" && sim>threshold){
          w[c1,c2]<-1
          w[c2,c1]<-1
        }
        if(mode=="dist"){
          w[c1,c2]<-sim
          w[c2,c1]<-sim
        }
      }
    }
  }
  return(w)
}
