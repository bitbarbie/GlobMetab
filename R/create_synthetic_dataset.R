#' @title Create a synthetic dataset based on KEGG IDs.
#' @description This method was written to create different semi-synthetic
#' dataset for validation purposes.
#' @param x is an vector of masses that can be named or a vector of KEGG ids.
#' @param input a character string defining the content of x. It can be 
#' either \"mass\" or \"id\".
#' @param mode a character string defining the values of the matrix. It can
#' be \"likelihood\" or \"value\". (see Details)
#' @param value a numeric value that should is set for each candidate
#' if mode was "value".
#' @author Sarah Scharfenberg
#' @details For each mass / id a compound search using keggFind is done and
#' each match is stored with a value. If likelihood was chosen a probability
#' is calculated using an error function taking the mean and standard 
#' deviation of the sample into account. If value was chosen each candidate
#' having a mass in a fix range gets the same value.
#' @return A CxM matrix containing a value for each candidate of each mass
#' and 0 else.
#' @import KEGGREST
#' @export
create_synthetic_dataset_kegg<-function(x, input=c("mass","id"), 
                                   mode=c("likelihood", "value"),
                                   value = 0.5)
{
  if(length(input) > 1 || !is.element(input, c("mass","id"))){
    stop("Input must be one of c(\"mass\",\"id\")")
  }
  
  if(length(mode) > 1 || !is.element(mode,c("likelihood", "value","mass"))){
    stop("Mode must be one of c(\"likelihood\", \"value\", \"mass\")")
  }
  
  M <-length(x)
  p<-matrix(0, ncol=M)
  m<-1
  delete_cols <- NULL
  
  massDB <- NULL
  massObs <- NULL
  
  if(input=="mass"){
  
    for(mass in x){
      # search kegg for candidates in range, sadly not 10ppm
      clist <- keggFind("compound", round(mass,6), "exact_mass")
      if(length(clist)>0){
        p0 <- matrix(0, ncol=M, nrow= length(clist))
        p0[,m] <- as.numeric(clist)
        
        if(mode=="likelihood"){
          massDB <- c(massDB,p0[,m])
          massObs <- c(massObs,rep(mass,times=length(clist)))
          p0[,m] <-abs(mass-p0[,m])/p0[,m]
        }
        
        if(mode=="value"){
          p0[,m]<-rep(value,times=length(clist))
        }
                
        # cut cpd:       
        rownames(p0)<-substr(names(clist),5,nchar(names(clist)))
        p <- rbind(p,p0)
      }else{
        delete_cols <- c(delete_cols,m)
      }
      m <- m+1
    }
    
    if(!is.null(names(x))){
      colnames(p)<-names(x)
    }
  }
  
  if(input=="id"){
    
    for(id in x){
      # search kegg for candidates in range, sadly not 10ppm
      mass <- as.numeric(keggGet(id)[[1]]$EXACT_MASS)
      clist <- keggFind("compound", round(mass,6), "exact_mass")
      if(length(clist)>0){
        p0 <- matrix(0, ncol=M, nrow= length(clist))
        p0[,m] <- as.numeric(clist)
        
        if(mode=="likelihood"){
          massDB <- c(massDB,p0[,m])
          massObs <- c(massObs,rep(mass,times=length(clist)))
          p0[,m] <-abs(mass-p0[,m])/p0[,m]
        }
        if(mode=="value"){
          p0[,m]<-rep(value,times=length(clist))
        }
        
        # cut cpd:       
        rownames(p0)<-substr(names(clist),5,nchar(names(clist)))
        p <- rbind(p,p0)
      }else{
        delete_cols <- c(delete_cols,m)
      }
      m <- m+1
    }
    colnames(p)<-x
  }
  
  p <- p[-1,]

  count <- 0
  for(c in delete_cols){
    p <- p[,-(c-count)]
    count <- count+1
  }
  
  if(mode == "likelihood"){
    centered <- abs(massObs-massDB)/massDB
    mean_centered <- mean(centered)
    sd_centered <- sd(centered)
    p <- as.numeric(p!=0) * 1 * pnorm(((p-mean_centered)/sd_centered) * sqrt(2), lower.tail = FALSE)
  }
  
  return(p)
}



#' @title Create a synthetic dataset based on PubChem IDs.
#' @description This method was written to create different semi-synthetic
#' dataset for validation purposes.
#' @param x is an vector of masses that can be named or a vector of 
#' PubChem ids.
#' @param input a character string defining the content of x. It can be 
#' either \"mass\" or \"id\".
#' @param mode a character string defining the values of the matrix. It can
#' be \"likelihood\" or \"value\" or \"mass\". (see Details)
#' @param value a numeric value that should is set for each candidate
#' if mode was "value".
#' @author Sarah Scharfenberg
#' @details For each mass / id a compound search using keggFind is done and
#' each match is stored with a value. If likelihood was chosen a probability
#' is calculated using an error function taking the mean and standard 
#' deviation of the sample into account. If value was chosen each candidate
#' having a mass in a fix range gets the same value. If mass was chosen 
#' each entry contains the mass of the candidate.
#' @return A CxM matrix containing a value for each candidate of each mass
#' and 0 else.
#' @import metfRag
#' @export
create_synthetic_dataset_pubchem<-function(x, input=c("mass","id"), 
                                        mode=c("likelihood", "value", "mass"),
                                        value = 0.5)
{
  
  if(length(input) > 1 || !is.element(input, c("mass","id"))){
    stop("Input must be one of c(\"mass\",\"id\")")
  }
  
  if(length(mode) > 1 || !is.element(mode,c("likelihood", "value","mass"))){
    stop("Mode must be one of c(\"likelihood\", \"value\", \"mass\")")
  }
  

  M <-length(x)
  p<-matrix(0, ncol=M)
  m<-1
  delete_cols <- NULL
  
  massDB <- NULL
  massObs <- NULL
  
  if(input=="mass"){
    
    for(mass in x){
      # search kegg for candidates in range, sadly not 10ppm
      params <- list(mass=mass, range=0.00005);
      clist <-db.pubchem.getId(params);

      if(length(clist)>0){
        p0 <- matrix(0, ncol=M, nrow= length(clist))
        
        # get masses
        i<-1
        for(id in clist){
  
          pubchem.container <- db.pubchem.getMoleculeContainer(id)
          this.mass <- as.numeric(get.properties(pubchem.container[[1]])$MONOISOTOPIC_WEIGHT)
          p0[i,m] <- this.mass
          i<-i+1
        }
              
        if(mode=="likelihood"){
          massDB <- c(massDB,p0[,m])
          massObs <- c(massObs,rep(mass,times=length(clist)))
          p0[,m] <-abs(mass-p0[,m])/p0[,m]
        }
        
        if(mode=="value"){
          p0[,m]<-rep(value,times=length(clist))
        }
        
        rownames(p0)<-unlist(clist)
        p <- rbind(p,p0)
      }else{
        delete_cols <- c(delete_cols,m)
      }
      m <- m+1
    }
    
    if(!is.null(names(x))){
      colnames(p)<-names(x)
    }
  }
#   
#   if(input=="id"){
#     
#     for(id in x){
#       # search kegg for candidates in range, sadly not 10ppm
#       mass <- as.numeric(keggGet(id)[[1]]$EXACT_MASS)
#       clist <- keggFind("compound", round(mass,6), "exact_mass")
#       if(length(clist)>0){
#         p0 <- matrix(0, ncol=M, nrow= length(clist))
#         p0[,m] <- as.numeric(clist)
#         
#         if(mode=="likelihood"){
#           massDB <- c(massDB,p0[,m])
#           massObs <- c(massObs,rep(mass,times=length(clist)))
#           p0[,m] <-abs(mass-p0[,m])/p0[,m]
#         }
#         if(mode=="value"){
#           p0[,m]<-rep(value,times=length(clist))
#         }
#         
#         # cut cpd:       
#         rownames(p0)<-substr(names(clist),5,nchar(names(clist)))
#         p <- rbind(p,p0)
#       }else{
#         delete_cols <- c(delete_cols,m)
#       }
#       m <- m+1
#     }
#     colnames(p)<-x
#   }
  
  p <- p[-1,]
  
  count <- 0
  for(c in delete_cols){
    p <- p[,-(c-count)]
    count <- count+1
  }
  
  if(mode == "likelihood"){
    centered <- abs(massObs-massDB)/massDB
    mean_centered <- mean(centered)
    sd_centered <- sd(centered)
    p <- as.numeric(p!=0) * 1 * pnorm(((p-mean_centered)/sd_centered) * sqrt(2), lower.tail = FALSE)
  }
  
  return(p)
}

