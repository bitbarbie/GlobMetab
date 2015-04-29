#'@title Recieve PubChem ID.
#'@description This function takes a list of KEGG Ids and tries to
#'convert them into PubChem IDs using the package KEGGREST. 
#'@param keggIDs a vector of character strings containing the KEGG ids
#'@details As in KEGG only the SID is stored the match to CID is done
#'using the additional package rpubchem.
#'Access via internet. Be careful, long input lists might overacces
#'the KEGG server.
#'@return A matrix with header keggID, SID, CID with one row foreach
#'input id storing the matched PubChem ids.
#'@author Sarah Scharfenberg
#'@import KEGGREST
#'@import rpubchem
#'@export
get_pubchemID_from_keggID <- function(keggIDs){
  
  res <- matrix("", nrow=length(keggIDs), ncol=3)
  colnames(res)<-c("keggID","SID","CID")
  res[,1]<-keggIDs
  
  k<-1
  K <- length(keggIDs)
  for(id in keggIDs){
    print (paste(k,"/",K))
    k<-k+1
    # keggGet to get SID
    links <- keggGet(id)[[1]]$DBLINKS
    for(link in links){
      if(grepl("PubChem", link)){
        res[which(res[,1]==id),2] <- strsplit(link,"PubChem: ")[[1]][2]
      }
    }
    
    #rpubchem get.cid.list 
    tmp <- unlist(rpubchem::get.cid.list(res[which(res[,1]==id),2])["CID"])
    tmp <- paste(tmp, collapse = ',')
    res[which(res[,1]==id),3] <- tmp
  }

return(res)
}