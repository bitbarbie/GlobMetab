#'@title A function to read list of sdf file into a matrix p
#'@param path a character string specifying the absolute path of
#'the folder containing the sdf files to be read.
#'@param id a character string identifying the property the 
#'compounds should identified by eg DatabaseID, Identifyer or
#'InChIKey1
#'@param add a numeric value that should be added to each input score
#'to distinguish a scroe of 0 from a candidate not appaering in the
#'candidate list at all.
#'@param filter_by_formula a logical indicating whether the canidate
#'lists should be cut to only those having the same molecular formula
#'@param formulas a named character vector specifying a formula for
#'each filename in the folder. If not provided when filter_by_formula
#'is TRUE then the corresponding list is not filtered.
#'@param filter_by_peaks a logical indicating whether the candidate
#'lists should be cut to only those explaining at least a given
#'number of peaks.
#'@param npe a character string defining the sdf property where the
#'number of explained peaks can be found
#'@param NumberPeaksExplained the number of peaks that should be at least
#'explained by the MetFrag-candidates peakspectrum
#'@param formula_by_SMILES a logical indicating whether the molecular 
#'formula for the candidates should be calculated using the SMILES 
#'property rather than the MolecularFormula, e.g. if MolecularFormula
#'is not provided.
#'@import iterators
#'@import rcdk
#'@export
read.sdf <- function(path, id = "DatabaseID",add=0.01, 
                     filter_by_formula=FALSE, formulas =NULL,
                     filter_by_peaks=FALSE,  npe = NULL, 
                     NumberPeaksExplained=NULL,
                     formula_by_SMILES=FALSE){
  
  folder <- list.files(path, all.files=FALSE, full.names = TRUE)
  
  # each query provides one sdf input file
  M <- length(folder)
  p <- matrix(0,ncol=M)
  current <- 0
  tmp_rownames <- vector("character")
  tmp_colnames <- vector("character")
  for(file in folder){
    
    base <- tail(strsplit(x = file,split='/')[[1]],n=1)
        
    moliter <- iload.molecules(file, type="sdf")
    
    while(hasNext(moliter)) {
      mol <- nextElem(moliter)
      
      if(filter_by_peaks){
        
        if(is.null(npe)||is.null(NumberPeaksExplained)){
          stop("npe and NumberPeaksExplained must be provided when filter_by_peaks is set.")
        }
        if(rcdk::get.property(mol,npe)>=NumberPeaksExplained){
          
          if(filter_by_formula ){
            if(is.element(base,names(formulas))){
              
              if(formula_by_SMILES){
                smiles <- rcdk::get.property(mol,"SMILES")
                molecule <- rcdk::parse.smiles(smiles)[[1]]
                formula <- rcdk::get.mol2formula(molecule,charge=0)@string
              }else{
                formula <- rcdk::get.property(mol,"MolecularFormula")
              }
              if( formula == formulas[base]){
                p <- rbind(p,c(rep(0,times=current),(as.numeric(rcdk::get.property(mol,"Score"))+add),rep(0, times=M-current-1)))
                tmp_rownames <- c(tmp_rownames, rcdk::get.property(mol,id))
              }
              
            }else{
              print(paste0("No formula provided for ",base,". Return unfiltered list here."))
            }
          }else{
            p <- rbind(p,c(rep(0,times=current),(as.numeric(rcdk::get.property(mol,"Score"))+add),rep(0, times=M-current-1)))
            tmp_rownames <- c(tmp_rownames, rcdk::get.property(mol,id))
          }
        }
        
      }else{ #filter not by peaks
        
        if(filter_by_formula ){
          if(is.element(base,names(formulas))){
            
            if(formula_by_SMILES){
              smiles <- rcdk::get.property(mol,"SMILES")
              molecule <- rcdk::parse.smiles(smiles)[[1]]
              formula <- rcdk::get.mol2formula(molecule,charge=0)@string
            }else{
              formula <- rcdk::get.property(mol,"MolecularFormula")
            }
            if( formula == formulas[base]){
              p <- rbind(p,c(rep(0,times=current),(as.numeric(rcdk::get.property(mol,"Score"))+add),rep(0, times=M-current-1)))
              tmp_rownames <- c(tmp_rownames, rcdk::get.property(mol,id))
            }
            
          }else{
            print(paste0("No formula provided for ",base,". Return unfiltered list here."))
          }
        }else{
          p <- rbind(p,c(rep(0,times=current),(as.numeric(rcdk::get.property(mol,"Score"))+add),rep(0, times=M-current-1)))
          tmp_rownames <- c(tmp_rownames, rcdk::get.property(mol,id))
        }
        
      }
    }
    current <- current+1
    
    # get colname from file name
    splits <- strsplit(file,"/|\\.")[[1]]
    tmp_colnames <- c(tmp_colnames,splits[length(splits)-1])
  }  
  p <- p[-1,]
  rownames(p) <- tmp_rownames
  colnames(p) <- tmp_colnames
  remove <- vector("numeric")
  for(i in 1:ncol(p)){
    if(sum(p[,i])==0){
      remove <- c(remove,i)
      print(paste0("remove ",colnames(p)[i]))
    }
  }
  if(length(remove)>0){
    p <- p[,-remove]
  }
  return(p)
}

#
# p is NOT sparse 
#
# read.sdf_sparse <- function(path, id = "DatabaseID",add=0.01, 
#                      filter_by_formula=FALSE, formulas =NULL,
#                      filter_by_peaks=FALSE,  npe = NULL, 
#                      NumberPeaksExplained=NULL,
#                      formula_by_SMILES=FALSE){
#   
#   folder <- list.files(path, all.files=FALSE, full.names = TRUE)
#   
#   # each query provides one sdf input file
#   M <- length(folder)
#   
#   #  p <- matrix(0,ncol=M)
#   pRows <- vector()
#   pCols <- vector()
#   pVals <- vector()
# 
#   current_col <- 1
#   current_row <- 1
#   set <- FALSE
#   
#   tmp_rownames <- vector("character")
#   tmp_colnames <- vector("character")
#   for(file in folder){
#     
#     base <- tail(strsplit(x = file,split='/')[[1]],n=1)
#     
#     moliter <- iload.molecules(file, type="sdf")
#     
#     while(hasNext(moliter)) {
#       mol <- nextElem(moliter)
#       seize <- FALSE
#       if(filter_by_peaks){
#         
#         if(is.null(npe)||is.null(NumberPeaksExplained)){
#           stop("npe and NumberPeaksExplained must be provided when filter_by_peaks is set.")
#         }
#         if(rcdk::get.property(mol,npe)>=NumberPeaksExplained){
#           
#           if(filter_by_formula ){
#             if(is.element(base,names(formulas))){
#               
#               if(formula_by_SMILES){
#                 smiles <- rcdk::get.property(mol,"SMILES")
#                 molecule <- rcdk::parse.smiles(smiles)[[1]]
#                 formula <- rcdk::get.mol2formula(molecule,charge=0)@string
#               }else{
#                 formula <- rcdk::get.property(mol,"MolecularFormula")
#               }
#               if( formula == formulas[base]){
#                 seize <- TRUE
#               }
#               
#             }else{
#               print(paste0("No formula provided for ",base,". Return unfiltered list here."))
#             }
#           }else{
#             seize <- TRUE
#           }
#         }
#         
#       }else{ #filter not by peaks
#         
#         if(filter_by_formula ){
#           if(is.element(base,names(formulas))){
#             
#             if(formula_by_SMILES){
#               smiles <- rcdk::get.property(mol,"SMILES")
#               molecule <- rcdk::parse.smiles(smiles)[[1]]
#               formula <- rcdk::get.mol2formula(molecule,charge=0)@string
#             }else{
#               formula <- rcdk::get.property(mol,"MolecularFormula")
#             }
#             if( formula == formulas[base]){
#               seize <- TRUE
#             }
#             
#           }else{
#             print(paste0("No formula provided for ",base,". Return unfiltered list here."))
#           }
#         }else{
#           seize <- TRUE
#         }
#         
#       }
#       
#       if(seize){
#         # p <- rbind(p,c(rep(0,times=current),(as.numeric(rcdk::get.property(mol,"Score"))+add),rep(0, times=M-current-1)))
#         rname <- rcdk::get.property(mol,id)
#         if(is.element(rname,tmp_rownames)){
#           pRows <- c(pRows,which(tmp_rownames==rname))
#         }else{
#           pRows <- c(pRows,current_row)
#           current_row <- current_row +1
#         }    
#         pCols <- c(pCols,current_col)
#         pVals <- c(pVals,as.numeric(rcdk::get.property(mol,"Score"))+add)
#         tmp_rownames <- c(tmp_rownames, rname)
#         set <- TRUE
#       }
#     }
#     if(set){
#       current_col <- current_col+1
#       # get colname from file name
#       splits <- strsplit(file,"/|\\.")[[1]]
#       tmp_colnames <- c(tmp_colnames,splits[length(splits)-1])
#       set <- FALSE
#     }
#   }  
#   #p <- p[-1,]
#   p <- sparseMatrix(i = pRows, j = pCols, x = pVals,
#                     symmetric = FALSE, index1 = TRUE)
#   rownames(p) <- tmp_rownames
#   colnames(p) <- tmp_colnames
# #  remove <- vector("numeric")
# #   for(i in 1:ncol(p)){
# #     if(sum(p[,i])==0){
# #       remove <- c(remove,i)
# #       print(paste0("remove ",colnames(p)[i]))
# #     }
# #   }
# #   if(length(remove)>0){
# #     p <- p[,-remove]
# #   }
#   return(p)
# }

#'@title Eliminate redundancies.
#'@description This function matchs redundand rows into one row.
#'@param p a matrix
#'@param log a character String specifying the path to a pubchem database
#'log file, containing information about the latest update, omitting
#'ids greater that the last id specified there
#'@return the reduced matrix.
#'@export  
eliminate_redundancies <- function(p, log = NULL){
  rows <- rownames(p)
  uniques <- unique(rows)
  res <- matrix(0.0, nrow = length(uniques),ncol=ncol(p))
  rownames(res)<-uniques
  colnames(res)<-colnames(p)
  for(row in 1:nrow(p)){
    index <- which(p[row,]>0)
    res[rows[row],index]<-p[row,index]
  }
  
  if(!is.null(log)){
    # omit not sortable ids
    #    if(!is.na(as.integer(rownames(p)[1]))){
    last_id <- as.integer(scan(log,what=character(), skip=1, 
                               nlines=1)[3])
    remove <- which(as.integer(rownames(res))>last_id)
    res <- res[-remove,]
    #    }
  }
  return(res)
}
