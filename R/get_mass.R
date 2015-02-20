#' @title Get mass from kegg entry by kegg-id.
#' @description This function searches KEGG by the given ids and returns
#' the exact mass for each request in a named vector.
#' @param ids a vector of KEGG ids as character strings.
#' @author Sarah Scharfenberg
#' @return A named vector storing the exact mass foreach id.
#' @import KEGGREST
#' @export
get_keggmass_by_id <- function(ids){
  masses <- list()
  for(id in ids){
    masses[[id]] <-as.numeric(keggGet(id)[[1]]$EXACT_MASS)
  }
  return(unlist(masses))
}




#' @title Get mass from kegg entry by molecular formula.
#' @description This function searches KEGG by the given formula and returns
#' the exact mass for each request in a named vector.
#' @param formulas a vector of molecular formulas as character strings.
#' @author Sarah Scharfenberg
#' @return A named vector storing the exact mass foreach formula.
#' @import KEGGREST
#' @export
get_keggmass_by_formula <- function(formulas){
  masses <- list()
  for(form in formulas){
    id <- names(keggFind("compound",form,"formula"))
    if(length(id)>1){
      id <- id[1]
    }  
    masses[[form]]<-as.numeric(keggGet(id)[[1]]$EXACT_MASS)
  }
  return(unlist(masses))
}


