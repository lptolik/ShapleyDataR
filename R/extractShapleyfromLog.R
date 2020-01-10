#' Extract Shapley data from temporary savings
#'
#' @param fname temporary save file to extract data from
#'
#' @return  data.frame having the same column names as properties of shapley class.
#' @export

extractShapley<-function(fname){
  e1 <- new.env(parent = baseenv())
load(fname,envir = e1)
phiLast<-e1$phi[[length(e1$phi)]]
idx<-match(e1$reId,e1$mod@react_id)
phiDF<-data.frame(id=e1$reId,
                  name=e1$mod@react_name[idx],
                  shapley=-1*phiLast,
                  sd=e1$sd,
                  err01=e1$e[,1],
                  err05=e1$e[,2],
                  err10=e1$e[,3])
return(phiDF)
}