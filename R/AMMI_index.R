#' AMMI Stability Function
#'
#' Performs AMMI stability analysis
#' @param data A numeric data frame with the columns in the following order: genotypes,
#' environments, replications and the response variable.
#'
#' @export
#' @importFrom agricolae AMMI index.AMMI


AMMI_index<-function(data){

  model<- with(data,AMMI(data[,2], data[,1], data[,3], data[,4]))
  ranking<-as.data.frame(cbind(levels(as.factor(data[,1])),round(index.AMMI(model)$means,3),round(index.AMMI(model)$ASV,3),index.AMMI(model)$rASV,
                               index.AMMI(model)$YSI,index.AMMI(model)$rYSI))

  colnames(ranking)<-c("gen","Mean","ASV","rASV","YSI","rYSI")
  return(ranking)
}
