#' AMMI Stability Function
#'
#' Performs AMMI stability analysis. Calculate AMMI stability value (ASV) and Yield stability index (YSI).
#'
#' @param data A numeric data frame with the columns in the following order: genotypes,
#' environments, replications and the response variable.
#'
#'
#' @return A data frame containing:\itemize{
#' \item{means}{average genotype by environment}
#'  \item{ASV}{AMMI stability value}
#'  \item{rASV}{Rank of AMMI stability value}
#'  \item{YSI}{Yield stability index}
#'  \item{rYSI}{Rank of yield stability index}
#'  }
#' @external
#' @export
#'
#' @importFrom agricolae AMMI index.AMMI


AMMI_index<-function(data){

  model<- with(data,AMMI(data[,2], data[,1], data[,3], data[,4]))
  ranking<-as.data.frame(cbind(levels(as.factor(data[,1])),round(index.AMMI(model)$means,3),round(index.AMMI(model)$ASV,3),index.AMMI(model)$rASV,
                               index.AMMI(model)$YSI,index.AMMI(model)$rYSI))

  colnames(ranking)<-c("gen","Mean","ASV","rASV","YSI","rYSI")
  return(ranking)
}
