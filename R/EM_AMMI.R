#' EM-AMMI imputation
#'
#' @param X a data frame or matrix that contains the genotypes in the rows and the
#' environments in the columns when there are no replications of the experiment.
#' @param PC.nb the number of principal components in the AMMI model that will be
#' used; the default value is 1. For PC.nb=0 only main effects are used to estimate cells in the
#' data table (the interaction is ignored). The number of principal components must not be
#' greater than min(number of rows in the X table, number of columns in the X table)–2.
#' @param initial.values (optional) the initial values of the missing cells. It can be a single value,
#' which then will be used for all empty cells, or a vector of length equal to the number of
#' missing cells (starting from the missing values in the first column). If omitted, the initial
#' values will be obtained by the main effects from the corresponding model, that is, by the
#' grand mean of the observed data increased (or decreased) by row and column main effects.
#' @param precision (optional) the algorithm converges if the maximal change in the values of the
#' missing cells in two subsequent steps is not greater than this value (the default is 0.01);
#' @param max.iter (optional) a maximum permissible number of iterations (that is, number of
#' repeats of the algorithm’s steps 2 through 5); the default value is 1000;
#' @param change.factor (optional) introduced by analogy to step size in gradient descent method,
#' this parameter that can shorten the time of executing the algorithm by decreasing the
#' number of iterations. The change.factor=1 (default) defines that the previous
#' approximation is changed with the new values of missing cells (standard EM-AMMI algorithm). However, when change.factor<1, then the new approximations are computed
#' and the values of missing cells are changed in the direction of this new approximation but
#' the change is smaller. It could be useful if the changes are cyclic and thus convergence could
#' not be reached. Usually, this argument should not affect the final outcome (that is, the
#' imputed values) as compared to the default value of change.factor=1.
#' @param simplified.model the AMMI model contains the general mean, effects of rows, columns
#'and interaction terms. So the EM-AMMI algorithm in step 2 calculates the current effects of
#'rows and columns; these effects change from iteration to iteration because the empty (at the
#' outset) cells in each iteration are filled with different values. In step 3 EM-AMMI uses those
#'effects to re-estimate cells marked as missed (as default, simplified.model=FALSE). It is,
#'however, possible that this procedure will not converge. Thus the user is offered a simplified
#'EM-AMMI procedure that calculates the general mean and effects of rows and columns only
#'in the first iteration and in next iterations uses these values (simplified.model=TRUE). In
#'this simplified procedure the initial values affect the outcome (whilst EM-AMMI results
#' usually do not depend on initial values). For the simplified procedure the number of
#'iterations to convergence is usually smaller and, furthermore, convergence will be reached
#'even in some cases where the regular procedure fails. If the regular procedure does not
#'converge for the standard initial values (see the description of the argument
#'initial.values), the simplified model can be used to determine a better set of initial values.
#'@references Paderewski, J. (2013). \emph{An R function for imputation of missing cells in two-way
#'initial.values), the simplified model can be used to determine a better set of initial values.
#'data sets by EM-AMMI algorithm.}. Communications in Biometry and Crop Science 8 (2), 60–69.
#'@return A list of class \code{EM_AMMI} containing:
#'\itemize{
#'\item  X: the imputed matrix (filled in with the missing values estimated by the EM-AMMI procedure);
#'\item PC.SS: the sum of squares representing variation explained by the principal components (the squares of eigenvalues of singular value decomposition);
#'\item iteration: the final number of iterations;
#'\item precision.final: the maximum change of the estimated values for missing cells in the last step of iteration (the precision of convergence). If the algorithm converged, this value is slightly smaller than the argument precision;
#'\item PC.nb.final: a number of principal components that were eventually used by the EM.AMMI() function. The function checks if there are too many missing cells to unambiguously compute the parameters by the SVD decomposition (Gauch and Zobel, 1990). In that case the final number of principal components used is smaller than that which was passed on to the function through the argument PC.nb;
#'\item convergence: the value TRUE means that the algorithm converged in the last iteration.
#' }
#'@keywords internal
#'

EM.AMMI<-function(X, PC.nb=1, initial.values=NA, precision=0.01,
                  max.iter=1000, change.factor=1, simplified.model=FALSE)
{
  X<-as.matrix(X)
  X.missing<-matrix(1,nrow(X),ncol(X))
  X.missing[is.na(X)]<-0
  max.IPC=min(c(rowSums(X.missing),colSums(X.missing)))-1
  if (max.IPC<PC.nb)
  {PC.nb.used<-max.IPC} else {PC.nb.used<-PC.nb}
  X.ini<-X
  if (length(initial.values)==1)
  {
    if (!is.na(initial.values))
    {
      initial.values<-matrix(initial.values,nrow(X),ncol(X))
      X.ini[is.na(X.ini)]<-initial.values[is.na(X)]
    } else {
      X.mean<-mean(c(X),na.rm = TRUE)
      row.m<-matrix(rowMeans(X,na.rm = TRUE),nrow(X),ncol(X))
      col.m<-t(matrix(colMeans(X,na.rm = TRUE),ncol(X),nrow(X)))
      estimated<-(-X.mean)+row.m+col.m
      X.ini[is.na(X.ini)]<-estimated[is.na(X)]
    }
  } else {
    X.ini[is.na(X.ini)]<-initial.values[is.na(X)]
  }
  iteration<-1
  X.new<-X.ini
  change<-precision+1
  while ((change>precision)&(iteration<max.iter))
  {
    if (iteration==1 | !simplified.model)
    {
      x.mean<-mean(X.new)
      X.new.Ie<-X.new-x.mean
      X.new.Ie<-scale(X.new.Ie,center = TRUE, scale = FALSE)
      x.col.center<-attr(X.new.Ie,"scaled:center")
      X.new.Ie<-t(scale(t(X.new.Ie),center = TRUE, scale = FALSE))
      x.row.center<-attr(X.new.Ie,"scaled:center")
    } else { X.new.Ie<-X.new-x.mean-row.eff-col.eff }
    if (PC.nb.used>=1)
    {
      SVD <- La.svd(X.new.Ie)
      SVD$d<-(SVD$d[1:PC.nb.used])
      SVD$u<-SVD$u[,1:PC.nb.used]
      SVD$v<-SVD$v[1:PC.nb.used,]
      diag.l<-diag(SVD$d,nrow=PC.nb.used)
      interaction.adj<-SVD$u%*%diag.l%*%SVD$v
    } else interaction.adj<-0
    if (iteration==1 | !simplified.model)
    {
      row.eff<-matrix(x.row.center,nrow(X),ncol(X))
      col.eff<-t(matrix(x.col.center,ncol(X),nrow(X)))
    }
    X.next<-X.new
    X.next[is.na(X)]<-(x.mean+row.eff+col.eff+interaction.adj)[is.na(X)]
    change<-max(abs(c( (X.next-X.new)[is.na(X)] )))
    iteration<-iteration+1
    X.new<-change.factor*X.next+(1-change.factor)*X.new
  }
  if (change<=precision) {state=TRUE} else {state=FALSE}
  if (PC.nb.used<PC.nb) {state=FALSE}
  if (PC.nb.used==0) {SVD<-list(d=0)}
  return(list(X=X.new,iteration=iteration,precision.final=change,PC.nb.final=PC.nb.used, convergence=state))
}
