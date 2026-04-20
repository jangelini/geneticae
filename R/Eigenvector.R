#' Eigenvector Imputation Function (Internal)
#' 
#' Internal function for GxE imputation using the Krzanowski (1988) eigenvector 
#' approach with a leave-one-out strategy.
#'
#' @param X A matrix with missing values (NAs).
#' @param f Number of components (rank) to use for the reconstruction.
#' @return A list with the number of iterations, convergence status, 
#' final rank used, and the imputed matrix.
#' @importFrom MASS ginv
#' @importFrom corpcor fast.svd
#' @importFrom stats var
#' @keywords internal

Eigenvectorfun <- function(X, f) {
  
  X 
  is.matrix(X) 
  nfilasX<-nrow(X) 
  ncolX<-ncol(X)
  nfilasX
  ncolX
  totalelementos<-nfilasX*ncolX
  
  # Mode <- function(x) {
  #   ux <- unique(x)
  #   ux[which.max(tabulate(match(x, ux)))]
  # }
  # bestrank=numeric(1000)
  # for(i in 1:1000){
  #   koptimalmatrix <- summary(cv.svd.wold(X))
  #   bestrank[i] <- koptimalmatrix$rank.best
  # }
  # #table(bestrank)
  # f<-Mode(bestrank)
  if (f==ncol(X)) {f<-f-1}
  
  

  XMISSING=X
  indicamissing<-is.na(XMISSING)
  indicadora<-indicamissing*1
  totalfaltantes<-sum(indicadora)
  medcol<-t(matrix(colMeans(XMISSING, na.rm=TRUE))) 

  incompleta<-XMISSING
  XMISSING
  posicfaltantes<-matrix(0,totalfaltantes,2)
  indi3<-0
  for (indi1 in 1:nfilasX){
    for(indi2 in 1:ncolX){
      if (indicadora[indi1,indi2]==1){
        indi3<-indi3+1
        posicfaltantes[indi3,1]=indi1
        posicfaltantes[indi3,2]=indi2
        XMISSING[indi1,indi2]=medcol[1,indi2]
        
      } 
    }
  } 
  
  posicfaltantes
  XMISSING
  completed_matrix<-XMISSING
  maxiter<-50
  iter<-0
  desvest<-sqrt(t(diag(var(XMISSING))))
  #print(desvest)
  Xmissingestand<-scale(XMISSING)
  epsilon<-1*(10**(-2))
  stabilitycrit<-1
  stabilityiter<-0
  
  while (stabilitycrit>epsilon & stabilityiter<=500){
    
    for (i in 1:totalfaltantes){
      
      oneout<-Xmissingestand[-posicfaltantes[i,1],]
      dvs_oneout<-fast.svd(oneout)
      
      if(f==1){
        u<-dvs_oneout$u[,f]
        d<-diag(dvs_oneout$d)
        v<-dvs_oneout$v[,f]
        T<-matrix(u*d[f,f])
        P<-matrix(v)
      }
      
      else{
        u<-dvs_oneout$u[,1:f]
        d<-diag(dvs_oneout$d)
        v<-dvs_oneout$v[,1:f]
        T<-u%*%d[1:f,1:f]
        P<-v
      }
      
      tj_T<-t(Xmissingestand[posicfaltantes[i,1],-posicfaltantes[i,2]])%*%P[-posicfaltantes[i,2],]%*%ginv(t(P[-posicfaltantes[i,2],])%*%P[-posicfaltantes[i,2],])
      x_ij<-tj_T%*%P[posicfaltantes[i,2],]
      completed_matrix[posicfaltantes[i,1],posicfaltantes[i,2]]<-medcol[1,posicfaltantes[i,2]]+(x_ij*desvest[1,posicfaltantes[i,2]])

      
    } 
    
    convergence.old<-stabilitycrit
    stabilitycrit1<-sqrt(  (1/totalfaltantes)*(sum((completed_matrix-XMISSING)**2))  )
    stabilitycrit2<-sqrt( (1/(totalelementos-totalfaltantes))*(sum((XMISSING*(1-indicadora))**2))  )
    stabilitycrit<-stabilitycrit1/stabilitycrit2
    stabilityiter<-stabilityiter+1
    XMISSING<-completed_matrix
    Xmissingestand<-scale(completed_matrix)
    medcol<-t(matrix(colMeans(completed_matrix)))
    desvest<-sqrt(t(diag(var(completed_matrix))))
    convergence.new<-stabilitycrit

    
    if(f!=1 & convergence.new>convergence.old){
      f<-f-1
      stabilitycrit<-1
      stabilityiter<-0
    }
    
    
  }
  
  list(EigenItera=stabilityiter, EigenConverg=stabilitycrit, EigenComp=f, EigenImputaciones=completed_matrix)
  
}