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
  
  #CREANDO UNA MATRIZ A PARTIR DE KRZANOWSKI(1988) MISSING
  ######OJO QUE ANTES DE CORRER ESTE PROGRAMA TOCA CAMBIAR EL N?MERO DE COLUMNAS SI SE CAMBIA EL DATASET
  X # Imprime la matriz de datos
  is.matrix(X) # Pregunta si los datos que se tiene son tipo matriz
  nfilasX<-nrow(X) 
  ncolX<-ncol(X)
  nfilasX
  ncolX
  totalelementos<-nfilasX*ncolX
  
  # ###Escolha do posto ?timo para a matriz incompleta por meio da VC
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
  
  
  #######Creando un vector de las posiciones faltantes y reemplazando temporalmente los faltantes por la media de columna
  
  XMISSING=X
  indicamissing<-is.na(XMISSING)# indica la ocurrencia de faltantes en la matriz
  indicadora<-indicamissing*1
  totalfaltantes<-sum(indicadora)
  medcol<-t(matrix(colMeans(XMISSING, na.rm=TRUE))) # Calcula media por columnas sin tener en cuenta datos faltantes
  # Inicial<-impute.svd(XMISSING,k=f, maxiter=5)
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
        # XMISSING[indi1,indi2]=Inicial$x[indi1,indi2]
      } #FIN DEL if (indicadora[indi1,indi2] 
    }# FIN DEL for(indi2 in 1:ncolX)
  } # FIN DEL for (indi1 in 1:nfilasX)
  
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
      
      #print (i)
      #print (P)
      
      tj_T<-t(Xmissingestand[posicfaltantes[i,1],-posicfaltantes[i,2]])%*%P[-posicfaltantes[i,2],]%*%ginv(t(P[-posicfaltantes[i,2],])%*%P[-posicfaltantes[i,2],])
      x_ij<-tj_T%*%P[posicfaltantes[i,2],]
      completed_matrix[posicfaltantes[i,1],posicfaltantes[i,2]]<-medcol[1,posicfaltantes[i,2]]+(x_ij*desvest[1,posicfaltantes[i,2]])
      
      #print(i)
      #print (tj_T)	
      #print (P)
      #print (x_ij)
      
    } ###Fin del for (i in 1:nfilasX){ 
    
    ##stabilitycrit<-(sum((completed_matrix-XMISSING)**2))/(sum((XMISSING*indicadora)**2))
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
    
    #print (f)
    #print(stabilityiter)
    #print(stabilitycrit)
    
    ###print(completed_matrix)
    
    if(f!=1 & convergence.new>convergence.old){
      f<-f-1
      stabilitycrit<-1
      stabilityiter<-0
    } ### Fin de if(f!=1 & convergence.new>convergence.old){
    
    
  }### FIN de while (stabilitycrit>epsilon & stabilityiter<=500){
  
  list(EigenItera=stabilityiter, EigenConverg=stabilitycrit, EigenComp=f, EigenImputaciones=completed_matrix)
  
}