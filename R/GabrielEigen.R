#' Gabriel Eigen imputation
#'
#' @param X a data frame or matrix that contains the genotypes in the rows and
#'   the environments in the columns when there are no replications of the
#'   experiment.
#' @references Arciniegas-Alarcón S., García-Peña M., Dias C.T.S., Krzanowski
#'   W.J.. (2010). \emph{An alternative methodology for imputing missing data in
#'   trials with genotype-by-environment interaction}. Biometrical Letters, 47,
#'   1–14.
#' @return Imputed data with Gabriel Eigein method
#' @keywords internal
#' @importFrom MASS ginv
#'


Gabriel.Calinski<-function(X){

  if (missing(X)) stop("Need to provide X data frame or matrix")
  stopifnot(
  class(X) %in% c("matrix", "data.frame")
  )


  nfilasX<-nrow(X)
  ncolX<-ncol(X)
  #nfilasX
  #ncolX
  if(ncolX>nfilasX) {X<-t(X)}
  ######CREANDO UNA MATRIZ A PARTIR DE KRZANOWSKI(1988)
  ######OJO QUE ANTES DE CORRER ESTE PROGRAMA TOCA CAMBIAR EL N?MERO DE COLUMNAS SI SE CAMBIA EL DATASET
  XMISSING<-X
  is.matrix(X) # Pregunta si los datos que se tiene son tipo matriz
  totalelementos<-nfilasX*ncolX

  ######CREACI?N DE MATRIZ INDICADORA, C?LCULO DEL TOTAL DE DATOS FALTANTES Y C?LCULO INICIAL DE LA MEDIA A PARTIR DE MATRIZ INCOMPLETA
  indicamissing<-is.na(XMISSING)# indica la ocurrencia de faltantes en la matriz
  indicadora<-indicamissing*1
  totalfaltantes<-sum(indicadora)
  medcol<-t(matrix(colMeans(XMISSING, na.rm=TRUE))) # Calcula media por columnas sin tener en cuenta datos faltantes
  #MEDIA_inicial<-mean(XMISSING, na.rm=TRUE)

  ######CREACI?N DE MATRIZ CON LAS POSICIONES DE LOS DATOS FALTANTES
  ######REEMPLAZANDO DATOS FALTANTES CON LA MEDIA COLUMNA INICIAL
  posicfaltantes<-matrix(0,totalfaltantes,2)
  indi3<-0
  for (indi1 in 1:nfilasX){
    for(indi2 in 1:ncolX){
      if (indicadora[indi1,indi2]==1){
        indi3<-indi3+1
        posicfaltantes[indi3,1]=indi1
        posicfaltantes[indi3,2]=indi2
        XMISSING[indi1,indi2]=medcol[1,indi2]
        #XMISSING[indi1,indi2]=MEDIA_inicial
      } #FIN DEL if (indicadora[indi1,indi2]
    }# FIN DEL for(indi2 in 1:ncolX)
  } # FIN DEL for (indi1 in 1:nfilasX)
  #posicfaltantes
  #XMISSING

  ######DECLARACI?N DE ELEMENTOS PREVIOS PARA EL ESQUEMA ITERATIVO
  completed_matrix<-XMISSING
  maxiter<-1000
  iter<-0
  desvest<-sqrt(t(diag(var(XMISSING))))
  #print(desvest)
  Xmissingestand<-scale(XMISSING)
  Prediction_Gabriel<-matrix(0,nfilasX,ncolX)
  epsilon<-1*(10**(-6))
  stabilitycrit<-1
  stabilityiter<-1


  ######ESQUEMA ITERATIVO SEG?N ARCINIEGAS ET AL (2010)

  while (stabilitycrit>epsilon & stabilityiter<=100){


    for (i1 in 1:totalfaltantes){

      posic_linha<-posicfaltantes[i1,1]
      posic_col<-posicfaltantes[i1,2]
      auxiliar1<-rbind(Xmissingestand[posic_linha,],Xmissingestand[-posic_linha,])
      auxiliar2<-cbind(auxiliar1[,posic_col],auxiliar1[,-posic_col])
      #print(i1)
      #print(auxiliar2)
      Y1puntoTRANS<-auxiliar2[1,2:ncolX]
      Ypunto1<-auxiliar2[2:nfilasX,1]
      Ybarra11<-auxiliar2[-1,-1]

      ######Criterio para seleccionar el n?mero de componentes en la imputaci?n
      Ybarra11_Estand<-scale(Ybarra11)
      dvs_Ybarra11_Estand<-svd(Ybarra11_Estand)
      autovalores<-matrix(dvs_Ybarra11_Estand$d)**2
      sum_autovalores<-sum(autovalores)
      CritCalinski1<-(autovalores/sum_autovalores)
      CritCalinski2<-matrix(cumsum(CritCalinski1))
      CritCalinski3<-matrix(1,ncol(Ybarra11_Estand),1)
      for (Calins in 1:ncol(Ybarra11_Estand)){
        if (CritCalinski2[Calins,1]>0.70){CritCalinski3[Calins,1]=0}
      } # FIN de for (Calins in 1:ncol(Ybarra11_Estand)){
      f<-sum(CritCalinski3)+1
      #print(i1)
      #print(f)
      ####### FIN de Criterio para seleccionar el n?mero de componentes en la imputaci?n


      #print(Y1puntoTRANS)
      #print(Ypunto1)
      #print(Ybarra11)
      dvs.Ybarra11<-svd(Ybarra11)
      u_Ybarra11<-dvs.Ybarra11$u
      v_Ybarra11<-dvs.Ybarra11$v
      d_Ybarra11<-diag(dvs.Ybarra11$d)

      pred_estand<-Y1puntoTRANS%*%v_Ybarra11[,1:f]%*%MASS::ginv(d_Ybarra11[1:f,1:f])%*%t(u_Ybarra11[,1:f])%*%Ypunto1
      completed_matrix[posic_linha,posic_col]<-medcol[1,posic_col]+(pred_estand*desvest[1,posic_col])



    } ###FIN de for (i1 in 1:totalfaltantes){

    stabilitycrit<-(sum((completed_matrix-XMISSING)**2))/(sum((XMISSING*indicadora)**2))
    stabilityiter<-stabilityiter+1
    XMISSING<-completed_matrix
    Xmissingestand<-scale(completed_matrix)
    medcol<-t(matrix(colMeans(completed_matrix)))
    desvest<-sqrt(t(diag(var(completed_matrix))))

    #print(stabilityiter)
    #print(stabilitycrit)
    ######print (f)
    #print(completed_matrix)



  }### FIN de while (stabilitycrit>epsilon & stabilityiter<=100){

  list(NumeroIterGabriel=stabilityiter,CritConvergGabriel=stabilitycrit,GabrielImput=completed_matrix)
}### Fin function Gabriel.DI
